# NMFF-AFM: normal mode flexible fitting of proteins conformations to AFM images
#
# Copyright (c) 2024 TamaLab
# Authors: Xuan Wu, Nagoya University
#          Osamu Miyashita, RIKEN
#          Florence Tama, Nagoya University/RIKEN
#
# References:
# Modeling Conformational Transitions of Biomolecules from Atomic Force
# Microscopy Images using Normal Mode Analysis
# Xuan Wu, Osamu Miyashita, and Florence Tama
# https://doi.org/10.1021/acs.jpcb.4c04189
#
#
import os
import time
import shutil
import tomlkit
import argparse
import pandas as pd
from pathlib import Path

from module.modes_generator import make_block, run_rtb2, move_files_into_folder, generate_deformed_conformation, deformation_rtb
from module.images_generator import afm_generation
from module.slope_of_gradient import find_largest_slopes, process_data_tsv
from module.stop_iter import cc_judgement_avrg, create_log, cc_judgement_single
from module.rms_calculation import rmsd_calculation, process_rmsd_dictionary, rmsd_to_initial_and_reference
from module.scoring_and_export import scoring_conformations
from module.check_list import check_folders_and_print, confirmation, filename_of_next_iter, folder_of_next_iter, prepare_next_iter, summary_files


if __name__ == "__main__":
    print("\n=== NMFF-AFM initializing... ===\n")
    parser = argparse.ArgumentParser(description='AFM fitting by NMA')
    parser.add_argument('parameter_file', help='parameter file')
    args = parser.parse_args()
    with open(args.parameter_file, "r") as f:
        params = tomlkit.load(f)

    upper_folder_path = os.getcwd()
    initial_conformation_name = params['original_conformation']
    target_figure_name = params['target_conformation']

    combined_amplitude = params['combined_amplitude']
    res_x = params['res_x']
    res_y = params['res_y']
    res_z = params['res_z']
    size_x = params['size_x']
    size_y = params['size_y']
    start_from_this_mode = params['first_mode']
    stop_at_this_mode = params['last_mode']
    mode_selection = params['mode_selection']
    num_of_threads = params['how_many_threads']
    iteration_to_stop = num_iterations = params['num_iterations']  # Define the maximum number of iterations can be performed
    whether_calculate_rmsd_reference = params['calculate_rmsd_to_reference']
    use_which_figure = params['which_type_of_target']
    if whether_calculate_rmsd_reference == "yes":
        reference_pdb_name = params['file_name_of_reference_pdb']
    else:
        reference_pdb_name = ""
    probe_radius = params['radius_of_probe']
    probe_angle = params['angle_of_probe']

    # stop_mode indicating the way to stop the iteration. If using 'average', the iteration will be stopped by finding
    # the turning point of the average value of the CC from last five iterations, if using 'single', the iteration
    # will be stopped by finding the turning point of the CC from current and previous iteration
    termination_criterion = "numeric"

    folder_name = f"{initial_conformation_name}#s0"
    upper_folder_name = Path(upper_folder_path).stem
    overall_log_file_path = Path(upper_folder_path) / f"{upper_folder_name}_log.xlsx"
    initial_conformation_iter = initial_conformation_name
    stop_bool = True
    deformed_times = 0

    check_folders_and_print(
        upper_folder_path,
        initial_conformation_name,
        use_which_figure,
        target_figure_name,
        reference_pdb_name,
        whether_calculate_rmsd_reference,
        folder_name,
        overall_log_file_path,
        iteration_to_stop,
        combined_amplitude,
        res_x, res_y, res_z,
        size_x, size_y,
        start_from_this_mode,
        stop_at_this_mode,
        num_of_threads,
        probe_radius, probe_angle
    )

    number_of_modes = stop_at_this_mode - start_from_this_mode + 1

    if confirmation():
        print("=== NMFF-AFM will be started soon... ===\n")
        # Create ../All_conformation folder under the upper folder
        if not os.path.exists(Path(upper_folder_path) / "All_conformation"):
            summary_folder = Path(upper_folder_path) / "All_conformation"
            Path(summary_folder).mkdir()
            print("Summary folder created.")

        # Making the overall logfile ready
        create_log(overall_log_file_path, whether_calculate_rmsd_reference)

        # Create #s0 folder
        work_directory_iter = Path(upper_folder_path) / folder_name
        if not os.path.exists(Path(work_directory_iter)):
            Path(work_directory_iter).mkdir()
            print("First iteration folder #s0 created.")
        else:
            raise FileExistsError(f"The folder {folder_name} already exists.")

        # Copy initial PDB into #S0 folder, reference no need to be moved
        shutil.copy(Path(upper_folder_path) / f"{initial_conformation_name}.pdb", Path(work_directory_iter) / f"{initial_conformation_name}.pdb")
        print("Initial PDB file copied into #s0 folder.")

        for i in range(num_iterations + 10):
            deformation_directory = Path(work_directory_iter) / "1st"
            log_of_step = Path(work_directory_iter) / f"{initial_conformation_iter}.xlsx"
            target_figure_path = Path(upper_folder_path) / f"{target_figure_name}.{use_which_figure}"

            # This is the PDB file of reference
            reference_pdb = Path(work_directory_iter) / f"{target_figure_name}.pdb"
            moved_target_pdb_path = Path(upper_folder_path) / "All_conformation" / f"{target_figure_name}.pdb"

            # From run_nma, perform NMA calculation
            t1 = time.perf_counter(), time.process_time()
            make_block(work_directory_iter, initial_conformation_iter)
            run_rtb2(work_directory_iter, stop_at_this_mode)
            t2 = time.perf_counter(), time.process_time()
            print(f"NMA calculation: Real time: {t2[0] - t1[0]:.2f} seconds, CPU time: {t2[1] - t1[1]:.2f} seconds")

            # From modes_generator_rtb_100, deform the initial conformation along each mode
            move_files_into_folder(work_directory_iter, deformation_directory)
            generate_deformed_conformation(
                deformation_directory,
                initial_conformation_iter,
                combined_amplitude,
                start_from_this_mode,
                stop_at_this_mode
            )

            # From image_generator, generate the simulated AFM images
            t1 = time.perf_counter(), time.process_time()
            afm_generation(
                deformation_directory,
                res_x, res_y, res_z,
                size_x, size_y,
                num_of_threads,
                probe_radius, probe_angle)
            t2 = time.perf_counter(), time.process_time()
            print(f"AFM image simulation: Real time: {t2[0] - t1[0]:.2f} seconds, CPU time: {t2[1] - t1[1]:.2f} seconds")

            df_slopes, df_pcc = process_data_tsv(
                log_of_step,
                target_figure_path,
                initial_conformation_iter,
                start_from_this_mode,
                stop_at_this_mode,
                use_which_figure,
                deformation_directory
            )
            print("Slopes")
            pd.options.display.float_format = "{:.6g}".format
            df_slopes['abs_slope'] = df_slopes['slope'].apply(abs)
            print(df_slopes.sort_values(by=['abs_slope'], ascending=False))
            largest_mode, first_slope, second_largest_mode, second_slope = (
                find_largest_slopes(df_slopes)
            )
            print("Max CC points")
            print(df_pcc.iloc[0:3])

            if mode_selection == "slope":
                deform_amplitude_1 = first_slope / abs(first_slope) * combined_amplitude
                num_of_largest_mode = largest_mode

            elif mode_selection == "maxcc":
                num_of_largest_mode = df_pcc["mode"].iloc[0]
                deform_amplitude_1 = df_pcc["dq"].iloc[0]
            elif mode_selection == "maxcc_force_move":
                df_not_zero = df_pcc[df_pcc["dq"] != 0]
                if df_pcc["dq"].iloc[0] == 0:
                    print("Max CC points among dq != 0")
                    print(df_not_zero.iloc[0:3])
                num_of_largest_mode = df_not_zero["mode"].iloc[0]
                deform_amplitude_1 = df_not_zero["dq"].iloc[0]
            else:
                raise ValueError(
                    "Possible mode selection algorithms are 'slope' and 'maxcc'"
                )

            # Whether to stop the loop
            if termination_criterion == "average":
                stop_bool = cc_judgement_avrg(deformation_directory, initial_conformation_iter, target_figure_name,
                                              num_of_largest_mode, deform_amplitude_1, overall_log_file_path, use_which_figure)
            elif termination_criterion == "numeric":
                if deformed_times == iteration_to_stop or deform_amplitude_1 == 0:
                    cc_judgement_avrg(deformation_directory, initial_conformation_iter, target_figure_name,
                                      num_of_largest_mode, deform_amplitude_1, overall_log_file_path, use_which_figure)
                    stop_bool = False
                else:
                    cc_judgement_avrg(deformation_directory, initial_conformation_iter, target_figure_name,
                                      num_of_largest_mode, deform_amplitude_1, overall_log_file_path, use_which_figure)
                    stop_bool = True
            else:
                stop_bool = cc_judgement_single(deformation_directory, initial_conformation_iter, target_figure_name,
                                                num_of_largest_mode, deform_amplitude_1, overall_log_file_path,
                                                use_which_figure)

            # Copy all the initial files of the iteration into /All_conformation
            summary_files(upper_folder_path, deformation_directory, initial_conformation_iter)

            if not stop_bool:
                break

            # Prepare next iteration
            conf_name_next_iter = filename_of_next_iter(initial_conformation_iter)
            dir_name_next_itr = folder_of_next_iter(work_directory_iter)
            prepare_next_iter(dir_name_next_itr)
            print(
                "Deformation by mode and amplitude: ",
                num_of_largest_mode,
                deform_amplitude_1,
            )
            deformation_rtb(
                deformation_directory,
                Path(initial_conformation_iter).with_suffix(".pdb"),
                num_of_largest_mode,
                deform_amplitude_1,
                Path(dir_name_next_itr) / f"{conf_name_next_iter}.pdb",
            )

            work_directory_iter = folder_of_next_iter(work_directory_iter)
            initial_conformation_iter = filename_of_next_iter(initial_conformation_iter)

            deformed_times = deformed_times + 1

        # Calculate RMSD value basing on input files
        rmsd_to_initial_and_reference(
            upper_folder_path, initial_conformation_name, reference_pdb_name, overall_log_file_path, whether_calculate_rmsd_reference
        )

        # Scoring all the conformations and export the best one, then written into the log file.
        result_step = scoring_conformations(upper_folder_path, overall_log_file_path, whether_calculate_rmsd_reference)

        print("\n=== NMFF-AFM flexible fit-in finished ===\n")
        print("Result:\n")
        print(f"Best fitted conformation at s{result_step}\nThe conformation saved as {initial_conformation_name}#s{result_step}.pdb.")
    else:
        print("\n=== NMFF-AFM flexible fit-in stopped ===\n")
