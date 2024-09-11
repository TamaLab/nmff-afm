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
import sys
import shutil
from pathlib import Path

nma_folder = os.environ.get('NMA_FOLDER')
afmize_path = os.environ.get('AFMIZE_PATH')
profit_path = os.environ.get('PRO_FIT_PATH')


def check_env_vars(env_vars):
    """
    This function will check if specified environment variables are set
    :param env_vars: List of environment variable names to check.
    """
    for var in env_vars:
        if os.environ.get(var) is None:
            print(f"Error: Environment variable {var} is not set.")
            sys.exit(1)  # Exit the script with a non-zero status

    print("Environment variables checked.")


def check_package_and_scripts():
    """
    This function will check the related packages and scripts exists
    :return: None
    """
    nma_related = ["makebloc.pl", "rtb2", "movemode.pl"]
    for nma_file in nma_related:
        path_for_check = Path(nma_folder) / nma_file
        if not os.path.exists(path_for_check):
            print(f"{nma_file} not found, please check the exist of it.")

    path_afmize = Path(afmize_path)
    if not os.path.exists(path_afmize):
        print(f"Package afmize not found, please check the exist of it.")

    path_profit = Path(profit_path)
    if not os.path.exists(path_profit):
        print(f"Package Pro-Fit not found, please check the exist of it.")

    print("All related packages and scripts checked.")


def check_folders_and_print(
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
):
    """
    This function will check the existence of specified folders and prints status messages
    :param upper_folder_path: Path to the folder contains all the files
    :param initial_conformation_name: Name of the initial conformation file
    :param use_which_figure: Identifier for which figure to use
    :param target_figure_name: Name of the target figure file
    :param reference_pdb_name: Name of the reference PDB file
    :param whether_calculate_rmsd_reference: Flag indicating whether to calculate RMSD to reference PDB
    :param folder_name: Name of the folder for first step of the iteration
    :param overall_log_file_path: Path to the overall log file
    :param iteration_to_stop: Iteration number at which to stop the process
    :param combined_amplitude: Combined amplitude value
    :param res_x: Resolution of the simulated AFM figures in x-axis
    :param res_y: Resolution of the simulated AFM figures in y-axis
    :param res_z: Resolution of the simulated AFM figures in z-axis
    :param size_x: Size in the X dimension
    :param size_y: Size in the Y dimension
    :param start_from_this_mode: Number of mode from which to start
    :param stop_at_this_mode: Number of mode at which to stop
    :param num_of_threads: Number of threads to use for processing
    :param probe_radius: The radius of the simulated probe
    :param probe_angle: The angle of the simulated probe
    :returns: None
    """

    env_vars_to_check = ['AFMIZE_PATH', 'NMA_FOLDER', 'PRO_FIT_PATH']
    check_env_vars(env_vars_to_check)

    # Check if the related packages and scripts exist
    check_package_and_scripts()

    # Check whether the folder exists
    if os.path.exists(Path(upper_folder_path)) and os.path.isdir(Path(upper_folder_path)):
        print(f"Folder '{upper_folder_path}' checked.")
    else:
        raise FileExistsError(f"The folder '{upper_folder_path}' does not exist, please check the input path.")

    # Check whether the Initial PDB exists
    if os.path.exists(Path(upper_folder_path) / f"{initial_conformation_name}.pdb"):
        print(f"Initial PDB {initial_conformation_name} checked.")
    else:
        raise FileExistsError(
            f"Initial PDB file {initial_conformation_name} does not exist, please check the input file.")

    # Check whether # in the name of Initial PDB
    if '#' in initial_conformation_name:
        raise ValueError(
            f"Initial PDB file '{initial_conformation_name}' contains a '#'. To avoid potential issues,"
            f"please rename the file without '#' character."
        )
    else:
        print(f"Initial PDB filename '{initial_conformation_name}' checked.")

    # Check whether the Target figure exists
    if not os.path.exists(Path(upper_folder_path) / f"{target_figure_name}.{use_which_figure}"):
        raise FileExistsError(
            f"The Target {use_which_figure} file {target_figure_name} does not exist, please check the input file.")
    else:
        print(f"Target figure {target_figure_name}.{use_which_figure} checked.")

    # Check whether the Reference PDB exists
    if whether_calculate_rmsd_reference == "yes":
        if not os.path.exists(Path(upper_folder_path) / f"{reference_pdb_name}.pdb"):
            raise FileExistsError(
                f"Reference PDB file {reference_pdb_name}.pdb does not exist, please check the input file.")
        else:
            print(f"Reference PDB {reference_pdb_name}.pdb checked.")

    # Check whether the /All_conformation exists
    if os.path.exists(Path(upper_folder_path) / "All_conformation"):
        raise FileExistsError("Folder /All_conformation already exists, please check if it is empty and remove it.")

    # Check whether the overall logfile exists
    if os.path.exists(Path(overall_log_file_path)):
        raise FileExistsError(f"Log {overall_log_file_path} already exists.")

    # Check whether #s0 folder exists
    work_directory_iter = Path(upper_folder_path) / folder_name
    if os.path.exists(Path(work_directory_iter)):
        raise FileExistsError(f"{folder_name} already exists.")

    check_mark = ' \u2713'
    underline = '\033[4m'  # ANSI escape code for underline
    reset = '\033[0m'  # ANSI escape code to reset formatting
    degree_sign = "\u00B0"  # degree sign
    print("\n=== Check List ===\n")
    print(f"Work directory: {underline}{upper_folder_path}{reset}", check_mark)
    print(f"Initial PDB file: {underline}{initial_conformation_name}{reset}", check_mark)
    print(f"Target figure: {underline}{target_figure_name}{reset}", check_mark)
    print(f"Use which format as target figure: {underline}{use_which_figure}{reset}")
    print(f"How many steps will be performed: {underline}{iteration_to_stop}{reset}")
    print(f"Overall amplitude used in deformation: {underline}{combined_amplitude}{reset}")
    print(f"Resolution in x-axis of the AFM images: {underline}{res_x}{reset}nm")
    print(f"Resolution in y-axis of the AFM images: {underline}{res_y}{reset}nm")
    print(f"Resolution in z-axis of the AFM images: {underline}{res_z}{reset}Angstrom")
    print(f"The size of the AFM image will be: {underline}{size_x * 2}{reset}nm by {underline}{size_y * 2}{reset}nm")
    print(
        f"Frequency from {underline}{start_from_this_mode}{reset} to {underline}{stop_at_this_mode}{reset} will be used in the calculation.")
    print(f"{underline}{num_of_threads}{reset} threads will be used for AFM image simulation.\n")
    print(f"Calculate the RMSD-Reference value: {underline}{whether_calculate_rmsd_reference}{reset}")
    if whether_calculate_rmsd_reference == "yes":
        print(f"PDB file named {underline}{reference_pdb_name}{reset} will be used to calculate RMSD-Reference.\n")
    print(f"Radius of the probe in the simulated figure: {underline}{probe_radius}{reset}nm")
    print(f"Angle of the probe in the simulated figure: {underline}{probe_angle}{reset}{degree_sign}")


def confirmation():
    """
    This function will confirm whether to start the AFM-NMA flexible fit-in
    :return: None
    """

    while True:
        user_input = input(
            f"Start flexible fit-in with above parameters? [yes/no]: ").strip().lower()

        if user_input == 'yes':
            print("Starting the flexible fit-in...")
            return True
        elif user_input == 'no':
            print("Program stopped.")
            return False
        else:
            print("Invalid input. Please enter 'yes' to continue or 'no' to stop.")


def filename_of_next_iter(previous_name):
    """
    This function will generate the **NAME** of the next iteration, usually for folder name and file name
    :param previous_name: The name of the current folder
    :return: The name for next iteration
    """
    if "#" in previous_name:
        temp = previous_name.split("#s")[0]
        postfix = int(previous_name.split("#s")[1]) + 1
        new_name = f"{temp}#s{str(postfix)}"
    else:
        new_name = f"{previous_name}#s1"
    return new_name


def folder_of_next_iter(previous_dir):
    """
    This function will prepare the **PATH** of the folder for next iteration
    :param previous_dir: The path to the current folder
    :return: The new path for next iteration
    """
    path_of_folder = Path(previous_dir).stem
    new_folder_name = filename_of_next_iter(str(path_of_folder))
    temp = Path(previous_dir).parent
    new_path = temp / new_folder_name
    return str(new_path)


def prepare_next_iter(next_iter_path):
    """
    This function will prepare a folder for the next iteration, creating and naming it, then put the PDB files into it
    :param next_iter_path: The path to the folder for the next step
    :return: None
    """
    Path(next_iter_path).mkdir()
    print("New folder has been created.")


def summary_files(upper_path, deform_path, initial_name):
    """
    This function will copy the generated figures into All_conformation folder
    :param upper_path: The path to the folder contains all the subfolders
    :param deform_path: The path to the folder contains the PDB files
    :param initial_name: The name of the initial conformation
    :return: None
    """
    # Copy the initial files of this iteration into ../All_conformation
    shutil.copy(Path(deform_path) / f"{initial_name}.svg",
                Path(upper_path) / "All_conformation" / f"{initial_name}.svg")
    shutil.copy(Path(deform_path) / f"{initial_name}.tsv",
                Path(upper_path) / "All_conformation" / f"{initial_name}.tsv")
    shutil.copy(Path(deform_path) / f"{initial_name}.pdb",
                Path(upper_path) / "All_conformation" / f"{initial_name}.pdb")
