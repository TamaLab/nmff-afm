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
import shutil
import subprocess
from tqdm import tqdm
from pathlib import Path

nma_folder = os.environ.get('NMA_FOLDER')


def make_block(work_path, initial_conformation_name):
    """
    This function creates blocks using the original conformation PDB file
    :param work_path: The path where the calculation will perform
    :param initial_conformation_name: The name of the initial conformation
    :return: None
    """
    # Step 1: Make Blocks
    cmd_mk_bloc = f"cd {work_path} && {nma_folder}/makebloc.pl {initial_conformation_name}.pdb > pdb"
    subprocess.run(cmd_mk_bloc, shell=True, check=True)

    print("makebloc performed")


def run_rtb2(work_path, num_to_compute):
    """
    Creates an input file for the RTB calculation and executes the RTB2 command
    :param work_path: The path to the working directory where the `rtb.inp` input file will be created and the
        RTB2 command will be executed
    :param num_to_compute: How many frequencies will be calculated
    :return: None
    """
    with open(Path(work_path) / "rtb.inp", "w") as file:
        file.write(" &inputs\n")
        file.write("   cutoff = 8.00,\n")
        file.write("   ncv = 60,\n")
        file.write("   tol = 1e-18\n")
        file.write("   nstep = 0,\n")
        file.write(f"   nvec = {num_to_compute}\n")
        file.write(" /&\n")
        file.write("\n")

    cmd_rtb = f"cd {work_path} && {nma_folder}/rtb2"
    subprocess.run(cmd_rtb, shell=True, check=True)
    print("rtb2 performed, modes generated")


def mode_formatting(mode_in_single):
    """
    This function is used to fit the new version of rtb, the file names changed from mov0.modX to mov000.modXXX
    :param mode_in_single: The number of the mode file will be used
    """
    formatted_mode = f"{int(mode_in_single):03d}"
    return formatted_mode


def deformation_rtb(work_folder, initial_pdb, frequency, amplitude, deformed_pdb):
    """
    This function will deform the PDB file along one mode and amplitude
    :param work_folder: The path to the folder where perform the deformation
    :param initial_pdb: ***PATH*** of the initial PDB file
    :param frequency: The number of the frequency used in the deformation
    :param amplitude: The value of the amplitude used in the deformation
    :param deformed_pdb: ***PATH*** of the deformed PDB file
    """
    pdb_org = Path(work_folder) / initial_pdb
    pdb_deformed = Path(work_folder) / deformed_pdb
    mode = Path(work_folder) / f"mov000.mod{mode_formatting(frequency)}"
    cmd = f"{nma_folder}/movemode.pl {pdb_org} {mode} {amplitude} > {pdb_deformed}"
    result = subprocess.run(cmd, shell=True, check=True)


def generate_deformed_conformation(work_directory, initial_conformation_name, amplitude, start_mode, end_mode):
    """
    This function will deform PDB along all the modes and amplitudes
    :param work_directory: The path where to perform the deformation
    :param initial_conformation_name: ***NAME*** of the initial conformation
    :param amplitude: The amplitude used in the deformation
    :param start_mode: Number of the mode from which used in the deformation
    :param end_mode: Number of the mode till which used in the deformation
    """
    assert amplitude > 0
    frequencies = []
    for num in range(start_mode, end_mode + 1, 1):
        frequencies.append(str(num))

    one_half = int(amplitude / 2)  # enforce to be integer for name formatting
    print('Generating deformed PDBs')
    for frequency in tqdm(frequencies):
        for dq in (-amplitude, -one_half, 0, one_half, amplitude):
            output_pdb = Path(work_directory) / f"{frequency}#{dq}.pdb"
            if os.path.exists(output_pdb):
                os.remove(output_pdb)
            initial_pdb = Path(initial_conformation_name).with_suffix(".pdb")
            deformation_rtb(work_directory, initial_pdb, frequency, dq, output_pdb)


def move_files_into_folder(work_path, sub_folder_path):
    """
    This function will move all the files after NMA into a folder named 1st, where the other steps will be performed
    :param work_path: The folder path contains all the generated files
    :param sub_folder_path: The folder newly created, named /1st
    """
    Path(Path(work_path) / Path(sub_folder_path)).mkdir()
    print("new folder created.")

    # Move files in the parent folder into /1st
    file_list = os.listdir(work_path)
    for file_name in file_list:
        if os.path.isfile(os.path.join(work_path, file_name)):
            shutil.move(os.path.join(work_path, file_name), os.path.join(sub_folder_path, file_name))
    print("files moved.")
