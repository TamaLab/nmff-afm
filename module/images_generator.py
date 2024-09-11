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
import glob
import subprocess
from pathlib import Path
from multiprocessing import Pool

afmize_path = os.environ.get('AFMIZE_PATH')


def text_generation(file_dir, res_x, res_y, res_z, range_x, range_y, probe_radius, probe_angle):
    """
    Generate a blank pattern for keywords replacement
    :param file_dir: The path where the pattern generated
    :param res_x: The resolution in x-axis of the simulated AFM figure
    :param res_y: The resolution in y-axis of the simulated AFM figure
    :param res_z: The resolution in z-axis of the simulated AFM figure
    :param range_x: The range of a half of the simulated AFM figure in x-axis
    :param range_y: The range of a half of the simulated AFM figure in y-axis
    :param probe_radius: The radius of the sampling probe
    :param probe_angle: The angle of the sampling probe
    :return: None
    """
    text = [
        'file.input           = "dbInput.pdb"',
        'file.output.basename = "dbOutput"',
        'file.output.formats  = ["tsv", "svg"]',
        'probe.size           = {radius = "' + f'{probe_radius}nm' + '", angle = ' + f'{probe_angle}' + '}',
        f'resolution.x         = "{res_x}nm"',
        f'resolution.y         = "{res_y}nm"',
        'resolution.z         = "' + f'{res_z}angstrom' + '"',
        f'range.x              = ["-{range_x}nm", "{range_x}nm"]',
        f'range.y              = ["-{range_y}nm", "{range_y}nm"]',
        'scale_bar.length     = "0.0nm"',
        "stage.align          = true",
        "stage.position       = 0.0",
        'noise                = "0.0nm"',
    ]
    f = open(Path(file_dir) / "blank_toml.txt", "w")
    for i in text:
        f.write(str(i) + "\n")
    f.close()


def toml_generation(file_dir, name_data, afmize_input_name):
    """
    Replacing the keywords in line by using dictionary
    :param file_dir: The path to the blank .toml exists
    :param name_data: The dictionary contain keywords and value of the input and output file name of simulated figures
    :param afmize_input_name: The name of the input name for afmize
    :return: None
    """
    f = open(file_dir / 'blank_toml.txt', 'r+')
    output = open(file_dir / afmize_input_name, 'w')
    for line in f:
        for f_key, f_value in name_data.items():
            if f_key in line:
                line = line.replace(f_key, f_value)
        output.write(str(line))
    output.close()


def filename_dict_generation(pdb_name):
    """
    Generate the dictionary for following replacement
    :param pdb_name: The full path without suffix of input PDB file the keywords for dictionary
    :return: Generated dictionary
    """
    list_key = ["dbInput", "dbOutput"]
    list_pdb_name = [pdb_name, pdb_name]
    dict_replace = dict(zip(list_key, list_pdb_name))
    return dict_replace


def run_afmize(pdb_full_name):
    """
    This function will generate the simulated AFM images
    :param pdb_full_name: The name with suffix of the PDB file will be used to generate simulated AFM images
    :return: None
    """
    print("running afmize for:", pdb_full_name)
    new_folder = Path(pdb_full_name).parent
    pdb_name = pdb_full_name.split(".")[0]
    dict_replace = filename_dict_generation(pdb_name)
    afmize_input_name = f"{Path(pdb_name).stem}_gen_image.toml"
    toml_generation(new_folder, dict_replace, afmize_input_name)
    command = f"cd {new_folder} && {afmize_path} {afmize_input_name}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print("=== Standard Output ===")
        print(result.stdout)
        print("=== Standard Error ===")
        print(result.stderr)
        raise RuntimeError("An error occurred during the command execution: {command}")


def afm_generation(figure_folder, resolution_x, resolution_y, resolution_z, range_x, range_y,
                   num_of_threads,
                   radius_of_probe, angle_of_probe):
    """
    This function will generate the input file of afmize and then generate simulated AFM images
    :param figure_folder: The path to the folder where to generate the simulated figures
    :param resolution_x: The resolution used for the simulated figures in x-axis
    :param resolution_y: The resolution used for the simulated figures in y-axis
    :param resolution_z: The resolution used for the simulated figures in z-axis
    :param range_x: The size of the simulated figure in x-axis
    :param range_y: The size of the simulated figure in y-axis
    :param num_of_threads: The number of the threads will be used
    :param radius_of_probe: The radius of the sampling probe
    :param angle_of_probe: The angle of the sampling probe
    :return: None
    """

    text_generation(figure_folder, resolution_x, resolution_y, resolution_z, range_x, range_y, radius_of_probe, angle_of_probe)
    files = glob.glob(os.path.join(figure_folder, "*.pdb"))
    num_processes = num_of_threads
    with Pool(processes=num_processes) as pool:
        pool.map(run_afmize, files)
