# afm-fitting-nma
# Usage
To perform a flexible fit-in, follow these steps:
#### 1. Create a Directory: Begin by creating a new directory where all necessary files will be stored.
#### 2. Prepare Required Files: Place the following files into the newly created directory:
* PDB File: A PDB file representing the initial conformation. Filenames must not contain #, spaces, or any of the following special characters: @, !, $, %, ^, &, *, (, ), +, =, {, }, [, ], |, \, :, ;, ", ', <, >, ,, ?, /, or ~.
* TSV File: The TSV file containing the target figure data.
* TOML File: An input TOML file with the required configuration settings (detailed below).
#### 3. Optional: Reference Conformation for RMSD Calculation:
An RMSD curve will be generated default, plotting the RMSD between the deformed conformation at each step and the initial conformation. If you want to calculate RMSD relative to a reference conformation:
* Place an additional PDB file (the reference conformation) in the same directory.
* Set `calculate_rmsd_to_reference = 'yes'` and input the name of the reference PDB in the TOML file.
* This will add an extra curve to the resulting figure, representing the RMSD relative to the reference conformation.
#### 4. In the above directory and activate the visual environment with the required modules (detailed at Requirements section), using:
```
$ python path/to/nmff_afm.py parameter.toml
```

# Configuration
### Example file
This is an example of the parameter file.
```
num_iterations = 70
original_conformation = '4ake' # No need to add the suffix
target_conformation = '1ake' # No need to add the suffix
which_type_of_target = 'tsv'
combined_amplitude = 25
res_x = 0.5 # in nanometer
res_y = 0.5 # in nanometer
res_z = 0.64 # in Angstrom
size_x = 6.0 # This range represents half the size of the generated image in x-axis
size_y = 6.0 # This range represents half the size of the generated image in y-axis
first_mode = 1
last_mode = 16
mode_selection = 'slope' # ['slope', 'maxcc', 'maxcc_force_move']
how_many_threads = 8

# Whether to calculate the RMSD to a reference PDB file
calculate_rmsd_to_reference = 'yes' # ['yes', 'no']
file_name_of_reference_pdb = '1ake' # No need to add the suffix

# Shape of the probe
radius_of_probe = 2.0 # in nanometer
angle_of_probe = 10 # in degree
```

### Parameters in the Input File
* `num_iterations`: Int\
A number indicating how many times the input PDB file will be deformed.
* `original_conformation`: String\
A name of initial PDB file. The suffix will be added.
* `target_conformation`: String\
A name of target figure. The suffix will be added.
* `combined_amplitude`: Int\
A nmber of the overall amplitude used in the deformation.
* `res_x`: Float\
A number of the resolution used in simulating AFM images in x-axis.\
It uses Nanometer as the unit of measurement.
* `res_y`: Float\
A number of the resolution used in simulating AFM images in y-axis.\
It uses Nanometer as the unit of measurement.
* `res_z`: Float\
A number of the resolution used in simulating AFM images in z-axis.\
It uses Angstrom as the unit of measurement.
* `size_x`: Float\
A number of the range of the simulated AFM images in x-axis. This range starts from the negative half-axis and ends at the positive half-axis.\
It uses nanometers as the unit of measurement.
* `size_y`: Float\
A number of the range of the simulated AFM images in y-axis. This range starts from the negative half-axis and ends at the positive half-axis.\
It uses nanometers as the unit of measurement.
* `first_mode`: Int\
A number of the lowest frequency used in deforming the initial PDB file.
* `last_mode`: Int\
A number of the highest frequency used in deforming the initial PDB file.
* `mode_selection`: String\
A string for selecting the mode for selecting the largest mode.
By deafult, `slope`.
* `how_many_threads`: Int\
A number indicating how many threads using in the flexible fit-in.
* `calculate_rmsd_to_reference`: String\
The string indicating whether to calculate the RMSD to a reference conformation.
* `file_name_of_reference_pdb`: String\
A name of a PDB file used as the reference for RMSD calculation. The suffix will be added.
* `radius_of_probe`: Float\
A number indicating the radius of the simulated probe used generating the pseudo-AFM images.\
It uses nanometers as the unit of measurement.
* `angle_of_probe`: Float\
A number indicating the angle of the simulated probe used generating the pseudo-AFM images.

# Installation
## Requirements
Before running the calculation, several softwares should be installed first:

&emsp;&bull;&emsp;Pro-Fit: Available at [Pro-Fit](http://www.bioinf.org.uk)

&emsp;&bull;&emsp;afmize: Available at [afmize](https://github.com/ToruNiina/afmize)

&emsp;&bull;&emsp;NMA scripts: Available at [NMA scripts](https://github.com/TamaLab/nma)

To install all required Python modules, run:
```
pip install -r requirements.txt
```

## Environment variable
Before running a flexible fit-in, add the following environment variables to your system:
```
# Path the excutable file of afmize
export AFMIZE_PATH=/to/your/afmize_directory/afmize

# Path of the folder containing NMA program and scripts. rtb2, makebloc.pl and movemode.pl should be saved in this folder.
export NMA_FOLDER=/to/your/nma_directory

# Path the excutable file of Pro-Fit
export PRO_FIT_PATH=/to/your/profit_directory/profit
```
