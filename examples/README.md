# Example Files
This example directory includes files for running a flexible fit-in simulation with or without a reference conformation.\
In `/example/with_reference`:
 * Initial conformation: `4ake_x150y330.pdb`
 * Target figure: `1ake_x150y330.tsv`
 * Reference conformation: `1ake_x150y330.pdb`
 * Parameter file: `parameter_ake_reference.toml`

In `/example/non_reference`:
 * Initial conformation: `4ake_x150y330.pdb`
 * Target figure: `1ake_x150y330.tsv`
 * Parameter file: `parameter_ake_non_reference.toml`

# Usage
Running the calculation, in the folder contains all the files, using:
```
$ python path/to/nmff_afm.py parameter_ake_reference.toml
or
$ python path/to/nmff_afm.py parameter_ake_non_reference.toml
```
