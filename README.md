# Desmond Builder
 Prepare a system for MD simulation in desmond starting from a .mae file

## Features

It allows to prepare a system for molecular dynamics simulation using the
Desmond engine. The input file it should be a .mae file that contains proteins
and/or molecules as ligands, cofactors, metals, etc.
the script receives a configuration file with the options to consider to build 
the system and to carry out the equilibration and MD protocols.

## Requirements

* Python 3.0 or higher versions.

* Desmond software. 
It does not require the paid version of Schrödinger but it can also be used.

## Usage

```
python run_desmond_builder.py -i config.dat
```

## Options

The configuration file must contain 3 headers:
* [settings]  
Contains global settings.

* [build_geometry]  
Contains the options for building the system

* [protocol]  
Contains options for molecular dynamics protocols.

[settings]
* workdir < Working directory name >
Acceptable values: any name
Default values: md_run
It could be a directory o a subfolder "directory1/directory2".

* file < Specify the path to the input .mae file >
Acceptable values: .mae file
Default values: None

* desmond_path < Specify the path to the schrodinger installation >
Acceptable values: any path
Default values: $SCHRODINGER

[build_geometry]
* counterions < Add counterions? >
Acceptable values: yes, true, on or no, false, off
Default values: no
If yes, the number of counterions will be calculated according to the
total charge of the system.
The default positive ion is Na and the default negative ion is Cl.
The default options could be changed in the configuration file as:
* counterions_positive_ion < Positive ion >
Acceptable values: Na, Li, K, Rb, Cs
* counterions_negative_ion < Negative ion >
Acceptable values: Cl, F, Br, I.
Note: The available ions to be used as conterions are different to the ions
to be used as Salt.

* ions_away < Add ions away from structure >
Acceptable values: yes, true, on or no, false, off
Default values: no
If yes, it should be added the options: ions_awaydistance and ions_awayfrom.

* ions_awaydistance < Minimum distance between ions and structure >
Default values: 5.0

* ions_awayfrom < Add ions away from structure >
Default values: protein
It should be an acceptable ASL selection.


## License

Licensed under the MIT license, see the separate LICENSE file.