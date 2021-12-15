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
python run_desmond_builder.py -top topology_file.xxx -traj trajectory_file.xxx -l "ligand_selection" -p "protein_selection" -o contact_surface_data
```