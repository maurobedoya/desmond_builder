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
* typing-extensions: ```pip install typing-extensions```
* Desmond software. 
It does not require the Schrödinger´s full version but it can also be used.

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

## [settings]
* ``workdir``: < Working directory name >  
Acceptable values: any name  
Default values: md_run  
It could be a directory or a subfolder "directory1/directory2".  

* ``file``: < Specify the path to the input .mae file >  
Acceptable values: .mae file  
Default values: None  

* ``desmond_path``: < Specify the path to the Desmond installation >  
Acceptable values: any path  
Default values: $SCHRODINGER  

## [build_geometry]
* ``counterions``: < Add counterions? >  
Acceptable values: yes, true, on or no, false, off  
Default values: no  
If yes, the number of counterions will be calculated according to the
total charge of the system.  
The default positive ion is 'Na' and the default negative ion is 'Cl'.  
The default options could be changed in the configuration file as:  
* ``counterions_positive_ion``: < Positive ion >  
Acceptable values: Na, Li, K, Rb, Cs  
* ``counterions_negative_ion``: < Negative ion >  
Acceptable values: Cl, F, Br, I.  
Note: The available ions to be used as conterions are different to the ions
to be used as Salt.  

* ``ions_away``: < Add ions away from structure >  
Acceptable values: yes, true, on or no, false, off  
Default values: no  
If yes, it should be added the options: ions_awaydistance and ions_awayfrom.  

* ``ions_awaydistance``: < Minimum distance between ions and structure >  
Default values: 5.0  

* ``ions_awayfrom``: < Add ions away from structure >  
Default values: protein  
It should be an acceptable ASL selection.  

* ``shape``: < Add box shape >  
Acceptable values: orthorhombic  
Default values: orthorhombic  
Different box shapes will be implemented in the future.  

* ``size``: < Add box size >  
Acceptable values: Three numbers separated by blank spaces.  
Default values: 10.0 10.0 10.0  

* ``size_type``: < Add box size type >  
Acceptable values: buffer, absolute  
Default values: buffer  

* ``salt``: < Add salt? >  
Acceptable values: yes, true, on or no, false, off  
Default values: no  
If yes, it should be added the option: concentration.  

* ``concentration``: < Salt concentration (M) >  
Default values: 0.15  

Optional settings:  
* ``positive_ion``: < Positive ion for salt >  
Acceptable values: Li, Na, K, Rb, Cs, Mg2, Ca2, Zn2, fe2, fe3  
Default values: Na  

* ``negative_ion``: < Negative ion for salt >  
Acceptable values: Cl, F, Br, I  
Default values: Cl  

* ``solvent``: < Solvent >  
Acceptable values: SPC, TIP3P, TIP4P, TIP4PEW, TIP5P, TIP4PD, DMSO, METHANOL, OCTANOL  
Default values: SPC  

## [protocol]  

The default protocols for desmond relaxation and the production are activated.  
To deactivate any step, the option should be set to 'no', 'off' or 'false'.  

* ``stage{x}``: < Relaxation stage for x=1,2,3,4,5 >  
Acceptable values: yes, true, on or no, false, off  
Default values: yes  

* ``stage{x}_time``: < Relaxation stage time (ps) for x=1,2,3,4,5 and production_time >  
Default values:  
``stage1_time``     = 100  
``stage2_time``     = 12  
``stage3_time``     = 12  
``stage4_time``     = 12  
``stage5_time``     = 24  
``production_time`` = 100000 ps.   

* ``stage{x}_temp``: < Relaxation stage temperature (K) for x=1,2,3,4,5 and production_temp >  
Default values:  
``stage1_temp``     = 10.0  
``stage2_temp``     = 10.0  
``stage3_temp``     = 10.0  
``stage4_temp``     = 300.0  
``stage5_temp``     = 300.0  
``production_temp`` = 300.0  

* ``stage{x}_ensemble``: < Relaxation stage ensemble for x=1,2,3,4,5 and production_ensemble >  
Acceptable values: NVE, NVT, NPT  
Note: NPAT and NPYT will be implemented in the future.  
Default values:  
``stage1_ensemble``     = NVT  
``stage2_ensemble``     = NVT  
``stage3_ensemble``     = NPT  
``stage4_ensemble``     = NPT  
``stage5_ensemble``     = NPT  
``production_ensemble`` = NPT  

* ``stage{x}_method``: < Relaxation stage method for x=1,2,3,4,5 and production_method >  
Acceptable values: NVT, NPT  
Default values:  
``stage1_method``     = Brownie  
``stage2_method``     = Berendsen  
``stage3_method``     = Berendsen  
``stage4_method``     = Berendsen  
``stage5_method``     = Berendsen  
``production_method`` = MTK  

* ``stage{x}_thermostat_tau``: < Relaxation stage thermostat tau (ps) for x=1,2,3,4,5 and production_thermostat_tau >  
Default values:  
``stage1_thermostat_tau`` = 0.1  
``stage2_thermostat_tau`` = 0.1  
``stage3_thermostat_tau`` = 0.1  
``stage4_thermostat_tau`` = 0.1  
``stage5_thermostat_tau`` = 0.1  
``production_thermostat_tau`` = 1.0  

* ``stage{x}_barostat_tau``: < Relaxation stage barostat tau (ps) for x=1,2,3,4,5 and production_barostat_tau >  
Default values:  
``stage1_barostat_tau``     = None  
``stage2_barostat_tau``     = None  
``stage3_barostat_tau``     = 50.0  
``stage4_barostat_tau``     = 20.0  
``stage5_barostat_tau``     = 2.0  
``production_barostat_tau`` = 2.0  

### Restraints

Note: For now there is only support for positional, distance, angle and improper harmonic restraints. Flat-bottomed restraints are not currently supported.  

Each type of restraint is defined with the preffix of "``stage{x}_restraints``" (for restraints from stage1 to stage5, x = 1,2,3,4,5) and for the restraints of additional stages the preffix "``additional_stage_restraints``" is used. The terminal of the name is the type of restraint: For "``positional``", "``distance``", "``angle``" or "``improper``" the terminals are "``pos``", "``dist``", "``ang``" and "``imp``" respectively.  

* ```stage{x}_restraints_number_{type}```: < number of restraints in stage "x" (x=1,2,3,4,5) and "type" (type = ``pos``, ``dist``, ``ang`` and ``imp``)>  
Default values:  
``stage{x}_restraints_number_pos``     = 1 (for x=1,2,3,4)  
``stage5_restraints_number_pos``       = 0  
``stage{x}_restraints_number_{type}``  = 0 (for x=1,2,3,4,5 and type = dist, ang and imp)  

* ``stage{x}_restraints_atoms_{type}``: < atoms (ASL) in stage "x" (x=1,2,3,4,5) and "type" (type = pos, dist, ang and imp) >  
Default values:  
``stage{x}_restraints_atoms_pos``       = solute_heavy_atom (for for x=1,2,3,4,5)  
``stage{x}_restraints_atoms_{type}``    = "None" (for x=1,2,3,4,5 and type = dist, ang and imp)  

* ``stage{x}_restraints_forces_{type}``: < Restraint_force (kcal·mol-1·Å-2) in stage "x" (x = 1,2,3,4,5) and "type" (type = pos, dist, ang and imp) >  
Default values:  
``stage{x}_restraint_forces_pos``      = 50.0 (for x=1,2,3,4)   
``stage5_restraints_forces_pos``       = 0.0  
``stage{x}_restraints_forces_{type}``  = "None" (for x=1,2,3,4,5 and type = dist, ang and imp)  

* ``stage{x}_restraints__{c}_{type}``: < Equilibrium constant in stage "x" (x = 1,2,3,4,5), "c" (r0, theta0, phi0) and "type" (type = pos, dist, ang and imp) >  
Default values:  
``stage{x}_restraint_r0_dist``       = "None" (for x=1,2,3,4,5)  
``stage{x}_restraints_theta0_ang``   = "None" (for x=1,2,3,4,5)  
``stage{x}_restraints_phi0_imp``     = "None" (for x=1,2,3,4,5)  

### Additional stages

* ``additional_stages``: < Number of additional stages to run >  
Acceptable values: integer number  
Default values: None  

* ``additional_stages_time``: < Time for additional stages (ps) >  
Acceptable values: integer number in ps.  
If only one value is given, it will be used for all stages, otherwise it should be a list of numbers separated by comma. The list of number should be the same length as the number of additional stages.
Default values: None  

* ``additional_stages_temp``: < Temperature for additional stages (K) >  
Acceptable values: integer number in K.  
If only one value is given, it will be used for all stages, otherwise it should be a list of numbers separated by comma. The list of number should be the same length as the number of additional stages.  
Default values: None  

* ``additional_stages_ensembles``: < Ensemble for additional stages >  
Acceptable values: NVT, NPT  
If only one value is given, it will be used for all stages, otherwise it should be a list of strings separated by comma. The list of strings should be the same length as the number of additional stages.  
Default values: None  

* ``additional_stages_methods``: < Method for additional stages >  
Acceptable values: Berendsen, Langevin  
If only one value is given, it will be used for all stages, otherwise it should be a list of strings separated by comma. The list of strings should be the same length as the number of additional stages.  
Default values: None  

* ``additional_stages_thermostat_tau``: < Thermostat tau for additional stages (ps) >  
Acceptable values: integer number in ps.  
If only one value is given, it will be used for all stages, otherwise it should be a list of numbers separated by comma. The list of number should be the same length as the number of additional stages.  
Default values: 0.1  

* ``additional_stages_barostat_tau``: < Barostat tau for additional stages (ps) >  
Acceptable values: integer number in ps.  
If only one value is given, it will be used for all stages, otherwise it should be a list of numbers separated by comma. The list of number should be the same length as the number of additional stages.  
Default values: 2.0  

### Additional stages restraints
Note: The notation for the additional stages restraints are similar to the stage{x} restraints. See above.
The difference with respect to the restraints of the stages is that in this case lists of restraints must be passed in the option: "``additional_stage_restraints_number_{type}``" according to the amount of "``additional_stages``" defined.
If three additional stages are defined and are required two positional restraints in each additional step and in the last step just one positional restraint, it should be added like this:  
``additional_stages`` = 3  
``additional_stage_restraints_number_pos`` = 2,2,1  

Additionally, it could be added the different type of restraints for each additional stage:  
e.g. ``additional_stage_restraints_number_dist`` = 1,2,1  

* ``additional_stages_restraints_number_{type}``: < Number of restraints for additional stages (int) for "type" (type = pos, dist, ang and imp) >  
Acceptable values: Any int number sepparated by comma according to number of "additional_stages".  
Default values:  
``additional_stages_restraints_number_{type}`` = "None"  

* ``additional_stages_restraints_forces_{type}``: < atoms (ASL) for additional stages (kcal·mol-1·Å-2) for "type" (type = pos, dist, ang and imp) >  
Acceptable values: Float number sepparated by comma according to number of "additional_stages".  
Note: The number of selections depends on the type of restraint. "pos" requires one atom selection, "dist" requires two atom selections, "ang" requires three atom selections and "imp" requires four atom selections.
Default values: 
``additional_stages_restraints_forces_{type}`` = "None"

* ``additional_stages_restraints_forces_{type}``: < Restraint forces for additional stages (kcal·mol-1·Å-2) for "type" (type = pos, dist, ang and imp) >  
Acceptable values: Float number sepparated by comma according to number of "additional_stages".  
Default values: 
``additional_stages_restraints_forces_{type}`` = "None"

* ``production_cutoff``: < Production cutoff (Å) >  
Default values: 9.0  

* ``production_timestep_bonded``: < Production timestep bonded (ps) >  
Default values: 0.002  

* ``production_timestep_near``: < Production timestep near (ps) >  
Default values: 0.002  

* ``production_timestep_far``: < Production timestep far (ps) >  
Default values: 0.006  

* ``production_traj_frames_per_file``: < Production: frames per-file >  
Default values: 250  

* ``production_traj_interval``: < Recording interval for trajectory (ps) >  
Default values: 50.0  

* ``production_pressure``: < Production pressure (atm) >  
Default values: 1.01325  

* ``run_preparation``: < Run preparation stage? >  
Acceptable values: True, yes, on, or False, no, off.  
Default values: False  

* ``run_protocol``: < Run MD protocols? >  
Acceptable values: True, yes, on, or False, no, off.  
Default values: False  

## License

Licensed under the MIT license, see the separate LICENSE file.
