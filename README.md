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
The default positive ion is 'Na' and the default negative ion is 'Cl'.  
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

* shape < Add box shape >
Acceptable values: orthorhombic  
Default values: orthorhombic  
Different box shapes will be implemented in the future.  

* size < Add box size >
Acceptable values: Three numbers separated by blank spaces.  
Default values: 10.0 10.0 10.0

* size_type < Add box size type >
Acceptable values: buffer, absolute  
Default values: buffer  

* salt < Add salt >
Acceptable values: yes, true, on or no, false, off  
Default values: no  
If yes, it should be added the option: concentration.  

* concentration < Salt concentration (M) >
Default values: 0.15  

Optional settings:  
* positive_ion < Positive ion > 
Acceptable values: Li, Na, K, Rb, Cs, Mg2, Ca2, Zn2, fe2, fe3  
Default values: Na  

* negative_ion < Negative ion >
Acceptable values: Cl, F, Br, I  
Default values: Cl  

* solvent < Solvent >
Acceptable values: SPC, TIP3P, TIP4P, TIP4PEW, TIP5P, TIP4PD, DMSO, METHANOL, OCTANOL  
Default values: SPC  

[protocol]

The default protocols for desmond relaxation and the production are activated.
To deactivate any step, the option should be set to 'no', 'off' or 'false'.  

* stage1 < Relaxation stage 1 >  
Acceptable values: yes, true, on or no, false, off  
Default values: yes  
The same procedure applies to stage2, stage3, stage4, and stage5.  

* stage1_time < Relaxation stage 1 time (ps) >  
Default values stage1: 100, stage2: 12, stage3: 12, stage4: 12, stage5: 24 and production_time: 100000ps = 100ns.  

* production_timestep_bonded: < Production timestep bonded (ps) >  
Default values: 0.002  

* production_timestep_near: < Production timestep near (ps) >  
Default values: 0.002  

* production_timestep_far: < Production timestep far (ps) >  
Default values: 0.006  

* stage1_temp < Relaxation stage 1 temperature (K) >  
Default values: 10.0  
The default values for the other stages are:  
stage1_temp: 10.0, stage2_temp: 10.0, stage3_temp: 10.0, stage4_temp: 300.0, stage5_temp: 300.0, production_temp: 10.0  

* stage1_ensemble < Relaxation stage 1 ensemble >  
Acceptable values: NVE, NVT, NPT  
Note, NPAT and NPYT will be implemented in the future.  
Default values: stage1_ensemble: NVT, stage2_ensemble: NVT, stage3_ensemble: NPT, stage4_ensemble: NPT, stage5_ensemble: NPT, production_ensemble: NPT.  

* stage1_method < Relaxation stage 1 method >  
Acceptable values: NVT, NPT  

## License

Licensed under the MIT license, see the separate LICENSE file.