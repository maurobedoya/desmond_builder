"""
Version: 1.0
"If this script is useful for you, consider giving acknowledgments comments in the publication."
Contact:
Mauricio Bedoya
maurobedoyat@gmail.com"
"""
from __future__ import print_function

import argparse
import os
from posix import sched_param
import subprocess
import sys
import glob
from typing_extensions import TypeAlias
import numpy as np
from os import path, PathLike, supports_fd, write

from dataclasses import dataclass, fields
from typing import Dict, List, Optional, Set, TextIO, Tuple
import configparser
import random


@dataclass
class Args:
    """Argument parser for overall settings."""

    input: str
    file: str
    desmond_path: str
    windows: str = "false"
    workdir: str = "md_run"


@dataclass
class BuilderOptions:
    """Variables for Builder settings."""

    basename: Optional[str] = "None"
    counterions: Optional[str] = None
    counterions_positive_ion = "Na"
    counterions_negative_ion = "Cl"
    ion: Optional[str] = None
    number: Optional[str] = "neutralize_system"

    shape: Optional[str] = "orthorhombic"
    size: Optional[str] = "10.0 10.0 10.0"
    size_type: Optional[str] = "buffer"

    ions_away: Optional[str] = None
    ion_awaydistance: Optional[str] = "5.0"
    ion_awayfrom: Optional[str] = "protein"

    override_forcefield: Optional[str] = "OPLS_2005"
    rezero_system: Optional[str] = "true"
    forcefield: Optional[str] = "OPLS_2005"

    salt: Optional[str] = "false"
    concentration: Optional[float] = 0.15
    positive_ion: Optional[str] = "Na"
    negative_ion: Optional[str] = "Cl"
    solvent: Optional[str] = "SPC"

    def __init__(
        self,
        opts: Dict[str, str],
        file: str,
        filename: str,
        desmond_path: str,
        file_path: str,
        windows: str,
    ) -> None:
        self.opts = opts
        self.file = file
        self.filename = filename
        self.basename = str(filename).split(".", maxsplit=1)[0]
        self.desmond_path = desmond_path
        self.windows = windows
        self.file_path = file_path
        for key in self.opts:
            setattr(self, key, self.opts[key])

    def __getattr__(self, item):
        if item not in self.opts:
            return None
        return self.opts[item]


@dataclass
class ProtocolOptions:
    """Variables for Protocol settings."""

    stage1: Optional[bool] = True
    stage2: Optional[bool] = True
    stage3: Optional[bool] = True
    stage4: Optional[bool] = True
    stage5: Optional[bool] = True
    production: Optional[bool] = True

    stage1_time: Optional[int] = 100
    stage2_time: Optional[int] = 12
    stage3_time: Optional[int] = 12
    stage4_time: Optional[int] = 12
    stage5_time: Optional[int] = 24
    production_time: Optional[int] = 100000

    stage1_timestep: Optional[str] = "0.001 0.001 0.003"
    stage2_timestep: Optional[str] = "0.001 0.001 0.003"
    production_timestep_bonded: Optional[float] = 0.002
    production_timestep_near: Optional[float] = 0.002
    production_timestep_far: Optional[float] = 0.006

    stage1_temp: Optional[float] = 10.0
    stage2_temp: Optional[float] = 10.0
    stage3_temp: Optional[float] = 10.0
    stage4_temp: Optional[float] = 300.0
    stage5_temp: Optional[float] = 300.0
    production_temp: Optional[float] = 300.0
    production_temp_group: Optional[str] = "0"

    stage1_ensemble: Optional[str] = "NVT"
    stage2_ensemble: Optional[str] = "NVT"
    stage3_ensemble: Optional[str] = "NPT"
    stage4_ensemble: Optional[str] = "NPT"
    stage5_ensemble: Optional[str] = "NPT"
    production_ensemble: Optional[str] = "NPT"

    stage1_method: Optional[str] = "Brownie"
    stage2_method: Optional[str] = "Berendsen"
    stage3_method: Optional[str] = "Berendsen"
    stage4_method: Optional[str] = "Berendsen"
    stage5_method: Optional[str] = "Berendsen"
    production_method: Optional[str] = "MTK"

    stage1_thermostat_tau: Optional[float] = None
    stage2_thermostat_tau: Optional[float] = 0.1
    stage3_thermostat_tau: Optional[float] = 0.1
    stage4_thermostat_tau: Optional[float] = 0.1
    stage5_thermostat_tau: Optional[float] = 0.1
    production_thermostat_tau: Optional[float] = 1.0

    stage1_barostat_tau: Optional[float] = None
    stage2_barostat_tau: Optional[float] = None
    stage3_barostat_tau: Optional[float] = 50.0
    stage4_barostat_tau: Optional[float] = 50.0
    stage5_barostat_tau: Optional[float] = 2.0
    production_barostat_tau: Optional[float] = 2.0

    stage1_traj_center: Optional[str] = "[]"
    stage2_traj_center: Optional[str] = "[]"
    stage3_traj_center: Optional[str] = "[]"
    stage4_traj_center: Optional[str] = "[]"
    stage5_traj_center: Optional[str] = "solute"

    stage1_restraints_number_pos: Optional[int] = 1
    stage1_restraints_atoms_pos: Optional[str] = "solute_heavy_atom"
    stage1_restraints_forces_pos: Optional[float] = 50.0

    stage1_restraints_number_dist: Optional[int] = 0
    stage1_restraints_atoms_dist: Optional[str] = None
    stage1_restraints_forces_dist: Optional[float] = None
    stage1_restraints_r0_dist: Optional[float] = None

    stage1_restraints_number_ang: Optional[int] = 0
    stage1_restraints_atoms_ang: Optional[str] = None
    stage1_restraints_forces_ang: Optional[float] = None
    stage1_restraints_theta0_ang: Optional[float] = None

    stage1_restraints_number_imp: Optional[int] = 0
    stage1_restraints_atoms_imp: Optional[str] = None
    stage1_restraints_forces_imp: Optional[float] = None
    stage1_restraints_phi0_imp: Optional[float] = None

    stage2_restraints_number_pos: Optional[int] = 1
    stage2_restraints_atoms_pos: Optional[str] = "solute_heavy_atom"
    stage2_restraints_forces_pos: Optional[float] = 50.0

    stage2_restraints_number_dist: Optional[int] = 0
    stage2_restraints_atoms_dist: Optional[str] = None
    stage2_restraints_forces_dist: Optional[float] = None
    stage2_restraints_r0_dist: Optional[float] = None

    stage2_restraints_number_ang: Optional[int] = 0
    stage2_restraints_atoms_ang: Optional[str] = None
    stage2_restraints_forces_ang: Optional[float] = None
    stage2_restraints_theta0_ang: Optional[float] = None

    stage2_restraints_number_imp: Optional[int] = 0
    stage2_restraints_atoms_imp: Optional[str] = None
    stage2_restraints_forces_imp: Optional[float] = None
    stage2_restraints_phi0_imp: Optional[float] = None

    stage3_restraints_number_pos: Optional[int] = 1
    stage3_restraints_atoms_pos: Optional[str] = "solute_heavy_atom"
    stage3_restraints_forces_pos: Optional[float] = 50.0

    stage3_restraints_number_dist: Optional[int] = 0
    stage3_restraints_atoms_dist: Optional[str] = None
    stage3_restraints_forces_dist: Optional[float] = None
    stage3_restraints_r0_dist: Optional[float] = None

    stage3_restraints_number_ang: Optional[int] = 0
    stage3_restraints_atoms_ang: Optional[str] = None
    stage3_restraints_forces_ang: Optional[float] = None
    stage3_restraints_theta0_ang: Optional[float] = None

    stage3_restraints_number_imp: Optional[int] = 0
    stage3_restraints_atoms_imp: Optional[str] = None
    stage3_restraints_forces_imp: Optional[float] = None
    stage3_restraints_phi0_imp: Optional[float] = None

    stage4_restraints_number_pos: Optional[int] = 1
    stage4_restraints_atoms_pos: Optional[str] = "solute_heavy_atom"
    stage4_restraints_forces_pos: Optional[float] = 50.0

    stage4_restraints_number_dist: Optional[int] = 0
    stage4_restraints_atoms_dist: Optional[str] = None
    stage4_restraints_forces_dist: Optional[float] = None
    stage4_restraints_r0_dist: Optional[float] = None

    stage4_restraints_number_ang: Optional[int] = 0
    stage4_restraints_atoms_ang: Optional[str] = None
    stage4_restraints_forces_ang: Optional[float] = None
    stage4_restraints_theta0_ang: Optional[float] = None

    stage4_restraints_number_imp: Optional[int] = 0
    stage4_restraints_atoms_imp: Optional[str] = None
    stage4_restraints_forces_imp: Optional[float] = None
    stage4_restraints_phi0_imp: Optional[float] = None

    stage5_restraints_number_pos: Optional[int] = 0
    stage5_restraints_atoms_pos: Optional[str] = "solute_heavy_atom"
    stage5_restraints_forces_pos: Optional[float] = None

    stage5_restraints_number_dist: Optional[int] = 0
    stage5_restraints_atoms_dist: Optional[str] = None
    stage5_restraints_forces_dist: Optional[float] = None
    stage5_restraints_r0_dist: Optional[float] = None

    stage5_restraints_number_ang: Optional[int] = 0
    stage5_restraints_atoms_ang: Optional[str] = None
    stage5_restraints_forces_ang: Optional[float] = None
    stage5_restraints_theta0_ang: Optional[float] = None

    stage5_restraints_number_imp: Optional[int] = 0
    stage5_restraints_atoms_imp: Optional[str] = None
    stage5_restraints_forces_imp: Optional[float] = None
    stage5_restraints_phi0_imp: Optional[float] = None

    production_bigger_rclone: Optional[str] = "false"
    production_checkpt_first: Optional[float] = 0.0
    production_checkpt_interval: Optional[float] = 240.06
    production_write_last_step: Optional[str] = "true"
    production_cutoff: Optional[float] = 9.0
    production_elapsed_time: Optional[float] = 0.0
    production_energy_group: Optional[str] = "false"
    production_eneseq_first: Optional[float] = 0.0
    production_eneseq_interval: Optional[float] = 1.2
    production_glue: Optional[str] = "solute"
    production_maeff_first: Optional[float] = 0.0
    production_maeff_interval: Optional[float] = 120.0
    production_maeff_periodicfix: Optional[str] = "true"
    production_meta: Optional[str] = "false"
    production_pressure: Optional[float] = 1.01325
    production_pressure_type: Optional[str] = "isotropic"
    production_randomize_vel_first: Optional[float] = 0.0
    production_randomize_vel_interval: Optional[str] = "inf"
    production_randomize_vel_seed: Optional[int] = random.randint(0000, 9999)
    production_simbox_first: Optional[float] = 0.0
    production_simbox_interval: Optional[float] = 1.2
    production_surface_tension: Optional[float] = 0.0
    production_taper: Optional[str] = "false"
    production_traj_center: Optional[str] = ""
    production_traj_first: Optional[float] = 0.0
    production_traj_format: Optional[str] = "dtr"
    production_traj_frames_per_file: Optional[int] = 250
    production_traj_interval: Optional[float] = 50.0
    production_traj_periodicfix: Optional[str] = "true"
    production_traj_write_velocity: Optional[str] = "false"

    production_restraints_number_pos: Optional[int] = 0
    production_restraints_atoms_pos: Optional[str] = None
    production_restraints_forces_pos: Optional[float] = None

    production_restraints_number_dist: Optional[int] = 0
    production_restraints_atoms_dist: Optional[str] = None
    production_restraints_forces_dist: Optional[float] = None
    production_restraints_r0_dist: Optional[float] = None

    production_restraints_number_ang: Optional[int] = 0
    production_restraints_atoms_ang: Optional[str] = None
    production_restraints_forces_ang: Optional[float] = None
    production_restraints_theta0_ang: Optional[float] = None

    production_restraints_number_imp: Optional[int] = 0
    production_restraints_atoms_imp: Optional[str] = None
    production_restraints_forces_imp: Optional[float] = None
    production_restraints_phi0_imp: Optional[float] = None

    name_pos: Optional[str] = "posre_harm"
    name_dist: Optional[str] = "stretch_harm"
    name_ang: Optional[str] = "angle_harm"
    name_imp: Optional[str] = "improper_harm"

    # Additional stages
    additional_stages: Optional[int] = 0
    additional_stage_times: Optional[str] = 0
    additional_stage_temps: Optional[str] = 0
    additional_stage_ensembles: Optional[str] = 0
    additional_stage_methods: Optional[str] = 0
    additional_stage_thermostat_tau: Optional[float] = 0.1
    additional_stage_barostat_tau: Optional[float] = 2.0

    additional_stage_restraints_number_pos: Optional[int] = 0
    additional_stage_restraints_atoms_pos: Optional[str] = None
    additional_stage_restraints_forces_pos: Optional[float] = None

    additional_stage_restraints_number_dist: Optional[int] = 0
    additional_stage_restraints_atoms_dist: Optional[str] = None
    additional_stage_restraints_forces_dist: Optional[float] = None
    additional_stage_restraints_r0_dist: Optional[float] = None

    additional_stage_restraints_number_ang: Optional[int] = 0
    additional_stage_restraints_atoms_ang: Optional[str] = None
    additional_stage_restraints_forces_ang: Optional[float] = None
    additional_stage_restraints_theta0_ang: Optional[float] = None

    additional_stage_restraints_number_imp: Optional[int] = 0
    additional_stage_restraints_atoms_imp: Optional[str] = None
    additional_stage_restraints_forces_imp: Optional[float] = None
    additional_stage_restraints_phi0_imp: Optional[float] = None

    additional_stage_traj_center: Optional[str] = "solute"
    # Run protocols
    run_preparation: Optional[str] = "false"
    run_protocols: Optional[str] = "false"

    def __init__(self, opts: Dict) -> None:
        self.opts = opts
        for key in self.opts:
            setattr(self, key, self.opts[key])

        input_names = [opt.name for opt in fields(ProtocolOptions)]
        for key in self.opts:
            try:
                if key not in input_names:
                    raise InputError(key)
            except InputError as e_rror:
                print(f"Error: {e_rror.args[0]}")
                print("Please check the input file.")
                sys.exit()
        # Define title for stages
        self.stage1_title: Optional[
            str] = f"{getattr(self, 'stage1_method')} {getattr(self, 'stage1_ensemble')}, T = {getattr(self, 'stage1_temp')} K, {getattr(self, 'stage1_time')}ps"
        self.stage2_title: Optional[
            str] = f"{getattr(self, 'stage2_method')} {getattr(self, 'stage2_ensemble')}, T = {getattr(self, 'stage2_temp')} K, {getattr(self, 'stage2_time')}ps"
        self.stage3_title: Optional[
            str] = f"{getattr(self, 'stage3_method')} {getattr(self, 'stage3_ensemble')}, T = {getattr(self, 'stage3_temp')} K, {getattr(self, 'stage3_time')}ps"
        self.stage4_title: Optional[
            str] = f"{getattr(self, 'stage4_method')} {getattr(self, 'stage4_ensemble')}, T = {getattr(self, 'stage4_temp')} K, {getattr(self, 'stage4_time')}ps"
        self.stage5_title: Optional[
            str] = f"{getattr(self, 'stage5_method')} {getattr(self, 'stage5_ensemble')}, T = {getattr(self, 'stage5_temp')} K, {getattr(self, 'stage5_time')}ps"
        self.production_title: Optional[
            str] = f"{getattr(self, 'production_method')} {getattr(self, 'production_ensemble')}, T = {getattr(self, 'production_temp')} K, {getattr(self, 'production_time')}ps"

        # TO:DO Check if all required options are set here.
        # To check type and number of additional positional restraints
        # TO:DO Add check for missing constants variables. r0, theta0, phi0.
        # TO:DO Add check for length of ensembles and methods, bug if not equal or long variables

        # To check type of stages 1-5 restraints
        for stage in range(1, 6):
            for restraint in ["pos", "dist", "ang", "imp"]:
                if getattr(self,
                           f"stage{stage}_restraints_number_{restraint}") != 0:
                    try:
                        value = getattr(
                            self,
                            f"stage{stage}_restraints_number_{restraint}")
                        # print(value)
                        value = int(value)
                    except ValueError as e_rror:
                        print(f"Error: {e_rror.args[0]}")
                        print(
                            f"stage{stage}_restraints_number_{restraint} must be an integer, not '{value}'"
                        )
                        print("Please check the input file.")
                        sys.exit()
        # To check length of stages 1-5 restraints atoms, forces and constants.
        for stage in range(1, 6):
            for restraint in ["pos", "dist", "ang", "imp"]:
                if restraint == "pos":
                    factor = 1
                if restraint == "dist":
                    constant = "r0"
                    factor = 2
                elif restraint == "ang":
                    constant = "theta0"
                    factor = 3
                elif restraint == "imp":
                    constant = "phi0"
                    factor = 4
                if int(
                        getattr(self,
                                f"stage{stage}_restraints_number_{restraint}")
                ) != 0:
                    try:
                        len1 = int(
                            getattr(
                                self,
                                f"stage{stage}_restraints_number_{restraint}"))

                        len2a = getattr(
                            self, f"stage{stage}_restraints_atoms_{restraint}")
                        len2 = len(list(str(len2a).split(",")))
                        len3a = getattr(
                            self,
                            f"stage{stage}_restraints_forces_{restraint}")
                        len3 = len(list(str(len3a).split(",")))
                        if (restraint == "dist" or restraint == "ang"
                                or restraint == "imp"):
                            len4a = getattr(
                                self,
                                f"stage{stage}_restraints_{constant}_{restraint}",
                            )
                            len4 = len(list(str(len4a).split(",")))
                            if len1 != len4:  # or len1 != len3 or len1 != len4:
                                raise LenError(
                                    f"stage{stage}_restraints_number_{restraint}",
                                    len1,
                                    f"stage{stage}_restraints_{constant}_{restraint}",
                                    len4,
                                    len4a,
                                )
                        if len1 != len2 / factor:  # or len1 != len3 or len1 != len4:
                            raise LenError2(
                                f"stage{stage}_restraints_number_{restraint}",
                                len1,
                                f"stage{stage}_restraints_atoms_{restraint}",
                                len2,
                                len2a,
                            )
                        if len1 != len3:
                            raise LenError(
                                f"stage{stage}_restraints_number_{restraint}",
                                len1,
                                f"stage{stage}_restraints_forces_{restraint}",
                                len3,
                                len3a,
                            )
                    except LenError as e_rror:
                        print(f"Error: {e_rror.args[0]}")
                        # print(
                        #    f"Length of stage{stage}_restraints_atoms_{restraint} and stage{stage}_restraints_forces_{restraint} and stage{stage}_restraints_constants_{restraint} must be equal to stage{stage}_restraints_number_{restraint}"
                        # )
                        print("Please check the input file.")
                        sys.exit()
                    except LenError2 as e_rror:
                        print(f"Error: {e_rror.args[0]}")
                        # print(
                        #    f"Length of stage{stage}_restraints_atoms_{restraint} and stage{stage}_restraints_forces_{restraint} and stage{stage}_restraints_constants_{restraint} must be equal to stage{stage}_restraints_number_{restraint}"
                        # )
                        print("Please check the input file.")
                        sys.exit()
        # To check lenght of stage_times, temps, ensembles, methods, barostat_tau and thermostat_tau.
        if self.additional_stages != 0:
            add_stages = int(self.additional_stages)
            try:
                len1 = int(self.additional_stages)
                len2 = len(list(str(self.additional_stage_times).split(",")))
                len3 = len(list(str(self.additional_stage_temps).split(",")))
                len4 = len(
                    list(str(self.additional_stage_ensembles).split(",")))
                len5 = len(list(str(self.additional_stage_methods).split(",")))
                len6 = len(
                    list(str(self.additional_stage_barostat_tau).split(",")))
                len7 = len(
                    list(str(self.additional_stage_thermostat_tau).split(",")))

                if len1 != len2 and len2 != 1:
                    raise LenError3(
                        "additional_stages",
                        len1,
                        "additional_stage_times",
                        len2,
                        self.additional_stage_times,
                    )
                if len1 != len3 and len3 != 1:
                    raise LenError3(
                        "additional_stages",
                        len1,
                        "additional_stage_temps",
                        len3,
                        self.additional_stage_temps,
                    )
                if len1 != len4 and len4 != 1:
                    raise LenError3(
                        "additional_stages",
                        len1,
                        "additional_stage_ensembles",
                        len4,
                        self.additional_stage_ensembles,
                    )
                if len1 != len5 and len5 != 1:
                    raise LenError3(
                        "additional_stages",
                        len1,
                        "additional_stage_methods",
                        len5,
                        self.additional_stage_methods,
                    )
                if len1 != len6 and len6 != 1:
                    raise LenError3(
                        "additional_stages",
                        len1,
                        "additional_stage_barostat_tau",
                        len6,
                        self.additional_stage_barostat_tau,
                    )
                if len1 != len7 and len7 != 1:
                    raise LenError3(
                        "additional_stages",
                        len1,
                        "additional_stage_thermostat_tau",
                        len7,
                        self.additional_stage_thermostat_tau,
                    )
            except LenError3 as e_rror:
                print(f"Error: {e_rror.args[0]}")
                print("Please check the input file.")
                sys.exit()
        # To check type and number of additional positional restraints
        if self.additional_stage_restraints_number_pos != 0:
            add_number = self.additional_stage_restraints_number_pos.split(",")
            add_stages = int(self.additional_stages)
            try:
                if len(add_number) != add_stages:
                    raise LenError(
                        "additional_stages",
                        add_stages,
                        "additional_stage_restraints_number_pos",
                        len(add_number),
                        add_number,
                    )
            except LenError as e_rror:
                print(f"Error: {e_rror.args[0]}")
                sys.exit()
            for value in self.additional_stage_restraints_number_pos.split(
                    ","):
                try:
                    value = int(value)
                except ValueError as e_rror:
                    print(f"Error: {e_rror.args[0]}")
                    print(
                        "Please check the values of 'additional_stage_restraints_number_pos'."
                    )
                    sys.exit()
        # To check type and number of additional distance restraints
        if self.additional_stage_restraints_number_dist != 0:
            add_number = self.additional_stage_restraints_number_dist.split(
                ",")
            add_stages = int(self.additional_stages)
            try:
                if len(add_number) != add_stages:
                    raise LenError(
                        "additional_stages",
                        add_stages,
                        "additional_stage_restraints_number_dist",
                        len(add_number),
                        add_number,
                    )
            except LenError as e_rror:
                print(f"Error: {e_rror.args[0]}")
                sys.exit()
            for value in self.additional_stage_restraints_number_dist.split(
                    ","):
                try:
                    value = int(value)
                except ValueError as e_rror:
                    print(f"Error: {e_rror.args[0]}")
                    print(
                        "Please check the values of 'additional_stage_restraints_number_dist'."
                    )
                    sys.exit()

        # To check type and number of additional angle restraints
        if self.additional_stage_restraints_number_ang != 0:
            add_number = self.additional_stage_restraints_number_ang.split(",")
            add_stages = int(self.additional_stages)
            try:
                if len(add_number) != add_stages:
                    raise LenError(
                        "additional_stages",
                        add_stages,
                        "additional_stage_restraints_number_ang",
                        len(add_number),
                        add_number,
                    )
            except LenError as e_rror:
                print(f"Error: {e_rror.args[0]}")
                sys.exit()
            for value in self.additional_stage_restraints_number_ang.split(
                    ","):
                try:
                    value = int(value)
                except ValueError as e_rror:
                    print(f"Error: {e_rror.args[0]}")
                    print(
                        "Please check the values of 'additional_stage_restraints_number_ang'."
                    )
                    sys.exit()
        # To check type and number of additional improper restraints
        if self.additional_stage_restraints_number_imp != 0:
            add_number = self.additional_stage_restraints_number_imp.split(",")
            add_stages = int(self.additional_stages)
            try:
                if len(add_number) != add_stages:
                    raise LenError(
                        "additional_stages",
                        add_stages,
                        "additional_stage_restraints_number_imp",
                        len(add_number),
                        add_number,
                    )
            except LenError as e_rror:
                print(f"Error: {e_rror.args[0]}")
                sys.exit()
            for value in self.additional_stage_restraints_number_imp.split(
                    ","):
                try:
                    value = int(value)
                except ValueError as e_rror:
                    print(f"Error: {e_rror.args[0]}")
                    print(
                        "Please check the values of 'additional_stage_restraints_number_imp'."
                    )
                    sys.exit()


class Error(Exception):
    """Base class for other exceptions"""

    pass


class InputError(Error):
    """
    Custom error class for wrong option in the input."""
    def __init__(self, key, message="Unknown option "):
        self.key = key
        self.message = f"{message}'{key}'"
        super().__init__(self.message)


class LenError(Error):
    """
    Custom error class for wrong number of values in the input."""
    def __init__(self,
                 var1,
                 len1,
                 var2,
                 len2,
                 passed,
                 err1="Wrong length of values."):
        self.var1 = var1
        self.var2 = var2
        self.len1 = len1
        self.len2 = len2
        self.passed = passed
        self.message1 = f"{err1}"
        self.message2 = f"'{var2}' option should have '{len1}' values as '{var1}' option, but it has '{len2}' values: '[{passed}]'."
        self.message = f"{self.message1} {self.message2}"
        super().__init__(self.message)


class LenError2(Error):
    """
    Custom error class for wrong number of values in the input."""
    def __init__(self,
                 var1,
                 len1,
                 var2,
                 len2,
                 passed,
                 err1="Wrong length of values."):
        self.var1 = var1
        self.var2 = var2
        self.len1 = len1
        self.len2 = len2
        self.passed = passed
        self.message1 = f"{err1}"
        self.message2 = f"'{var2}' option should have '{len1}' pairs of values as ['{var1}']x2, but it has '{len2}' values: '[{passed}]'."
        self.message = f"{self.message1} {self.message2}"
        super().__init__(self.message)


class LenError3(Error):
    """
    Custom error class for wrong number of values in the input. When the accepted values are 1 or the length of the list."""
    def __init__(self,
                 var1,
                 len1,
                 var2,
                 len2,
                 passed,
                 err1="Wrong length of values."):
        self.var1 = var1
        self.var2 = var2
        self.len1 = len1
        self.len2 = len2
        self.passed = passed
        self.message1 = f"{err1}"
        self.message2 = f"'{var2}' option should have '{len1}' values as '{var1}' option or only one value, but it has '{len2}' values: '[{passed}]'."
        self.message = f"{self.message1} {self.message2}"
        super().__init__(self.message)


def identation(indentvar: int = 0) -> Tuple[str, str]:
    indent = indentvar
    if indentvar == 0:
        outer_space = " " * (indent - 1)
        inner_space = " " * (indent + 4 - 1)
    elif indentvar == 1:
        outer_space = " " * (indent + 2)
        inner_space = " " * (indent + 6)
    elif indentvar == 2:
        outer_space = " " * (indent + 5)
        inner_space = " " * (indent + 9)
    else:
        print("Error: indentation is not 0 or 1")
    return outer_space, inner_space


class Protocol:
    def __init__(
        self,
        file: str = None,
        builder_opts: BuilderOptions = None,
        protocol_opts: ProtocolOptions = None,
    ) -> None:
        self.builder_opts = builder_opts
        self.p_opts = protocol_opts
        self.stage1 = protocol_opts.stage1
        self.stage2 = self.p_opts.stage2
        self.stage3 = self.p_opts.stage3
        self.stage4 = self.p_opts.stage4
        self.stage5 = self.p_opts.stage5
        self.production = self.p_opts.production
        self.file = file
        self.basename = builder_opts.basename
        # self.outputname = self.builder_opts.outputname

    def write(self) -> None:
        path_preparation = str(self.basename + "_md.msj")
        outer_space, inner_space = identation(0)
        eq = "= "
        q = '"'
        with open(path_preparation, "w", encoding="utf8") as fd:
            print("Preparing input files for MD protocol...")
            print("# Desmond protocol", file=fd)
            print("# Time units are in ps", file=fd)
            print("# Energy units are in kcal/mol", file=fd)
            print(file=fd)
            print(f"{outer_space}task {'{'}", file=fd)
            # outer_space, inner_space = identation(1)
            print(f"{inner_space} {'task':<16}{eq}{q}{'desmond:auto'}{q}",
                  file=fd)
            print(f"{inner_space} {'set_family':<16}{eq}{'{'}", file=fd)

            outer_space, inner_space = identation(2)
            print(f"{outer_space} {'desmond':<12}{eq}{'{'}", file=fd)
            print(f"{inner_space} {'checkpt.write_last_step':<16} {eq}{'no'}",
                  file=fd)
            print(f"{outer_space} {'}'}", file=fd)
            outer_space, inner_space = identation(0)
            print(f"{inner_space} {'}'}", file=fd)
            print(f"{outer_space}{'}'}", file=fd)
            print(file=fd)
            # Stage 1 block
            if str(self.p_opts.stage1).lower() in ["yes", "on", "true"]:
                outer_space, inner_space = identation(0)
                print(f"{outer_space}{'simulate':<20}{'{'}", file=fd)
                print(
                    f"{inner_space} {'title':<16}{eq}{q}{self.p_opts.stage1_title}{q}",
                    file=fd,
                )
                print(f"{inner_space} {'annealing':<16}{eq}{'off'}", file=fd)
                print(
                    f"{inner_space} {'time':<16}{eq}{self.p_opts.stage1_time}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'timestep':<16}{eq}{'['}{self.p_opts.stage1_timestep}{']'}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'temperature':<16}{eq}{self.p_opts.stage1_temp}",
                    file=fd,
                )
                outer_space, inner_space = identation(1)
                print(f"{outer_space} {'ensemble':<16}{eq}{'{'}", file=fd)
                print(
                    f"{inner_space} {'class':<11} {eq}{self.p_opts.stage1_ensemble}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'method':<11} {eq}{self.p_opts.stage1_method}",
                    file=fd,
                )
                outer_space, inner_space = identation(2)
                print(f"{outer_space} {'brownie':<17}{eq}{'{'}", file=fd)
                print(f"{inner_space} {'delta_max':<12} {eq}{'0.1'}", file=fd)
                print(f"{outer_space} {'}'}", file=fd)
                outer_space, inner_space = identation(1)
                print(f"{outer_space} {'}'}", file=fd)
                # ==============================================================
                # Restraints block
                # ==============================================================
                # Block for stage1 positional-restraints
                if (int(self.p_opts.stage1_restraints_number_pos) != 0
                        or int(self.p_opts.stage1_restraints_number_dist) != 0
                        or int(self.p_opts.stage1_restraints_number_ang) != 0
                        or int(self.p_opts.stage1_restraints_number_imp) != 0):
                    outer_space, inner_space = identation(1)
                    print(f"{outer_space} {'restraints.new':<16}{eq}{'['}",
                          file=fd)
                    if int(self.p_opts.stage1_restraints_number_pos) != 0:
                        atoms, forces = self.set_restraint(
                            "stage1",
                            self.p_opts.stage1_restraints_number_pos,
                            self.p_opts.stage1_restraints_atoms_pos,
                            self.p_opts.stage1_restraints_forces_pos,
                            "positional",
                            None,
                        )
                        for i in range(
                                int(self.p_opts.stage1_restraints_number_pos)):
                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_pos}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[i]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]} {forces[i]} {forces[i]}]",
                                file=fd,
                            )
                            print(f"{outer_space} {'}'}", file=fd)

                    if int(self.p_opts.stage1_restraints_number_dist) != 0:
                        # Block for stage1 distance-restraints
                        outer_space, inner_space = identation(1)
                        atoms, forces, constants = self.set_restraint(
                            "stage1",
                            self.p_opts.stage1_restraints_number_dist,
                            self.p_opts.stage1_restraints_atoms_dist,
                            self.p_opts.stage1_restraints_forces_dist,
                            "distance",
                            self.p_opts.stage1_restraints_r0_dist,
                        )
                        j = 0
                        for i in range(
                                int(self.p_opts.stage1_restraints_number_dist)
                        ):

                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_dist}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'r0':<11} {eq}{constants[i]}",
                                file=fd)
                            print(f"{outer_space} {'}'}", file=fd)
                            j += 2
                    if int(self.p_opts.stage1_restraints_number_ang) != 0:
                        # Block for stage1 angle-restraints
                        outer_space, inner_space = identation(1)
                        atoms, forces, constants = self.set_restraint(
                            "stage1",
                            self.p_opts.stage1_restraints_number_ang,
                            self.p_opts.stage1_restraints_atoms_ang,
                            self.p_opts.stage1_restraints_forces_ang,
                            "angle",
                            self.p_opts.stage1_restraints_theta0_ang,
                        )
                        j = 0
                        for i in range(
                                int(self.p_opts.stage1_restraints_number_ang)):
                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_ang}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q} {q}{atoms[j+2]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'theta0':<11} {eq}{constants[i]}",
                                file=fd,
                            )
                            print(f"{outer_space} {'}'}", file=fd)
                            j += 3
                    if int(self.p_opts.stage1_restraints_number_imp) != 0:
                        # Block for stage1 improper-restraints
                        outer_space, inner_space = identation(1)
                        atoms, forces, constants = self.set_restraint(
                            "stage1",
                            self.p_opts.stage1_restraints_number_imp,
                            self.p_opts.stage1_restraints_atoms_imp,
                            self.p_opts.stage1_restraints_forces_imp,
                            "improper",
                            self.p_opts.stage1_restraints_phi0_imp,
                        )
                        j = 0
                        for i in range(
                                int(self.p_opts.stage1_restraints_number_imp)):
                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_imp}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q} {q}{atoms[j+2]}{q} {q}{atoms[j+3]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'phi0':<11} {eq}{constants[i]}",
                                file=fd,
                            )
                            print(f"{outer_space} {'}'}", file=fd)
                            j += 4
                    outer_space, inner_space = identation(1)
                    print(f"{outer_space} {']'}", file=fd)
                outer_space, inner_space = identation(0)
                print(f"{outer_space}{'}'}", file=fd)
                print(file=fd)
            # ==============================================================
            # Restraints block
            # ==============================================================
            # Stage 2 block
            if self.p_opts.stage2.lower() in ["yes", "on", "true"]:
                outer_space, inner_space = identation(0)
                print(f"{outer_space}{'simulate':<20}{'{'}", file=fd)
                print(
                    f"{inner_space} {'title':<16}{eq}{q}{self.p_opts.stage2_title}{q}",
                    file=fd,
                )
                gpu_text = '[["==" "-gpu" "@*.*.jlaunch_opt[-1]"] \'ensemble.method = Langevin\']'
                print(f"{inner_space} {'effect_if':<16}{eq}{gpu_text}",
                      file=fd)
                print(f"{inner_space} {'annealing':<16}{eq}{'off'}", file=fd)
                print(
                    f"{inner_space} {'time':<16}{eq}{self.p_opts.stage2_time}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'timestep':<16}{eq}{'['}{self.p_opts.stage2_timestep}{']'}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'temperature':<16}{eq}{self.p_opts.stage2_temp}",
                    file=fd,
                )
                outer_space, inner_space = identation(1)
                print(f"{outer_space} {'ensemble':<16}{eq}{'{'}", file=fd)
                print(
                    f"{inner_space} {'class':<11} {eq}{self.p_opts.stage2_ensemble}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'method':<11} {eq}{self.p_opts.stage2_method}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'thermostat.tau':<11} {eq}{self.p_opts.stage2_thermostat_tau}",
                    file=fd,
                )
                outer_space, inner_space = identation(1)
                print(f"{outer_space} {'}'}", file=fd)
                # Restraints block
                # ==============================================================
                # Restraints block
                # ==============================================================
                # Block for stage2 positional-restraints
                if (int(self.p_opts.stage2_restraints_number_pos) != 0
                        or int(self.p_opts.stage2_restraints_number_dist) != 0
                        or int(self.p_opts.stage2_restraints_number_ang) != 0
                        or int(self.p_opts.stage2_restraints_number_imp) != 0):
                    outer_space, inner_space = identation(1)
                    print(f"{outer_space} {'restraints.new':<16}{eq}{'['}",
                          file=fd)
                    if int(self.p_opts.stage2_restraints_number_pos) != 0:
                        atoms, forces = self.set_restraint(
                            "stage2",
                            self.p_opts.stage2_restraints_number_pos,
                            self.p_opts.stage2_restraints_atoms_pos,
                            self.p_opts.stage2_restraints_forces_pos,
                            "positional",
                            None,
                        )
                        for i in range(
                                int(self.p_opts.stage2_restraints_number_pos)):
                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_pos}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[i]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]} {forces[i]} {forces[i]}]",
                                file=fd,
                            )
                            print(f"{outer_space} {'}'}", file=fd)

                    if int(self.p_opts.stage2_restraints_number_dist) != 0:
                        # Block for stage2 distance-restraints
                        outer_space, inner_space = identation(1)
                        atoms, forces, constants = self.set_restraint(
                            "stage2",
                            self.p_opts.stage2_restraints_number_dist,
                            self.p_opts.stage2_restraints_atoms_dist,
                            self.p_opts.stage2_restraints_forces_dist,
                            "distance",
                            self.p_opts.stage2_restraints_r0_dist,
                        )
                        j = 0
                        for i in range(
                                int(self.p_opts.stage2_restraints_number_dist)
                        ):

                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_dist}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'r0':<11} {eq}{constants[i]}",
                                file=fd)
                            print(f"{outer_space} {'}'}", file=fd)
                            j += 2
                    if int(self.p_opts.stage2_restraints_number_ang) != 0:
                        # Block for stage2 angle-restraints
                        outer_space, inner_space = identation(1)
                        atoms, forces, constants = self.set_restraint(
                            "stage2",
                            self.p_opts.stage2_restraints_number_ang,
                            self.p_opts.stage2_restraints_atoms_ang,
                            self.p_opts.stage2_restraints_forces_ang,
                            "angle",
                            self.p_opts.stage2_restraints_theta0_ang,
                        )
                        j = 0
                        for i in range(
                                int(self.p_opts.stage2_restraints_number_ang)):
                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_ang}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q} {q}{atoms[j+2]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'theta0':<11} {eq}{constants[i]}",
                                file=fd,
                            )
                            print(f"{outer_space} {'}'}", file=fd)
                            j += 3
                    if int(self.p_opts.stage2_restraints_number_imp) != 0:
                        # Block for stage2 improper-restraints
                        outer_space, inner_space = identation(1)
                        atoms, forces, constants = self.set_restraint(
                            "stage2",
                            self.p_opts.stage2_restraints_number_imp,
                            self.p_opts.stage2_restraints_atoms_imp,
                            self.p_opts.stage2_restraints_forces_imp,
                            "improper",
                            self.p_opts.stage2_restraints_phi0_imp,
                        )
                        j = 0
                        for i in range(
                                int(self.p_opts.stage2_restraints_number_imp)):
                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_imp}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q} {q}{atoms[j+2]}{q} {q}{atoms[j+3]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'phi0':<11} {eq}{constants[i]}",
                                file=fd,
                            )
                            print(f"{outer_space} {'}'}", file=fd)
                            j += 4
                    outer_space, inner_space = identation(1)
                    print(f"{outer_space} {']'}", file=fd)
                    # End of restrains block
                print(file=fd)
                # ==============================================================
                # Restraints block END
                # ==============================================================
                outer_space, inner_space = identation(0)
                print(
                    f"{inner_space} {'randomize_velocity.interval':<29} {eq}{'1.0'}",
                    file=fd,
                )
                print(f"{inner_space} {'eneseq.interval':<29} {eq}{'0.3'}",
                      file=fd)

                print(
                    f"{inner_space} {'trajectory.center':<29} {eq}{self.p_opts.stage2_traj_center}",
                    file=fd)

                print(f"{outer_space}{'}'}", file=fd)
                print(file=fd)
            # Stage 3 block
            if self.p_opts.stage3.lower() in ["yes", "on", "true"]:
                outer_space, inner_space = identation(0)
                print(f"{outer_space}{'simulate':<20}{'{'}", file=fd)
                if self.p_opts.stage3_ensemble != "NVT":
                    print(
                        f"{inner_space} {'title':<16}{eq}{q}{self.p_opts.stage3_title}{q}",
                        file=fd,
                    )
                else:
                    self.p_opts.stage3_title = f"{self.p_opts.stage3_method} {self.p_opts.stage3_ensemble}, {self.p_opts.stage3_time}ps"
                    print(
                        f"{inner_space} {'title':<16}{eq}{q}{self.p_opts.stage3_title}{q}",
                        file=fd,
                    )

                if self.p_opts.stage3_ensemble != "NVT":
                    gpu_text = '[["==" "-gpu" "@*.*.jlaunch_opt[-1]"] \'ensemble.method = Langevin\']'
                    print(f"{inner_space} {'effect_if':<16}{eq}{gpu_text}",
                          file=fd)
                    print(f"{inner_space} {'annealing':<16}{eq}{'off'}",
                          file=fd)
                    print(
                        f"{inner_space} {'temperature':<16}{eq}{self.p_opts.stage3_temp}",
                        file=fd,
                    )
                else:
                    gpu_text1 = '[["@*.*.annealing"] \'annealing = off temperature = "@*.*.temperature[0][0]"\''
                    gpu_text2 = '["==" "-gpu" "@*.*.jlaunch_opt[-1]"] \'ensemble.method = Langevin\']'
                    print(f"{inner_space} {'effect_if':<16}{eq}{gpu_text1}",
                          file=fd)
                    print(f"{inner_space} {' ':<19}{gpu_text2}", file=fd)

                print(
                    f"{inner_space} {'time':<16}{eq}{self.p_opts.stage3_time}",
                    file=fd,
                )

                outer_space, inner_space = identation(1)
                print(f"{outer_space} {'ensemble':<16}{eq}{'{'}", file=fd)
                print(
                    f"{inner_space} {'class':<11} {eq}{self.p_opts.stage3_ensemble}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'method':<11} {eq}{self.p_opts.stage3_method}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'thermostat.tau':<11} {eq}{self.p_opts.stage3_thermostat_tau}",
                    file=fd,
                )
                if self.p_opts.stage3_ensemble != "NVT":
                    print(
                        f"{inner_space} {'barostat.tau':<11} {eq}{self.p_opts.stage3_barostat_tau}",
                        file=fd,
                    )
                outer_space, inner_space = identation(1)
                print(f"{outer_space} {'}'}", file=fd)
                # ==============================================================
                # Restraints block
                # ==============================================================
                # Block for stage3 positional-restraints
                if (int(self.p_opts.stage3_restraints_number_pos) != 0
                        or int(self.p_opts.stage3_restraints_number_dist) != 0
                        or int(self.p_opts.stage3_restraints_number_ang) != 0
                        or int(self.p_opts.stage3_restraints_number_imp) != 0):
                    outer_space, inner_space = identation(1)
                    print(f"{outer_space} {'restraints.new':<16}{eq}{'['}",
                          file=fd)
                    if int(self.p_opts.stage3_restraints_number_pos) != 0:
                        atoms, forces = self.set_restraint(
                            "stage3",
                            self.p_opts.stage3_restraints_number_pos,
                            self.p_opts.stage3_restraints_atoms_pos,
                            self.p_opts.stage3_restraints_forces_pos,
                            "positional",
                            None,
                        )
                        for i in range(
                                int(self.p_opts.stage3_restraints_number_pos)):
                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_pos}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[i]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]} {forces[i]} {forces[i]}]",
                                file=fd,
                            )
                            print(f"{outer_space} {'}'}", file=fd)

                    if int(self.p_opts.stage3_restraints_number_dist) != 0:
                        # Block for stage3 distance-restraints
                        outer_space, inner_space = identation(1)
                        atoms, forces, constants = self.set_restraint(
                            "stage3",
                            self.p_opts.stage3_restraints_number_dist,
                            self.p_opts.stage3_restraints_atoms_dist,
                            self.p_opts.stage3_restraints_forces_dist,
                            "distance",
                            self.p_opts.stage3_restraints_r0_dist,
                        )
                        j = 0
                        for i in range(
                                int(self.p_opts.stage3_restraints_number_dist)
                        ):

                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_dist}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'r0':<11} {eq}{constants[i]}",
                                file=fd)
                            print(f"{outer_space} {'}'}", file=fd)
                            j += 2
                    if int(self.p_opts.stage3_restraints_number_ang) != 0:
                        # Block for stage3 angle-restraints
                        outer_space, inner_space = identation(1)
                        atoms, forces, constants = self.set_restraint(
                            "stage3",
                            self.p_opts.stage3_restraints_number_ang,
                            self.p_opts.stage3_restraints_atoms_ang,
                            self.p_opts.stage3_restraints_forces_ang,
                            "angle",
                            self.p_opts.stage3_restraints_theta0_ang,
                        )
                        j = 0
                        for i in range(
                                int(self.p_opts.stage3_restraints_number_ang)):
                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_ang}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q} {q}{atoms[j+2]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'theta0':<11} {eq}{constants[i]}",
                                file=fd,
                            )
                            print(f"{outer_space} {'}'}", file=fd)
                            j += 3
                    if int(self.p_opts.stage3_restraints_number_imp) != 0:
                        # Block for stage3 improper-restraints
                        outer_space, inner_space = identation(1)
                        atoms, forces, constants = self.set_restraint(
                            "stage3",
                            self.p_opts.stage3_restraints_number_imp,
                            self.p_opts.stage3_restraints_atoms_imp,
                            self.p_opts.stage3_restraints_forces_imp,
                            "improper",
                            self.p_opts.stage3_restraints_phi0_imp,
                        )
                        j = 0
                        for i in range(
                                int(self.p_opts.stage3_restraints_number_imp)):
                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_imp}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q} {q}{atoms[j+2]}{q} {q}{atoms[j+3]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'phi0':<11} {eq}{constants[i]}",
                                file=fd,
                            )
                            print(f"{outer_space} {'}'}", file=fd)
                            j += 4
                    outer_space, inner_space = identation(1)
                    print(f"{outer_space} {']'}", file=fd)
                print(file=fd)
                # ==============================================================
                # Restraints block end
                # ==============================================================
                outer_space, inner_space = identation(0)

                if self.p_opts.stage3_ensemble != "NVT":
                    print(
                        f"{inner_space} {'randomize_velocity.interval':<29} {eq}{'1.0'}",
                        file=fd,
                    )
                print(f"{inner_space} {'eneseq.interval':<29} {eq}{'0.3'}",
                      file=fd)
                print(
                    f"{inner_space} {'trajectory.center':<29} {eq}{self.p_opts.stage3_traj_center}",
                    file=fd)
                print(f"{outer_space}{'}'}", file=fd)
                # solvate pocket block
                ##outer_space, inner_space = identation(0)
                # print(file=fd)
                ##print(f"{outer_space}{'solvate_pocket':<16}{'{'}", file=fd)
                # print(
                ##    f"{inner_space} {'should_skip':<11} {eq}{'true'}",
                # file=fd,
                # )
                # print(
                ##    f"{inner_space} {'ligand_file':<11} {eq}{'?'}",
                # file=fd,
                # )
                ##print(f"{outer_space}{'}'}", file=fd)
                print(file=fd)
            # Stage 4 block
            if self.p_opts.stage4.lower() in ["yes", "on", "true"]:
                outer_space, inner_space = identation(0)
                print(f"{outer_space}{'simulate':<20}{'{'}", file=fd)
                print(
                    f"{inner_space} {'title':<16}{eq}{q}{self.p_opts.stage4_title}{q}",
                    file=fd,
                )
                gpu_text1 = '[["@*.*.annealing"] \'annealing = off temperature = "@*.*.temperature[0][0]"\''
                gpu_text2 = '["==" "-gpu" "@*.*.jlaunch_opt[-1]"] \'ensemble.method = Langevin\']'
                print(f"{inner_space} {'effect_if':<16}{eq}{gpu_text1}",
                      file=fd)
                print(f"{inner_space} {' ':<19}{gpu_text2}", file=fd)

                print(
                    f"{inner_space} {'time':<16}{eq}{self.p_opts.stage4_time}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'temperature':<16}{eq}{self.p_opts.stage4_temp}",
                    file=fd,
                )
                outer_space, inner_space = identation(1)
                print(f"{outer_space} {'ensemble':<16}{eq}{'{'}", file=fd)
                print(
                    f"{inner_space} {'class':<11} {eq}{self.p_opts.stage4_ensemble}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'method':<11} {eq}{self.p_opts.stage4_method}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'thermostat.tau':<11} {eq}{self.p_opts.stage4_thermostat_tau}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'barostat.tau':<11} {eq}{self.p_opts.stage4_barostat_tau}",
                    file=fd,
                )
                outer_space, inner_space = identation(1)
                print(f"{outer_space} {'}'}", file=fd)
                # ==============================================================
                # Restraints block
                # ==============================================================
                # Block for stage4 positional-restraints
                if (int(self.p_opts.stage4_restraints_number_pos) != 0
                        or int(self.p_opts.stage4_restraints_number_dist) != 0
                        or int(self.p_opts.stage4_restraints_number_ang) != 0
                        or int(self.p_opts.stage4_restraints_number_imp) != 0):
                    outer_space, inner_space = identation(1)
                    print(f"{outer_space} {'restraints.new':<16}{eq}{'['}",
                          file=fd)
                    if int(self.p_opts.stage4_restraints_number_pos) != 0:
                        atoms, forces = self.set_restraint(
                            "stage4",
                            self.p_opts.stage4_restraints_number_pos,
                            self.p_opts.stage4_restraints_atoms_pos,
                            self.p_opts.stage4_restraints_forces_pos,
                            "positional",
                            None,
                        )
                        for i in range(
                                int(self.p_opts.stage4_restraints_number_pos)):
                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_pos}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[i]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]} {forces[i]} {forces[i]}]",
                                file=fd,
                            )
                            print(f"{outer_space} {'}'}", file=fd)

                    if int(self.p_opts.stage4_restraints_number_dist) != 0:
                        # Block for stage4 distance-restraints
                        outer_space, inner_space = identation(1)
                        atoms, forces, constants = self.set_restraint(
                            "stage4",
                            self.p_opts.stage4_restraints_number_dist,
                            self.p_opts.stage4_restraints_atoms_dist,
                            self.p_opts.stage4_restraints_forces_dist,
                            "distance",
                            self.p_opts.stage4_restraints_r0_dist,
                        )
                        j = 0
                        for i in range(
                                int(self.p_opts.stage4_restraints_number_dist)
                        ):

                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_dist}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'r0':<11} {eq}{constants[i]}",
                                file=fd)
                            print(f"{outer_space} {'}'}", file=fd)
                            j += 2
                    if int(self.p_opts.stage4_restraints_number_ang) != 0:
                        # Block for stage4 angle-restraints
                        outer_space, inner_space = identation(1)
                        atoms, forces, constants = self.set_restraint(
                            "stage4",
                            self.p_opts.stage4_restraints_number_ang,
                            self.p_opts.stage4_restraints_atoms_ang,
                            self.p_opts.stage4_restraints_forces_ang,
                            "angle",
                            self.p_opts.stage4_restraints_theta0_ang,
                        )
                        j = 0
                        for i in range(
                                int(self.p_opts.stage4_restraints_number_ang)):
                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_ang}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q} {q}{atoms[j+2]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'theta0':<11} {eq}{constants[i]}",
                                file=fd,
                            )
                            print(f"{outer_space} {'}'}", file=fd)
                            j += 3
                    if int(self.p_opts.stage4_restraints_number_imp) != 0:
                        # Block for stage4 improper-restraints
                        outer_space, inner_space = identation(1)
                        atoms, forces, constants = self.set_restraint(
                            "stage4",
                            self.p_opts.stage4_restraints_number_imp,
                            self.p_opts.stage4_restraints_atoms_imp,
                            self.p_opts.stage4_restraints_forces_imp,
                            "improper",
                            self.p_opts.stage4_restraints_phi0_imp,
                        )
                        j = 0
                        for i in range(
                                int(self.p_opts.stage4_restraints_number_imp)):
                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_imp}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q} {q}{atoms[j+2]}{q} {q}{atoms[j+3]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'phi0':<11} {eq}{constants[i]}",
                                file=fd,
                            )
                            print(f"{outer_space} {'}'}", file=fd)
                            j += 4
                    outer_space, inner_space = identation(1)
                    print(f"{outer_space} {']'}", file=fd)
                print(file=fd)
                # ==============================================================
                # Restraints block end
                # ==============================================================
                outer_space, inner_space = identation(0)
                print(
                    f"{inner_space} {'randomize_velocity.interval':<29} {eq}{'1.0'}",
                    file=fd,
                )
                print(f"{inner_space} {'eneseq.interval':<29} {eq}{'0.3'}",
                      file=fd)
                print(
                    f"{inner_space} {'trajectory.center':<29} {eq}{self.p_opts.stage4_traj_center}",
                    file=fd)
                print(f"{outer_space}{'}'}", file=fd)
                print(file=fd)
            # Stage 5 block
            if self.p_opts.stage5.lower() in ["yes", "on", "true"]:
                outer_space, inner_space = identation(0)
                print(f"{outer_space}{'simulate':<20}{'{'}", file=fd)
                print(
                    f"{inner_space} {'title':<16}{eq}{q}{self.p_opts.stage5_title}{q}",
                    file=fd,
                )
                gpu_text1 = '[["@*.*.annealing"] \'annealing = off temperature = "@*.*.temperature[0][0]"\''
                gpu_text2 = '["==" "-gpu" "@*.*.jlaunch_opt[-1]"] \'ensemble.method = Langevin\']'
                print(f"{inner_space} {'effect_if':<16}{eq}{gpu_text1}",
                      file=fd)
                print(f"{inner_space} {' ':<19}{gpu_text2}", file=fd)
                print(
                    f"{inner_space} {'time':<16}{eq}{self.p_opts.stage5_time}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'temperature':<16}{eq}{self.p_opts.stage5_temp}",
                    file=fd,
                )
                outer_space, inner_space = identation(1)
                print(f"{outer_space} {'ensemble':<16}{eq}{'{'}", file=fd)
                print(
                    f"{inner_space} {'class':<11} {eq}{self.p_opts.stage5_ensemble}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'method':<11} {eq}{self.p_opts.stage5_method}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'thermostat.tau':<11} {eq}{self.p_opts.stage5_thermostat_tau}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'barostat.tau':<11} {eq}{self.p_opts.stage5_barostat_tau}",
                    file=fd,
                )
                outer_space, inner_space = identation(1)
                print(f"{outer_space} {'}'}", file=fd)
                ### Restraints block START ###
                # Block for stage5 positional-restraints
                if (int(self.p_opts.stage5_restraints_number_pos) != 0
                        or int(self.p_opts.stage5_restraints_number_dist) != 0
                        or int(self.p_opts.stage5_restraints_number_ang) != 0
                        or int(self.p_opts.stage5_restraints_number_imp) != 0):
                    outer_space, inner_space = identation(1)
                    print(f"{outer_space} {'restraints.new':<16}{eq}{'['}",
                          file=fd)
                    if int(self.p_opts.stage5_restraints_number_pos) != 0:
                        atoms, forces = self.set_restraint(
                            "stage5",
                            self.p_opts.stage5_restraints_number_pos,
                            self.p_opts.stage5_restraints_atoms_pos,
                            self.p_opts.stage5_restraints_forces_pos,
                            "positional",
                            None,
                        )
                        for i in range(
                                int(self.p_opts.stage5_restraints_number_pos)):
                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_pos}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[i]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]} {forces[i]} {forces[i]}]",
                                file=fd,
                            )
                            print(f"{outer_space} {'}'}", file=fd)
                    if int(self.p_opts.stage5_restraints_number_dist) != 0:
                        # Block for stage5 distance-restraints
                        outer_space, inner_space = identation(1)
                        atoms, forces, constants = self.set_restraint(
                            "stage5",
                            self.p_opts.stage5_restraints_number_dist,
                            self.p_opts.stage5_restraints_atoms_dist,
                            self.p_opts.stage5_restraints_forces_dist,
                            "distance",
                            self.p_opts.stage5_restraints_r0_dist,
                        )
                        j = 0
                        for i in range(
                                int(self.p_opts.stage5_restraints_number_dist)
                        ):

                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_dist}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'r0':<11} {eq}{constants[i]}",
                                file=fd)
                            print(f"{outer_space} {'}'}", file=fd)
                            j += 2
                    if int(self.p_opts.stage5_restraints_number_ang) != 0:
                        # Block for stage5 angle-restraints
                        outer_space, inner_space = identation(1)
                        atoms, forces, constants = self.set_restraint(
                            "stage5",
                            self.p_opts.stage5_restraints_number_ang,
                            self.p_opts.stage5_restraints_atoms_ang,
                            self.p_opts.stage5_restraints_forces_ang,
                            "angle",
                            self.p_opts.stage5_restraints_theta0_ang,
                        )
                        j = 0
                        for i in range(
                                int(self.p_opts.stage5_restraints_number_ang)):
                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_ang}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q} {q}{atoms[j+2]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'theta0':<11} {eq}{constants[i]}",
                                file=fd,
                            )
                            print(f"{outer_space} {'}'}", file=fd)
                            j += 3
                    if int(self.p_opts.stage5_restraints_number_imp) != 0:
                        # Block for stage5 improper-restraints
                        outer_space, inner_space = identation(1)
                        atoms, forces, constants = self.set_restraint(
                            "stage5",
                            self.p_opts.stage5_restraints_number_imp,
                            self.p_opts.stage5_restraints_atoms_imp,
                            self.p_opts.stage5_restraints_forces_imp,
                            "improper",
                            self.p_opts.stage5_restraints_phi0_imp,
                        )
                        j = 0
                        for i in range(
                                int(self.p_opts.stage5_restraints_number_imp)):
                            outer_space, inner_space = identation(2)
                            print(f"{outer_space} {'{'}", file=fd)
                            print(
                                f"{inner_space} {'name':<11} {eq}{self.p_opts.name_imp}",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q} {q}{atoms[j+2]}{q} {q}{atoms[j+3]}{q}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'force_constants':<11} {eq}[{forces[i]}]",
                                file=fd,
                            )
                            print(
                                f"{inner_space} {'phi0':<11} {eq}{constants[i]}",
                                file=fd,
                            )
                            print(f"{outer_space} {'}'}", file=fd)
                            j += 4
                    outer_space, inner_space = identation(1)
                    print(f"{outer_space} {']'}", file=fd)
                print(file=fd)
                ### Restraints block END ###
                outer_space, inner_space = identation(0)
                print(f"{inner_space} {'eneseq.interval':<29} {eq}{'0.3'}",
                      file=fd)
                print(
                    f"{inner_space} {'trajectory.center':<29} {eq}{self.p_opts.stage5_traj_center}",
                    file=fd)
                print(f"{outer_space}{'}'}", file=fd)
                print(file=fd)
            # Additional stages block
            if self.p_opts.additional_stages != None:
                stages = int(self.p_opts.additional_stages)
                stages_time = list(
                    str(self.p_opts.additional_stage_times).split(","))
                stages_temp = list(
                    str(self.p_opts.additional_stage_temps).split(","))
                stages_ensemble = list(
                    str(self.p_opts.additional_stage_ensembles).split(","))
                stages_method = list(
                    str(self.p_opts.additional_stage_methods).split(","))
                stages_thermostat_tau = list(
                    str(self.p_opts.additional_stage_thermostat_tau).split(
                        ","))
                stages_barostat_tau = list(
                    str(self.p_opts.additional_stage_barostat_tau).split(","))

                for stage in range(1, stages + 1):
                    outer_space, inner_space = identation(0)
                    print(f"{outer_space}{'simulate':<20}{'{'}", file=fd)
                    print(
                        f"{inner_space} {'title':<16}{eq}{q}Additional stage = {stage}{q}",
                        file=fd,
                    )
                    gpu_text1 = '[["@*.*.annealing"] \'annealing = off temperature = "@*.*.temperature[0][0]"\''
                    gpu_text2 = '["==" "-gpu" "@*.*.jlaunch_opt[-1]"] \'ensemble.method = Langevin\']'
                    print(f"{inner_space} {'effect_if':<16}{eq}{gpu_text1}",
                          file=fd)
                    print(f"{inner_space} {' ':<19}{gpu_text2}", file=fd)
                    if (self.p_opts.additional_stage_times != None
                            and len(stages_time) == 1):
                        print(
                            f"{inner_space} {'time':<16}{eq}{self.p_opts.additional_stage_times}",
                            file=fd,
                        )
                    elif (self.p_opts.additional_stage_times != None
                          and len(stages_time) == stages):
                        print(
                            f"{inner_space} {'time':<16}{eq}{stages_time[stage-1]}",
                            file=fd,
                        )
                    if (self.p_opts.additional_stage_temps != None
                            and len(stages_time) == 1):
                        print(
                            f"{inner_space} {'temperature':<16}{eq}{self.p_opts.additional_stage_temps}",
                            file=fd,
                        )
                    elif (self.p_opts.additional_stage_temps != None
                          and len(stages_temp) == stages):
                        print(
                            f"{inner_space} {'temperature':<16}{eq}{stages_temp[stage-1]}",
                            file=fd,
                        )
                    outer_space, inner_space = identation(1)

                    print(f"{outer_space} {'ensemble':<16}{eq}{'{'}", file=fd)
                    # Ensemble block
                    if (self.p_opts.additional_stage_ensembles != None
                            and len(stages_ensemble) == 1):
                        print(
                            f"{inner_space} {'class':<11}{eq}{self.p_opts.additional_stage_ensembles}",
                            file=fd,
                        )
                    elif (self.p_opts.additional_stage_ensembles != None
                          and len(stages_ensemble) == stages):
                        print(
                            f"{inner_space} {'class':<11}{eq}{stages_ensemble[stage-1]}",
                            file=fd,
                        )
                    # Method block
                    if (self.p_opts.additional_stage_methods != None
                            and len(stages_method) == 1):
                        print(
                            f"{inner_space} {'method':<11}{eq}{self.p_opts.additional_stage_methods}",
                            file=fd,
                        )
                    elif (self.p_opts.additional_stage_methods != None
                          and len(stages_method) == stages):
                        print(
                            f"{inner_space} {'method':<11}{eq}{stages_method[stage-1]}",
                            file=fd,
                        )
                    # Thermostat block
                    if (self.p_opts.additional_stage_thermostat_tau != None
                            and len(stages_thermostat_tau) == 1):
                        print(
                            f"{inner_space} {'thermostat.tau':<15}{eq}{self.p_opts.additional_stage_thermostat_tau}",
                            file=fd,
                        )
                    elif (self.p_opts.additional_stage_thermostat_tau != None
                          and len(stages_thermostat_tau) == stages):
                        print(
                            f"{inner_space} {'thermostat.tau':<15}{eq}{stages_thermostat_tau[stage-1]}",
                            file=fd,
                        )
                    # Barostat block
                    if (self.p_opts.additional_stage_barostat_tau != None
                            and len(stages_barostat_tau) == 1):
                        print(
                            f"{inner_space} {'barostat.tau':<15}{eq}{self.p_opts.additional_stage_barostat_tau}",
                            file=fd,
                        )
                    elif (self.p_opts.additional_stage_barostat_tau != None
                          and len(stages_barostat_tau) == stages):
                        print(
                            f"{inner_space} {'barostat.tau':<15}{eq}{stages_barostat_tau[stage-1]}",
                            file=fd,
                        )

                    outer_space, inner_space = identation(1)
                    print(f"{outer_space} {'}'}", file=fd)
                    ### Restraints block START ###
                    if (self.p_opts.additional_stage_restraints_number_pos != 0
                            or
                            self.p_opts.additional_stage_restraints_number_dist
                            != 0 or
                            self.p_opts.additional_stage_restraints_number_ang
                            != 0 or
                            self.p_opts.additional_stage_restraints_number_imp
                            != 0):
                        outer_space, inner_space = identation(1)
                        stage = stage - 1
                        header_rest = True
                        if self.p_opts.additional_stage_restraints_number_pos != 0:
                            pos_object = self.set_restraint_multi(
                                stage=stage,
                                stage_name="additional_stage",
                                stage_restraints_number=self.p_opts.
                                additional_stage_restraints_number_pos,
                                stage_restraints_atoms=self.p_opts.
                                additional_stage_restraints_atoms_pos,
                                stage_restraints_forces=self.p_opts.
                                additional_stage_restraints_forces_pos,
                                rest_type="positional",
                                stage_restraints_constants=None,
                            )
                            # Positional restraints variables
                            long_pos = int(pos_object.number)
                            atoms_slice_pos = pos_object.atoms
                            forces_slice_pos = pos_object.forces
                            if long_pos != 0:
                                if header_rest:
                                    print(
                                        f"{outer_space} {'restraints.new':<16}{eq}{'['}",
                                        file=fd,
                                    )
                                    header_rest = False
                                for r in range(long_pos):
                                    print(f"{inner_space} {'{'}", file=fd)
                                    print(
                                        f"{inner_space} {'name':<16}{eq}{self.p_opts.name_pos}",
                                        file=fd,
                                    )
                                    print(
                                        f"{inner_space} {'atoms':<16}{eq}[{q}{atoms_slice_pos[r]}{q}]",
                                        file=fd,
                                    )
                                    print(
                                        f"{inner_space} {'force_constants':<16}{eq}[{forces_slice_pos[r]} {forces_slice_pos[r]} {forces_slice_pos[r]}]",
                                        file=fd,
                                    )
                                    print(f"{inner_space} {'}'}", file=fd)

                        if self.p_opts.additional_stage_restraints_number_dist != 0:
                            dist_object = self.set_restraint_multi(
                                stage=stage,
                                stage_name="additional_stage",
                                stage_restraints_number=self.p_opts.
                                additional_stage_restraints_number_dist,
                                stage_restraints_atoms=self.p_opts.
                                additional_stage_restraints_atoms_dist,
                                stage_restraints_forces=self.p_opts.
                                additional_stage_restraints_forces_dist,
                                rest_type="distance",
                                stage_restraints_constants=self.p_opts.
                                additional_stage_restraints_r0_dist,
                            )
                            # Distance restraints variables
                            long_dist = int(dist_object.number)
                            atoms_slice_dist = dist_object.atoms
                            forces_slice_dist = dist_object.forces
                            constants_slice_dist = dist_object.constants
                            if long_dist != 0:
                                if header_rest:
                                    print(
                                        f"{outer_space} {'restraints.new':<16}{eq}{'['}",
                                        file=fd,
                                    )
                                    header_rest = False
                                for r in range(long_dist):
                                    print(f"{inner_space} {'{'}", file=fd)
                                    print(
                                        f"{inner_space} {'name':<16}{eq}{self.p_opts.name_dist}",
                                        file=fd,
                                    )
                                    print(
                                        f"{inner_space} {'atoms':<16}{eq}[{q}{atoms_slice_dist[r][0]}{q} {q}{atoms_slice_dist[r][1]}{q}]",
                                        file=fd,
                                    )
                                    print(
                                        f"{inner_space} {'force_constants':<16}{eq}[{forces_slice_dist[r]}]",
                                        file=fd,
                                    )
                                    print(
                                        f"{inner_space} {'r0':<16}{eq}{constants_slice_dist[r]}",
                                        file=fd,
                                    )
                                    print(f"{inner_space} {'}'}", file=fd)

                        if self.p_opts.additional_stage_restraints_number_ang != 0:
                            ang_object = self.set_restraint_multi(
                                stage=stage,
                                stage_name="additional_stage",
                                stage_restraints_number=self.p_opts.
                                additional_stage_restraints_number_ang,
                                stage_restraints_atoms=self.p_opts.
                                additional_stage_restraints_atoms_ang,
                                stage_restraints_forces=self.p_opts.
                                additional_stage_restraints_forces_ang,
                                rest_type="angle",
                                stage_restraints_constants=self.p_opts.
                                additional_stage_restraints_theta0_ang,
                            )
                            # Distance restraints variables
                            long_ang = int(ang_object.number)
                            atoms_slice_ang = ang_object.atoms
                            forces_slice_ang = ang_object.forces
                            constants_slice_ang = ang_object.constants
                            if long_ang != 0:
                                if header_rest:
                                    print(
                                        f"{outer_space} {'restraints.new':<16}{eq}{'['}",
                                        file=fd,
                                    )
                                    header_rest = False
                                for r in range(long_ang):
                                    print(f"{inner_space} {'{'}", file=fd)
                                    print(
                                        f"{inner_space} {'name':<16}{eq}{self.p_opts.name_ang}",
                                        file=fd,
                                    )
                                    print(
                                        f"{inner_space} {'atoms':<16}{eq}[{q}{atoms_slice_ang[r][0]}{q} {q}{atoms_slice_ang[r][1]}{q} {q}{atoms_slice_ang[r][2]}{q}]",
                                        file=fd,
                                    )
                                    print(
                                        f"{inner_space} {'force_constants':<16}{eq}[{forces_slice_ang[r]}]",
                                        file=fd,
                                    )
                                    print(
                                        f"{inner_space} {'theta0':<16}{eq}{constants_slice_ang[r]}",
                                        file=fd,
                                    )
                                    print(f"{inner_space} {'}'}", file=fd)

                        if self.p_opts.additional_stage_restraints_number_imp != 0:
                            imp_object = self.set_restraint_multi(
                                stage=stage,
                                stage_name="additional_stage",
                                stage_restraints_number=self.p_opts.
                                additional_stage_restraints_number_imp,
                                stage_restraints_atoms=self.p_opts.
                                additional_stage_restraints_atoms_imp,
                                stage_restraints_forces=self.p_opts.
                                additional_stage_restraints_forces_imp,
                                rest_type="improper",
                                stage_restraints_constants=self.p_opts.
                                additional_stage_restraints_phi0_imp,
                            )
                            # Distance restraints variables
                            long_imp = int(imp_object.number)
                            atoms_slice_imp = imp_object.atoms
                            forces_slice_imp = imp_object.forces
                            constants_slice_imp = imp_object.constants
                            if long_imp != 0:
                                if header_rest:
                                    print(
                                        f"{outer_space} {'restraints.new':<16}{eq}{'['}",
                                        file=fd,
                                    )
                                    header_rest = False
                                for r in range(long_imp):
                                    print(f"{inner_space} {'{'}", file=fd)
                                    print(
                                        f"{inner_space} {'name':<16}{eq}{self.p_opts.name_imp}",
                                        file=fd,
                                    )
                                    print(
                                        f"{inner_space} {'atoms':<16}{eq}[{q}{atoms_slice_imp[r][0]}{q} {q}{atoms_slice_imp[r][1]}{q} {q}{atoms_slice_imp[r][2]}{q} {q}{atoms_slice_imp[r][3]}{q}]",
                                        file=fd,
                                    )
                                    print(
                                        f"{inner_space} {'force_constants':<16}{eq}[{forces_slice_imp[r]}]",
                                        file=fd,
                                    )
                                    print(
                                        f"{inner_space} {'theta0':<16}{eq}{constants_slice_imp[r]}",
                                        file=fd,
                                    )
                                    print(f"{inner_space} {'}'}", file=fd)

                        if header_rest == False:
                            print(f"{outer_space} {']'}", file=fd)
                        #### Restraints block end ####
                    outer_space, inner_space = identation(0)
                    print(file=fd)
                    print(f"{inner_space} {'eneseq.interval':<29} {eq}{'0.3'}",
                          file=fd)
                    print(
                        f"{inner_space} {'trajectory.center':<29} {eq}{self.p_opts.additional_stage_traj_center}",
                        file=fd,
                    )
                    print(f"{outer_space}{'}'}", file=fd)
                    print(file=fd)
            # Stage 6 (Production) block
            if self.p_opts.production.lower() in ["yes", "on", "true"]:
                outer_space, inner_space = identation(0)
                print(f"{outer_space}{'simulate':<20}{'{'}", file=fd)
                print(
                    f"{inner_space} {'cfg_file':<16}{eq}{q}{self.basename}_md.cfg{q}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'jobname':<16}{eq}{q}{'$MASTERJOBNAME'}{q}",
                    file=fd,
                )
                print(f"{inner_space} {'dir':<16}{eq}{q}{'.'}{q}", file=fd)
                print(f"{inner_space} {'compress':<16}{eq}{q}{''}{q}", file=fd)
                print(f"{outer_space}{'}'}", file=fd)
                print(file=fd)
                self.write_cfg_file()

    def set_restraint_multi(
        self,
        stage_name: str,
        stage_restraints_number: str,
        stage_restraints_atoms: str,
        stage_restraints_forces: str,
        rest_type: str,
        constant: Optional[str],
    ) -> None:
        stage_name = stage_name
        number = list(str(stage_restraints_number).split(","))
        list_atoms = list(str(stage_restraints_atoms).split(","))
        list_forces = list(str(stage_restraints_forces).split(","))
        rest_type = rest_type
        list_constants = list(str(constant).split(",")) if constant != None else None
        list_number_long = sum(map(int, number))
        # Block to check the number of restraints for positional restraints
        if rest_type == "positional":
            if len(list_forces) != list_number_long:
                raise ValueError(
                    f"The number of forces in the restraint ({stage_name}_restraints_forces_pos) must be equal to the number of {stage_name}_restraints_number_pos"
                )
            elif len(list_atoms) != list_number_long:
                raise ValueError(
                    f"The number of ASL in the restraint ({stage_name}_restraints_atoms_pos) must be equal to the number of {stage_name}_restraints_number_pos"
                )

            atoms = []
            forces = []

            for i in range(len(number)):
                for j in range(int(number[i])):
                    atoms.append(list_atoms[0])
                    forces.append(list_forces[0])
                    del list_atoms[0]
                    del list_forces[0]
            return atoms, forces, number
        # Block to check the number of restraints for distance restraints
        if rest_type == "distance":
            if len(list_forces) != list_number_long:
                raise ValueError(
                    f"The number of forces in the restraint ({stage_name}_restraints_forces_dist) must be equal to the number of {stage_name}_restraints_number_dist"
                )
            elif len(list_atoms) != list_number_long * 2:
                raise ValueError(
                    f"The number of ASL in the restraint ({stage_name}_restraints_atoms_dist) must be twice the number of {stage_name}_restraints_number_dist"
                )
            elif len(list_constants) != list_number_long:
                raise ValueError(
                    f"The number of constants in the restraint ({stage_name}_restraints_r0_dist) must be equal to the number of {stage_name}_restraints_number_dist"
                )

            atoms = []
            forces = []
            constants = []

            for i in range(len(number)):
                for j in range(int(number[i])):
                    atoms.append(list_atoms[0:2])
                    forces.append(list_forces[0])
                    constants.append(list_constants[0])
                    del list_atoms[0:2]
                    del list_forces[0]
                    del list_constants[0]
            return atoms, forces, number, constants

        # Block to check the number of restraints for angle restraints
        if rest_type == "angle":
            if len(list_forces) != list_number_long:
                raise ValueError(
                    f"The number of forces in the restraint ({stage_name}_restraints_forces_ang) must be equal to the number of {stage_name}_restraints_number_ang"
                )
            elif len(list_atoms) != list_number_long * 3:
                raise ValueError(
                    f"The number of ASL in the restraint ({stage_name}_restraints_atoms_ang) must be three times the number of {stage_name}_restraints_number_ang"
                )
            elif len(list_constants) != list_number_long:
                raise ValueError(
                    f"The number of constants in the restraint ({stage_name}_restraints_r0_ang) must be equal to the number of {stage_name}_restraints_number_ang"
                )

            atoms = []
            forces = []
            constants = []

            for i in range(len(number)):
                for j in range(int(number[i])):
                    atoms.append(list_atoms[0:3])
                    forces.append(list_forces[0])
                    constants.append(list_constants[0])
                    del list_atoms[0:3]
                    del list_forces[0]
                    del list_constants[0]
            return atoms, forces, number, constants

        # Block to check the number of restraints for improper restraints
        if rest_type == "improper":
            if len(list_forces) != list_number_long:
                raise ValueError(
                    f"The number of forces in the restraint ({stage_name}_restraints_forces_imp) must be equal to the number of {stage_name}_restraints_number_imp"
                )
            elif len(list_atoms) != list_number_long * 4:
                raise ValueError(
                    f"The number of ASL in the restraint ({stage_name}_restraints_atoms_imp) must be four times the number of {stage_name}_restraints_number_imp"
                )
            elif len(list_constants) != list_number_long:
                raise ValueError(
                    f"The number of constants in the restraint ({stage_name}_restraints_r0_imp) must be equal to the number of {stage_name}_restraints_number_imp"
                )

            atoms = []
            forces = []
            constants = []

            for i in range(len(number)):
                for j in range(int(number[i])):
                    atoms.append(list_atoms[0:4])
                    forces.append(list_forces[0])
                    constants.append(list_constants[0])
                    del list_atoms[0:4]
                    del list_forces[0]
                    del list_constants[0]
            return atoms, forces, number, constants

    def set_restraint(
        self,
        stage_name: str,
        stage_restraints_number: str,
        stage_restraints_atoms: str,
        stage_restraints_forces: str,
        rest_type: str,
        constant: Optional[str],
    ) -> None:
        stage_name = stage_name
        stage_rest_number = int(stage_restraints_number)
        stage_rest_atoms = list(str(stage_restraints_atoms).split(","))
        stage_rest_forces = list(str(stage_restraints_forces).split(","))
        rest_type = rest_type
        constant = list(str(constant).split(",")) if constant != None else None
        eq = "= "
        q = '"'

        # Block to check the number of restraints for positional restraints
        if rest_type == "positional":
            if len(stage_rest_forces) >= 2 and stage_rest_number != len(
                stage_rest_forces
            ):
                raise ValueError(
                    f"The number of forces in the restraint ({stage_name}_restraints_forces_pos) must be one value or equal to the number of {stage_name}_restraints_number_pos"
                )
            elif len(stage_rest_atoms) >= 2 and stage_rest_number != len(
                stage_rest_atoms
            ):
                raise ValueError(
                    f"The number of ASL in the restraint ({stage_name}_restraints_atoms_pos) must be one value or equal to the number of {stage_name}_restraints_number_pos"
                )
            # Block to set the restraints
            atoms = []
            forces = []
            if stage_rest_number != 0 and len(stage_rest_atoms) != 1:
                for rest in range(1, stage_rest_number + 1):
                    atoms.append(stage_rest_atoms[rest - 1])

            if stage_rest_number != 0 and len(stage_rest_forces) != 1:
                for rest in range(1, stage_rest_number + 1):
                    forces.append(stage_rest_forces[rest - 1])

            if stage_rest_number != 0 and len(stage_rest_atoms) == 1:
                for rest in range(1, stage_rest_number + 1):
                    atoms.append(stage_rest_atoms[0])

            if stage_rest_number != 0 and len(stage_rest_forces) == 1:
                for rest in range(1, stage_rest_number + 1):
                    forces.append(stage_rest_forces[0])

            return atoms, forces

        # Block to check the number of restraints for distance restraints
        if rest_type == "distance":
            if len(stage_rest_atoms) / 2 != stage_rest_number:
                raise ValueError(
                    f"The number of selections ({stage_name}_restraints_atoms_dist) should be twice the number of restraints. Check {stage_name}_restraints_number_dist and {stage_name}_restraints_atoms_dist"
                )
            elif len(stage_rest_forces) >= 2 and stage_rest_number != len(
                stage_rest_forces
            ):
                raise ValueError(
                    f"The number of forces in the restraint ({stage_name}_restraints_forces_dist) must be one value or equal to the number of {stage_name}_restraints_number_dist"
                )
            elif len(constant) >= 2 and stage_rest_number != len(constant):
                raise ValueError(
                    f"The number of constants in the restraint ({stage_name}_restraints_r0_dist) must be one value or equal to the number of {stage_name}_restraints_number_dist"
                )
            # Block to set the restraints
            atoms = []
            forces = []
            constants = []
            if stage_rest_number != 0 and len(stage_rest_atoms) != 1:
                for rest in range(1, stage_rest_number * 2 + 1):
                    atoms.append(stage_rest_atoms[rest - 1])

            if stage_rest_number != 0 and len(stage_rest_forces) != 1:
                for rest in range(1, stage_rest_number + 1):
                    forces.append(stage_rest_forces[rest - 1])

            if stage_rest_number != 0 and len(constant) != 1:
                for rest in range(1, stage_rest_number + 1):
                    constants.append(constant[rest - 1])

            if stage_rest_number != 0 and len(stage_rest_atoms) == 1:
                for rest in range(1, stage_rest_number * 2 + 1):
                    atoms.append(stage_rest_atoms[0])

            if stage_rest_number != 0 and len(stage_rest_forces) == 1:
                for rest in range(1, stage_rest_number + 1):
                    forces.append(stage_rest_forces[0])

            if stage_rest_number != 0 and len(constant) == 1:
                for rest in range(1, stage_rest_number + 1):
                    constants.append(constant[0])

            return atoms, forces, constants
        # Block to check the number of restraints for angle restraints
        if rest_type == "angle":
            if len(stage_rest_atoms) / 3 != stage_rest_number:
                raise ValueError(
                    f"The number of selections ({stage_name}_restraints_atoms_ang) should be three times the number of restraints. Check {stage_name}_restraints_number_ang and {stage_name}_restraints_atoms_ang"
                )
            elif len(stage_rest_forces) >= 2 and stage_rest_number != len(
                stage_rest_forces
            ):
                raise ValueError(
                    f"The number of forces in the restraint ({stage_name}_restraints_forces_ang) must be one value or equal to the number of {stage_name}_restraints_number_ang"
                )
            elif len(constant) >= 2 and stage_rest_number != len(constant):
                raise ValueError(
                    f"The number of constants in the restraint ({stage_name}_restraints_theta0_ang) must be one value or equal to the number of {stage_name}_restraints_number_ang"
                )
            # Block to set the restraints
            atoms = []
            forces = []
            constants = []
            if stage_rest_number != 0 and len(stage_rest_atoms) != 1:
                for rest in range(1, stage_rest_number * 3 + 1):
                    atoms.append(stage_rest_atoms[rest - 1])

            if stage_rest_number != 0 and len(stage_rest_forces) != 1:
                for rest in range(1, stage_rest_number + 1):
                    forces.append(stage_rest_forces[rest - 1])

            if stage_rest_number != 0 and len(constant) != 1:
                for rest in range(1, stage_rest_number + 1):
                    constants.append(constant[rest - 1])

            if stage_rest_number != 0 and len(stage_rest_atoms) == 1:
                for rest in range(1, stage_rest_number * 3 + 1):
                    atoms.append(stage_rest_atoms[0])

            if stage_rest_number != 0 and len(stage_rest_forces) == 1:
                for rest in range(1, stage_rest_number + 1):
                    forces.append(stage_rest_forces[0])

            if stage_rest_number != 0 and len(constant) == 1:
                for rest in range(1, stage_rest_number + 1):
                    constants.append(constant[0])

            return atoms, forces, constants
        # Block to check the number of restraints for improper restraints
        if rest_type == "improper":
            if len(stage_rest_atoms) / 4 != stage_rest_number:
                raise ValueError(
                    f"The number of selections ({stage_name}_restraints_atoms_imp) should be three times the number of restraints. Check {stage_name}_restraints_number_imp and {stage_name}_restraints_atoms_imp"
                )
            elif len(stage_rest_forces) >= 2 and stage_rest_number != len(
                stage_rest_forces
            ):
                raise ValueError(
                    f"The number of forces in the restraint ({stage_name}_restraints_forces_imp) must be one value or equal to the number of {stage_name}_restraints_number_imp"
                )
            elif len(constant) >= 2 and stage_rest_number != len(constant):
                raise ValueError(
                    f"The number of constants in the restraint ({stage_name}_restraints_theta0_imp) must be one value or equal to the number of {stage_name}_restraints_number_imp"
                )
            # Block to set the restraints
            atoms = []
            forces = []
            constants = []
            if stage_rest_number != 0 and len(stage_rest_atoms) != 1:
                for rest in range(1, stage_rest_number * 4 + 1):
                    atoms.append(stage_rest_atoms[rest - 1])

            if stage_rest_number != 0 and len(stage_rest_forces) != 1:
                for rest in range(1, stage_rest_number + 1):
                    forces.append(stage_rest_forces[rest - 1])

            if stage_rest_number != 0 and len(constant) != 1:
                for rest in range(1, stage_rest_number + 1):
                    constants.append(constant[rest - 1])

            if stage_rest_number != 0 and len(stage_rest_atoms) == 1:
                for rest in range(1, stage_rest_number * 4 + 1):
                    atoms.append(stage_rest_atoms[0])

            if stage_rest_number != 0 and len(stage_rest_forces) == 1:
                for rest in range(1, stage_rest_number + 1):
                    forces.append(stage_rest_forces[0])

            if stage_rest_number != 0 and len(constant) == 1:
                for rest in range(1, stage_rest_number + 1):
                    constants.append(constant[0])

            return atoms, forces, constants

    def write_cfg_file(self) -> None:
        eq = "= "
        q = '"'
        path_preparation = str(self.basename + "_md.cfg")
        with open(path_preparation, "w", encoding="utf8") as fd:
            outer_space, inner_space = identation(0)
            print(f"{outer_space}{'annealing':<20}{eq}{'false'}", file=fd)
            print(f"{outer_space}{'backend':<20}{eq}{'{'}", file=fd)
            print(f"{outer_space}{'}'}", file=fd)
            print(
                f"{outer_space}{'bigger_rclone':<20}{eq}{self.p_opts.production_bigger_rclone}",
                file=fd,
            )
            print(f"{outer_space}{'checkpt':<20}{eq}{'{'}", file=fd)
            print(
                f"{inner_space} {'first':<16}{eq}{self.p_opts.production_checkpt_first}",
                file=fd,
            )
            print(
                f"{inner_space} {'interval':<16}{eq}{self.p_opts.production_checkpt_interval}",
                file=fd,
            )
            print(f"{inner_space} {'name':<16}{eq}{q}{'$JOBNAME.cpt'}{q}", file=fd)
            print(
                f"{inner_space} {'write_last_step':<16}{eq}{self.p_opts.production_write_last_step}",
                file=fd,
            )
            print(f"{outer_space}{'}'}", file=fd)
            print(f"{outer_space}{'cpu':<20}{eq}{'1'}", file=fd)
            print(
                f"{outer_space}{'cutoff_radius':<20}{eq}{self.p_opts.production_cutoff}",
                file=fd,
            )
            print(
                f"{outer_space}{'elapsed_time':<20}{eq}{self.p_opts.production_elapsed_time}",
                file=fd,
            )
            print(
                f"{outer_space}{'energy_group':<20}{eq}{self.p_opts.production_energy_group}",
                file=fd,
            )
            print(f"{outer_space}{'eneseq':<20}{eq}{'{'}", file=fd)
            print(
                f"{inner_space} {'first':<16}{eq}{self.p_opts.production_eneseq_first}",
                file=fd,
            )
            print(
                f"{inner_space} {'interval':<16}{eq}{self.p_opts.production_eneseq_interval}",
                file=fd,
            )
            print(f"{outer_space}{'}'}", file=fd)
            print(f"{outer_space}{'ensemble':<20}{eq}{'{'}", file=fd)
            outer_space, inner_space = identation(1)
            print(
                f"{outer_space} {'class':<16}{eq}{self.p_opts.production_ensemble}",
                file=fd,
            )
            print(
                f"{outer_space} {'method':<16}{eq}{self.p_opts.production_method}",
                file=fd,
            )
            print(f"{outer_space} {'barostat':<16}{eq}{'{'}", file=fd)
            print(
                f"{inner_space} {'tau':<12}{eq}{self.p_opts.production_barostat_tau}",
                file=fd,
            )
            print(f"{outer_space} {'}'}", file=fd)
            print(f"{outer_space} {'thermostat':<16}{eq}{'{'}", file=fd)
            print(
                f"{inner_space} {'tau':<12}{eq}{self.p_opts.production_thermostat_tau}",
                file=fd,
            )
            print(f"{outer_space} {'}'}", file=fd)
            outer_space, inner_space = identation(0)
            print(f"{outer_space}{'}'}", file=fd)
            print(
                f"{outer_space}{'glue':<20}{eq}{self.p_opts.production_glue}", file=fd
            )
            print(f"{outer_space}{'maeff_output':<20}{eq}{'{'}", file=fd)
            print(
                f"{inner_space} {'first':<16}{eq}{self.p_opts.production_maeff_first}",
                file=fd,
            )
            print(
                f"{inner_space} {'interval':<16}{eq}{self.p_opts.production_maeff_interval}",
                file=fd,
            )
            name_maeff = '"$JOBNAME$[_replica$REPLICA$]-out.cms"'
            name_trjdir = '"$JOBNAME$[_replica$REPLICA$]_trj"'
            print(f"{inner_space} {'name':<16}{eq}{name_maeff}", file=fd)
            print(
                f"{inner_space} {'periodicfix':<16}{eq}{self.p_opts.production_maeff_periodicfix}",
                file=fd,
            )
            print(f"{inner_space} {'trjdir':<16}{eq}{name_trjdir}", file=fd)
            print(f"{outer_space}{'}'}", file=fd)
            print(
                f"{outer_space}{'meta':<20}{eq}{self.p_opts.production_meta}", file=fd
            )
            print(f"{outer_space}{'meta_file':<20}{eq}{'?'}", file=fd)
            print(
                f"{outer_space}{'pressure':<20}{eq}[{self.p_opts.production_pressure} {self.p_opts.production_pressure_type}]",
                file=fd,
            )
            print(f"{outer_space}{'randomize_velocity':<20}{eq}{'{'}", file=fd)
            print(
                f"{inner_space} {'first':<16}{eq}{self.p_opts.production_randomize_vel_first}",
                file=fd,
            )
            print(
                f"{inner_space} {'interval':<16}{eq}{self.p_opts.production_randomize_vel_interval}",
                file=fd,
            )
            print(
                f"{inner_space} {'seed':<16}{eq}{self.p_opts.production_randomize_vel_seed}",
                file=fd,
            )
            print(
                f"{inner_space} {'temperature':<16}{eq}{q}{'@*.temperature'}{q}",
                file=fd,
            )
            print(f"{outer_space}{'}'}", file=fd)
            ### Restraints block START ###
            # Block for production positional-restraints
            if (
                int(self.p_opts.production_restraints_number_pos) != 0
                or int(self.p_opts.production_restraints_number_dist) != 0
                or int(self.p_opts.production_restraints_number_ang) != 0
                or int(self.p_opts.production_restraints_number_imp) != 0
            ):
                outer_space, inner_space = identation(0)
                print(f"{outer_space}{'restraints.new':<16}{eq}{'['}", file=fd)
                if int(self.p_opts.production_restraints_number_pos) != 0:
                    atoms, forces = self.set_restraint(
                        "production",
                        self.p_opts.production_restraints_number_pos,
                        self.p_opts.production_restraints_atoms_pos,
                        self.p_opts.production_restraints_forces_pos,
                        "positional",
                        None,
                    )
                    for i in range(int(self.p_opts.production_restraints_number_pos)):
                        outer_space, inner_space = identation(1)
                        print(f"{outer_space}{'{'}", file=fd)
                        print(
                            f"{inner_space}{'name':<11} {eq}{self.p_opts.name_pos}",
                            file=fd,
                        )
                        print(
                            f"{inner_space}{'atoms':<11} {eq}[{q}{atoms[i]}{q}]",
                            file=fd,
                        )
                        print(
                            f"{inner_space}{'force_constants':<11} {eq}[{forces[i]} {forces[i]} {forces[i]}]",
                            file=fd,
                        )
                        print(f"{outer_space} {'}'}", file=fd)
                if int(self.p_opts.production_restraints_number_dist) != 0:
                    # Block for production distance-restraints
                    outer_space, inner_space = identation(0)
                    atoms, forces, constants = self.set_restraint(
                        "production",
                        self.p_opts.production_restraints_number_dist,
                        self.p_opts.production_restraints_atoms_dist,
                        self.p_opts.production_restraints_forces_dist,
                        "distance",
                        self.p_opts.production_restraints_r0_dist,
                    )
                    j = 0
                    for i in range(int(self.p_opts.production_restraints_number_dist)):
                        outer_space, inner_space = identation(1)
                        print(f"{outer_space}{'{'}", file=fd)
                        print(
                            f"{inner_space}{'name':<11} {eq}{self.p_opts.name_dist}",
                            file=fd,
                        )
                        print(
                            f"{inner_space}{'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q}]",
                            file=fd,
                        )
                        print(
                            f"{inner_space}{'force_constants':<11} {eq}[{forces[i]}]",
                            file=fd,
                        )
                        print(f"{inner_space}{'r0':<11} {eq}{constants[i]}", file=fd)
                        print(f"{outer_space}{'}'}", file=fd)
                        j += 2
                if int(self.p_opts.production_restraints_number_ang) != 0:
                    # Block for production angle-restraints
                    outer_space, inner_space = identation(0)
                    atoms, forces, constants = self.set_restraint(
                        "production",
                        self.p_opts.production_restraints_number_ang,
                        self.p_opts.production_restraints_atoms_ang,
                        self.p_opts.production_restraints_forces_ang,
                        "angle",
                        self.p_opts.production_restraints_theta0_ang,
                    )
                    j = 0
                    for i in range(int(self.p_opts.production_restraints_number_ang)):
                        outer_space, inner_space = identation(1)
                        print(f"{outer_space}{'{'}", file=fd)
                        print(
                            f"{inner_space}{'name':<11} {eq}{self.p_opts.name_ang}",
                            file=fd,
                        )
                        print(
                            f"{inner_space}{'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q} {q}{atoms[j+2]}{q}]",
                            file=fd,
                        )
                        print(
                            f"{inner_space}{'force_constants':<11} {eq}[{forces[i]}]",
                            file=fd,
                        )
                        print(
                            f"{inner_space}{'theta0':<11} {eq}{constants[i]}",
                            file=fd,
                        )
                        print(f"{outer_space}{'}'}", file=fd)
                        j += 3
                if int(self.p_opts.production_restraints_number_imp) != 0:
                    # Block for production improper-restraints
                    outer_space, inner_space = identation(0)
                    atoms, forces, constants = self.set_restraint(
                        "production",
                        self.p_opts.production_restraints_number_imp,
                        self.p_opts.production_restraints_atoms_imp,
                        self.p_opts.production_restraints_forces_imp,
                        "improper",
                        self.p_opts.production_restraints_phi0_imp,
                    )
                    j = 0
                    for i in range(int(self.p_opts.production_restraints_number_imp)):
                        outer_space, inner_space = identation(1)
                        print(f"{outer_space}{'{'}", file=fd)
                        print(
                            f"{inner_space}{'name':<11} {eq}{self.p_opts.name_imp}",
                            file=fd,
                        )
                        print(
                            f"{inner_space}{'atoms':<11} {eq}[{q}{atoms[j]}{q} {q}{atoms[j+1]}{q} {q}{atoms[j+2]}{q} {q}{atoms[j+3]}{q}]",
                            file=fd,
                        )
                        print(
                            f"{inner_space}{'force_constants':<11} {eq}[{forces[i]}]",
                            file=fd,
                        )
                        print(
                            f"{inner_space}{'phi0':<11} {eq}{constants[i]}",
                            file=fd,
                        )
                        print(f"{outer_space}{'}'}", file=fd)
                        j += 4
                print(f"{outer_space}{']'}", file=fd)
            print(file=fd)
            ### Restraints block END ###

            ##print(f"{outer_space}{'restrain':<20}{eq}{'{'}", file=fd)
            ##print(
            ##    f"{inner_space} {'atom':<15} {eq}{q}{self.p_opts.production_restraint}{q}",
            ##    file=fd,
            ##)
            ##if self.p_opts.production_restraint_force is None:
            ##    raise ValueError(
            ##        "You must set the force for the restraint")
            ##print(
            ##    f"{inner_space} {'force_constant':<15} {eq}{self.p_opts.production_restraint_force}",
            ##    file=fd,
            ##)
            ##print(f"{outer_space}{'}'}", file=fd)
            outer_space, inner_space = identation(0)
            print(f"{outer_space}{'simbox':<20}{eq}{'{'}", file=fd)
            print(
                f"{inner_space} {'first':<16}{eq}{self.p_opts.production_simbox_first}",
                file=fd,
            )
            print(
                f"{inner_space} {'interval':<16}{eq}{self.p_opts.production_simbox_interval}",
                file=fd,
            )
            name_simbox = '"$JOBNAME$[_replica$REPLICA$]_simbox.dat"'
            print(f"{inner_space} {'name':<16}{eq}{name_simbox}", file=fd)
            print(f"{outer_space}{'}'}", file=fd)
            print(
                f"{outer_space}{'surface_tension':<20}{eq}{self.p_opts.production_surface_tension}",
                file=fd,
            )
            print(
                f"{outer_space}{'taper':<20}{eq}{self.p_opts.production_taper}",
                file=fd,
            )
            print(
                f"{outer_space}{'temperature':<20}{eq}[ [{self.p_opts.production_temp} {self.p_opts.production_temp_group}] ]",
                file=fd,
            )
            print(
                f"{outer_space}{'time':<20}{eq}{self.p_opts.production_time}",
                file=fd,
            )
            print(
                f"{outer_space}{'timestep':<20}{eq}[{self.p_opts.production_timestep_bonded} {self.p_opts.production_timestep_near} {self.p_opts.production_timestep_far}]",
                file=fd,
            )
            print(f"{outer_space}{'trajectory':<20}{eq}{'{'}", file=fd)
            print(
                f"{inner_space} {'center':<16}{eq}[{self.p_opts.production_traj_center}]",
                file=fd,
            )
            print(
                f"{inner_space} {'first':<16}{eq}{self.p_opts.production_traj_first}",
                file=fd,
            )
            print(
                f"{inner_space} {'format':<16}{eq}{self.p_opts.production_traj_format}",
                file=fd,
            )
            print(
                f"{inner_space} {'frames_per_file':<16}{eq}{self.p_opts.production_traj_frames_per_file}",
                file=fd,
            )
            print(
                f"{inner_space} {'interval':<16}{eq}{self.p_opts.production_traj_interval}",
                file=fd,
            )
            name_traj = '"$JOBNAME$[_replica$REPLICA$]_trj"'
            print(f"{inner_space} {'name':<16}{eq}{name_traj}", file=fd)
            print(
                f"{inner_space} {'periodicfix':<16}{eq}{self.p_opts.production_traj_periodicfix}",
                file=fd,
            )
            print(
                f"{inner_space} {'write_velocity':<16}{eq}{self.p_opts.production_traj_write_velocity}",
                file=fd,
            )
            print(f"{outer_space}{'}'}", file=fd)

    def write_protocol_sh(self) -> None:
        if self.builder_opts.windows.lower() in [
                "yes",
                "on",
                "true",
        ]:
            executable = os.path.join(self.builder_opts.desmond_path,
                                      "utilities/multisim.exe")
        else:
            executable = os.path.join(self.builder_opts.desmond_path,
                                      "utilities/multisim")
        path_preparation_sh = str(self.basename + "_md.sh")
        input_msj = str(self.basename + "_md.msj")
        input_cms = self.basename + "_preparation" + "-out.cms"
        input_cfg = self.basename + "_md.cfg"
        gpu_opts = 'stage[1].set_family.md.jlaunch_opt=["-gpu"]'
        args1 = f"-HOST localhost -JOBNAME {self.basename}_md -maxjob 1 -cpu 1 -m {input_msj} -c {input_cfg} {input_cms}"
        args2 = f"-mode umbrella -set '{gpu_opts}' -o {self.basename}_md-out.cms -LOCAL"
        with open(path_preparation_sh, "w", encoding="utf8") as fd:
            print(executable, args1, args2, file=fd)

    def run_protocol(self) -> None:
        path_preparation_sh = str(self.basename + "_md.sh")
        subprocess.run(["bash", path_preparation_sh])


def parse_args(argv):
    conf_parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
    )
    conf_parser.add_argument(
        "-i", "--input", help="Specify a configuration file", metavar="FILE"
    )
    args, remaining_argv = conf_parser.parse_known_args()

    defaults = {"desmond_path": "$SCHRODINGER"}

    print(
        """
    █▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀█
    █   ┌┬┐┌─┐┌─┐┌┬┐┌─┐┌┐┌┌┬┐  ┌┐ ┬ ┬┬┬  ┌┬┐┌─┐┬─┐       █
    █    ││├┤ └─┐││││ ││││ ││  ├┴┐│ │││   ││├┤ ├┬┘       █
    █   ─┴┘└─┘└─┘┴ ┴└─┘┘└┘─┴┘  └─┘└─┘┴┴─┘─┴┘└─┘┴└─ v1.0  █
    █▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄█

    Prepare a system for MD simulation in desmond starting from a .mae file.
    1. Pass a configuration file to the script.
    2. Execute as: python run_schrod_mds.py -i options_file.dat
    3. The full list of configuration file options can be found at: 
    https://github.com/maurobedoya/desmond_builder

    "If this script is useful for you, consider giving acknowledgments comments in the publication."
    Contact: 
    Mauricio Bedoya
    maurobedoyat@gmail.com"
    
    """
    )

    try:
        if args.input:
            config = configparser.ConfigParser()
            config.read([args.input])
            defaults.update(dict(config.items("settings")))
            # defaults.update(dict(config.items("build_geometry")))
    except:
        raise ValueError("You must specify a configuration file")

    # Parse rest of arguments
    # Don't suppress add_help here so it will handle -h
    parser = argparse.ArgumentParser(
        # Inherit options from config_parser
        parents=[conf_parser]
    )
    parser.set_defaults(**defaults)
    args = parser.parse_args(remaining_argv)
    build_opts = {}
    build_opts.update(dict(config.items("build_geometry")))
    filename = vars(args)["file"].split("/")[-1]
    file = vars(args)["file"]
    desmond_path = vars(args)["desmond_path"]
    file_path = os.path.abspath(vars(args)["file"])
    if "windows" in vars(args):
        windows = vars(args)["windows"]
    else:
        windows = "false"

    protocol_opts = {}
    protocol_opts.update(dict(config.items("protocol")))

    return (
        Args(**vars(args)),
        BuilderOptions(build_opts, file, filename, desmond_path, file_path,
                       windows),
        file_path,
        ProtocolOptions(protocol_opts),
    )


@dataclass
class maefile:
    file: str
    charge: int
    desmond_path: str


class ReadMaefile:
    def __init__(self, file: str, desmond_path, windows: str) -> None:
        self.file = os.path.relpath(file)
        self.charge: int = 0
        self.desmond_path = desmond_path
        self.windows = windows

    def get_charge(self):
        if self.windows.lower() in [
                "yes",
                "on",
                "true",
        ]:
            executable = os.path.join(self.desmond_path, "run.exe")
        else:
            executable = os.path.join(self.desmond_path, "run")
        command = "schrod_script.py"
        arg1 = "-i"
        arg2 = "-get"
        option = "charge"
        charge = subprocess.check_output(
            [executable, command, arg1, self.file, arg2, option],
            universal_newlines=True,
        )[:-1]
        return charge

    def get_atoms_number(self, file, asl):
        if self.windows.lower() in [
                "yes",
                "on",
                "true",
        ]:
            executable = os.path.join(self.desmond_path, "run.exe")
        else:
            executable = os.path.join(self.desmond_path, "run")
        command = "schrod_script.py"
        arg1 = "-i"
        arg2 = "-get"
        option = "atoms_number"
        arg3 = "-asl"
        arg4 = asl
        atoms = subprocess.check_output(
            [executable, command, arg1, self.file, arg2, option, arg3, arg4],
            universal_newlines=True,
        )[:-1]
        return atoms.replace(",", "")


def check_folder_analysis(folder_name: str):
    if path.isdir(folder_name):
        raise ValueError(
            f"Folder '{folder_name}' exists, remove it before to continue")
    else:
        os.makedirs(folder_name)
        os.chdir(folder_name)


class Builder:
    def __init__(
        self,
        options: BuilderOptions,
        charge: int,
        atoms_number: Optional[str] = None,
    ) -> None:
        self.options = options
        self.counterions = options.counterions
        self.charge = charge
        self.basename = options.basename
        self.file = options.file
        self.desmond_path = options.desmond_path
        self.ions_away = options.ions_away
        self.atoms_number = atoms_number

    def write_input(self) -> None:
        path_preparation = str(self.basename + "_preparation.msj")
        outer_space, inner_space = identation(0)
        eq = "= "
        with open(path_preparation, "w", encoding="utf8") as fd:
            print("Preparing input files for system building...")
            print(f"{outer_space}task {'{'}", file=fd)
            print(f"{inner_space} task {eq}", '"desmond:auto"', file=fd)
            print(f"{outer_space}{'}'}", file=fd)
            print(file=fd)
            print(f"{outer_space}build_geometry {'{'}", file=fd)
            # add_counterions block
            outer_space, inner_space = identation(1)
            if (
                self.counterions.lower() in ["yes", "on", "true"]
                and int(self.charge) != 0
            ):
                print(f"{outer_space} {'add_counterion = {'}", file=fd)
                if int(self.charge) < 0:
                    self.options.ion = self.options.counterions_positive_ion
                elif int(self.charge) > 0:
                    self.options.ion = self.options.counterions_negative_ion
                print(f"{inner_space} {'ion':<16}{eq}{self.options.ion}", file=fd)
                print(f"{inner_space} {'number':<16}{eq}{self.options.number}", file=fd)
                print(f"{outer_space} {'}'}", file=fd)
            # box block (mandatory)
            print(f"{outer_space} {'box = {'}", file=fd)
            print(f"{inner_space} {'shape':<16}{eq}{self.options.shape}", file=fd)
            print(f"{inner_space} {'size':<16}{eq}[{self.options.size}]", file=fd)
            print(
                f"{inner_space} {'size_type':<16}{eq}{self.options.size_type}", file=fd
            )
            print(f"{outer_space} {'}'}", file=fd)
            outer_space, inner_space = identation(0)
            # ions_away block
            if self.ions_away.lower() in ["yes", "on", "true"]:
                print(
                    f"{inner_space} {'ion_awaydistance':<20}{eq}{self.options.ion_awaydistance}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'ion_awayfrom':<20}{eq}{self.atoms_number}",
                    file=fd,
                )
            # override_forcefield block
            print(
                f"{inner_space} {'override_forcefield':<20}{eq}{self.options.override_forcefield}",
                file=fd,
            )
            # rezero_system block
            print(
                f"{inner_space} {'rezero_system':<20}{eq}{self.options.rezero_system}",
                file=fd,
            )
            # salt block
            outer_space, inner_space = identation(1)
            if self.counterions.lower() in ["yes", "on", "true"]:
                print(f"{outer_space} {'salt = {'}", file=fd)
                print(
                    f"{inner_space} {'concentration':<16}{eq}{self.options.concentration}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'negative_ion':<16}{eq}{self.options.negative_ion}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'positive_ion':<16}{eq}{self.options.positive_ion}",
                    file=fd,
                )
                print(f"{outer_space} {'}'}", file=fd)
            # solvent block
            outer_space, inner_space = identation(0)
            print(f"{inner_space} {'solvent':<20}{eq}{self.options.solvent}", file=fd)
            print(f"{outer_space}{'}'}", file=fd)
            print(file=fd)
            # assign_forcefield block
            outer_space, inner_space = identation(0)
            print(f"{outer_space}{'assign_forcefield {'}", file=fd)
            print(
                f"{inner_space} {'forcefield':<20}{eq}{self.options.forcefield}",
                file=fd,
            )
            print(f"{outer_space}{'}'}", file=fd)

    def write_preparation_sh(self) -> None:
        q = '"'
        self.suffix_prep = "_preparation"
        if self.options.windows.lower() in [
                "yes",
                "on",
                "true",
        ]:
            executable = os.path.join(self.desmond_path,
                                      "utilities/multisim.exe")
        else:
            executable = os.path.join(self.desmond_path, "utilities/multisim")
        path_preparation_sh = str(self.basename + f"{self.suffix_prep}.sh")
        path_preparation = str(self.basename + f"{self.suffix_prep}.msj")
        input_mae = self.options.file_path
        os.system(f"cp {input_mae} ./")
        actual_mae = self.options.basename + ".mae"
        basename = self.options.basename
        args = f"-HOST {q}localhost{q} -maxjob 1 -JOBNAME preparation -m {path_preparation} {actual_mae} -o {basename}{self.suffix_prep}-out.cms -WAIT"
        with open(path_preparation_sh, "w", encoding="utf8") as fd:
            print(executable, args, file=fd)

    def run_preparation(self) -> None:
        path_preparation_sh = str(self.basename + f"{self.suffix_prep}.sh")
        print("Preparing system...")
        subprocess.run(["bash", path_preparation_sh])


def write_schrod_script() -> None:
    suffix = "schrod_script"
    path_script = f"{suffix}.py"
    script = """

import argparse
import sys

from schrodinger.structure import StructureReader
from schrodinger.application.jaguar.utils import get_total_charge
from schrodinger.structutils.analyze import evaluate_asl


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", help="Input file.")
    parser.add_argument("-get", help="Option.")
    parser.add_argument("-asl", help="ASL.")
    opts = parser.parse_args(argv)
    return vars(opts)


def get_total_charge_function(reader):
    for st in reader:
        charge = get_total_charge(st)
        return charge


def get_atoms_number_function(reader, asl):
    for st in reader:
        idxs = evaluate_asl(st, asl)
        return idxs


def main(argv):
    opts = parse_args(argv)
    maefile = opts["input"]
    reader = StructureReader(maefile)
    if opts["get"] == "charge":
        charge = get_total_charge_function(reader)
        print(charge)
    elif opts["get"] == "atoms_number":
        atoms = get_atoms_number_function(reader, opts["asl"])
        print(atoms)


if __name__ == "__main__":
    main(sys.argv[1:])
"""
    with open(path_script, "w", encoding="utf8") as fd:
        print(script, file=fd)


def main(argv):
    opts, build_opts, file_path, protocol_opts = parse_args(argv)
    # Prepare the system
    check_folder_analysis(opts.workdir)
    file = file_path
    basename = os.path.basename(file).split(".")[0]
    system = ReadMaefile(file, opts.desmond_path, opts.windows)
    write_schrod_script()
    charge = system.get_charge()
    if build_opts.ions_away.lower() in ["yes", "on", "true"]:
        atoms_number = system.get_atoms_number(file, build_opts.ion_awayfrom)
        builder = Builder(build_opts, charge, atoms_number)
        builder.write_input()
        builder.write_preparation_sh()
    else:
        builder = Builder(build_opts, charge)
        builder.write_input()
        builder.write_preparation_sh()
    # Run the preparation
    if str(protocol_opts.run_preparation).lower() in ["yes", "on", "true"]:
        builder.run_preparation()
    output_name_builder = os.path.join(opts.workdir, basename + "_system-out.cms")
    # Simulation protocol
    protocol = Protocol(output_name_builder, build_opts, protocol_opts)
    protocol.write()
    protocol.write_protocol_sh()
    # Run the simulation protocol
    if str(protocol_opts.run_protocols).lower() in ["yes", "on", "true"]:
        protocol.run_protocol()


if __name__ == "__main__":
    main(sys.argv[1:])
