#!/usr/bin/env python3
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
import numpy as np
from os import path, PathLike, supports_fd, write

from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple
import configparser
import random


@dataclass
class Args:
    workdir: str
    input: str
    file: str
    schrod_path: str
    # file_path: os.PathLike
    # output: str


@dataclass
class BuilderOptions:
    basename: Optional[str] = "None"
    counterions: Optional[str] = None
    ion: Optional[str] = None
    number: Optional[str] = "neutralize_system"

    shape: Optional[str] = "orthorhombic"
    size: Optional[str] = "10.0 10.0 10.0"
    size_type: Optional[str] = "buffer"

    # box: Optional[str] = True

    ions_away: Optional[str] = None
    ion_awaydistance: Optional[str] = "5.0"
    ion_awayfrom: Optional[str] = "protein"

    override_forcefield: Optional[str] = "OPLS_2005"
    rezero_system: Optional[str] = "true"

    salt: Optional[str] = "false"
    concentration: Optional[float] = 0.15
    negative_ion: Optional[str] = "Cl"
    positive_ion: Optional[str] = "Na"

    solvent: Optional[str] = "SPC"

    forcefield: Optional[str] = "OPLS_2005"

    def __init__(
        self,
        opts: Dict[str, str],
        file: str,
        filename: str,
        schrod_path: str,
        file_path: str,
    ) -> None:
        self.opts = opts
        self.file = file
        self.filename = filename
        self.basename = str(filename).split(".")[0]
        self.schrod_path = schrod_path
        self.file_path = file_path
        for key in self.opts:
            setattr(self, key, self.opts[key])

    def __getattr__(self, item):
        if item not in self.opts:
            return None
        return self.opts[item]


# @dataclass
class ProtocolOptions:
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
    stage4_temp: Optional[float] = 300
    stage5_temp: Optional[float] = 300
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

    stage1_restraint: Optional[str] = "solute_heavy_atom"
    stage2_restraint: Optional[str] = "solute_heavy_atom"
    stage3_restraint: Optional[str] = "solute_heavy_atom"
    stage4_restraint: Optional[str] = "solute_heavy_atom"
    stage5_restraint: Optional[str] = None
    production_restraint: Optional[str] = None

    stage1_restraint_force: Optional[float] = 50.0
    stage2_restraint_force: Optional[float] = 50.0
    stage3_restraint_force: Optional[float] = 50.0
    stage4_restraint_force: Optional[float] = 50.0
    stage5_restraint_force: Optional[float] = None
    production_restraint_force: Optional[float] = None

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

    # Additional stages
    additional_stages: Optional[int] = None
    additional_stage_times: Optional[str] = None
    additional_stage_temps: Optional[str] = None
    additional_stage_ensembles: Optional[str] = None
    additional_stage_methods: Optional[str] = None
    additional_stage_thermostat_tau: Optional[float] = 0.1
    additional_stage_barostat_tau: Optional[float] = 2.0
    additional_stage_restraints: Optional[str] = None
    additional_stage_restraints_forces: Optional[str] = None

    # Run protocols
    run_preparation: Optional[str] = "false"
    run_protocols: Optional[str] = "false"

    stage1_title: Optional[
        str
    ] = f"{stage1_method} {stage1_ensemble}, T = {stage1_temp} K, restraints on {stage1_restraint} of {stage1_restraint_force} {stage1_time}ps"
    stage2_title: Optional[
        str
    ] = f"{stage2_method} {stage2_ensemble}, T = {stage2_temp} K, restraints on {stage2_restraint} of {stage2_restraint_force}, {stage2_time}ps"
    stage3_title: Optional[
        str
    ] = f"{stage3_method} {stage3_ensemble}, T = {stage3_temp} K, restraints on {stage3_restraint} of {stage3_restraint_force}, {stage3_time}ps"
    stage4_title: Optional[
        str
    ] = f"{stage4_method} {stage4_ensemble}, T = {stage4_temp} K, restraints on {stage4_restraint} of {stage4_restraint_force}, {stage4_time}ps"
    stage5_title: Optional[
        str
    ] = f"{stage5_method} {stage5_ensemble}, T = {stage5_temp} K, restraints on {stage5_restraint} of {stage5_restraint_force}, {stage5_time}ps"
    production_title: Optional[
        str
    ] = f"Production run: {production_method} {production_ensemble}, T = {production_temp} K, restraints on {production_restraint} of {production_restraint_force}, {production_time}ps"

    def __init__(self, opts: Dict) -> None:
        self.opts = opts
        for key in self.opts:
            setattr(self, key, self.opts[key])

    def __getattr__(self, item):
        if item not in self.opts:
            return None
        return self.opts[item]


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
            print(f"{inner_space} {'task':<16}{eq}{q}{'desmond:auto'}{q}", file=fd)
            print(f"{inner_space} {'set_family':<16}{eq}{'{'}", file=fd)

            outer_space, inner_space = identation(2)
            print(f"{outer_space} {'desmond':<12}{eq}{'{'}", file=fd)
            print(f"{inner_space} {'checkpt.write_last_step':<16}{eq}{'no'}", file=fd)
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
                # Restraints block
            if self.p_opts.stage1_restraint.lower() is not None:
                outer_space, inner_space = identation(1)
                print(f"{outer_space} {'restrain':<16}{eq}{'{'}", file=fd)
                print(
                    f"{inner_space} {'atom':<16} {eq}{q}{self.p_opts.stage1_restraint}{q}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'force_constant':<16} {eq}{self.p_opts.stage1_restraint_force}",
                    file=fd,
                )
                print(f"{outer_space} {'}'}", file=fd)
                outer_space, inner_space = identation(0)
                print(f"{outer_space}{'}'}", file=fd)
                print(file=fd)
            # Stage 2 block
            if self.p_opts.stage2.lower() in ["yes", "on", "true"]:
                outer_space, inner_space = identation(0)
                print(f"{outer_space}{'simulate':<20}{'{'}", file=fd)
                print(
                    f"{inner_space} {'title':<16}{eq}{q}{self.p_opts.stage2_title}{q}",
                    file=fd,
                )
                gpu_text = '[["==" "-gpu" "@*.*.jlaunch_opt[-1]"] \'ensemble.method = Langevin\']'
                print(f"{inner_space} {'effect_if':<16}{eq}{gpu_text}", file=fd)
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
            if self.p_opts.stage2_restraint.lower() is not None:
                outer_space, inner_space = identation(1)
                print(f"{outer_space} {'restrain':<16}{eq}{'{'}", file=fd)
                print(
                    f"{inner_space} {'atom':<16} {eq}{q}{self.p_opts.stage2_restraint}{q}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'force_constant':<16} {eq}{self.p_opts.stage2_restraint_force}",
                    file=fd,
                )
                print(f"{outer_space} {'}'}", file=fd)
                outer_space, inner_space = identation(0)
                print(file=fd)
                print(
                    f"{inner_space} {'randomize_velocity.interval':<29} {eq}{'1.0'}",
                    file=fd,
                )
                print(f"{inner_space} {'eneseq.interval':<29} {eq}{'0.3'}", file=fd)
                print(f"{inner_space} {'trajectory.center':<29} {eq}{'[]'}", file=fd)
                print(f"{outer_space}{'}'}", file=fd)
                print(file=fd)
            # Stage 3 block
            if self.p_opts.stage3.lower() in ["yes", "on", "true"]:
                outer_space, inner_space = identation(0)
                print(f"{outer_space}{'simulate':<20}{'{'}", file=fd)
                print(
                    f"{inner_space} {'title':<16}{eq}{q}{self.p_opts.stage3_title}{q}",
                    file=fd,
                )
                gpu_text = '[["==" "-gpu" "@*.*.jlaunch_opt[-1]"] \'ensemble.method = Langevin\']'
                print(f"{inner_space} {'effect_if':<16}{eq}{gpu_text}", file=fd)
                print(f"{inner_space} {'annealing':<16}{eq}{'off'}", file=fd)
                print(
                    f"{inner_space} {'time':<16}{eq}{self.p_opts.stage3_time}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'temperature':<16}{eq}{self.p_opts.stage3_temp}",
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
                print(
                    f"{inner_space} {'barostat.tau':<11} {eq}{self.p_opts.stage3_barostat_tau}",
                    file=fd,
                )
                outer_space, inner_space = identation(1)
                print(f"{outer_space} {'}'}", file=fd)
                # Restraints block
                if self.p_opts.stage3_restraint.lower() is not None:
                    outer_space, inner_space = identation(1)
                    print(f"{outer_space} {'restrain':<16}{eq}{'{'}", file=fd)
                    print(
                        f"{inner_space} {'atom':<16} {eq}{q}{self.p_opts.stage3_restraint}{q}",
                        file=fd,
                    )
                    print(
                        f"{inner_space} {'force_constant':<16} {eq}{self.p_opts.stage3_restraint_force}",
                        file=fd,
                    )
                    print(f"{outer_space} {'}'}", file=fd)
                    outer_space, inner_space = identation(0)
                    print(file=fd)
                    print(
                        f"{inner_space} {'randomize_velocity.interval':<29} {eq}{'1.0'}",
                        file=fd,
                    )
                    print(f"{inner_space} {'eneseq.interval':<29} {eq}{'0.3'}", file=fd)
                    print(
                        f"{inner_space} {'trajectory.center':<29} {eq}{'[]'}", file=fd
                    )
                    print(f"{outer_space}{'}'}", file=fd)
                # solvate pocket block
                outer_space, inner_space = identation(0)
                print(file=fd)
                print(f"{outer_space}{'solvate_pocket':<16}{'{'}", file=fd)
                print(
                    f"{inner_space} {'should_skip':<11} {eq}{'true'}",
                    file=fd,
                )
                print(
                    f"{inner_space} {'ligand_file':<11} {eq}{'?'}",
                    file=fd,
                )
                print(f"{outer_space}{'}'}", file=fd)
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
                print(f"{inner_space} {'effect_if':<16}{eq}{gpu_text1}", file=fd)
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
                # Restraints block
                if self.p_opts.stage4_restraint.lower() is not None:
                    outer_space, inner_space = identation(1)
                    print(f"{outer_space} {'restrain':<16}{eq}{'{'}", file=fd)
                    print(
                        f"{inner_space} {'atom':<16} {eq}{q}{self.p_opts.stage4_restraint}{q}",
                        file=fd,
                    )
                    print(
                        f"{inner_space} {'force_constant':<16} {eq}{self.p_opts.stage4_restraint_force}",
                        file=fd,
                    )
                    print(f"{outer_space} {'}'}", file=fd)
                outer_space, inner_space = identation(0)
                print(file=fd)
                print(
                    f"{inner_space} {'randomize_velocity.interval':<29} {eq}{'1.0'}",
                    file=fd,
                )
                print(f"{inner_space} {'eneseq.interval':<29} {eq}{'0.3'}", file=fd)
                print(f"{inner_space} {'trajectory.center':<29} {eq}{'[]'}", file=fd)
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
                print(f"{inner_space} {'effect_if':<16}{eq}{gpu_text1}", file=fd)
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
                # Restraints block
                if str(self.p_opts.stage5_restraint).lower() != "none":
                    outer_space, inner_space = identation(1)
                    print(f"{outer_space} {'restrain':<16}{eq}{'{'}", file=fd)
                    print(
                        f"{inner_space} {'atom':<16} {eq}{q}{self.p_opts.stage5_restraint}{q}",
                        file=fd,
                    )
                    print(
                        f"{inner_space} {'force_constant':<16} {eq}{self.p_opts.stage5_restraint_force}",
                        file=fd,
                    )
                    print(f"{outer_space} {'}'}", file=fd)
                outer_space, inner_space = identation(0)
                print(file=fd)
                print(f"{inner_space} {'eneseq.interval':<29} {eq}{'0.3'}", file=fd)
                print(
                    f"{inner_space} {'trajectory.center':<29} {eq}{'solute'}", file=fd
                )
                print(f"{outer_space}{'}'}", file=fd)
                print(file=fd)
            # Additional stages block
            if self.p_opts.additional_stages != None:
                stages = int(self.p_opts.additional_stages)
                stages_time = list(str(self.p_opts.additional_stage_times).split(","))
                stages_temp = list(str(self.p_opts.additional_stage_temps).split(","))
                stages_ensemble = list(
                    str(self.p_opts.additional_stage_ensembles).split(",")
                )
                stages_method = list(
                    str(self.p_opts.additional_stage_methods).split(",")
                )
                stages_thermostat_tau = list(
                    str(self.p_opts.additional_stage_thermostat_tau).split(",")
                )
                stages_barostat_tau = list(
                    str(self.p_opts.additional_stage_barostat_tau).split(",")
                )
                stages_restraint = list(
                    str(self.p_opts.additional_stage_restraints).split(",")
                )
                stages_restraint_force = list(
                    str(self.p_opts.additional_stage_restraints_forces).split(",")
                )
                if len(stages_time) >= 2 and stages != len(stages_time):
                    raise ValueError(
                        "The number of times in additional_stage_times must be one value or equal to the number of stages"
                    )
                if len(stages_temp) >= 2 and stages != len(stages_temp):
                    raise ValueError(
                        "The number of temperatures in additional_stage_temps must be one value or equal to the number of stages"
                    )
                if len(stages_ensemble) >= 2 and stages != len(stages_ensemble):
                    raise ValueError(
                        "The number of ensembles in additional_stage_ensembles must be one value or equal to the number of stages"
                    )
                if len(stages_method) >= 2 and stages != len(stages_method):
                    raise ValueError(
                        "The number of methods in additional_stage_methods must be one value or equal to the number of stages"
                    )
                if len(stages_restraint) >= 2 and stages != len(stages_restraint):
                    raise ValueError(
                        "The number of restraints in additional_stage_restraints must be one value or equal to the number of stages"
                    )
                if len(stages_restraint_force) >= 2 and stages != len(
                    stages_restraint_force
                ):
                    raise ValueError(
                        "The number of restraint_forces in additional_stage_restraints_forces must be one value or equal to the number of stages"
                    )
                print(file=fd)
                print("Additional stages =", stages)
                print("Additional stage times =", stages_time)
                print("Additional stage temperatures =", stages_temp)
                print("Additional stage ensembles =", stages_ensemble)
                print("Additional stage methods =", stages_method)
                print("Additional stage thermostat_tau =", stages_thermostat_tau)
                print("Additional stage barostat_tau =", stages_barostat_tau)
                print("Additional stage restraints =", stages_restraint)
                print("Additional stage restraint_forces =", stages_restraint_force)

                for stage in range(1, stages + 1):
                    outer_space, inner_space = identation(0)
                    print(f"{outer_space}{'simulate':<20}{'{'}", file=fd)
                    print(
                        f"{inner_space} {'title':<16}{eq}{q}Additional stage = {stage}, {q}",
                        file=fd,
                    )
                    gpu_text1 = '[["@*.*.annealing"] \'annealing = off temperature = "@*.*.temperature[0][0]"\''
                    gpu_text2 = '["==" "-gpu" "@*.*.jlaunch_opt[-1]"] \'ensemble.method = Langevin\']'
                    print(f"{inner_space} {'effect_if':<16}{eq}{gpu_text1}", file=fd)
                    print(f"{inner_space} {' ':<19}{gpu_text2}", file=fd)
                    if (
                        self.p_opts.additional_stage_times != None
                        and len(stages_time) == 1
                    ):
                        print(
                            f"{inner_space} {'time':<16}{eq}{self.p_opts.additional_stage_times}",
                            file=fd,
                        )
                    elif (
                        self.p_opts.additional_stage_times != None
                        and len(stages_time) == stages
                    ):
                        print(
                            f"{inner_space} {'time':<16}{eq}{stages_time[stage-1]}",
                            file=fd,
                        )
                    if (
                        self.p_opts.additional_stage_temps != None
                        and len(stages_time) == 1
                    ):
                        print(
                            f"{inner_space} {'temperature':<16}{eq}{self.p_opts.additional_stage_temps}",
                            file=fd,
                        )
                    elif (
                        self.p_opts.additional_stage_temps != None
                        and len(stages_temp) == stages
                    ):
                        print(
                            f"{inner_space} {'temperature':<16}{eq}{stages_temp[stage-1]}",
                            file=fd,
                        )
                    outer_space, inner_space = identation(1)

                    print(f"{outer_space} {'ensemble':<16}{eq}{'{'}", file=fd)
                    # Ensemble block
                    if (
                        self.p_opts.additional_stage_ensembles != None
                        and len(stages_ensemble) == 1
                    ):
                        print(
                            f"{inner_space} {'class':<11}{eq}{self.p_opts.additional_stage_ensembles}",
                            file=fd,
                        )
                    elif (
                        self.p_opts.additional_stage_ensembles != None
                        and len(stages_ensemble) == stages
                    ):
                        print(
                            f"{inner_space} {'class':<11}{eq}{stages_ensemble[stage-1]}",
                            file=fd,
                        )
                    # Method block
                    if (
                        self.p_opts.additional_stage_methods != None
                        and len(stages_method) == 1
                    ):
                        print(
                            f"{inner_space} {'method':<11}{eq}{self.p_opts.additional_stage_methods}",
                            file=fd,
                        )
                    elif (
                        self.p_opts.additional_stage_methods != None
                        and len(stages_method) == stages
                    ):
                        print(
                            f"{inner_space} {'method':<11}{eq}{stages_method[stage-1]}",
                            file=fd,
                        )
                    # Thermostat block
                    if (
                        self.p_opts.additional_stage_thermostat_tau != None
                        and len(stages_thermostat_tau) == 1
                    ):
                        print(
                            f"{inner_space} {'thermostat.tau':<15}{eq}{self.p_opts.additional_stage_thermostat_tau}",
                            file=fd,
                        )
                    elif (
                        self.p_opts.additional_stage_thermostat_tau != None
                        and len(stages_thermostat_tau) == stages
                    ):
                        print(
                            f"{inner_space} {'thermostat.tau':<15}{eq}{stages_thermostat_tau[stage-1]}",
                            file=fd,
                        )
                    # Barostat block
                    if (
                        self.p_opts.additional_stage_barostat_tau != None
                        and len(stages_barostat_tau) == 1
                    ):
                        print(
                            f"{inner_space} {'barostat.tau':<15}{eq}{self.p_opts.additional_stage_barostat_tau}",
                            file=fd,
                        )
                    elif (
                        self.p_opts.additional_stage_barostat_tau != None
                        and len(stages_barostat_tau) == stages
                    ):
                        print(
                            f"{inner_space} {'barostat.tau':<15}{eq}{stages_barostat_tau[stage-1]}",
                            file=fd,
                        )

                    outer_space, inner_space = identation(1)
                    print(f"{outer_space} {'}'}", file=fd)
                    # Restraints block
                    if (
                        str(self.p_opts.additional_stage_restraints_forces).lower()
                        != "none"
                        and float(stages_restraint_force[stage - 1]) != 0.0
                    ):
                        outer_space, inner_space = identation(1)
                        print(f"{outer_space} {'restrain':<16}{eq}{'{'}", file=fd)
                        # Atom block
                        if (
                            self.p_opts.additional_stage_restraints != None
                            and len(stages_restraint) == 1
                        ):
                            print(
                                f"{inner_space} {'atom':<16} {eq}{q}{self.p_opts.additional_stage_restraints}{q}",
                                file=fd,
                            )
                        elif (
                            self.p_opts.additional_stage_restraints != None
                            and len(stages_restraint) == stages
                        ):
                            print(
                                f"{inner_space} {'atom':<16} {eq}{q}{stages_restraint[stage-1]}{q}",
                                file=fd,
                            )
                        # Forces block
                        if (
                            self.p_opts.additional_stage_restraints_forces != None
                            and len(stages_restraint_force) == 1
                        ):
                            print(
                                f"{inner_space} {'force_constant':<16} {eq}{self.p_opts.additional_stage_restraints_forces}",
                                file=fd,
                            )
                        elif (
                            self.p_opts.additional_stage_restraints_forces != None
                            and len(stages_restraint_force) == stages
                        ):
                            print(
                                f"{inner_space} {'force_constant':<16} {eq}{stages_restraint_force[stage-1]}",
                                file=fd,
                            )
                        print(f"{outer_space} {'}'}", file=fd)
                    outer_space, inner_space = identation(0)
                    print(file=fd)
                    print(f"{inner_space} {'eneseq.interval':<29} {eq}{'0.3'}", file=fd)
                    print(
                        f"{inner_space} {'trajectory.center':<29} {eq}{'solute'}",
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
            if str(self.p_opts.production_restraint).lower() != "none":
                outer_space, inner_space = identation(0)
                print(f"{outer_space}{'restrain':<20}{eq}{'{'}", file=fd)
                print(
                    f"{inner_space} {'atom':<15} {eq}{q}{self.p_opts.production_restraint}{q}",
                    file=fd,
                )
                if self.p_opts.production_restraint_force is None:
                    raise ValueError("You must set the force for the restraint")
                print(
                    f"{inner_space} {'force_constant':<15} {eq}{self.p_opts.production_restraint_force}",
                    file=fd,
                )
                print(f"{outer_space}{'}'}", file=fd)
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
        executable = os.path.join(self.builder_opts.schrod_path, "utilities/multisim")
        path_preparation_sh = str(self.basename + "_md.sh")
        input_msj = str(self.basename + "_md.msj")
        input_cms = self.basename + "_preparation" + "-out.cms"
        input_cfg = self.basename + "_md.cfg"
        gpu_opts = 'stage[1].set_family.md.jlaunch_opt=["-gpu"]'
        args1 = f"-HOST localhost -JOBNAME {self.basename}_md -maxjob 1 -cpu 1 -m {input_msj} -c {input_cfg} {input_cms}"
        args2 = f"-mode umbrella -set '{gpu_opts}' -o {self.basename}_md-out.cms"
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

    defaults = {"schrod_path": "$SCHRODINGER"}

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
    schrod_path = vars(args)["schrod_path"]
    file_path = os.path.abspath(vars(args)["file"])

    protocol_opts = {}
    protocol_opts.update(dict(config.items("protocol")))

    return (
        Args(**vars(args)),
        BuilderOptions(build_opts, file, filename, schrod_path, file_path),
        file_path,
        ProtocolOptions(protocol_opts),
    )


@dataclass
class maefile:
    file: str
    charge: int
    schrod_path: str


class ReadMaefile:
    def __init__(self, file: str, schrod_path) -> None:
        self.file = os.path.relpath(file)
        self.charge: int = 0
        self.schrod_path = schrod_path

    def get_charge(self):
        executable = os.path.join(self.schrod_path, "run")
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
        executable = os.path.join(self.schrod_path, "run")
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
        raise ValueError(f"Folder '{folder_name}' exists, remove it before to continue")
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
        self.schrod_path = options.schrod_path
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
                    self.options.ion = self.options.positive_ion
                elif int(self.charge) > 0:
                    self.options.ion = self.options.negative_ion
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
        executable = os.path.join(self.schrod_path, "utilities/multisim")
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
    system = ReadMaefile(file, opts.schrod_path)
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
