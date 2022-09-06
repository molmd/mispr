# coding: utf-8


# Implements a core class MaestroRunner for assigning force field
# parameters on a molecule using Maestro software.

import itertools
import logging
import os
import re
import subprocess

from configparser import ConfigParser

import numpy as np
import pandas as pd

from fireworks.fw_config import CONFIG_FILE_DIR
from pymatgen.core.structure import Molecule

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Aug 2022"
__version__ = "0.0.2"

logger = logging.getLogger(__name__)

NONBONDED_HEADER = "OPLSAA FORCE FIELD TYPE ASSIGNED\n"
BONDS_HEADER = " Stretch            k            r0    quality         bt        comment\n"
ANGLES_HEADER = " Bending                      k       theta0    quality   at  comment\n"
DIHEDRALS_HEADER = " proper Torsion                     V1      V2      V3      V4    quality  tt  comment\n"
IMPROPER_HEADER = " improper Torsion                   V2    quality  comment\n"


class MaestroRunner:
    """
    Wrapper for assigning OPLS_2005 force field parameters for a
    molecule using Maestro software of Schrondinger. The OPLS_2005
    parameters are described in:

    Banks, J.L.; Beard, H.S.; Cao, Y.; Cho, A.E.; Damm, W.; Farid, R.;
    Felts, A.K.; Halgren, T.A.; Mainz, D.T.; Maple, J.R.; Murphy, R.;
    Philipp, D.M.; Repasky, M.P.; Zhang, L.Y.; Berne, B.J.; Friesner,
    R.A.; Gallicchio, E.; Levy. R.M. Integrated Modeling Program,
    Applied Chemical Theory (IMPACT). J. Comp. Chem. 2005, 26, 1752.

    The OPLS_2005 parameters are located in (you need to replace the
    vversion in the path with the correct Schrodinger version
    installed):

    $SCHRODINGER/mmshare-vversion/data/f14/

    To use this class, you need to do the following:

    1. Download the free version of Maestro
    (https://www.schrodinger.com/freemaestro)
    2. Install Maestro and set the environment variable $SCHRODINGER
    (e.g. 'export SCHRODINGER=/opt/schrodinger/suites2022-2'.
    Check https://www.schrodinger.com/kb/1842 for more details.

    """
    # TODO: process improper torsions

    def __init__(self, name, input_file, mae_cmd=None, ffld_cmd=None, working_dir=None):
        """
        Args:
            name (str): name of the molecule to use when saving the mae
                and log files generated by Maestro
            input_file (str): molecule/structure file to use for
                converting to a mae file. For a full list of supported
                formats, check https://www.schrodinger.com/kb/1278
            mae_cmd (str): full command to use for converting the
                structure format to mae file; if not specified, will
                try to parse the command from the MISPR configuration
                file
            ffld_cmd (str): full command to use for generating the force
                field parameters from the mae file; if not specified,
                will try to parse the command from the MISPR
                configuration file
            working_dir (str): working directory where the input_file
                is located and Maestro files to be generated
        """
        self.name = name
        self.input_file = input_file
        self.mae_cmd = mae_cmd
        self.ffld_cmd = ffld_cmd
        self.working_dir = working_dir or os.getcwd()

        if not os.path.isabs(self.input_file):
            self.input_file = f"{self.working_dir}/{self.input_file}"

    def get_opls_params(self, cleanup=True):
        """
        Wrapper function that converts input structure to mae format,
        uses it to generate a log file with the force field parameters
        of the molecule, parses the log file to generate a
        dictionary of force parameters in a format compatible with the
        lammps workflow in MISPR, and optionally cleans up the working
        directory by removing the intermediate files generated in the
        process.

        Args:
            cleanup (bool): whether to clean up the working directory

        Returns:
            dict: dictionary with force field parameters
        """
        if not self.input_file.endswith(".mae"):
            self.convert_mol_to_mae()
        self.get_ff_log()
        cosolvent_param = self.convert_ff_log_to_json()
        if cleanup:
            self.cleanup()
        return cosolvent_param

    @staticmethod
    def _get_num_lines(file_):
        """
        Returns the number of lines in a file.

        Args:
            file_ (str): file to read

        Returns:
            int
        """
        return sum(1 for line in open(file_))

    @staticmethod
    def _get_footer(file_, footer_str):
        """
        Returns the number of lines remaining in the file after a footer
        string.

        Args:
            file_ (str): file to read
            footer_str (str): line to read

        Returns:
            int
        """
        with open(file_) as f:
            g = itertools.dropwhile(lambda x: x != footer_str, f)
            footer_len = len([i for i, _ in enumerate(g)])
        return footer_len

    @staticmethod
    def _get_atom_type(atom_string):
        """
        Finds the weight of an atom from a string representation.

        Args:
            atom_string (str): atom representation, e.g. Cl1, H10, etc.
        Returns:
            float
        """
        match = re.match(r"([a-z]+)([0-9]+)", atom_string, re.I)
        if match:
            items = match.groups()
            mol = Molecule.from_str(f"[{items[0]}]", fmt="smi")
            return mol[0].species.weight
        else:
            raise ValueError(f"Atom type cannot be determined from {atom_string}")

    def convert_mol_to_mae(self):
        """
        Converts a structure input file to a mae file.
        """
        if not self.mae_cmd:
            cfg = ConfigParser()
            cfg.read(CONFIG_FILE_DIR + "/config.ini")
            self.mae_cmd = cfg["MaestroCalc"]["mae_cmd"]
        cmd = self.mae_cmd.replace("$input_file$", self.input_file).replace(
            "$output_file$", f"{self.working_dir}/{self.name}.mae"
        )
        logger.info("Running command: {}".format(cmd))
        return_code = subprocess.call(cmd, env=os.environ, shell=True)
        logger.info("Finished running with return code: {}".format(return_code))

    def get_ff_log(self):
        """
        Generates a log file with force field parameters using a mae
        file.
        """
        if not self.ffld_cmd:
            cfg = ConfigParser()
            cfg.read(CONFIG_FILE_DIR + "/config.ini")
            self.ffld_cmd = cfg["MaestroCalc"]["ffld_cmd"]
        cmd = self.ffld_cmd.replace(
            "$input_file$", f"{self.working_dir}/{self.name}.mae"
        ).replace("$output_file$", f"{self.working_dir}/{self.name}.log")
        logger.info("Running command: {}".format(cmd))
        return_code = subprocess.call(cmd, env=os.environ, shell=True)
        logger.info("Finished running with return code: {}".format(return_code))

    def convert_ff_log_to_json(self):
        """
        Converts the force field parameters in the log file to a
        dictionary format that is compatible with MISPR.
        """
        log_file = f"{self.working_dir}/{self.name}.log"
        num_lines = self._get_num_lines(log_file)

        # non bonded
        row_skips = num_lines - self._get_footer(log_file, NONBONDED_HEADER) + 3
        footer_skips = (
            self._get_footer(
                log_file,
                BONDS_HEADER,
            )
            + 3
        )
        nonbonded_df = pd.read_csv(
            log_file,
            skiprows=row_skips,
            skipfooter=footer_skips,
            delimiter=r"\s\s+",
            engine="python",
        ).reset_index()

        df_cols = []
        for i in nonbonded_df.columns[:-1]:
            if nonbonded_df[i].dtype == "object":
                df_cols.append(nonbonded_df[i].str.split(expand=True))
            else:
                df_cols.append(nonbonded_df[i])
        nonbonded_df = pd.concat(df_cols, axis=1, ignore_index=True)
        nonbonded_df = nonbonded_df[[0, 3, 4, 5, 6]]
        nonbonded_df.columns = ["Atom", "Type", "Charge", "Sigma", "Epsilon"]
        labels = nonbonded_df["Type"].to_list()
        charges = nonbonded_df["Charge"].to_list()
        nonbonded_df_unique = nonbonded_df.drop_duplicates(subset=["Type"]).reset_index(
            drop=True
        )
        atoms_dict = {
            i: self._get_atom_type(j)
            for i, j in zip(nonbonded_df_unique["Type"], nonbonded_df_unique["Atom"])
        }
        non_bonded = [
            [i, j] for i, j in zip(nonbonded_df_unique["Epsilon"], nonbonded_df_unique["Sigma"])
        ]

        # bonds
        row_skips = num_lines - self._get_footer(
            log_file,
            BONDS_HEADER,
        )
        footer_skips = self._get_footer(
            log_file,
            ANGLES_HEADER,
        )
        bonds_df = pd.read_csv(
            log_file,
            skiprows=row_skips,
            skipfooter=footer_skips,
            delimiter=r"\s+",
            engine="python",
        ).reset_index()

        bond_data = []
        if not bonds_df.empty:
            bonds_df = bonds_df[["quality", "comment", "level_2", "level_3"]]
            df = bonds_df[["quality", "level_2", "level_3"]]
            df.columns = ["comment", "level_2", "level_3"]
            del bonds_df["quality"]
            bonds_df = bonds_df.append(df)
            bonds_df.columns = ["Bond", "K", "r"]
            bonds_df = bonds_df.join(
                bonds_df["Bond"]
                .str.split("-", expand=True)
                .rename(columns={0: "t1", 1: "t2"})
            )
            bonds_df["Bond"] = bonds_df.apply(
                lambda x: "-".join(
                    (
                        [x["t1"], x["t2"]]
                        if not (
                            x["t1"] > x["t2"] or (x["t1"] == x["t2"] and x["t1"] > x["t2"])
                        )
                        else [x["t2"], x["t1"]]
                    )
                ),
                axis=1,
            )
            bonds_df.drop(["t1", "t2"], axis=1, inplace=True)
            bonds_df = bonds_df.drop_duplicates(subset=["Bond"]).reset_index(drop=True)

            for index, row in bonds_df.iterrows():
                bond_data.append(
                    {"coeffs": [row[1], row[2]], "types": [tuple(row[0].split("-"))]}
                )

        # angles
        row_skips = num_lines - self._get_footer(
            log_file,
            ANGLES_HEADER,
        )
        footer_skips = self._get_footer(
            log_file,
            DIHEDRALS_HEADER,
        )
        angles_df = pd.read_csv(
            log_file,
            skiprows=row_skips,
            skipfooter=footer_skips,
            delimiter=r"\s+",
            engine="python",
        ).reset_index()

        angle_data = []
        if not angles_df.empty:
            angles_df = angles_df[["quality", "comment", "level_3", "Bending"]]
            df = angles_df[["quality", "level_3", "Bending"]]
            df.columns = ["comment", "level_3", "Bending"]
            del angles_df["quality"]
            angles_df = angles_df.append(df)
            angles_df.columns = ["Angle", "K", "Theta"]
            angles_df = angles_df.join(
                angles_df["Angle"]
                .str.split("-", expand=True)
                .rename(columns={0: "t1", 1: "t2", 2: "t3"})
            )
            angles_df["Angle"] = angles_df.apply(
                lambda x: "-".join(
                    (
                        [x["t1"], x["t2"], x["t3"]]
                        if not (
                            x["t1"] > x["t3"]
                            or (x["t1"] == x["t3"] and x["t1"] > x["t3"])
                        )
                        else [x["t3"], x["t2"], x["t1"]]
                    )
                ),
                axis=1,
            )
            angles_df.drop(["t1", "t2", "t3"], axis=1, inplace=True)
            angles_df = angles_df.drop_duplicates(subset=["Angle"]).reset_index(
                drop=True
            )

            for index, row in angles_df.iterrows():
                angle_data.append(
                    {"coeffs": [row[1], row[2]], "types": [tuple(row[0].split("-"))]}
                )

        # dihedrals
        row_skips = num_lines - self._get_footer(
            log_file,
            DIHEDRALS_HEADER,
        )
        footer_skips = self._get_footer(
            log_file,
            IMPROPER_HEADER,
        )
        dihedrals_df = pd.read_csv(
            log_file,
            skiprows=row_skips,
            skipfooter=footer_skips,
            delimiter=r"\s+",
            engine="python",
        ).reset_index()

        dihedral_data = []
        if not dihedrals_df.empty:
            dihedrals_df = dihedrals_df[
                ["quality", "comment", "proper", "Torsion", "V1", "V2"]
            ]
            df = dihedrals_df[["quality", "proper", "Torsion", "V1", "V2"]]
            df.columns = ["comment", "proper", "Torsion", "V1", "V2"]
            del dihedrals_df["quality"]
            dihedrals_df = dihedrals_df.append(df)
            dihedrals_df.columns = ["Dihedral", "V1", "V2", "V3", "V4"]
            dihedrals_df = dihedrals_df.join(
                dihedrals_df["Dihedral"]
                .str.split("-", expand=True)
                .rename(columns={0: "t1", 1: "t2", 2: "t3", 3: "t4"})
            )
            dihedrals_df["Dihedral"] = dihedrals_df.apply(
                lambda x: "-".join(
                    (
                        [x["t1"], x["t2"], x["t3"], x["t4"]]
                        if not (
                            x["t1"] > x["t4"]
                            or (x["t1"] == x["t4"] and x["t1"] > x["t4"])
                        )
                        else [x["t4"], x["t3"], x["t2"], x["t1"]]
                    )
                ),
                axis=1,
            )
            dihedrals_df.drop(["t1", "t2", "t3", "t4"], axis=1, inplace=True)
            dihedrals_df[["V1", "V2", "V3", "V4"]] = (
                dihedrals_df[["V1", "V2", "V3", "V4"]] / 2
            )
            dihedrals_df = dihedrals_df.drop_duplicates(
                subset=["Dihedral"]
            ).reset_index(drop=True)
            dihedrals_df = dihedrals_df[
                ~dihedrals_df["Dihedral"].str.contains("?", regex=False)
            ]

            for index, row in dihedrals_df.iterrows():
                dihedral_data.append(
                    {
                        "coeffs": [
                            row[i] for i in range(len(dihedrals_df.columns) - 1)
                        ],
                        "types": [tuple(row[0].split("-"))],
                    }
                )
            mapping = {0: [1, 0], 1: [2, 180], 2: [3, 0]}
            for i in range(len(dihedral_data)):
                list_ = dihedral_data[i]["coeffs"][1:]
                data = [i for i, e in enumerate(list_) if e != 0]
                if not data:
                    data = [0]
                dihedrals = [[list_[i]] + mapping[i] for i in data]
                dihedrals = [item for sublist in dihedrals for item in sublist]
                dihedral_data[i]["coeffs"] = ["fourier", len(data)] + dihedrals

        # impropers
        row_skips = num_lines - self._get_footer(log_file, IMPROPER_HEADER)
        improper_df = pd.read_csv(
            log_file, skiprows=row_skips, delimiter=r"\s+", engine="python"
        ).reset_index()

        improper_data = []
        if not improper_df.empty:
            improper_df = improper_df[
                ["level_0", "level_1", "level_2", "improper", "Torsion"]
            ].copy()
            improper_df.replace(
                nonbonded_df[['Atom', 'Type']].set_index('Atom').squeeze().to_dict(),
                inplace=True)
            for index, row in improper_df.iterrows():
                improper_data.append(
                    {
                        "coeffs": ["cvff", row[-1]/2, -1, 2],
                        "types": [row[i] for i in range(len(improper_df.columns) - 1)],
                    }
                )

        logger.info("Finished converting OPLS parameters to MISPR format")

        ff_params = {
            "Labels": labels,
            "Masses": atoms_dict,
            "Nonbond": non_bonded,
            "Bonds": bond_data,
            "Angles": angle_data,
            "Dihedrals": dihedral_data,
            "Impropers": improper_data,
            "Improper Topologies": None,
            "Charges": charges,
        }
        return ff_params

    def cleanup(self):
        """
        Deletes the log and mae files created by Maestro
        """
        for ext in ["log", "mae"]:
            os.remove(f"{self.working_dir}/{self.name}.{ext}")
