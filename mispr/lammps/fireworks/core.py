# coding: utf-8


# Defines core fireworks for ambertools and lammps.

import os
import logging

from fireworks import Firework

from mdproptools.structural.rdf_cn import calc_atomic_rdf, calc_molecular_rdf

from mispr.gaussian.workflows.base.core import _process_mol_check
from mispr.gaussian.firetasks.geo_transformation import ProcessMoleculeInput
from mispr.lammps.firetasks.run import (
    RunTleap,
    RunLammpsDirect,
    RunParmchk,
    RunAntechamber,
    RunMaestro,
)
from mispr.gaussian.utilities.metadata import get_chem_schema
from mispr.lammps.firetasks.write_inputs import (
    LabelFFDict,
    WriteControlFile,
    WriteTleapScript,
    LabelFFDictFromDB,
)
from mispr.lammps.firetasks.parse_outputs import GetRDF, CalcDiff, ProcessPrmtop

__author__ = "Matthew Bliss"
__maintainer__ = "Matthew Bliss"
__email__ = "matthew.bliss@stonybrook.edu"
__status__ = "Development"
__date__ = "Apr 2020"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)

FIREWORK_KWARGS = Firework.__init__.__code__.co_varnames


def ambertools_tasks(**kwargs):
    """
    Define a list of common tasks for generating GAFF parameters for a 
    molecule. This is a helper function for the GetFFDictFW Firework.

    Args:
        kwargs: other kwargs that are passed to:

            1. ``mispr.lammps.firetasks.run.RunAntechamber``
            2. ``mispr.lammps.firetasks.run.RunParmchk``
            3. ``mispr.lammps.firetasks.write_inputs.WriteTleapScript``
            4. ``mispr.lammps.firetasks.run.RunTleap``
            5. ``mispr.lammps.firetasks.parse_outputs.ProcessPrmtop``
    
    Returns:
        List of Firetasks.
    """
    common_t = [
        RunAntechamber(
            **{
                i: j
                for i, j in kwargs.items()
                if i in RunAntechamber.required_params + RunAntechamber.optional_params
            }
        ),
        RunParmchk(
            **{
                i: j
                for i, j in kwargs.items()
                if i in RunParmchk.required_params + RunParmchk.optional_params
            }
        ),
        WriteTleapScript(
            **{
                i: j
                for i, j in kwargs.items()
                if i
                in WriteTleapScript.required_params + WriteTleapScript.optional_params
            }
        ),
        RunTleap(
            **{
                i: j
                for i, j in kwargs.items()
                if i in RunTleap.required_params + RunTleap.optional_params
            }
        ),
        ProcessPrmtop(
            **{
                i: j
                for i, j in kwargs.items()
                if i in ProcessPrmtop.required_params + ProcessPrmtop.optional_params
            }
        ),
    ]
    return common_t


class GetFFDictFW(Firework):
    """
    Generate the ff parameters for a molecule.
    """

    def __init__(
        self,
        mol,
        mol_operation_type,
        data,
        operation_type="get_from_esp",
        label="",
        name="get_ff_dict",
        parents=None,
        working_dir=None,
        db=None,
        save_ff_to_db=False,
        save_ff_to_file=True,
        ff_filename="ff.json",
        tag="unknown",
        **kwargs
    ):
        """
        Args:
            mol (Molecule, GaussianOutput, str, dict): Source of the 
                molecule to be processed. Should match the 
                ``mol_operation_type``.
            mol_operation_type (str): The type of molecule operation. 
                See process_mol defined in 
                ``mispr/gaussian/utilities/mol.py`` for supported 
                operations.
            data (str, dict): Data to be processed, e.g. path to the 
                esp file if operation_type is 'get_from_esp'; path to 
                the prmtop file if operation_type is 'get_from_prmtop';
                etc.
            operation_type (str, optional): The operation to perform on 
                the data to read or generate the force field parameters. 
                Defaults to ``get_from_esp``. Supported commands:

                1. 'get_from_esp': If the input is an ESP file.
                2. 'get_from_prmtop': If the input is a prmtop file.
                3. 'get_from_opls': If the input is a molecule file to 
                    be used to generate opls ff parameters.
                4. 'get_from_dict': If the input is a dictionary of
                    force field parameters. e.g.

                    .. code-block:: python

                        {"Molecule": pymatgen_molecule,
                         "Labels": atom_labels,
                         "Masses": atomtype_masses,
                         "Nonbond": nonbond_params,
                         "Bonds": bond_params,
                         "Angles": angle_params,
                         "Dihedrals": dihedral_params,
                         "Impropers": improper_params,
                         "Improper Topologies": improper_topologies,
                         "Charges": charges}

                5. 'get_from_file': If the input is a json file of the
                    force field parameters.
                6. 'get_from_db': If the input is a filter for the
                    database to search for the force field parameters.
                
            label (str, optional): Label for the molecule. This should 
                be unique for each different molecular species in the 
                system. Defaults to an empty string. In this case, the 
                label will be obtained based on the molecular formula.
            name (str, optional): Name of the Firework. Defaults to 
                "get_ff_dict".
            parents (Firework or [Firework], optional): List of parent 
                Fireworks that this Firework depends on. Defaults to 
                ``None``.
            working_dir (str, optional): Directory to run the Firework 
                in. Defaults to the current working directory.
            db (str or dict, optional): Database credentials. Could be 
                provided as the path to the db.json file or in the form 
                of a dict. If none is provided, attempts to read it 
                from the configuration files to save the Firework to.
            save_ff_to_db (bool, optional): Whether to save the force 
                field to the database. Defaults to False.
            save_ff_to_file (bool, optional): Whether to save the force 
                field to a file. Defaults to True.
            ff_filename (str, optional): Filename to save the force 
                field to. Defaults to "ff.json".
            tag (str): Tag for the Firework. The provided tag will be 
                stored in the db documents for easy retrieval. Defaults 
                to "unknown".
            kwargs: Other kwargs that are passed to:

                1. Firework.__init__
                2. ``mispr.gaussian.firetasks.geo_transformation.ProcessMoleculeInput``
                3. ``mispr.lammps.firetasks.run.RunAntechamber``
                4. ``mispr.lammps.firetasks.run.RunParmchk``
                5. ``mispr.lammps.firetasks.write_inputs.WriteTleapScript``
                6. ``mispr.lammps.firetasks.run.RunTleap``
                7. ``mispr.lammps.firetasks.parse_outputs.ProcessPrmtop``
                8. ``mispr.lammps.firetasks.run.RunMaestro``
                9. ``mispr.lammps.firetasks.write_inputs.LabelFFDict``
                10. ``mispr.lammps.firetasks.write_inputs.LabelFFDictFromDB``

        """
        tasks = []
        working_dir = working_dir or os.getcwd()

        # TODO: add logic for ensuring that directory name is legal or at \
        #  least usable
        os.makedirs(working_dir, exist_ok=True)

        tasks.append(
            ProcessMoleculeInput(
                mol=mol,
                operation_type=mol_operation_type,
                db=db,
                working_dir=working_dir,
                **{
                    i: j
                    for i, j in kwargs.items()
                    if i
                    in ProcessMoleculeInput.required_params
                    + ProcessMoleculeInput.optional_params
                }
            )
        )

        if operation_type == "get_from_esp":
            if not isinstance(data, str):
                raise TypeError(
                    "data must be a str of the path to the esp "
                    'file for operation_type="get_from_esp"'
                )
            tasks += ambertools_tasks(
                db=db,
                working_dir=working_dir,
                input_filename_a=data,
                unique_molecule_name=label,
                save_ff_to_db=save_ff_to_db,
                save_ff_to_file=save_ff_to_file,
                ff_filename=ff_filename,
                **{
                    i: j
                    for i, j in kwargs.items()
                    if i
                    in RunAntechamber.required_params
                    + RunAntechamber.optional_params
                    + RunParmchk.required_params
                    + RunParmchk.optional_params
                    + WriteTleapScript.required_params
                    + WriteTleapScript.optional_params
                    + RunTleap.required_params
                    + RunTleap.optional_params
                    + ProcessPrmtop.required_params
                    + ProcessPrmtop.optional_params
                }
            )

        elif operation_type == "get_from_prmtop":
            if not isinstance(data, str):
                raise TypeError(
                    "data must be a str of the path to the prmtop "
                    'file for operation_type="get_from_prmtop"'
                )
            tasks.append(
                ProcessPrmtop(
                    working_dir=working_dir,
                    db=db,
                    prmtop_path=data,
                    unique_molecule_name=label,
                    save_ff_to_db=save_ff_to_db,
                    save_ff_to_file=save_ff_to_file,
                    ff_filename=ff_filename,
                    **{
                        i: j
                        for i, j in kwargs.items()
                        if i
                        in ProcessPrmtop.required_params + ProcessPrmtop.optional_params
                    }
                )
            )

        elif operation_type == "get_from_opls":
            if not isinstance(data, str):
                raise TypeError(
                    "data must be a str of the path to the molecule "
                    "file for operation_type='get_from_opls'. "
                    "Check https://www.schrodinger.com/kb/1278 for "
                    "supported structure file formats"
                )
            tasks.append(
                RunMaestro(
                    input_file=data,
                    label=label,
                    db=db,
                    working_dir=working_dir,
                    save_ff_to_db=save_ff_to_db,
                    save_ff_to_file=save_ff_to_file,
                    ff_filename=ff_filename,
                    **{
                        i: j
                        for i, j in kwargs.items()
                        if i in RunMaestro.required_params + RunMaestro.optional_params
                    }
                )
            )

        elif operation_type == "get_from_dict":
            if not isinstance(data, dict):
                raise TypeError(
                    "data must be a dict for " 'operation_type="get_from_dict"'
                )
            tasks.append(
                LabelFFDict(
                    unlabeled_dict=data,
                    working_dir=working_dir,
                    label=label,
                    **{
                        i: j
                        for i, j in kwargs.items()
                        if i
                        in LabelFFDict.required_params + LabelFFDict.optional_params
                    }
                )
            )

        elif operation_type == "get_from_file":
            if not isinstance(data, str):
                raise TypeError(
                    "data must be a str of the path to the ff "
                    'json file for operation_type="get_from_file"'
                )
            tasks.append(
                LabelFFDict(
                    ff_file=data,
                    working_dir=working_dir,
                    label=label,
                    **{
                        i: j
                        for i, j in kwargs.items()
                        if i
                        in LabelFFDict.required_params + LabelFFDict.optional_params
                    }
                )
            )

        elif operation_type == "get_from_db":
            tasks.append(
                LabelFFDictFromDB(
                    filter=data,
                    db=db,
                    working_dir=working_dir,
                    label=label,
                    **{
                        i: j
                        for i, j in kwargs.items()
                        if i
                        in LabelFFDictFromDB.required_params
                        + LabelFFDictFromDB.optional_params
                    }
                )
            )

        spec = kwargs.pop("spec", {})
        spec.update({"tag": tag, "_launch_dir": working_dir})
        super(GetFFDictFW, self).__init__(
            tasks,
            parents=parents,
            name=name,
            spec=spec,
            **{i: j for i, j in kwargs.items() if i in FIREWORK_KWARGS}
        )


class RunLammpsFW(Firework):
    """
    Run a lammps simulation.
    """
    def __init__(
        self,
        control_file=None,
        db=None,
        name="run_lammps",
        parents=None,
        working_dir=None,
        save_run_to_db=True,
        save_run_to_file=False,
        tag="unknown",
        **kwargs
    ):
        """
        Args:
            control_file (str, optional): Path to the control file. If 
                not provided, the control file will be generated. 
                Defaults to None.
            db (str or dict, optional): Database credentials. Could be 
                provided as the path to the db.json file or in the form 
                of a dict. If none is provided, attempts to read it 
                from the configuration files to save the Firework to.
            name (str, optional): Name of the Firework. Defaults to 
                "run_lammps".
            parents (Firework or [Firework], optional): List of parent 
                Fireworks that this Firework depends on. Defaults to 
                ``None``.
            working_dir (str, optional): Directory to run the Firework 
                in. Defaults to the current working directory.
            save_run_to_db (bool, optional): Whether to save the run to 
                the database. Defaults to True.
            save_run_to_file (bool, optional): Whether to save the run 
                to a file. Defaults to False.
            tag (str): Tag for the Firework. The provided tag will be 
                stored in the db documents for easy retrieval. Defaults 
                to "unknown".
            kwargs: Other kwargs that are passed to:

                1. Firework.__init__
                2. ``mispr.lammps.firetasks.write_inputs.WriteControlFile``
                3. ``mispr.lammps.firetasks.run.RunLammpsDirect``
        """
        tasks = []
        working_dir = working_dir or os.getcwd()

        if not control_file:
            tasks.append(
                WriteControlFile(
                    working_dir=working_dir,
                    db=db,
                    save_runs_to_db=save_run_to_db,
                    save_runs_to_file=save_run_to_file,
                    **{
                        i: j
                        for i, j in kwargs.items()
                        if i
                        in WriteControlFile.required_params
                        + WriteControlFile.optional_params
                    }
                )
            )

        tasks.append(
            RunLammpsDirect(
                working_dir=working_dir,
                **{
                    i: j
                    for i, j in kwargs.items()
                    if i
                    in RunLammpsDirect.required_params + RunLammpsDirect.optional_params
                }
            )
        )

        spec = kwargs.pop("spec", {})
        spec.update({"tag": tag, "_launch_dir": working_dir})

        super(RunLammpsFW, self).__init__(
            tasks,
            parents=parents,
            name=name,
            spec=spec,
            **{i: j for i, j in kwargs.items() if i in FIREWORK_KWARGS}
        )


class RunAnalysisFW(Firework):
    """
    Run an analysis on a LAMMPS trajectory.
    """
    def __init__(
        self,
        md_property,
        name="run_analysis",
        parents=None,
        working_dir=None,
        tag="unknown",
        **kwargs
    ):
        """
        Args:
            md_property (str): The property to calculate. Supported 
                properties:
                
                1. 'diffusion': Calculate the diffusion coefficient.
                2. 'rdf': Calculate the radial distribution function.
            name (str, optional): Name of the Firework. Defaults to 
                "run_analysis".
            parents (Firework or [Firework], optional): List of parent 
                Fireworks that this Firework depends on. Defaults to 
                ``None``.
            working_dir (str, optional): Directory to run the Firework 
                in. Defaults to the current working directory.
            tag (str): Tag for the Firework. The provided tag will be 
                stored in the db documents for easy retrieval. Defaults 
                to "unknown".
            kwargs: Other kwargs that are passed to:

                1. Firework.__init__
                2. ``mispr.lammps.firetasks.parse_outputs.CalcDiff``
                3. ``mispr.lammps.firetasks.parse_outputs.GetRDF``
                4. ``mispr.lammps.firetasks.parse_outputs.CalcDiff``
        """
        tasks = []
        working_dir = working_dir or os.getcwd()

        if md_property == "diffusion":
            diff_kwargs = {
                i: j
                for i, j in kwargs.items()
                if i in CalcDiff.required_params + CalcDiff.optional_params
            }
            diff_kwargs.update({"working_dir": working_dir})
            tasks.append(CalcDiff(**diff_kwargs))

        elif md_property == "rdf":
            rdf_kwargs = {
                i: j
                for i, j in kwargs.items()
                if i
                in inspect.getfullargspec(calc_atomic_rdf).args
                + inspect.getfullargspec(calc_molecular_rdf).args
                + GetRDF.required_params
                + GetRDF.optional_params
            }
            rdf_kwargs.update({"working_dir": working_dir})
            tasks.append(GetRDF(**rdf_kwargs))
        spec = kwargs.pop("spec", {})
        spec.update({"tag": tag, "_launch_dir": working_dir})
        super(RunAnalysisFW, self).__init__(
            tasks,
            parents=parents,
            name=name,
            spec=spec,
            **{i: j for i, j in kwargs.items() if i in FIREWORK_KWARGS}
        )
