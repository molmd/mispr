import os
import ntpath
import shutil

from copy import deepcopy
from fireworks.core.firework import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from mispr.gaussian.workflows.base.nmr import get_nmr_tensors
from mispr.gaussian.utilities.fw_utilities import run_fake_gaussian


@explicit_serialize
class NMRFromMD(FiretaskBase):
    required_params = []
    optional_params = [
        "working_dir",
        "db",
        "opt_gaussian_inputs",
        "freq_gaussian_inputs",
        "nmr_gaussian_inputs",
        "solvent_gaussian_inputs",
        "solvent_properties",
        "cart_coords",
        "oxidation_states",
        "additional_kwargs",
    ]

    def run_task(self, fw_spec):
        nmr_wfs = []
        working_dir = self.get("working_dir", os.getcwd())
        additional_kwargs = self.get("additional_kwargs", {})
        for key in ["mol_name", "skips", "process_mol_func", "charge"]:
            additional_kwargs.pop(key, None)
        top_config_files = sorted(fw_spec.get("top_config_files"))

        # used only if running fake gaussian calculations
        fake_gaussian_kwargs = additional_kwargs.get("fake_gaussian_kwargs", {})
        ref_dirs = fake_gaussian_kwargs.get("ref_dirs", [])
        input_files = fake_gaussian_kwargs.get(
            "input_files", ["mol.com"] * 3 * len(top_config_files)
        )

        for ind, file in enumerate(top_config_files):
            config_file = ntpath.basename(file)
            shutil.copy(file, f"{working_dir}/{config_file}")
            nmr_wf = get_nmr_tensors(
                mol_operation_type="get_from_file",
                mol=config_file,
                db=self.get("db"),
                working_dir=working_dir,
                opt_gaussian_inputs=deepcopy(self.get("opt_gaussian_inputs")),
                freq_gaussian_inputs=deepcopy(self.get("freq_gaussian_inputs")),
                nmr_gaussian_inputs=deepcopy(self.get("nmr_gaussian_inputs")),
                solvent_gaussian_inputs=self.get("solvent_gaussian_inputs"),
                solvent_properties=self.get("solvent_properties"),
                cart_coords=self.get("cart_coords"),
                oxidation_states=self.get("oxidation_states"),
                skips=None,
                process_mol_func=False,
                mol_name=config_file.strip(".xyz"),
                **additional_kwargs,
            )

            # added just for testing purposes
            if fake_gaussian_kwargs:
                nmr_wf = run_fake_gaussian(
                    nmr_wf,
                    ref_dirs=ref_dirs[ind * 3: ind * 3 + 3],
                    input_files=input_files[ind * 3: ind * 3 + 3],
                    tolerance=fake_gaussian_kwargs.get("tolerance")
                )
            nmr_wfs.append(nmr_wf)
        return FWAction(detours=nmr_wfs)
