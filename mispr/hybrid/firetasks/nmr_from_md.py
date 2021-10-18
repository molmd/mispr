import os
import shutil
import ntpath
from fireworks.core.firework import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize
from mispr.gaussian.workflows.base.nmr import get_nmr_tensors


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
        top_config_files = fw_spec.get("top_config_files")
        for file in top_config_files:
            config_file = ntpath.basename(file)
            shutil.copy(file, f"{working_dir}/{config_file}")
            nmr_wf = get_nmr_tensors(
                    mol_operation_type="get_from_file",
                    mol=config_file,
                    db=self.get("db"),
                    working_dir=working_dir,
                    opt_gaussian_inputs=self.get("opt_gaussian_inputs"),
                    freq_gaussian_inputs=self.get("freq_gaussian_inputs"),
                    nmr_gaussian_inputs=self.get("nmr_gaussian_inputs"),
                    solvent_gaussian_inputs=self.get("solvent_gaussian_inputs"),
                    solvent_properties=self.get("solvent_properties"),
                    cart_coords=self.get("cart_coords"),
                    oxidation_states=self.get("oxidation_states"),
                    skips=None,
                    process_mol_func=False,
                    mol_name=config_file.strip(".xyz"),
                    **additional_kwargs,
                )
            nmr_wfs.append(nmr_wf)
        return FWAction(detours=nmr_wfs)
