from fireworks import LaunchPad

from mispr.gaussian.workflows.base.esp import get_esp_charges

lpad = LaunchPad.auto_load()

wf, _ = get_esp_charges(
    mol_operation_type="get_from_pubchem",  # (1)!
    mol="monoglyme",
    format_chk=True,
    save_to_db=True,
    save_to_file=True,
    additional_prop_doc_fields={"name": "monoglyme"},
    tag="mispr_tutorial",
)
lpad.add_wf(wf)  # (2)!
