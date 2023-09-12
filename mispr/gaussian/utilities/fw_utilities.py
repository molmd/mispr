"""Define utility functions for modifying workflow settings. Based on atomate powerups."""

import logging

from mispr.gaussian.firetasks.run_calc import (
    RunGaussianFake,
    RunGaussianDirect,
    RunGaussianCustodian,
)

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.4"
__credits__ = "Anubhav Jain, Kiran Mathew"

logger = logging.getLogger(__name__)


def get_list_fireworks_and_tasks(
    workflow, firework_substring=None, task_substring=None
):
    """
    Return a list of (firework_index, task_index) tuples for all fireworks and tasks in
    a workflow.

    Args:
        workflow (Workflow): The workflow to search.
        firework_substring (str, optional): A substring to search for in the Firework
            names to exclude certain fireworks.
        task_substring (str, optional): A substring to search for in the Firetask
            names to exclude certain Firetasks.

    Returns:
        list: A list of (firework_index, task_index) tuples.
    """
    list_fireworks_and_tasks = []
    for i_firework, firework in enumerate(workflow.fws):
        if not firework_substring or any(
            n in firework.name for n in firework_substring
        ):
            for i_task, task in enumerate(firework.tasks):
                if not task_substring or any(n in str(task) for n in task_substring):
                    list_fireworks_and_tasks.append((i_firework, i_task))
    return list_fireworks_and_tasks


def control_worker(
    workflow, firework_substring=None, task_substring=None, fworker=None, category=None
):
    """
    Modify the Firework's fworker name and category in a workflow. Can be used when
    running workflows on multiple workers at the same time to specify which
    worker/machine to use.

    Args:
        workflow (Workflow): The workflow to control.
        firework_substring (str, optional): A substring to search for in the Firework
            names to exclude certain fireworks.
        task_substring (str, optional): A substring to search for in the Firetask
            names to exclude certain Firetasks.
        fworker (str, optional): The name of the fworker to use for the Firework;
            should be consistent with the one specified in the FireWorker
            (my_fworker.yaml file).
        category (str, optional): The category to be assigned for the Firework; should
            be consistent with the one specified in the FireWorker (my_fworker.yaml file).

    Returns:
        Workflow: The modified workflow with the specified fworker and/or category.
    """
    list_fireworks_and_tasks = get_list_fireworks_and_tasks(
        workflow, firework_substring=firework_substring, task_substring=task_substring
    )
    for i_firework, i_task in list_fireworks_and_tasks:
        if fworker:
            workflow.fws[i_firework].spec["_fworker"] = fworker
        if category:
            workflow.fws[i_firework].spec["_category"] = category
    return workflow


def modify_queue_parameters(
    workflow,
    ntasks_per_node=None,
    walltime=None,
    queue=None,
    pre_rocket=None,
    other_parameters=None,
    firework_substring=None,
    task_substring=None,
):
    """
    Modify the default Firework's queue parameters in a workflow. Default ones are
    specified in the my_qadapter.yaml file. Helpful when different workflows requires
    different computational resources (e.g. number of CPUs, memory, etc.).

    Args:
        workflow (Workflow): The workflow to modify.
        ntasks_per_node (int, optional): The number of tasks to run on each node.
        walltime (str, optional): The walltime for the job.
        queue (str, optional): The queue/partition to run the job on.
        pre_rocket (str, optional): The pre-rocket command to run before the job.
        other_parameters (dict, optional): Other parameters to be added to the
            queueadapter.
        firework_substring (str, optional): A substring to search for in the Firework
            names to exclude certain fireworks.
        task_substring (str, optional): A substring to search for in the Firetask
            names to exclude certain Firetasks.

    Returns:
        Workflow: The modified workflow with the specified queue parameters.
    """
    queue_parameters = {}
    if ntasks_per_node:
        queue_parameters.update({"ntasks_per_node": ntasks_per_node})
    if walltime:
        queue_parameters.update({"walltime": walltime})
    if queue:
        queue_parameters.update({"queue": queue})
    if pre_rocket:
        queue_parameters.update({"pre_rocket": pre_rocket})
    if other_parameters:
        queue_parameters.update(other_parameters)

    list_fireworks_and_tasks = get_list_fireworks_and_tasks(
        workflow, firework_substring=firework_substring, task_substring=task_substring
    )
    for i_firework, i_task in list_fireworks_and_tasks:
        workflow.fws[i_firework].spec.update({"_queueadapter": queue_parameters})
    return workflow


def replace_runtask(
    workflow,
    firework_substring=None,
    operation="remove_custodian",
    additional_params=None,
):
    """
    Replace all tasks with ``RunGaussian`` (e.g. RunGaussianDirect) with
    RunGaussianCustodian or vice versa.

    Args:
        workflow (Workflow): The workflow to modify.
        firework_substring (str, optional): A substring to search for in the Firework
            names to exclude certain fireworks.
        operation (str, optional): The operation to perform on the Firetask; supported
            ones are ``remove_custodian`` and ``use_custodian``.
        additional_params (dict, optional): Additional parameters to be added to the
            new Firetask that are not included in the original Firetask; refer to the
            corresponding Firetask documentation for supported parameters.

    Returns:
        Workflow: The workflow with the replaced run Firetasks.
    """
    if additional_params is None:
        additional_params = {}
    list_fireworks_and_tasks = get_list_fireworks_and_tasks(
        workflow, firework_substring, task_substring=["RunGaussian"]
    )

    for i_firework, i_task in list_fireworks_and_tasks:
        params = workflow.fws[i_firework].tasks[i_task].as_dict()
        task_params = {**params, **additional_params}
        if operation == "remove_custodian":
            firetask = RunGaussianDirect
        elif operation == "use_custodian":
            firetask = RunGaussianCustodian
        else:
            raise ValueError(f"Unsupported operation type: {operation}")
        workflow.fws[i_firework].tasks[i_task] = firetask(
            **{
                i: j
                for i, j in task_params.items()
                if i in firetask.required_params + firetask.optional_params
            }
        )
    return workflow


def run_fake_gaussian(workflow, ref_dirs, input_files=None, tolerance=None):
    """
    Replace all tasks with ``RunGaussian`` (i.e. RunGaussianDirect, RunGaussianCustodian)
    with RunGaussianFake that runs a fake Gaussian job. We do not actually run Gaussian
    but copy existing inputs and outputs. Useful for testing purposes.

    Args:
        workflow (Workflow): The workflow to modify.
        ref_dirs (list): A list of directories containing the reference calculations
            for the fake Gaussian job (e.g. ['home/opt', 'home/freq']).
        input_files (list, optional): A list of input files for the fake Gaussian job;
            order should match that in ref_dirs; e.g. ["opt.com", "freq.com"].
        tolerance (float, optional): The tolerance for the comparison of the provided
            input file with the existing one.

    Returns:
        Workflow: The workflow with the replaced run Firetasks.
    """
    list_fireworks_and_tasks = get_list_fireworks_and_tasks(
        workflow, task_substring=["RunGaussian"]
    )
    if not input_files:
        input_files = ["mol.com"] * len(ref_dirs)

    for ind, (i_firework, i_task) in enumerate(list_fireworks_and_tasks):
        workflow.fws[i_firework].tasks[i_task] = RunGaussianFake(
            ref_dir=ref_dirs[ind],
            input_file=input_files[ind],
            tolerance=tolerance,
        )
    return workflow


def add_common_mods(workflow, fw_mods=None):
    """
    Wrapper function to add common modifications to a workflow.

    Args:
        workflow (Workflow): The workflow to modify.
        fw_mods (dict, optional): A dictionary of modifications to be applied to the
            workflow; supported ones are ``CONTROL_WORKER``, ``MODIFY_QUEUE_PARAMETERS``,
            ``REPLACE_RUNTASK``, and ``RUN_FAKE_GAUSSIAN`` (see the docstring of each
            function for more details); values of the dictionary are the inputs to the
            corresponding function.

    Returns:
        Workflow: The modified workflow.
    """
    fw_mods = fw_mods or {}

    if fw_mods.get("CONTROL_WORKER"):
        workflow = control_worker(workflow, **fw_mods.get("CONTROL_WORKER", {}))

    if fw_mods.get("MODIFY_QUEUE_PARAMETERS"):
        workflow = modify_queue_parameters(
            workflow, **fw_mods.get("MODIFY_QUEUE_PARAMETERS", {})
        )

    if fw_mods.get("REPLACE_RUNTASK"):
        workflow = replace_runtask(workflow, **fw_mods.get("REPLACE_RUNTASK", {}))

    if fw_mods.get("RUN_FAKE_GAUSSIAN"):
        workflow = run_fake_gaussian(workflow, **fw_mods.get("RUN_FAKE_GAUSSIAN", {}))

    return workflow
