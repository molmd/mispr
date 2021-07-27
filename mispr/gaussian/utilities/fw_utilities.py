# coding: utf-8


# Contains utility functions for modifying workflow settings. 

import logging

from mispr.gaussian.firetasks.run_calc import RunGaussianDirect, \
    RunGaussianCustodian

__author__ = 'Rasha Atwi'
__maintainer__ = 'Rasha Atwi'
__email__ = 'rasha.atwi@stonybrook.edu'
__status__ = 'Development'
__date__ = 'Jan 2021'
__version__ = 0.2
__credits__ = 'Anubhav Jain, Kiran Mathew'

logger = logging.getLogger(__name__)


def get_list_fireworks_and_tasks(workflow, firework_substring=None,
                                 task_substring=None):
    list_fireworks_and_tasks = []
    for i_firework, firework in enumerate(workflow.fws):
        if not firework_substring or \
                any(n in firework.name for n in firework_substring):
            for i_task, task in enumerate(firework.tasks):
                if not task_substring or \
                        any(n in str(task) for n in task_substring):
                    list_fireworks_and_tasks.append((i_firework, i_task))
    return list_fireworks_and_tasks


def control_worker(workflow, firework_substring=None, task_substring=None,
                   fworker=None, category=None):
    list_fireworks_and_tasks = \
        get_list_fireworks_and_tasks(workflow,
                                     firework_substring=firework_substring,
                                     task_substring=task_substring)
    for i_firework, i_task in list_fireworks_and_tasks:
        if fworker:
            workflow.fws[i_firework].spec['_fworker'] = fworker
        if category:
            workflow.fws[i_firework].spec['_category'] = category
    return workflow


def modify_queue_parameters(workflow, ntasks_per_node=None, walltime=None,
                            queue=None, pre_rocket=None, other_parameters=None,
                            firework_substring=None, task_substring=None):
    queue_parameters = {}
    if ntasks_per_node:
        queue_parameters.update({'ntasks_per_node': ntasks_per_node})
    if walltime:
        queue_parameters.update({'walltime': walltime})
    if queue:
        queue_parameters.update({'queue': queue})
    if pre_rocket:
        queue_parameters.update({'pre_rocket': pre_rocket})
    if other_parameters:
        queue_parameters.update(other_parameters)

    list_fireworks_and_tasks = \
        get_list_fireworks_and_tasks(workflow,
                                     firework_substring=firework_substring,
                                     task_substring=task_substring)
    for i_firework, i_task in list_fireworks_and_tasks:
        workflow.fws[i_firework].spec.update(
            {'_queueadapter': queue_parameters})
    return workflow


def replace_runtask(workflow, firework_substring=None,
                    operation='remove_custodian', additional_params=None):
    if additional_params is None:
        additional_params = {}
    list_fireworks_and_tasks = \
        get_list_fireworks_and_tasks(workflow, firework_substring,
                                     task_substring=['RunGaussian'])

    for i_firework, i_task in list_fireworks_and_tasks:
        params = workflow.fws[i_firework].tasks[i_task].as_dict()
        task_params = {**params, **additional_params}
        if operation == 'remove_custodian':
            firetask = RunGaussianDirect
        elif operation == 'use_custodian':
            firetask = RunGaussianCustodian
        else:
            raise ValueError(f'Unsupported operation type: {operation}')
        workflow.fws[i_firework].tasks[i_task] = \
            firetask(**{i: j for i, j in task_params.items()
                        if i in firetask.required_params +
                        firetask.optional_params})
    return workflow
