# coding: utf-8


# Contains utility functions for modifying workflow settings. 

import logging

__author__ = 'Rasha Atwi'
__maintainer__ = 'Rasha Atwi'
__email__ = 'rasha.atwi@stonybrook.edu'
__status__ = 'Development'
__date__ = 'Jan 2021'
__version__ = 0.2

logger = logging.getLogger(__name__)


def control_worker(workflow, fworker=None, category=None):
    list_fireworks_and_tasks = get_list_fireworks_and_tasks(workflow)
    for i_firework, i_task in list_fireworks_and_tasks:
        if fworker:
            workflow.fws[i_firework].spec['_fworker'] = fworker
        if category:
            workflow.fws[i_firework].spec['_category'] = category
    return workflow


def get_list_fireworks_and_tasks(workflow):
    list_fireworks_and_tasks = []
    for i_firework, firework in enumerate(workflow.fws):
        for i_task, task in enumerate(firework.tasks):
            list_fireworks_and_tasks.append((i_firework, i_task))
    return list_fireworks_and_tasks


def modify_queue_parameters(workflow, ntasks_per_node=None, walltime=None,
                            queue=None, pre_rocket=None, other_parameters=None):
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

    list_fireworks_and_tasks = get_list_fireworks_and_tasks(workflow)
    for i_firework, i_task in list_fireworks_and_tasks:
        workflow.fws[i_firework].spec.update(
            {'_queueadapter': queue_parameters})
    return workflow
