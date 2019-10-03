from fireworks import LaunchPad, Firework, FiretaskBase, FWAction, Workflow
from fireworks.core.rocket_launcher import rapidfire, launch_rocket
from gaussian_input_task import ConvertToMoleculesTask, ConvertMoleculeToGaussianInputTask, RunGaussianDirect
from gaussian_output_task import ParseGaussianOutputFile

if __name__ == "__main__":
    # set up the LaunchPad and reset it
    launchpad = LaunchPad(host="", uri_mode=True)

    # define tasks
    task1 = ConvertToMoleculesTask()

    # assemble FireWork from tasks and give the FireWork a unique id
    fire_work1 = Firework(task1, spec={'filesDir': '/Users/rashaatwi/Desktop'}, name='ConvertToMol', fw_id=1)

    task2 = ConvertMoleculeToGaussianInputTask()
    fire_work2 = Firework(task2, spec={'filesDir': '/Users/rashaatwi/Desktop','functional': 'B3LYP',
                                    'basis_set': 'TZVP', 'route_parameters': {'Opt': None, "NoSymmetry": None},
                                    'link0_parameters': {"%mem": "2GB", "%NProcShared": "8", "%Chk": "Optimization.chk"},
                                    'input_parameters': {}, 'var': 'molecules'}, name='ConvertToInput', fw_id=2)

    task3 = ParseGaussianOutputFile()
    fire_work3 = Firework(task3, spec={'filesDir': '/Users/rashaatwi/Desktop'}, name='ParseOutput', fw_id=3)

    task4 = ConvertMoleculeToGaussianInputTask()
    fire_work4 = Firework(task4, spec={'filesDir': '/Users/rashaatwi/Desktop', 'functional': 'B3LYP',
                                    'basis_set': 'TZVP', 'route_parameters':
                                    {'Freq': None, "NoSymmetry": None},
                                    'link0_parameters': {"%mem": "2GB", "%NProcShared": "8", "%Chk": "Frequency.chk"},
                                    'input_parameters': {}, 'var': 'outputs'}, name='ConvertOutputToInput', fw_id=4)

    # assemble Workflow from FireWorks and their connections by id
    wf = Workflow([fire_work1, fire_work2, fire_work3, fire_work4], links_dict={1: 2, 2: 3, 3: 4})

    # store workflow and launch it
    launchpad.add_wf(wf)
    rapidfire(launchpad)
