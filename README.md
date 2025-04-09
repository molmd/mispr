# <img alt="mispr" src="https://raw.githubusercontent.com/molmd/mispr/master/docs/source/_static/logo.png" width="500">

[![Downloads][downloads-badge]][downloads-link]
[![Monthly Downloads][monthly-downloads-badge]][downloads-link]
[![GitHub tag][github-tag-badge]][github-tag-link]
[![Codacy Badge][codacy-badge]][codacy-link]
[![Website][website-badge]][website-link]
[![GitHub commit activity][commit-badge]][commit-link]

## Overview

MISPR is a free, open-source software for executing, managing, and storing computational materials science simulations. It can:

- Predict various materials properties through pre-defined workflows
- Integrate density functional theory (DFT) and molecular dynamics (MD) simulations seamlessly
- Process, analyze, and visualize simulation results with built-in utilities
- Handle errors automatically using recipe-like fixes
- Build and manage computational databases in MongoDB backend for easy querying, analysis, and sharing
- Support multiple queuing systems via FireWorks including SLURM, PBS, SGE, etc.
- Scale architecture to handle computations from individual materials to thousands of compounds
- Modify and chain workflows together flexibly
- Leverage open-source libraries including [pymatgen][pymatgen], [custodian][custodian], and [FireWorks][fireworks]

## Workflows

MISPR provides several pre-defined workflows for both Gaussian (DFT) and LAMMPS (MD) calculations:

### DFT (Gaussian) Workflows

- **ESP**: Electrostatic partial charges calculations for charge fitting
- **NMR**: Nuclear magnetic resonance chemical shift predictions
- **BDE**: Bond dissociation energy calculations
- **BE**: Binding energy calculations
- **IP_EA**: Redox potential calculations; methods supported include HOMO/LUMO, vertical IP/EA, adiabatic IP/EA, and sequential PCET; electron transfer calculations can be performed via single-step or multi-step pathways

### MD (LAMMPS) Workflows

A standard workflow for performing classical MD simulations of liquid solutions and subsequently deriving various structural and dynamical properties. The default operations run as follows:

1. Run ESP workflow on all species
2. Fit RESP charges and extract GAFF parameters (OPLS2005 and user-defined parameters are supported)
3. Build initial system configuration
4. Create LAMMPS data file containing initial atomic coordinations, molecular topology, and force field parameters
5. Run two step energy minimization, NPT equilibration run, melting and quenching, and an NVT production run
6. Perform analysis using the generated log and trajectory files (radial distribution function, coordination number, cluster analysis, mean squared displacement and diffusion) using [MDPropTools][mdproptools]

### Hybrid Workflows

- **General**: A workflow that combines the ESP DFT workflow with the MD workflow described above
- **NMR**: A workflow that first runs the general hybrid workflow, then analyzes solvation shells from MD trajectories, and finally runs the NMR DFT workflow on representative solvation structures to predict chemical shifts

All these workflows are easily customizable and can be modified to suit specific research needs.

## Installation

You can either download the source from GitHub and compile yourself, or install directly using pip.
Please see the [Installation][install-docs] page for detailed instructions.

## Useful Links

- [MISPR Website][mispr-website]: Visit this site to get an overview of MISPR, check the installation instructions, and follow MISPR tutorials
- [MISPR API Reference][api-docs]
- [Resources][resources]

## How to cite

Please include the following two citations ([paper1][paper1] and [paper2][paper2]) if MISPR and/or MDPropTools were used for an academic study:

```bib
@article{atwi2022mispr,
  title={MISPR: an open-source package for high-throughput multiscale molecular simulations},
  author={Atwi, Rasha and Bliss, Matthew and Makeev, Maxim and Rajput, Nav Nidhi},
  journal={Scientific Reports},
  volume={12},
  number={1},
  pages={15760},
  year={2022},
  publisher={Nature Publishing Group UK London}
}
```

```bib
@article{atwi2022automated,
  title={An automated framework for high-throughput predictions of NMR chemical shifts within liquid solutions},
  author={Atwi, Rasha and Chen, Ying and Han, Kee Sung and Mueller, Karl T and Murugesan, Vijayakumar and Rajput, Nav Nidhi},
  journal={Nature Computational Science},
  volume={2},
  number={2},
  pages={112--122},
  year={2022},
  publisher={Nature Publishing Group US New York}
}
```

## License Information

MISPR is a free, open-source software package (distributed under the [MIT license][license]).

[downloads-badge]: https://static.pepy.tech/badge/mispr
[monthly-downloads-badge]: https://static.pepy.tech/badge/mispr/month
[downloads-link]: https://pepy.tech/project/mispr
[github-tag-badge]: https://img.shields.io/github/tag/molmd/mispr
[github-tag-link]: https://GitHub.com/molmd/mispr/tags/
[codacy-badge]: https://app.codacy.com/project/badge/Grade/8c047110974a42af9baed409664d2547
[codacy-link]: https://www.codacy.com/gh/molmd/mispr/dashboard?utm_source=github.com&utm_medium=referral&utm_content=molmd/mispr&utm_campaign=Badge_Grade
[website-badge]: https://img.shields.io/website?down_message=down&label=mispr%20website&up_message=up&url=https%3A%2F%2Fmolmd.github.io%2Fmispr%2F
[website-link]: https://molmd.github.io/mispr/
[commit-badge]: https://img.shields.io/github/commit-activity/m/molmd/mispr
[commit-link]: https://github.com/molmd/mispr/commits/master
[mdproptools]: https://github.com/molmd/mdproptools
[fireworks]: https://materialsproject.github.io/fireworks/
[pymatgen]: https://pymatgen.org
[custodian]: https://materialsproject.github.io/custodian/
[install-docs]: https://molmd.github.io/mispr/html/installation/index.html
[mispr-website]: https://molmd.github.io/mispr/
[api-docs]: https://molmd.github.io/mispr/html/py-modindex.html
[resources]: https://molmd.github.io/mispr/html/resources/resources.html
[paper1]: https://www.nature.com/articles/s41598-022-20009-w
[paper2]: https://doi.org/10.1038/s43588-022-00200-9
[license]: https://github.com/molmd/mispr/blob/master/LICENSE
