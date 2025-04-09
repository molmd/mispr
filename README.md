# <img alt="mispr" src="https://raw.githubusercontent.com/molmd/mispr/master/docs/source/_static/logo.png" width="500">

[![Downloads][downloads-badge]][downloads-link]
[![Monthly Downloads][monthly-downloads-badge]][downloads-link]
[![GitHub tag][github-tag-badge]][github-tag-link]
[![Codacy Badge][codacy-badge]][codacy-link]
[![Website][website-badge]][website-link]
[![GitHub commit activity][commit-badge]][commit-link]

## Overview

MISPR is a software that executes, manages, and stores computational materials science
simulations. It contains pre-defined density functional theory (DFT) and molecular dynamics (MD) workflows to calculate and analyze different
properties of materials. MISPR uses [MDPropTools][mdproptools] to perform MD analysis.

## Installation

You can either download the source from GitHub and compile yourself, or install directly using pip.
Please see the [Installation][install-docs] page for detailed instructions.

## Useful Links

- [MISPR Website][mispr-website]: Visit this site to get an overview of MISPR, check the installation instructions, and follow MISPR tutorials
- [MISPR API Reference][api-docs]
- [Resources][resources]

## How to cite

Please include the following two citations if MISPR and/or MDPropTools were used for an academic study:

- Atwi, R., Bliss, M., Makeev, M., & Rajput, N. N. (2022). [MISPR: An automated infrastructure for high-throughput DFT and MD simulations][paper1]. Scientific Reports, 12(1), 1-16.
- Atwi, R., Chen, Y., Han, K. S., Mueller, K. T., Murugesan, V., & Rajput, N. N. (2022).
  [An automated framework for high-throughput predictions of NMR chemical shifts within liquid solutions][paper2].
  Nature Computational Science, 2(2), 112-122.

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
[install-docs]: https://molmd.github.io/mispr/html/installation/index.html
[mispr-website]: https://molmd.github.io/mispr/
[api-docs]: https://molmd.github.io/mispr/html/py-modindex.html
[resources]: https://molmd.github.io/mispr/html/resources/resources.html
[paper1]: https://www.nature.com/articles/s41598-022-20009-w
[paper2]: https://doi.org/10.1038/s43588-022-00200-9
[license]: https://github.com/molmd/mispr/blob/master/LICENSE
