# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- Add docstrings for lammps workflows.

### Fixed

- Update installation to use python 3.9 or higher.
- Update the "cluster_analysis" firetask inputs to match recent changes in the MDPropTools package.

## [0.0.4] - 2023-07-15

### Added

- Add support for retrieving molecules from pubchem directly and using them in the workflows.
- Add support for OPLS 2005 ff by running Maestro in the backend.
- Add an option for charge scaling of ionic species when preparing the lammps data file.
- Allow identification of system element types when creating dump files in lammps simulations.

## [0.0.2] - 2022-08-29
