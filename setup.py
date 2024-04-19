import os

import setuptools

module_dir = os.path.dirname(os.path.abspath(__file__))

setuptools.setup(
    name="mispr",
    version="0.0.4",
    author="Rasha Atwi, Matthew Bliss",
    author_email="rasha.atwi@stonybrook.edu, matthew.bliss@stonybrook.edu",
    description="mispr contains FireWorks workflows for Materials Science",
    long_description=open(os.path.join(module_dir, "README.md")).read(),
    url="https://github.com/molmd/mispr",
    install_requires=[
        "numpy >= 1.21.1",
        # "pymongo >= 3.11.0",
        # "pymongo <= 3.12.0",
        "pymongo>=3.3.0,<=3.12.0",
        "matplotlib >= 3.3.1",
        "networkx >= 2.5",
        "FireWorks >= 2.0.3",
        # "monty >= 4.0.0",
        "scipy >= 1.5.2",
        "pandas >= 1.1.2",
        "pubchempy",
        "parmed",
        "mdproptools",
        "dnspython",
        "ruamel.yaml>=0.15.35,<=0.17.40",  # should be removed in the future, but is needed since FireWorks has not updated their function calls to match newer versions of ruamel.yaml
    ],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Information Technology",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
    ],
    python_requires=">=3.9",
    package_data={"": ["gaussian/data/*.bib", "lammps/data/*.json"]},
)
