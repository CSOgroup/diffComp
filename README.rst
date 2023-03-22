.. These are examples of badges you might want to add to your README:
   please update the URLs accordingly

    .. image:: https://api.cirrus-ci.com/github/<USER>/diffComp.svg?branch=main
        :alt: Built Status
        :target: https://cirrus-ci.com/github/<USER>/diffComp
    .. image:: https://img.shields.io/coveralls/github/<USER>/diffComp/main.svg
        :alt: Coveralls
        :target: https://coveralls.io/r/<USER>/diffComp
    .. image:: https://img.shields.io/conda/vn/conda-forge/diffComp.svg
        :alt: Conda-Forge
        :target: https://anaconda.org/conda-forge/diffComp
    .. image:: https://pepy.tech/badge/diffComp/month
        :alt: Monthly Downloads
        :target: https://pepy.tech/project/diffComp


.. image:: https://readthedocs.org/projects/diffComp/badge/?version=latest
        :alt: ReadTheDocs
        :target: https://diffComp.readthedocs.io/en/stable/

.. image:: https://img.shields.io/pypi/v/diffComp.svg
        :alt: PyPI-Server
        :target: https://pypi.org/project/diffComp/


========
diffComp
========


    Python package for the analysis of differential Hi-C compartments


``diffComp`` is a software for the detection of differential chromatin compartmentalization from Hi-C data. It takes as input two compartment segmentation files produced by the `CALDER software <https://github.com/CSOgroup/CALDER>`_ [Liu2021]_.


Setup
=====

To install the library, first clone the repository::

    git clone https://github.com/lucananni93/diffComp.git

enter in the library directory::

    cd diffComp

and install the package by creating a new conda environment::

    conda env create -f environment.yml

Usage
=====

Currently, the package offers the following command line tools:

Detection of compartment repositioning events
---------------------------------------------
Compartment repositioning events (CoREs) are detected with the ``cores`` command::

    usage: cores [-h] [--algo ALGO] [--min_std MIN_STD] [--control1_path CONTROL1_PATH] [--control2_path CONTROL2_PATH] [--signal_path SIGNAL_PATH]
                 [--bed_path BED_PATH] [--coordinates COORDINATES] [--genome GENOME] [--chromosomes CHROMOSOMES] [--verbose] [--very-verbose] [--version]
                 sample1_path sample2_path binsize output_path

    Identifying Compartment Repositioning Events from Calder genomic segmentations. Given sample1 and sample2 Calder segmentations, it identifies regions
    undergoing statistically significant compartment repositioning in sample2 in comparison to sample1. Significance of the repositioning is determined
    using paired control samples (control2 vs control1), which are provided by the user. Usually, replicates of the same experiments are used to model the
    intrinsic biological noise in the Hi-C compartment calls.

    positional arguments:
      sample1_path          Path to the Calder segmentation of sample 1
      sample2_path          Path to the Calder segmentation of sample 2
      binsize               Resolution to use in the analysis
      output_path           Path where to store the identified regions

    optional arguments:
      -h, --help            show this help message and exit
      --algo ALGO           Which algorithm to use for finding CoREs
      --min_std MIN_STD     Maximum standard deviation allowed for segmented regions
      --control1_path CONTROL1_PATH
                            Path(s) to the Calder segmentation(s) to use to use as control 1 (comma-separated)
      --control2_path CONTROL2_PATH
                            Path(s) to the Calder segmentation(s) to use to use as control 2 (comma-separated)
      --signal_path SIGNAL_PATH
                            Path where to store the binned differential signal
      --bed_path BED_PATH   Path where to store the identified regions in BED format
      --coordinates COORDINATES
                            Coordinate system of the input files (zero-based / one-based)
      --genome GENOME       Genome (Default: hg19)
      --chromosomes CHROMOSOMES
                            List of chromosomes to perform the analysis on (Default: all chromosomes, comma-separated)
      --verbose             Set loglevel to INFO
      --very-verbose        Set loglevel to DEBUG
      --version             show program's version number and exit


References
==========

.. [Liu2021] Liu, Y., Nanni, L., Sungalee, S. et al. Systematic inference and comparison of multi-scale chromatin sub-compartments connects spatial organization to cell phenotypes. Nat Commun 12, 2439 (2021). https://doi.org/10.1038/s41467-021-22666-3


Note
====

.. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
    :alt: Project generated with PyScaffold
    :target: https://pyscaffold.org/

This project has been set up using PyScaffold 4.1.2. For details and usage
information on PyScaffold see https://pyscaffold.org/.
