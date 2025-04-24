########
ALPAR
########

ALPAR, short for Automated Learning Pipeline for Antimicrobial Resistance, serves as a tool crafted to generate antimicrobial resistance predicting machine learning models. Licensed under the MIT license, it is open source and conveniently accessible on
`GitHub <https://github.com/kalininalab/ALPAR>`_. Installation is made simple through
`conda <https://anaconda.org/kalininalab/ALPAR>`_.

Install
#######

ALPAR is installable from conda:.

.. note::
    It is recommended to use `mamba <https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html>`_
    for the installation because conda might not be able to resolve the dependencies of ALPAR successfully.

.. raw:: html
    :file: install.html

To install it into the new environment:

mamba create -n alpar -c conda-forge -c kalininalab -c bioconda -c etetoolkit alpar
conda activate alpar
pip install panacota
Or to install it into the already existing environment:

mamba install -c conda-forge -c kalininalab -c bioconda -c etetoolkit alpar
pip install panacota
