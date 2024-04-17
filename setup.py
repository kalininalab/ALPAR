from setuptools import setup, find_packages
from sr_amr.version import __version__

with open("README.md", "r") as desc_file:
    long_description = desc_file.read()

setup(
    name="SR-AMR",
    version=__version__,
    description="Pipeline for generating binary matrices, from single reference, suitable for machine learning from genomic fasta data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT',
    author="Alper Yurtseven",
    maintainer="Alper Yurtseven",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3.9",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    packages=find_packages(),
    # setup_requires=['setuptools_scm'],
    include_package_data=False,
    install_requires=[],
    package_data={},
    python_requires="==3.9.*",
    keywords="bioinformatics",
)