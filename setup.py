from setuptools import setup, find_packages
from scripts.version import __version__

with open("README.md", "r") as desc_file:
    long_description = desc_file.read()

setup(
    name="SR-AMR",
    version=__version__,
    description="Generating binary matrices suitable for machine learning from genomic fasta data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT',
    author="Alper Yurtseven",
    maintainer="Alper Yurtseven",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    packages=find_packages(),
    setup_requires=['setuptools_scm'],
    include_package_data=True,
    install_requires=[],
    package_data={},
    python_requires=">=3.8, <4.0.0",
    keywords="bioinformatics",
)