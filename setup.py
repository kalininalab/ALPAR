from setuptools import setup, find_packages
from sr_amr.version import __version__

with open("README.md", "r") as desc_file:
    long_description = desc_file.read()

setup(
    name="ALPAR",
    version=__version__,
    description="Automated Learning Pipeline for Antimicrobial Resistance",
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
    python_requires=">=3.9.*",
    keywords="bioinformatics",
)