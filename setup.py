from setuptools import setup, find_packages
import re

with open("README.md", "r") as desc_file:
    long_description = desc_file.read()

with open("sr_amr/version.py", "r") as version_file:
    version_content = version_file.read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_content, re.M)
    if version_match:
        __version__ = version_match.group(1)
    else:
        raise RuntimeError("Unable to find version string.")

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
        "Programming Language :: Python :: 3.12",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    packages=find_packages(),
    # setup_requires=['setuptools_scm'],
    include_package_data=False,
    install_requires=[],
    package_data={},
    python_requires=">=3.12",
    keywords="bioinformatics",
    entry_points={
        'console_scripts': [
            'alpar=sr_amr.amr:main',
        ],
    },
    )