import setuptools
from setuptools import setup
import os
import sys
import subprocess

long_description = "kmermaid: Ultrafast functional annotations of shotgut metagenomic sequencing into protein clusters based on K-mer frequencies"

requirements = [
    'numpy>=1.22.3',
    'pandas>=1.3.5',
    'argparse',
]


setup(
    name="kmermaid",
    version="1",
    author="Add list",
    author_email="Anastasia.Lucas@pennmedicine.upenn.edu, nauslander@wistar.org",
    description="Fast functional classification of metagenomic reads",
    long_description=long_description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    install_requires=requirements,
    test_suite='nose.collector',
    tests_require=['nose'],
    entry_points = {
        'console_scripts': [
            'kmermaid=kmermaid.command_line:kmermaid_predict',
        ],
    }
)
