from setuptools import setup, find_packages
import os

# Read the README file if available
long_description = open('README.md').read() if os.path.exists('README.md') else ''

setup(
    name='gene_copy_number_workflow',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'matplotlib',
        'seaborn',
        'scipy',
        'pysam',
        'gffutils'
    ],
    entry_points={
        'console_scripts': [
            'gene-copy-number=gene_copy_number.main:main',
        ],
    },
    author='Carl Malina',
    description='End-to-end gene copy number estimation workflow using BAM and GFF3 files.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.10',
)