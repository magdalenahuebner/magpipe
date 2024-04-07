from setuptools import find_packages, setup

setup(
    name='magpipe',
    packages=find_packages(where='src'),  # Finds packages in src
    package_dir={'': 'src'},  # Tells distutils packages are under src
    version='0.1.0',
    description='A bioinformatics pipeline for performing differential phosphosite occupancy analysis on phosphoproteomics data.',
    author='Magdalena Huebner',
    license='MIT',
)
