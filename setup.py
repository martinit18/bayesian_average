from setuptools import find_packages, setup
from bayesian_average._version import __version__ 

setup(
    name='bayesian_average',
    packages=find_packages(include = ['bayesian_average']),
    version=__version__,
    description='Python package to calculate the weighted average with Bayesian methods',
    long_description='The long description will be added in the future.',
    author='Marleen Maxton',
    install_requires=[],
    license = 'X11'
)
