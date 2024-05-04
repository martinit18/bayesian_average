import os
from setuptools import find_packages, setup
from bayesian_average._version import __version__ 
from bayesian_average.__init__ import __doc__ 


def read_file(filename):
    with open(os.path.join(os.path.dirname(__file__), filename)) as file:
        return file.read()

setup(
    name='bayesian_average',
    packages=find_packages(include = ['bayesian_average']),
    version=__version__,
    description='__doc__',
    long_description=read_file('README.md'),
    long_description_content_type='text/markdown',
    author='Marleen Maxton, Martino Trassinelli',
    install_requires=[],
    license = 'X11'
)
