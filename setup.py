import os
from setuptools import find_packages, setup


version_namespace = {}
with open(os.path.join(os.path.dirname(__file__), 'bayesian_average', '_version.py')) as file:
    exec(file.read(), version_namespace)


def read_file(filename):
    with open(os.path.join(os.path.dirname(__file__), filename)) as file:
        return file.read()

setup(
    name='bayesian_average',
    packages=find_packages(include = ['bayesian_average']),
    version=version_namespace['__version__'],
    description='Bayesian weighted averages for inconsistent data sets',
    long_description=read_file('README.md'),
    long_description_content_type='text/markdown',
    author='Marleen Maxton, Martino Trassinelli',
    install_requires=['numpy', 'sympy', 'scipy', 'matplotlib'],
    license = 'X11'
)
