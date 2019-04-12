from setuptools import setup, find_packages

setup(
    name='project',
    version='0.0.1',
    description='Generic PDE structured solver using firedrake',
    author='Michael Wathen',
    author_email='michael.wathen@stfc.ac.uk',
    packages=find_packages(exclude=('examples'))
)
