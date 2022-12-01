from setuptools import find_packages, setup

with open('deimos/__init__.py') as f:
    exec([x for x in f.readlines() if '__version__' in x][0])

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

with open('requirements.txt') as f:
    requirements = f.read()

pkgs = find_packages(exclude=('test'))

setup(
    name='uimfpy',
    version=__version__,
    description='A Python library for universal ion mobility format (UIMF) files',
    long_description=readme,
    url='https://github.com/PNNL-Comp-Mass-Spec/uimfpy',
    install_requires=requirements,
    python_requires='>=3.7',
    license=license,
    packages=pkgs
)
