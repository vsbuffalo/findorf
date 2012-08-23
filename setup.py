#!/usr/bin/env python
try:
    from setuptools import setup, find_packages
except ImportError:
    sys.exit("error: install setuptools")

setup(
    name='findorf',
    version=1.02,
    author='Vince Buffalo',
    author_email='vsbuffalo@gmail.com',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    entry_points = {
        'console_scripts': [
            'findorf = findorf.findorf:main']
        },
    url='http://github.com/vsbuffalo/findorf',
    license='GPL 2.0',
    description='An ORF finder and annotator.',
    requires=["BioPython"]
    )   
