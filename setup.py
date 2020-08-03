#!/usr/bin/env python
import os

from setuptools import setup

here = os.path.abspath(os.path.dirname(__file__))

# To update the package version number, edit dynamics/__version__.py
version = {}
with open(os.path.join(here, 'dynamics', '__version__.py')) as f:
    exec(f.read(), version)

with open('README.rst') as readme_file:
    readme = readme_file.read()

setup(
    name='dynamics',
    version=version['__version__'],
    description="Molecular dynamics using interpolate potential energy surfaces",
    long_description=readme + '\n\n',
    author="Felipe Zapata",
    author_email='f.zapata@esciencecenter.nl',
    url='https://github.com/resmol/dynamics',
    packages=[
        'dynamics',
    ],
    include_package_data=True,
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='dynamics',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    install_requires=[
        'mendeleev', 'numpy', 'pyparsing', 'scipy'],

    extras_require={
        'test': ['coverage', 'pycodestyle', 'pytest>=3.9', 'pytest-cov',
                 'pytest-mock'],
        'doc': ['sphinx', 'sphinx-autodoc-typehints', 'sphinx_rtd_theme']
    }
)
