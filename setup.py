#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import setup, find_packages
import sys, os

version = '0.1'

setup(name='iRNA',
      version=version,
      description="iRNA is an integrative pipeline for genome wide detection of sRNA targets",
      long_description="""\
iRNA provide an analysis based sequence analysis and enrichment analysis. iRNA result files can be visualized through iRNA_vis software.""",
      classifiers=["Development Status :: 4 - Beta","Topic :: Scientific/Engineering :: Bio-Informatics","License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)"], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='sRNA targets, enrichment analysis, intaRNA, bioinformatics',
      author='Amine Ghozlane',
      author_email='amineghozlane@gmail.com',
      url='http://www.cbib.u-bordeaux2.fr/fr/content/iRNA',
      license='LGPL',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=True,
      install_requires=["mpi4py","scipy","numpy","rpy2","sqlitebck","lxml","suds"],
      entry_points={
        'console_scripts':[
        'iRNA_stat = irna.iRNA_stat.iRNA_stat:main',
        'iRNA_seq = irna.iRNA_seq.iRNA_seq:main',
        'iRNA_pred = irna.iRNA_pred.iRNA_pred:main',
        'David2tulip = irna.David2tulip.David2tulip:main']},
      )
