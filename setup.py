#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import io
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

from setuptools import find_packages
from setuptools import setup


def read(*names, **kwargs):
    with io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ) as fh:
        return fh.read()


setup(
    name='ciftostr',
    version='0.2.0',
    license='Apache-2.0',
    description='A small utility to convert crystal structures in CIF format into TOPAS STR format.',
    long_description='%s\n%s' % (
        re.compile('^.. start-badges.*^.. end-badges', re.M | re.S).sub('', read('README.rst')),
        re.sub(':[a-z]+:`~?(.*?)`', r'``\1``', read('CHANGELOG.rst'))
    ),
    author='Matthew Rowles',
    author_email='rowlesmr@gmail.com',
    url='https://github.com/rowlesmr/ciftostr',
    package_dir={'': 'src'},
	packages=find_packages('src', exclude=['exe*','data*']),
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=False,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 4 - Beta',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
        # uncomment if you test on these interpreters:
        # 'Programming Language :: Python :: Implementation :: IronPython',
        # 'Programming Language :: Python :: Implementation :: Jython',
        # 'Programming Language :: Python :: Implementation :: Stackless',
        'Topic :: Utilities',
        'Topic :: Scientific/Engineering',
    ],
    project_urls={
        'Documentation': 'https://ciftostr.readthedocs.io/',
        'Changelog': 'https://ciftostr.readthedocs.io/en/latest/changelog.html',
        'Issue Tracker': 'https://github.com/rowlesmr/ciftostr/issues',
    },
    keywords=[
        'diffraction', 'X-ray', 'neutron',
        'synchrotron', 'powder diffraction', 'CIF',
        'crystallography', 'Rietveld', 'TOPAS',
    ],
    python_requires='>=3.6',
    install_requires=[
        # eg: 'aspectlib==1.1.1', 'six>=1.7',
	'PyCifRW>=4.4.3',
	'PySimpleGUI>=4.39.1',
    ],
    extras_require={
        # eg:
        #   'rst': ['docutils>=0.11'],
        #   ':python_version=="2.6"': ['argparse'],
    },
    entry_points={
        'console_scripts': [
            'ciftostr = ciftostr.cli:main',
        ]
    },
)
