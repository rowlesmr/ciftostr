========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |requires|
        |
    * - package
      - | |version| |wheel| |supported-versions| |supported-implementations|
        | |commits-since|
.. |docs| image:: https://readthedocs.org/projects/ciftostr/badge/?style=flat
    :target: https://ciftostr.readthedocs.io/
    :alt: Documentation Status

.. |requires| image:: https://requires.io/github/rowlesmr/ciftostr/requirements.svg?branch=master
    :alt: Requirements Status
    :target: https://requires.io/github/rowlesmr/ciftostr/requirements/?branch=master

.. |version| image:: https://img.shields.io/pypi/v/ciftostr.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/ciftostr

.. |wheel| image:: https://img.shields.io/pypi/wheel/ciftostr.svg
    :alt: PyPI Wheel
    :target: https://pypi.org/project/ciftostr

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/ciftostr.svg
    :alt: Supported versions
    :target: https://pypi.org/project/ciftostr

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/ciftostr.svg
    :alt: Supported implementations
    :target: https://pypi.org/project/ciftostr

.. |commits-since| image:: https://img.shields.io/github/commits-since/rowlesmr/ciftostr/v0.0.0.svg
    :alt: Commits since latest release
    :target: https://github.com/rowlesmr/ciftostr/compare/v0.0.0...master



.. end-badges

This program converts an arbitrary number of CIF (Crystallographic Information Format) files, each containing an arbitrary number of crystal structures, into a number of individual STR files that are compatible with the Rietveld analysis software, TOPAS. 

If the program is run with zero command line arguments, then a GUI is launched, otherwise it is command-line driven.

* Free software: Apache Software License 2.0

Installation
============

To install the release version from PyPi::

    pip install ciftostr

You can also install the in-development version with::

    pip install https://github.com/rowlesmr/ciftostr/archive/master.zip


How to run
==========

If you would like to run it in the command line, you need to provide some command-line arguments::

	ciftostr *.cif
	ciftostr C:\user\username\some\directory\name.cif D:\data\file.cif
	ciftostr /nfs/an/disks/jj/home/dir/file.cif home/folder/nest/deeper/important.cif

will convert all CIF files in the current directory or the two CIF files specifically - both Windows and Linux filepaths are accepter. The resulting STR files will be saved in the same path as their parent CIF file. Relative or absolute file paths can be used.

To launch the GUI::

	python ciftostr.pyz


Bonus Windows executable behaviour: You can drag a CIF file onto the program icon, and it will convert the CIF for you.

To see some other helpful information, try::

	python ciftostr.pyz --help
	python ciftostr.pyz --info


jEdit integration
=================

`Link text <link URL>`_

`John Evans <http://topas.dur.ac.uk>`_ has made available `a number of macros <http://topas.dur.ac.uk/topaswiki/doku.php?id=jedi>`_ for the text editor `jEdit <http://www.jedit.org/>`_ that make it integrate very well with TOPAS when operating in launch mode.

As a part of these macros, there is the ability to insert CIF files in STR format. If you wish to use ``ciftostr`` to generate the STR information, then you need to make the following changes:

0. ``pip install ciftostr`` (or use ``ciftostr.exe``)
1. Replace your copy of ``TAInsertCIF.bsh`` with `this one <TAInsertCIF.bsh>`_. (Your file is probably found in C:\\Users\\?????\\AppData\\Roaming\\jEdit\\macros)

Now, when you choose "``Insert CIFs in INP format``" from the plugin in jEdit, ``ciftostr`` will be used in the background to generate that format.


Documentation
=============


https://ciftostr.readthedocs.io/

This program is designed to reformat data from a CIF format into the STR format suitable for use by the Rietveld analysis software, TOPAS. It doesn't carry out any validation checks to ensure the data is consistent, but simply transcribes the data as given in the CIF.
    
Where similar or identical data could be given in several places, the values are taken in a specific order of precedence, as detailed in the sections below. In general, if a value exists in an earlier place, the later places are not looked at.

This program uses the `PyCifRW <https://bitbucket.org/jamesrhester/pycifrw/src/development>`_ library, written by James Hester, to parse the CIF files into a format easily used to remix the underlying data. The GUI was created using `PySimpleGUI <https://pysimplegui.readthedocs.io/en/latest>`_. 
    
The STR's phase\name is taken from ``\chemical\name\mineral``, ``\chemical\name\common``, ``\chemical\name\systematic``, or ``\chemical\name\structure\type``, in that order, appended with the value of the ``data`` block. If none of these keys are available, then the name of the ``data\`` block is used. This is also used as the filename of the STR file.
    
The unit cell parameters are taken from ``\cell\length\a``, ``\cell\length\b``, ``\cell\length\c``, ``\cell\angle\alpha``, ``\cell\angle\beta``, and ``\cell\angle\gamma``. Some comparisons are made to enable some standard macros to be used (eg Cubic, Tetragonal...). In any of these fail, then all parameters are given as a fail safe.

The space\group is taken from ``\symmetry\space\group\name\H-M``, ``\space\group\name\H-M\alt``, ``\symmetry\Int\Tables\number``, or ``\space\group\IT\number``, in that order. If none of these keys are available, then an empty string is written as the space group. The value of the space group is as exactly given in the CIF. No validation or editing is done.

The atomic sites are constructed as follows: The atom labels are taken from ``\atom\site\label``, with the fractional x, y, and z coordinates given by ``\atom\site\fract\x``, ``\y``, and ``\z``. If the decimal values of the fractional coordinates are consistent with the fractions 1/6, 1/3, 2/3, or 5/6, then the decimal value is replaced by the fractions. The site occupancy is given by ``\atom\site\occupancy``, or by ``1``, if that key is not given. The atom type is given by ``\atom\site\type\symbol``, where available, or by the first one or two characters of the site label. If these characters match an element symbol, then that is used, otherwise, the label is used in it's entirety, and the user must decide the correct atom type to use. If the site label starts with ``Wat``, (and no atom type is given) it is assumed that oxygen (from water) is correct. An attempt is also made to reorder the charge given on an atom, to ensure it is compatible with TOPAS ordering, eg ``Fe+2``, not ``Fe2+``.

Isotropic Atomic Displacement Parameters (ADPs; Biso), are taken from ``\atom\site\B\iso\or\equiv``, or from ``\atom\site\U\iso\or\equiv``, where they are then multiplied by 8*Pi^2 to give B values. If anisotropic ADPs are given, then ``\atom\site\aniso\B\11``, ``\atom\site\aniso\B\22``, and ``\atom\site\aniso\B\33`` are averaged to give an equivalent Biso value. Alternatively, the equivalent U values are used to calculate Biso. As anisotropic values could be given for a subset of the atoms in the structure, the labels given by ``\atom\site\label`` and ``\atom\site\aniso\label`` are matched, and if an atom doesn``t have an anisotropic value, it takes its isotropic value, or is assigned a value of ``1``.

The atomic site is also given a ``num\posns 0`` entry, which will update with the multiplicity of the site following a refinement. This will allow the user to compare this value with the CIF or Vol A to help ensure that the correct symmetry is being applied.

Finally, the STR is given a fixed Lorentzian crystallite size of 200 nm, and a refinable scale factor of 0.0001 to allow for an easy start to a refinement. All other values given in the STR are fixed, and require active intervention to name, refine, constrain, or restrain them.

If you have any feedback, please contact me. If you find any bugs, please provide the CIF which caused the error, a description of the error, and a description of how you believe the program should work in that instance.


Development
===========

Come and talk to me!