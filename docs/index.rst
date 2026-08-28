ocal documentation
==================

``ocal`` is a program for the calculation of the orbital tail fraction (OTF) of
organic semiconductor molecules.

Overview
--------

Starting from a molecular or crystal structure file, ``ocal`` generates a Gaussian
input, runs a geometry optimization, converts the checkpoint file, and generates
molecular-orbital cube files for the four frontier orbitals (NHOMO, HOMO, LUMO,
NLUMO). For each orbital it integrates the probability density over the voxels that
lie **outside** the van der Waals surface of the molecule, and reports that fraction
as the OTF.

The OTF quantifies how much of a frontier orbital spills out of the molecular van der
Waals volume, i.e. how much of it is available for intermolecular overlap.

.. raw:: html

   <p align="center">
     <img src="_static/OTF.gif" alt="Orbital tail fraction (OTF) illustration" width="400">
   </p>

.. note::

   The molecule is always treated as a neutral closed-shell singlet (charge/spin fixed
   to ``0 1``), and the functional and basis set are fixed to **B3LYP/6-31G(d,p)** so
   that OTF values remain comparable with previously published results.

Requirements
------------

* Python 3.11 or newer
* NumPy 2.0.2 or newer
* Pandas 2.3.3 or newer
* Gaussian 16 (the ``g16``, ``formchk``, and ``cubegen`` executables must be on ``$PATH``)

Installation
------------

``ocal`` is distributed on PyPI under the name ``yu-ocal`` (the import name and the
command name are both ``ocal``):

.. code-block:: bash

   pip install yu-ocal

NumPy and Pandas are installed automatically; Gaussian 16 is not included and must be
installed separately.

Verify the installation with:

.. code-block:: bash

   ocal --help

Basic usage
-----------

.. code-block:: bash

   ocal <filename> [options]

``ocal`` accepts three classes of input and enters the pipeline at the corresponding
stage:

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Input
     - Extensions
     - Pipeline stages executed
   * - Structure file
     - ``.gjf``, ``.com``, ``.xyz``, ``.mol``, ``.mol2``, ``.cif``
     - gjf generation → Gaussian (Opt=Tight) → formchk → cubegen → OTF
   * - Formatted checkpoint
     - ``.fchk``
     - cubegen → OTF
   * - Cube file
     - ``.cube``
     - OTF only

For ``.fchk`` and ``.cube`` input, ``-s, --skip-gaussian`` is required; conversely it
is rejected for structure-file input.

Options
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Default
     - Description
   * - ``-c, --cpu <number>``
     - ``4``
     - Number of CPUs used by Gaussian and ``cubegen``.
   * - ``-m, --mem <memory>``
     - ``10``
     - Amount of memory in GB used by Gaussian.
   * - ``-s, --skip-gaussian``
     - off
     - Skip the Gaussian calculation and reuse existing results.

Examples
~~~~~~~~

.. code-block:: bash

   # Full run from a molecular structure
   ocal xxx.xyz

   # Use 8 CPUs and 16 GB memory
   ocal xxx.mol -c 8 -m 16

   # Start from an existing Gaussian formatted checkpoint file
   ocal xxx.fchk -s

   # Compute OTF for a single cube file that already exists
   ocal xxx_HOMO.cube -s

The full usage manual, including the description of the generated files and
troubleshooting notes, is available in the
`README <https://github.com/matsui-lab-yamagata-univ/ocal>`_.

API reference
-------------

.. toctree::
   :maxdepth: 2

   ocal

Authors
-------

`Matsui Laboratory, Research Center for Organic Electronics (ROEL), Yamagata University
<https://matsui-lab.yz.yamagata-u.ac.jp/index-e.html>`_

Tomoharu Okada, Koki Ozawa, Hiroyuki Matsui

Email: h-matsui[at]yz.yamagata-u.ac.jp (please replace [at] with @)

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
