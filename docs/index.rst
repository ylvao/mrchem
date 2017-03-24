.. MRChem documentation master file, created by
   sphinx-quickstart on Tue Jan 26 15:03:29 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==================================
Welcome to MRChem's documentation!
==================================

MRChem is a numerical real-space code for molecular
electronic structure calculations within the self-consistent field (SCF)
approximations of quantum chemistry (Hartree-Fock and Density Functional
Theory). The code is divided in two main parts: the MultiResolution
Computation Program Package (MRCPP), which is a general purpose numerical
mathematics library based on multiresolution analysis and the multiwavelet
basis which provide low-scaling algorithms as well as rigorous error control
in numerical computations, and the MultiResolution Chemistry (MRChem) program
that uses the functionalities of MRCPP for computational chemistry applications.

The code is being developed at the `Centre for Theoretical and Computational
Chemistry <http://www.ctcc.no/>`_ (CTCC) at
`UiT - The Arctic University of Norway <http://en.uit.no>`_.

--------------------------------------------------------------------------------

We are currently in the process of rewriting the code and making it publicly
available, and the latest version (with limited functionality) can be found on
`GitHub <https://github.com/MRChemSoft/mrchem>`_. This is **not** a stable
version, expect major changes in the future.

Features as of March 2017:

* Restricted (closed-shell) Hartree-Fock and Kohn-Sham DFT (LDA, GGA and hybrid)
* Unrestricted Hartree-Fock (no spin-DFT)
* Properties:
    + Ground state energy
    + Dipole moment
* Localized orbitals
* Shared memory parallelization (OpenMP)
* Current limitations on a single medium-memory compute node (128GB):
    + nano-Hartree accuracy: ~10 orbitals
    + micro-Hartree accuracy: ~50 orbitals
    + milli-Hartree accuracy: ~100 orbitals

Upcoming features:

* Spin-polarized DFT
* Linear response properties (electric, magnetic)
* Distributed memory parallelization (MPI + OpenMP)


.. toctree::
   :maxdepth: 2

   installation
   mrcpp_api
   mrchem_manual

