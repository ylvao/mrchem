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
Theory). The code is divided in two main parts: the `MultiResolution
Computation Program Package <https://mrcpp.readthedocs.io/en/latest>`_ (MRCPP),
which is a general purpose numerical
mathematics library based on multiresolution analysis and the multiwavelet
basis which provide low-scaling algorithms as well as rigorous error control
in numerical computations, and the MultiResolution Chemistry (MRChem) program
that uses the functionalities of MRCPP for computational chemistry applications.

The code is being developed at the `Hylleraas Centre for Quantum Molecular
Sciences <http://www.ctcc.no/>`_ at
`UiT - The Arctic University of Norway <http://en.uit.no>`_.

--------------------------------------------------------------------------------

We are currently in the process of rewriting the code and making it publicly
available, and the latest version (with limited functionality) can be found on
`GitHub <https://github.com/MRChemSoft/mrchem>`_.

Features in MRChem-0.2.0:
-------------------------

* Wave functions:
    + Kohn-Sham DFT
        - Spin-polarized
        - Spin-unpolarized
        - LDA, GGA and hybrid functionals
    + Hartree-Fock
        - Restricted closed-shell
        - Unrestricted
    + Explicit external fields
        - Electric field
* Properties:
    + Ground state energy
    + Dipole moment
* Parallel implementation:
    + Shared memory (OpenMP): ~20 cores
    + Distributed memory (MPI): ~100 procs
    + Hybrid scheme (MPI + OpenMP): ~1000 cores
* Current size limitations:
    + ~200 orbitals on ~50 medium-memory (64GB) compute nodes
    + ~100 orbitals on a single high-memory (1TB) compute node

Upcoming features:
------------------

* Wave functions:
    + Meta-GGAs
* Properties:
    + Quadrupole moment
    + Polarizability
    + Hyperpolarizability
    + Optical rotation
    + Magnetizability
    + NMR shielding constant
    + Spin-spin coupling constant
    + Hyperfine coupling constant
    + Magnetically induced currents
    + Geometry optimization


.. toctree::
   :maxdepth: 2

   installation
   mrchem_manual
   input_reference
   code_reference/classes-and-functions
