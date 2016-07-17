===================
MRChem: User manual
===================

---------------------
The mrchem input file
---------------------

The input file is organized in sections and keywords that can be of different
type

.. code-block:: bash

     Section {
        keyword_1 = 1                       # int
        keyword_2 = 3.14                    # double
        keyword_3 = [1, 2, 3]               # int array
        keyword_4 = foo                     # string
        keyword_5 = true                    # bool
     }

Top section
-----------

The main input section contain two important keywords that specify the
polynomial order :math:`k` of the multiwavelet basis set, and the relative
precision :math:`\epsilon_{rel}` that will be guaranteed in the calculation.
The main input section is not specified by name, just write the keywords
directly, e.g

.. code-block:: bash

    order = 7                               # Polynomial order of MW basis
    rel_prec = 1.0e-5                       # Overall relative precision

Increased precision requires higher polynomial order (use e.g :math:`k = 5`
for :math:`\epsilon_{rel} = 10^{-3}`, and :math:`k = 13` for
:math:`\epsilon_{rel} = 10^{-9}`, and interpolate in between). If the ``order``
keyword is left out it will be set automatically according to

.. math:: k=-1.5*log_{10}(\epsilon_{rel})

The relative precision sets an upper limit for the number of correct digits
you are expected to get out of the computation (note that
:math:`\epsilon_{rel}=10^{-6}` yields :math:`\mu` Ha accuracy for the hydrogen
molecule, but only mHa accuracy for benzene). It is also possible to specify
an `absolute` precision for the molecular energy by replacing ``rel_prec``
with ``abs_prec``. This will provide e.g. mHa precision in the total energy
regardless of the molecular size (this might get `very` expensive for large
systems). In this case the magnitude of the energy is estimated as

.. math:: E = \sum_i^{nuc} Z_i^{5/2} 

and the relative precision is set as

.. math:: \epsilon_{rel} = \frac{\epsilon_{abs}}{E}

With the ``abs_prec`` keyword one can also use kcal/mol or kJ/mol as energy
unit instead of Hartree by setting the ``energy_unit`` keyword. The following
will set ``rel_prec`` sufficiently high so that the energy can be computed
within a kcal/mol

.. code-block:: bash

    abs_prec = 1.0
    energy_unit = kcal

Note that ``order`` and ``rel_prec`` are the fundamental input parameters that
are finally passed to the program, so if any of these are set explicitly in the
input file, they will always have precedence.

A final keyword in the main input section that is usually used for debugging is
the ``printlevel`` (should be zero for production calculations).

World
-----

This section will specify the computational domain (defaults shown)

.. code-block:: bash

     World {
        scale = 0                           # Size of each root box 2^{-scale}
        boxes = [ 1, 1, 1 ]                 # Number of root boxes
        corner = [ 0, 0, 0]                 # Translation of first root box
        gauge_origin = [0.0, 0.0, 0.0]      # Origin used in molecular properties
    }

The scale and translation of the boxes are absolute, which means that the only
way to get a symmetric world around the origin is to use two root boxes in each
direction and set corner at -1 (if this does not fit well with your molecular
geometry, use a larger box or translate your molecular coordinates). The
computational world should be large enough so
that the electron density vanishes at the boundaries. The ``gauge_origin`` can
also be specified (relevant for molecular properties), otherwise it will be the
molecular center of mass. The default computational domain displayed above
corresponds to the unit cube (in bohr).

Molecule
--------

This input section specifies the geometry, charge and spin multiplicity of the
molecule, e.g. for water (coords must be specified, otherwise
defaults are shown)

.. code-block:: bash

    Molecule {
        charge = 0                          # total charge of molecule
        multiplicity = 1                    # spin multiplicity
        angstrom = false                    # geometry given in angstrom
        $coords
        O   0.0000     0.0000     0.0000
        H   0.0000     1.4375     1.1500
        H   0.0000    -1.4375     1.1500
        $end
    }

WaveFunction
------------

Here we give the wavefunction method and whether we run spin restricted (alpha
and beta spins are forced to occupy the same spatial orbitals) or not (method
must be specified, otherwise defaults are shown) 

.. code-block:: bash

    WaveFunction {
        method = <wavefunction_method>      # Core, Hartree, HF or DFT
        restricted = true                   # Spin restricted/unrestricted
    }

There are currently four methods available: Core Hamiltonian, Hartree,
Hartree-Fock (HF) and Density Functional Theory (DFT). When running DFT the
functional(s) must be specified in a separate DFT section (see below)

DFT
---
 
This section specifies the exchange-correlation functional used in DFT. For HF
this section should be omitted (functionals must be specified, otherwise
defaults are shown)

.. code-block:: bash

    DFT {
        spin_polarized = false              # Use spin-polarized functionals
        exact_exchange = 0.0                # Amount of exact HF exchange
        density_cutoff = 0.0                # Cutoff to set XC potential to zero
        $functionals
        <func1>     1.0                     # Functional name and coefficient
        <func2>     1.0
        $end
    }

You can specify as many functionals as you want, and they will be added on top
of each other with the given coefficient. Both exchange and correlation
functinals must be set explicitly e.g. ``SLATERX`` and ``VWN5C`` for the
standard LDA functional. For hybrid functionals you must
specify the amount of exact Hartree-Fock exchange that should be used (0.2 for
B3LYP and 0.25 for PBE0 etc.). Option to use spin-polarized functionals (for
open-shell systems). XC functionals are provided by the `XCFun 
<https://github.com/dftlibs/xcfun>`_ library.

Properties
----------

Specify which properties to compute. Currently the following are available
(defaults shown)

.. code-block:: bash

    Properties {
        total_energy = false                # Compute total energy
        dipole_moment = false               # Compute dipole moment
    }

SCF
---

Specify the parameters for the SCF optimization of the ground state wave
function (defaults shown)

.. code-block:: bash

    SCF {
        run = true                          # Run SCF optimization
        orbital_thrs = 1.0                  # Convergence threshold orbitals
        property_thrs = 1.0                 # Convergence threshold energy
        orbital_prec = [1.0e-4, -1.0]       # Initial and final relative precision in SCF
        history = 0                         # Length of KAIN iterative subspace
        rotation = 0                        # Iterations between each localization/diagonalization
        max_iter = -1                       # Maximum number of SCF iterations
        localize = false                    # Use localized or canonical orbitals
        write_orbitals = false              # Write final orbitals to disk
        initial_guess = none                # Type of inital guess (none, gto, mw)
    }

With ``run=false`` no SCF optimization is performed, and the requested molecular
properties are computed directly from the initial guess wave function.

Here we specify the convergence thresholds for the orbitals
(:math:`\|\Delta \phi_i \|`) and the property (:math:`\Delta E`).
Notice that these corresponds to two separate optimizations: first the orbitals
are converged within ``orbital_thrs`` using a KAIN optimization that yields
energy accuracy that is linear in the orbital errors. Then a separate algorithm
that is quadratic in the orbital error is used (one that avoids the use of the
kinetic energy operator) to converge the energy within ``property_thrs``. This
algorithm does not use KAIN, and is thus not efficient for converging the 
orbitals. If one is `not` interested in the total energy to high precision, one
can avoid the second optimization by setting ``property_thrs = -1.0``, and
simply converge the orbitals to the desired precision. This should yield similar
accuracy for all properties. Notice also that even if the energy error is
quadratic using the second algorithm, it is still limited by the overall
precision ``rel_prec``. For instance, the following should yield 5 digits in
the total energy and three digits in other properties

.. code-block:: bash

    rel_prec = 1.0e-5

    SCF {
        orbital_thrs = 1.0e-3
        property_thrs = 1.0e-6
    }

To get 5 digits in all properties, choose the following (always keep at least
a factor of 10 between ``rel_prec`` and ``orbital_thrs`` to avoid numerical
instabilities)

.. code-block:: bash

    rel_prec = 1.0e-6

    SCF {
        orbital_thrs = 1.0e-5
        property_thrs = -1.0
    }


If these thresholds are not set explicitly they will be computed from the top
level ``rel_prec`` as

.. math:: \Delta E < \epsilon_{rel}/10
.. math:: \|\Delta \phi_i \| < \sqrt{\epsilon_{rel}/10}

The ``orbital_prec=[init,final]`` keyword controls the dynamic precision used
in the SCF iterations. To improve efficiency, the first iterations are done
with reduced precision, staring at ``init`` and gradually increased
to ``final``. The initial precision should not be set lower than
``init=1.0e-3``, and the final precision should not exceed the top level
``rel_prec``. Negative values sets them equal to ``rel_prec``. 

The ``rotation`` keyword says how often the Fock matrix should be
diagonalized/localized (for iterations in between a Löwdin orthonormalization
is used). Option to use Foster-Boys localization or Fock matrix diagonalization
in these rotations. Note that the KAIN history is cleared every time this
rotation is employed to avoid mixing of orbtials in the history, so
``rotation=1`` effectively cancels the KAIN accelerator. The default
``rotation=0`` will localize/diagonalize the first two iterations and then
perform Lövdin orthonormalizations (this is usually the way to go).

The ``history`` keyword sets the size of the iterative subspace that is used
in the KAIN accelerator for the orbital optimization.

You also need to specify which initial guess to use, "none" means starting
from hydrogen solutions, "gto" means start with a Gaussian type calculation
(basis and MO matrix input files must be provided) and "mw" means starting
from a previous MRChem calculation (compatible orbitals must have been written
to disk using the ``write_orbitals`` keyword).

