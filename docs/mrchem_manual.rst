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
        keyword_1 = 1
        keyword_2 = 3.14
        keyword_3 = [1, 2, 3]
        keyword_4 = "foo"
        keyword_5 = True
     }


    
The main input section contain two important keywords that specify the
polynomial order of the multiwavelet basis set, and the relative precision that
will be guaranteed in the calculation. The main input section is not specified
by name, just write the keywords directly, e.g

.. code-block:: bash

    order = 7 
    rel_prec = 1.0e-5

Increased precision requires higher polynomial order (use e.g order = 5 for
rel_prec = 1.0e-3, and order = 13 for rel_prec = 1.0e-9, and interpolate in
between).


WaveFunction
------------

Here we give the wavefunction method (HF or DFT) and whether we run
spin restricted (alpha and beta spins are forced to occupy the same spatial 
orbitals) or not. When running DFT we must also specify the functional to be 
used in a separate DFT section (for HF this section should be omitted)

.. code-block:: bash

    WaveFunction {
        method = <wavefunction_method>
        restricted = True
    }

    DFT {
        spin = False
        exact_exchange = 0.0
        $functionals
        <func1>     <coef1>
        <func2>     <coef2>
        $end
    }

You can specify as many functionals as you want, and they will be added on top
of each other with the given coefficient. For hybrid functionals you must 
specify the amount of exact Hartree-Fock
exchange that should be used (0.2 for B3LYP and 0.25 for PBE0 etc.). Option to
use spin-density functional theory (for open-shell systems).

World
-----

This section will specify the computational domain

.. code-block:: bash

     World {
        scale = -5
        world_origin = [0.0, 0.0, 0.0]
        gauge_origin = [0.0, 0.0, 0.0]
    }

where scale gives the size of the domain as :math:`2^{-scale}`. This will 
be symmetric around zero, so the above will define a computational domain of 
:math:`[-16,16]^3`. The computational World should be large enough so that the 
electron density vanishes at the boundaries. The world origin can be used to 
translate the domain away from the symmetric box around the origin. The gauge 
origin can also be specified (relevant for magnetic properties).

Molecule
--------

This input section specifies the geometry, charge and spin multiplicity of the 
molecule, e.g. for the water molecule
   
.. code-block:: bash

    Molecule {
        charge = 0
        multiplicity = 1
        angstrom = False
        $coords
        O   0.0000     0.0000     0.0000
        H   0.0000     1.4375     1.1500
        H   0.0000    -1.4375     1.1500
        $end
    }

LSDalton
--------

MRChem can use the LSDalton program to obtain an initial guess for the orbitals,
using a small Gaussian basis set, which is specified in this section
    
.. code-block:: bash

    LSDalton {
        run = True
        method = <wavefunction_method>
        basis = <basis_set>
    }
Currently, only HF (Hartree-Fock) and LDA can be used as 
``<wavefunction_method>``, and the ``<basis_set>`` must be quite small, as 
MRChem can only read s- p- and (uncontracted) d-functions. Option to run 
LSDalton or not.

Properties
----------

Specify which properties to compute. Currently the following are available

.. code-block:: bash

    Properties {
        ground_state = True
        dipole_moment = True
        quadrupole_moment = True
        polarizability = True
        magnetizability = True
        optrot_electric = True
        optrot_magnetic = True
        nmr_shielding = True
        nmr_nuclei = [<nuc1>, <nuc2>, ...]
        frequencies = [<omega1>, <omega2>, ...]
    }

Optical rotation can be computed using either electric or magnetic response.
When computing NMR shielding constants you can specify which atom(s) you want to
compute (the default is [-1] which computes for all nuclei). Here you also
specify the frequencies of the perturbing laser field (for dynamic properties),
default frequency is 0.0 (static field). Several properties can be computed at
once, and magnetic properties are always static, while the frequencies applies 
to polarizability and optical rotation.

SCF
---

Specify the parameters for the SCF optimization of the ground state wave 
function

.. code-block:: bash
 
    SCF {
        property_thrs = 1.0e-4
        orbital_thrs = 1.0e-3
        history = 4
        rotation = 50
        localize = False
        write_orbitals = False
        initial_guess = <guess>
    }

Here we specify the convergence thresholds for the orbitals and the property 
(total energy). The rotation keyword says how often the Fock matrix should be
diagonalized/localized. Option to use localized molecular orbitals, and whether
the final orbitals should be written to disk. You can set the length of the
iterative history that is used in the KAIN accelerator. You also need to specify 
which initial guess to use, "gto" means start with an LSDalton calculation, "mw" 
means that we start from a previous MRChem calculation (final orbitals must have 
been written).

Response
--------

Specify the parameters for the SCF optimization of the linear response wave 
function. This section must be included if any linear response properties 
are computed.

.. code-block:: bash
   
    Response {
        property_thrs = 1.0e-4
        orbital_thrs = 1.0e-3
        history = 6
        localize = False
    }

Convergence thresholds are specified for the molecular propery and the perturbed
orbitals. Option to use localized orbitals in the response solver (independent
of the localize option for the ground state calculation). You can also set the 
length of the iterative history that is used in the KAIN accelerator in the 
response solver. 

