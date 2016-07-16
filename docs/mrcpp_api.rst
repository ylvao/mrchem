====================================
MRCPP: Application Program Interface
====================================

------------
Introduction
------------

The main features of MRCPP are the numerical multiwavelet (MW) representations
of functions and operators. Two integral convolution operators are implemented
(the Poisson and Helmholtz operators), as well as the
partial derivative and arithmetic operators. In addition
to the numerical representations there are a limited number of analytic
functions that are usually used as starting point for the numerical
computations.

------------------
Analytic functions
------------------

The general way of defining an analytic function is to use lambdas
(or std::function), which provide lightweight functions that can be used on the
fly. However, some analytic functions, like Gaussians, are of special
importance and have been explicitly implemented with additional functionality.

In order to be accepted by the ``MWProjector`` (see below) the lambda needs to
have the following signature

.. code-block:: cpp

    auto f = [] (const double *r) -> double;

e.i. it must take a ``double`` pointer defining a set of Cartesian coordinates,
and return a ``double``. For instance, the electrostatic potential from a point
nuclear charge :math:`Z` (in atomic units) is

.. math:: f(r) = \frac{Z}{r}

which can be written as the lambda function

.. code-block:: cpp

    double Z = 1.0;                                         // Hydrogen nuclear charge
    auto f = [Z] (const double *r) -> double {
        double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return Z/R;
    }

Note that the function signature must be exactly as given above, which means
that any additional arguments (such as :math:`Z` in this case) must be given in
the capture list (square brackets), see e.g. `cppreference.com 
<http://en.cppreference.com/w/cpp/language/lambda>`_ for more
details on lambda functions and how to use the capture list.

------------------------------
MultiResolution Analysis (MRA)
------------------------------

In order to combine different functions and operators in mathematical
operations they must be compatible, that is, they must be
defined on the same computational domain and constructed using the same
polynomial basis (order and type). This information constitutes an MRA,
which needs to be defined and passed as argument to all function and operator
constructors, and only functions and operators with compatible MRAs can be
combined in subsequent calculations. An MRA is defined in two steps, first the
computational domain is given by a ``BoundingBox`` (D is the dimension)

.. code-block:: cpp

    int n;                                                  // Root scale defines box size 2^{-n}
    int l[D];                                               // Translation of first box
    int nb[D];                                              // Number of boxes
    NodeIndex<D> idx(n, l);                                 // Index defining the first box
    BoundingBox<D> world(idx, nb);

which is combined with a ``ScalingBasis`` to give an MRA

.. code-block:: cpp

    int k;                                                  // Polynomial order
    ScalingBasis basis(k);                                  // Legendre or Interpolating basis
    MultiResolutionAnalysis<D> MRA(world, basis);

Two types of ``ScalingBasis`` are supported (``LegendreBasis`` and
``InterpolatingBasis``), and they are both available at orders
:math:`k=1,2,\dots,40` (note that some operators are constructed using
intermediates of order :math:`2k`, so in that case the maximum order available
is :math:`k=20`).

------------------------
Function representations
------------------------

MW representations of functions are called ``FunctionTrees``, and are in
principle available in any dimension using the template parameter D (in
practice D=1,2,3). There are three different ways of constructing MW functions
(computing the expansion coefficients in the MW basis)

* Projection of analytic function
* Arithmetic operations
* Application of MW operator

and these will be described in the following sections.

The interface for constructing MW representations has a dual focus: on the one
hand we want a simple, intuitive way of producing adaptive numerical
approximations with guaranteed precision that does not require detailed
knowledge of the internals of the MW code and with a minimal number of
parameters that have to be set. On
the other hand we want the possibility for more detailed control of the
construction and refinement of the numerical grid where such control is
possible and even necessary. In the latter case it is important to be able to
reuse the existing grids in e.g. iterative algorithms without excessive
allocation/deallocation of memory.

FunctionTree
------------

Constructing a full grown ``FunctionTree`` involves a number of steps, including
setting up a memory allocator, constructing root nodes according to the given
MRA, building a tree structure and computing MW coefficients. For this reason
the ``FunctionTree`` constructor is made protected, and all construction is done
indirectly through ``TreeBuilder`` objects

.. code-block:: cpp

    FunctionTree<D> *tree = builder();

where ``builder`` is any of the ``TreeBuilders`` presented below which may or
may not take any arguments for the construction. Details on how the tree
structure is built and how the MW coefficients are computed are specified in
each particular ``TreeBuilder``. Since ``FunctionTrees`` always appear as
pointers, we will in the following use pointer notation for all trees.

Integrals are computed very efficiently in the orthonormal MW basis, and among
the important methods of the ``FunctionTree`` are estimating the error in the
representation (based on the wavelet norm), obtaining the squared
:math:`L^2`-norm of the function, as well as its integral and dot product with
another ``FunctionTree`` (both over the full computational domain)

.. code-block:: cpp

    double error = f_tree->estimateError();
    double sq_norm = f_tree->getSquareNorm();
    double integral = f_tree->integrate();
    double dot_prod = f_tree->dot(*g_tree);

FunctionTreeVector
------------------

The ``FunctionTreeVector`` is a convenience class for a collection of
``FunctionTrees`` which basically consists of two STL vectors, one containing
pointers to ``FunctionTrees`` and one with corresponding numerical coefficients.
Elements can be appended to the vector

.. code-block:: cpp
    
    FunctionTreeVector<D> tree_vec;
    tree_vec.push_back(2.0, tree_a);                            // Push back pointer to FunctionTree
    tree_vec.push_back(tree_b);                                 // Push back pointer to FunctionTree
    tree_vec.clear(false);                                      // Bool argument for tree destruction

where ``tree_b`` will be appended with a default coefficient of 1.0. Clearing
the vector means removing all its elements, and the ``bool`` argument tells if
the elements should be properly deallocated (default ``false``).

TreeBuilder
-----------

This is the class that is responsible for the construction of
``FunctionTrees``, which involves allocating memory, growing a tree structure
and calculating MW coefficients. The ``TreeBuilder`` has two important members:
a ``TreeCalculator`` that defines how the MW coefficients are computed, and a
``TreeAdaptor`` that defines how the tree structure is grown. There are four
different ways of computing MW coefficients (projection, addition,
multiplication and operator application), and we have the corresponding
``TreeBuilders`` (the MW prefix indicates that they compute MW coefficients)

* MWProjector
* MWAdder
* MWMultiplier
* MWOperator

Each of these is a specialization of the ``TreeBuilder`` class that differs in
the type of ``TreeCalculator``. They all contain a ``TreeAdaptor`` that
controls the accuracy of the function representations they build.
All ``TreeBuilders`` have the same fundamental building algorithm:

1. Start with an initial guess for the grid
2. Use the ``TreeCalculator`` to compute the output function on the current grid
3. Use the ``TreeAdaptor`` to refine the grid where needed
4. Iterate points 2 and 3 until the grid is converged

All builders are constructed using an MRA and one or more precision
parameters as arguments

.. code-block:: cpp

    double prec;
    MultiResolutionAnalysis<D> MRA;
    TreeBuilder<D> builder(MRA, prec);

Where the ``MRA`` defines the computational domain and type of MW basis and
``prec`` defines the relative precision used by the ``TreeAdaptor`` for the
contruction of the output tree structure. The ``MWOperator`` might take a
second precision parameter related to the construction of the operator itself.
The precision parameters have negative defaults, which means that no refinement
is made beyond the initial grid.

The interface for the ``TreeBuilders`` is mainly the ``operator()``, which
comes in two versions

.. code-block:: cpp

    out = builder(inp);
    builder(out, inp, max_iter);

where the former is a constructor that returns a pointer to a new
``FunctionTree``, while the latter will work on an already existing tree. The
main difference between the two is the choice of initial grid: the former will
always start at the root nodes; the latter will start at whatever grid is
already present in the output tree structure. The latter allows for more
detailed control for the user, however, the grids needs to be prepared in
advance, either using a ``GridGenerator`` to construct an empty grid or a
``GridCleaner`` to clear an existing ``FunctionTree`` (see advanced
initialization below). The final argument ``max_iter`` is used to stop the
building algorithm after a certain number of iterations beyond the initial
grid, even if the accuracy criterion is not met. This will of course not
guarantee the accuracy of the representation, but is useful in certain
situations, e.g. when you want to work on fixed grid sizes.

MWProjector
-----------

The ``MWProjector`` takes an analytic D-dimensional scalar function (which can
be defined as a lambda function or one of the explicitly implemented sub-classes
of the ``RepresentableFunction`` base class) and projects it
onto the MRA to the given precision. E.g. a unit charge Gaussian is
projected in the following way (the MRA must be initialized as above)

.. code-block:: cpp

    double beta = 10.0;                                     // Gaussian exponent
    double alpha = pow(beta/pi, 3.0/2.0);                   // Unit charge coefficient
    auto f = [alpha, beta] (const double *r) -> double {
        double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return alpha*exp(-beta*R*R);
    }

    double prec = 1.0e-5;
    MWProjector<3> Q(MRA, prec);
    FunctionTree<3> *f_tree = Q(f);

The projector will construct an initial grid containing only the root nodes of
the MRA and follow the builder algorithm (see above) to adaptively construct the
grid necessary to represent the function to the given precision (based on the
wavelet norm of the representation). Note that with a negative precision (which
is the default) the grid will not be refined beyond the initial grid, which
contains only root nodes in this case.

Arithmetic operations
---------------------

Arithmetic operations in the MW representation are performed using the
``FunctionTreeVector``, and the general sum :math:`g = \sum_i c_i f_i(x)` and
product :math:`h = \prod_i c_i f_i(x)` are done in the following way

.. code-block:: cpp

    FunctionTreeVector<D> inp_vec;
    inp_vec.push_back(c_1, f_tree_1);
    inp_vec.push_back(c_2, f_tree_2);
    inp_vec.push_back(c_3, f_tree_3);

    MWAdder<D> add(MRA, prec);
    FunctionTree<D> *g_tree = add(inp_vec);

    MWMultiplier<D> mult(MRA, prec);
    FunctionTree<D> *h_tree = mult(inp_vec);

The default initial grid is again only the root nodes, and a positive ``prec``
is required to build an adaptive tree structure for the result. The special
case of adding/multiplying two functions can be done directly without
initializing a ``FunctionTreeVector``

.. code-block:: cpp

    MWAdder<D> add(MRA, prec);
    FunctionTree<D> *g_tree = add(c_1, *f_tree_1, c_2, *f_tree_2);

    MWMultiplier<D> mult(MRA, prec);
    FunctionTree<D> *h_tree = mult(c_1*c_2, *f_tree_1, *f_tree_2);

Operator application
--------------------

Two types of operators are currently implemented in MRCPP:
the Cartesian derivative

.. math:: g(x) = \partial_x f(x)

and integral convolution

.. math::  g(r) = \int G(r-r') f(r') dr'

The syntax for construction and application follows closely the other
``TreeBuilders`` presented above.

### Derivative operator

The derivative operator is initialized with two parameters :math:`a` and
:math:`b` accounting for the boundary conditions between adjacent nodes 
(see `Alpert etal.
<http://www.sciencedirect.com/science/article/pii/S0021999102971603>`_,
in practice :math:`a=b=0` is the best choice),
and the Cartesian direction of application must be specified in advance
(otherwise it is applied in all directions consecutively, corresponding in 3D
to the :math:`\partial_{xyz}` operator)

.. code-block:: cpp

    double prec;                                            // Precision of operator application
    double a = 0.0, b = 0.0;                                // Boundary conditions for operator
    DerivativeOperator<3> D(MRA, prec, a, b);
    D.setApplyDir(1);                                       // Differentiate in y direction
    FunctionTree<3> *g_tree = D(*f_tree);                   // Build result adaptively

As for all ``TreeBuilders``, this operator will start at the root nodes and
build adaptively according to ``prec``. The derivative is usually applied
directly on the grid of the input function, without further refinement (see
advanced initialization below)

.. code-block:: cpp

    GridGenerator<3> G(MRA);                                // TreeBuilder that copy grids
    FunctionTree<3> *g_tree = G(*f_tree);                   // Copy grid from density function

    DerivativeOperator<3> D(MRA);                           // Default parameters prec = -1, a=b=0
    D.setApplyDir(1);                                       // Differentiate in y direction    
    D(*g_tree, *f_tree, 0);                                 // Compute derivative on given grid

### Poisson operator

The electrostatic potential :math:`g` arising from a charge distribution
:math:`f` are related through the Poisson equation

.. math:: -\nabla^2 g(r) = f(r)

This equation can be solved with respect to the potential by inverting the
differential operator into the Green's function integral convolution operator

.. math:: g(r) =  \int \frac{1}{4\pi\|r-r'\|} f(r') dr'

This operator is available in the MW representation, and can be solved with
arbitrary (finite) precision in linear complexity with respect to system size.
Given an arbitrary charge dirtribution ``f_tree`` in the MW representation, the
potential is computed in the following way

.. code-block:: cpp

    double apply_prec;                                      // Precision defining the operator application
    double build_prec;                                      // Precision defining the operator construction
    PoissonOperator P(MRA, apply_prec, build_prec);
    FunctionTree<3> *g_tree = P(*f_tree);                   // Apply operator adaptively

The Coulomb self-interaction energy can now be computed as the dot product

.. code-block:: cpp

    double E = g_tree->dot(*f_tree);

### Helmholtz operator

The Helmholtz operator is a generalization of the Poisson operator and is given
as the integral convolution

.. math:: g(r) =  \int \frac{e^{-\mu\|r-r'\|}}{4\pi\|r-r'\|} f(r') dr'

The operator is the inverse of the shifted Laplacian

.. math:: \big[-\nabla^2 + \mu^2 \big] g(r) = f(r)

and appears e.g. when solving the SCF equations. The construction and
application is similar to the Poisson operator, with an extra argument for the
:math:`\mu` parameter

.. code-block:: cpp

    double mu;                                              // Must be a positive real number
    double apply_prec;                                      // Precision defining the operator application
    double build_prec;                                      // Precision defining the operator construction
    HelmholtzOperator H(MRA, mu, apply_prec, build_prec);
    FunctionTree<3> *g_tree = H(*f_tree);                   // Apply operator adaptively

-----------------------
Advanced initialization
-----------------------

The ``TreeBuilders``, as presented above, have a clear and limited interface,
but there are two important drawbacks: every operation require the construction
of a new ``FunctionTree`` from scratch (including extensive memory allocation),
and the tree building algorithm always starts from a root node initial grid.
In many practical applications however (e.g. iterative algorithms), we are
recalculating the same functions over and over, where the requirements on the
numerical grids change only little between each iteration. In such situations it
will be beneficial to be able to reuse the existing grids without reallocating
the memory and recomputing all the coarse scale nodes in the building process.
For this purpose we have the following additional ``TreeBuilders`` (the Grid
prefix indicates that they do not compute MW coefficients):

* GridGenerator
* GridCleaner

where the former constructs empty grids from scratch and the latter clears the
MW coefficients on existing ``FunctionTrees``. The end result is in both cases
an empty tree skeleton with no MW coefficients (undefined function).

GridGenerator
-------------

Sometimes it is useful to construct an empty grid based on some available
information of the function that is about to be represented. This can be e.g.
that you want to copy the grid of an existing ``FunctionTree`` or that an
analytic function has more or less known grid requirements (like Gaussians).
Sometimes it is even necessary to force the grid refinement beyond the coarsest
scales in order for the ``TreeAdaptor`` to detect a wavelet "signal" that allows
it to do its job properly (this happens for narrow Gaussians where non of the
initial quadrature points hits a function value significantly different from
zero). In such cases we use a ``GridGenerator`` to build an initial tree
structure.

A special case of the ``GridGenerator`` (with no argument) corresponds to the
default constructor of the ``FunctionTree``

.. code-block:: cpp

    GridGenerator<D> G(MRA);
    FunctionTree<D> *f_tree = G();

which will construct a new ``FunctionTree`` with empty nodes (undefined
function with no MW coefficients), containing only the root nodes in the given
MRA. Passing an analytic function as argument to the generator will use a
``TreeAdaptor`` to build a grid based on some predefined information of the
function (if there are any, otherwise it is identical to the default
constructor)

.. code-block:: cpp

    FunctionTree<D> *f_tree = G(f_func);

The lambda analytic functions do `not` provide such information, this must be
explicitly implemented as a ``RepresentableFunction`` sub-class (see MRCPP
programmer's guide for details). Passing a ``FunctionTree`` to the generator
will copy its grid

.. code-block:: cpp

    FunctionTree<D> *f_tree = G(*g_tree);

Both of these will produce a skeleton ``FunctionTree`` with empty nodes. In
order to define a function in the new tree it is passed as the first argument
to the regular ``TreeBuilders`` presented above, e.g for projection

.. code-block:: cpp

    GridGenerator<D> G(MRA);
    MWProjector<D> Q(MRA, prec);
    FunctionTree<D> *f_tree = G(f_func);                    // Empty grid from analytic function
    Q(*f_tree, f_func, max_iter);                           // Starts projecting from given grid

This will first produce an empty grid suited for representing the analytic
function ``f_func`` (this is meant as a way to make sure that the projection
starts on a grid where the function is actually visible, as for very narrow
Gaussians, it's `not` meant to be a good approximation of the final grid) and
then perform the projection on the given numerical grid. With a negative
``prec`` (or ``max_iter = 0``) the projection will be performed strictly on the
given initial grid, with no further refinements.

Similar notation applies for all ``TreeBuilders``: if a ``FunctionTree`` is
given as the first argument, it will not construct a new tree but perform the
operation on the one given. E.g. the grid copy can be done in two steps as

.. code-block:: cpp

    f_tree = G();                                           // Construct empty grid of root nodes
    G(*f_tree, *g_tree);                                    // Extend grid with missing nodes relative to g

Actually, the effect of the ``GridGenerator`` is to *extend* the existing grid
with any missing nodes relative to the input. This means that we can build the
union of two grids by successive application of the generator

.. code-block:: cpp

    f_tree = G();                                           // Construct empty grid of root nodes
    G(*f_tree, *g_tree);                                    // Extend f with missing nodes relative to g
    G(*f_tree, *h_tree);                                    // Extend f with missing nodes relative to h

and one can make the grids of two functions equal to their union

.. code-block:: cpp

    G(*f_tree, *g_tree);                                    // Extend f with missing nodes relative to g
    G(*g_tree, *f_tree);                                    // Extend g with missing nodes relative to f

The union grid of several trees can be constructed in one go using a
``FunctionTreeVector``

.. code-block:: cpp

    FunctionTreeVector<D> inp_vec;
    inp_vec.push_back(tree_1);
    inp_vec.push_back(tree_2);
    inp_vec.push_back(tree_3);
    FunctionTree<D> *f_tree = G(inp_vec);

Addition of two functions is usually done on their union grid

.. code-block:: cpp

    MWAdder<D> add(MRA);                                    // Default negative precision
    GridGenerator<D> G(MRA);

    FunctionTree<D> *f_tree = G();                          // Construct empty root grid
    G(*f_tree, *g_tree);                                    // Copy grid of g
    G(*f_tree, *h_tree);                                    // Copy grid of h
    add(*f_tree, 1.0, *g_tree, 1.0, *h_tree);               // Add functions on union grid

Note that in the case of addition there is no extra information to be gained
by going beyond the finest refinement levels of the input functions, so the
union grid summation is simply the best you can do, and adding a positive
``prec`` will not make a difference. There are situations where you want to
use a `smaller` grid, though, e.g. when performing a unitary transformation
among a set of ``FunctionTrees``. In this case you usually don't want to
construct `all` the output functions on the union grid of `all` the input
functions, and this can be done by adding the functions adaptively starting
from root nodes.

For multiplications, however, there might be a loss of accuracy if
the product is restricted to the union grid. The reason for this is that the
product will contain signals of higher frequency than each of the input
functions, which require a higher grid refinement for accurate representation.
By specifying a positive ``prec`` you will allow the grid to adapt to the higher
frequencies, but it is usually a good idea to restrict to one extra refinement
level beyond the union grid (by setting ``max_iter=1``) as the grids are not
guaranteed to converge for such local operations (like arithmetics, derivatives
and function mappings)

.. code-block:: cpp

    double prec;
    MWMultiplier<D> mult(MRA, prec);
    GridGenerator<D> G(MRA);

    FunctionTree<D> *f_tree = G();                          // Construct empty root grid
    G(*f_tree, *g_tree);                                    // Copy grid of g
    G(*f_tree, *h_tree);                                    // Copy grid of h
    mult(*f_tree, 1.0, *g_tree, *h_tree, 1);                // Allow 1 extra refinement

If you have a summation over several functions but want to perform the
addition on the grid given by the `first` input function, you first copy the
wanted grid and then perform the operation on that grid

.. code-block:: cpp

    FunctionTreeVector<D> inp_vec;
    inp_vec.push_back(coef_1, tree_1);
    inp_vec.push_back(coef_2, tree_2);
    inp_vec.push_back(coef_3, tree_3);

    MWAdder<D> add(MRA);                                    // Default negative precision
    GridGenerator<D> G(MRA);

    FunctionTree<D> *f_tree = G(tree_1);                    // Copy grid of first input function
    add(*f_tree, inp_vec);                                  // Perform addition on given grid

Here you can of course also add a positive ``prec`` to the ``MWAdder``
and the resulting function will be built adaptively starting from the given
initial grid.

GridCleaner
-----------

Given a ``FunctionTree`` that is a valid function representation we can clear
its MW expansion coefficients (while keeping the grid refinement) with the
``GridCleaner`` (unlike the other ``TreeBuilders``, the ``GridCleaner`` will
not return a ``FunctionTree`` pointer, as it would always be the same as the
argument)

.. code-block:: cpp

    GridCleaner<D> C(MRA);
    C(f_tree);

This action will leave the ``FunctionTree`` in the same state as the
``GridGenerator`` (uninitialized function), and its coefficients can now be
re-computed.

In certain situations it might be desireable to separate the actions of the
``TreeCalculator`` and the ``TreeAdaptor``. For this we can combine the
``GridCleaner`` with the ``TreeAdaptor``, which will adaptively refine the
grid one level (based on the wavelet norm and the given precision) `before` it
is cleared

.. code-block:: cpp

    double prec;
    GridCleaner<D> C(MRA, prec);
    C(f_tree);

One example where this might be
useful is in iterative algorithms where you want to fix the grid size for
all calculations within one cycle and then relax the grid in the end in
preparation for the next iteration. The following is equivalent to the adaptive
projection above (the cleaner returns the number of new nodes that were created
in the process)

.. code-block:: cpp

    double prec;
    GridCleaner<D> C(MRA, prec);                            // The precision parameter is passed as
    MWProjector<D> Q(MRA);                                  // argument to the cleaner, not the projector

    int n_nodes = 1;
    while (n_nodes > 0) {
        Q(*f_tree, f);                                      // Project f on given grid
        n_nodes = C(*f_tree);                               // Refine grid and clear coefficients
    }
    Q(*f_tree, f);                                          // Project f on final converged grid


