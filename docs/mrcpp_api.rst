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

Note that all these analytic functions needs to be projected onto the MW basis
before they can be used in numerical computations (see below).


GaussFunc
---------

A very important analytic function is the Cartesian Gaussian

.. math:: f(r) = \alpha (x-x_0)^a (y-y_0)^b (z-z_0)^c e^{-\beta \|r-r_0\|^2}

which is initialized in the following way (in 3D)

.. code-block:: cpp

    double alpha, beta;
    double pos[3] = {x_0, y_0, z_0};
    int pow[3] = {a, b, c};
    GaussFunc<3> f_func(beta, alpha, pos, pow);

GaussPoly
---------

GaussPoly is a generalization of the GaussFunc, where there is an arbitrary
polynomial in front of the exponential

.. math:: f(r) = \alpha P(r-r_0) e^{-\beta \|r-r_0\|^2}

GaussExp
--------

A GaussExp is a collection of Gaussian in the form

.. math:: G(r) = \sum_i c_i g_i(r)

where :math:`g_i` can be either GaussFunc or GaussPoly

.. math:: g_i(r) =  \alpha_i P_i(r-r_i)e^{-\beta_i\|r-r_i\|^2}

PositionFunction
----------------

There is also a very simple analytic function for the Cartesian coordinate in D
dimensions relative to an origin, e.g. the z coordinate in 3D

.. math:: f(x, y, z) = z - z_0

This PositionFunction is available as

.. code-block:: cpp

    double r_0[3] = {0.0, 0.0, 0.0};
    PositionFunction<3> r_func(r_0);
    r_func.setDir(2);

------------------------------
MultiResolution Analysis (MRA)
------------------------------

In order to combine different functions and operators in mathematical operations
they must be compatible, that is, they have to be
defined on the same computational domain and constructed using the same
polynomial basis (order and type). This information constitutes an MRA, which
needs to be defined and passed as argument to all function and operator
constructors, and only functions and operators with compatible MRAs can be
combined in subsequent calculations. An MRA is defined in two steps, first the
computational domain is given by a BoundingBox (D is the dimension)

.. code-block:: cpp

    int n;                          // Root scale
    int l[D];                       // Translation of first box
    int nb[D];                      // Number of boxes
    NodeIndex<D> idx(n, l);         // Index defining the first box
    BoundingBox<D> world(idx, nb);

which is combined with a ScalingBasis to give an MRA

.. code-block:: cpp

    int k;                          // Polynomial order
    int type;                       // Polynomial type
    ScalingBasis basis(k, type);
    MultiResolutionAnalysis<D> MRA(world, basis);

Two polynomial types are supported (Legendre and Interpol), and they are
both available at orders :math:`k=1,2,\dots,40` (note that some operators are
constructed using intermediates of order :math:`2k`, so in that case the maximum
order available is :math:`k=20`).

------------------------
Function representations
------------------------

MW representations of functions are called ``FunctionTrees``, and are in principle
available in any dimension using the template parameter D (in practice D=1,2,3).
There are three different ways of constructing MW functions (computing the
expansion coefficients in the MW basis)

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
indirectly through ``TreeBuilder`` objects:

.. code-block:: cpp

    FunctionTree<D> *tree = builder();

where ``builder`` is any of the ``TreeBuilders`` presented below which may or may not
take any arguments for the construction. Details on how the tree structure is
built and how the MW coefficients are computed are specified in each particular
``TreeBuilder``.

Integrals are computed very efficiently in the orthonormal MW basis, and among
the important methods of the ``FunctionTree`` are estimating the error in the
representation (based on the wavelet norm), obtaining the squared
:math:`L^2`-norm of the function, as well as its integral and dot product with
another ``FunctionTree`` (both over the full computational domain)

.. code-block:: cpp

    double error = f_tree.estimateError();
    double sq_norm = f_tree.getSquareNorm();
    double integral = f_tree.integrate();
    double dot_prod = f_tree.dot(g_tree);

FunctionTreeVector
------------------

The ``FunctionTreeVector`` is a convenience class for a collection of ``FunctionTrees``
which basically consists of two STL vectors, one containing pointers to
``FunctionTrees`` and one with corresponding numerical coefficients.
Elements can be appended to the vector

.. code-block:: cpp

    FunctionTreeVector<D> tree_vec;
    tree_vec.push_back(2.0, tree_a);
    tree_vec.push_back(tree_b);
    tree_vec.clear();

where ``tree_b`` will be appended with a default coefficient of 1.0. Clearing the
vector means removing all its elements. Note that the vector holds only pointers
to ``FunctionTrees``, so clearing its content will `not` delete the ``FunctionTrees``.
This means that a particular tree can be part of several vectors, and the
construction/destruction of the tree must be handled outside.

TreeBuilder
-----------

This is the class that is responsible for the construction of
``FunctionTrees``, which involves allocating memory, growing a tree structure and
calculating MW coefficients. The ``TreeBuilder`` has two important members: a
``TreeCalculator`` that defines how the MW coefficients are computed, and a
``TreeAdaptor`` that defines how the tree structure is grown. There are four
different ways of computing MW coefficients (projection, addition,
multiplication and operator application), and we have the corresponding
``TreeBuilders`` (the MW prefix indicates that they compute MW coefficients)

* MWProjector
* MWAdder
* MWMultiplier
* MWOperator

Each of these is a specialization of the ``TreeBuilder`` class that differs in the
type of ``TreeCalculator``, and they all contain a ``WaveletAdaptor`` that controls the
accuracy of the function representations they build. All ``TreeBuilders`` have the
same fundamental building algorithm:

    1. Start with an initial guess for the grid
    2. Compute the output function on the current grid
    3. Use the ``WaveletAdaptor`` to refine the grid where needed
    4. Iterate point 2 and 3 until the grid is converged

All builders are constructed using the same arguments

.. code-block:: cpp

    int max_iter;
    double prec;
    MultiResolutionAnalysis<D> MRA;
    TreeBuilder<D> builder(MRA, prec, max_iter);

Where the ``MRA`` defines the computational domain and type of MW basis,
``prec`` defines the relative precision used by the ``WaveletAdaptor`` and
``max_iter`` will stop
the grid refinement even if the accuracy criterion is not met. The last two
arguments have negative defaults, which for ``prec`` means that no
refinement will take place beyond the initial grid, and for ``max_iter`` it
means that the number of iterations is unbounded.

The interface for the ``TreeBuilders`` is mainly the ``operator()``, which comes in
two versions

.. code-block:: cpp

    out = builder(inp);
    builder(out, inp);

where the former is a constructor that returns a pointer to a new
``FunctionTree``,
while the latter will work on an already existing tree. The main difference
between the two is the choice of initial grid: the former will automatically
construct a default grid based on the operation that is taking place (e.g.
arithmetic operators will copy the grids of the input functions); the latter
will use whatever grid is already present in the output tree structure which
allows for more detailed control for the user, however the grids needs to be
prepared in advance, either using a ``GridGenerator`` to construct an empty grid or
a ``GridCleaner`` to clear an existing function representation (see advanced
initialization below).

TreeCalculator
--------------

This class operates on the node level, computing MW coefficients based on the
proper input data (analytic functions in the case of projection,
``FunctionTrees`` in the case of operators). The ``TreeCalculator`` is hidden within the
``TreeBuilder``, and is not part of its interface. There is one calculator for each
of the MW-types of ``TreeBuilder``:

* ProjectionCalculator
* AdditionCalculator
* MultiplicationCalculator
* OperationCalculator

TreeAdaptor
-----------

Like the ``TreeCalculator``, this class operates on the node level, but instead of
computing coefficients, it decides whether each node needs to be split into
:math:`2^D` children nodes. There can be different reasons for splitting nodes,
the most important being to reduce the wavelet norm of the representation.
There are three different ``TreeAdaptors``:

* WaveletAdaptor
* AnalyticAdaptor
* CopyAdaptor

where the ``WaveletAdaptor`` tests the wavelet norm, the
``AnalyticAdaptor`` use some known information of an analytic function, and the
``CopyAdaptor`` will copy the node structure of another tree.

MWProjector
-----------

Given an analytic D-dimensional function ``f_func``, we can obtain its
numerical MW representation by projecting it onto the MW basis. For this we
have the ``MWProjector``

.. code-block:: cpp

    MWProjector<D> Q(MRA, prec, max_iter);
    FunctionTree<D> *f_tree = Q(f_func);

The projector will construct an
initial grid based on the analytic function (this is meant as a way to make sure
that the projection starts on a grid where the function is actually visible,
as for very narrow Gaussians, it's `not` meant to be a good approximation of the
final grid). The default projector (``prec`` and ``max_iter`` negative) will
simply do the projection on the initial grid with no further grid refinement.
By specifying a positive ``prec`` the grid will automatically be adapted to
represent the function to the given precision, based on the wavelet norm of
the representation. You can also allow the grid to be refined only a certain number
of iterations beyond the initial grid by specifying a positive ``max_iter``
(this will of course not guarantee the accuracy of the representation).

Arithmetic operations
---------------------

Arithmetic operations in the MW representation are performed using the
``FunctionTreeVector``, and the general sum :math:`g = \sum_i c_i f_i(x)` and
product :math:`h = \prod_i c_i f_i(x)` are done in the following way

.. code-block:: cpp

    FunctionTreeVector inp_vec;
    inp_vec.push_back(c_1, f_tree_1);
    inp_vec.push_back(c_2, f_tree_2);
    inp_vec.push_back(c_3, f_tree_3);

    MWAdder<D> add(MRA, prec, max_iter);
    FunctionTree<D> *g_tree = add(inp_vec);

    MWMultiplier<D> mult(MRA, prec, max_iter);
    FunctionTree<D> *h_tree = mult(inp_vec);

The default initial grid for the arithmetic operators is the union of the grids
of the input functions.

Note that in the case of addition there is no extra information to be gained
by going beyond the finest refinement levels of the input functions, so the
union grid summation is simply the best you can do, and adding a positive
``prec`` will not make a difference (there are situations where
you want to use a `smaller` grid, and this is discussed below under advanced
initialization).

In the case of multiplication, however, there might be a loss of accuracy if
the product is restricted to the union grid. The reason for this is that the
product will contain signals of higher frequency than each of the input
functions, which require a higher grid refinement for accurate representation.
By specifying a positive ``prec`` you will allow the grid to adapt to the higher
frequencies, but it is usually a good idea to restrict to one extra refinement
level (by setting ``max_iter=1``) as the grids are not guaranteed to
converge for such local operations (like arithmetics, derivatives and
function mappings).

-----------------------
Advanced initialization
-----------------------

The ``TreeBuilders``, as presented above, have a clear and limited interface, but
there is one important drawback: every operation require the construction
of a new ``FunctionTree`` from scratch (including extensive memory allocation).
In many practical applications however (e.g. iterative algorithms), we are
recalculating the same functions over and over, where the requirements on the
numerical grids change only little between each iteration. In such situations it
will be beneficial to be able to reuse the existing grids without reallocating
the memory. For this purpose we have the following additional ``TreeBuilder``
sub-classes (the Grid prefix indicates that they do not compute MW
coefficients):

* GridGenerator
* GridCleaner

where the former constructs empty grids from scratch and the latter clears the
MW coefficients on an existing ``FunctionTree``. The end result is in both cases an
empty tree skeleton with no MW coefficients (undefined function).

GridGenerator
-------------

Sometimes it is useful to construct an empty grid based on some available
information of the function that is about to be represented. This can be e.g.
that you want to copy the grid of an existing ``FunctionTree`` or that an analytic
function has more or less known grid requirements (like Gaussians). Sometimes it
is even necessary to force the grid refinement beyond the coarsest scales in
order for the ``WaveletAdaptor`` to detect a wavelet "signal" that allows it to do
its job properly (this happens for narrow Gaussians where non of the initial
quadrature points hits a function value significantly different from zero).
In such cases we use a ``GridGenerator`` to build the initial tree structure.

A special case of the ``GridGenerator`` (with no argument) corresponds to the
default constructor of the ``FunctionTree``

.. code-block:: cpp

    GridGenerator<D> G(MRA);
    FunctionTree<D> *f_tree = G();

which will construct a new ``FunctionTree`` with empty nodes (undefined
function with no MW coefs), containing only the root nodes of the given MRA.
Passing an analytic function as argument to the generator will use an
``AnalyticAdaptor`` to build a grid
based on some predefined knowledge of the function (if there are any, otherwise
it is identical to the default constructor)

.. code-block:: cpp

    FunctionTree<D> *f_tree = G(f_func);

while passing a ``FunctionTree`` to the generator will copy its grid using
a ``CopyAdaptor``

.. code-block:: cpp

    FunctionTree<D> *f_tree = G(g_tree);

Both of these will produce a skeleton ``FunctionTree`` with empty nodes. In order
to define a function in the new tree it is passed as the first argument to the
regular ``TreeBuilders`` presented above, e.g for projection

.. code-block:: cpp

    GridGenerator<D> G(MRA);
    MWProjector<D> Q(MRA, prec, max_iter);
    FunctionTree<D> *f_tree = G(f_func);
    Q(*f_tree, f_func);

This will first produce an empty grid suited for representing the analytic
function ``f_func`` and then perform the projection on the given numerical grid.
This will in fact be identical to the default projection

.. code-block:: cpp

    MWProjector<D> Q(MRA, prec, max_iter);
    FunctionTree<D> *f_tree = Q(f_func);

as the same ``GridGenerator`` and ``AnalyticAdaptor`` is used under the hood to
construct the default
initial grid. Similar notation applies for all ``TreeBuilders``, if an undefined
``FunctionTree`` is given as first argument, it will not construct a new tree
but perform the operation on the one given (the given tree is used as starting
point for the ``TreeBuilder``, and further grid refinements can occur if a
``TreeAdaptor`` is present), e.g. the grid copy can be done in two steps as

.. code-block:: cpp

    FunctionTree<D> *f_tree = G();
    G(*f_tree, g_tree);

Actually, the effect of the ``GridGenerator`` is to `extend` the existing grid
with any missing nodes relative to the input. This means that we can build the
union of two grids by successive application of the generator

.. code-block:: cpp

    G(f_tree, g_tree);
    G(f_tree, h_tree);

and one can make the grids of two functions equal to their union

.. code-block:: cpp

    G(f_tree, g_tree);
    G(g_tree, f_tree);

The union grid of several trees can be constructed using a ``FunctionTreeVector``

.. code-block:: cpp

    FunctionTreeVector<D> inp_vec;
    inp_vec.push_back(tree_1);
    inp_vec.push_back(tree_2);
    inp_vec.push_back(tree_3);
    FunctionTree<D> *f_tree = G(inp_vec);

If you don't want to use the default starting grids (e.g. union grid for
arithmetic operations) you can use the ``GridGenerator`` to construct the grid
you want. To multiply functions adaptively based on the wavelet norm,
you specify the precision of the ``MWMultiplier`` and start from a default
grid of root nodes

.. code-block:: cpp

    double prec;
    MWMultiplier<D> mult(MRA, prec);
    GridGenerator<D> G(MRA);

    FunctionTree<D> *f_tree = G();  // Construct empty grid of root nodes
    mult(*f_tree, inp_vec);         // Build result based on wavelet norm

If you have a summation over several functions but want to perform the
addition on the grid given by the `first` input function, you first copy the
wanted grid and then perform the operation on that grid

.. code-block:: cpp

    FunctionTreeVector<D> inp_vec;
    inp_vec.push_back(coef_1, tree_1);
    inp_vec.push_back(coef_2, tree_2);
    inp_vec.push_back(coef_3, tree_3);

    MWAdder<D> add(MRA);
    GridGenerator<D> G(MRA);

    FunctionTree<D> *f_tree = G(tree_1); // Copy grid of first input function
    add(*f_tree, inp_vec);               // Perform addition on given grid

In the last example you can of course also add a ``prec`` to the ``MWAdder``,
and the resulting function will be built adaptively starting from the given
initial grid. This in useful e.g. when performing a uniform transformation
among a set of ``FunctionTrees`` where you usually don't want to construct `all`
the output functions on the union grid of `all` the input functions. To allow for
individual grid adaptivity for the output functions you start all additions from
the root node and let the ``WaveletAdaptor`` build customized grids for each
function based on the precision

.. code-block:: cpp

    MWAdder<D> add(MRA, prec);
    GridGenerator<D> G(MRA);

    FunctionTreeVector<D> inp_vec;
    inp_vec.push_back(coef_1, inp_tree_1);
    inp_vec.push_back(coef_2, inp_tree_2);
    inp_vec.push_back(coef_3, inp_tree_3);

    FunctionTree<D> *out_tree = G(); // Start from simple grid of root nodes
    add(*out_tree, inp_vec);         // Refine adaptively based on precision


GridCleaner
-----------

Given a ``FunctionTree`` that is a valid function representation we can clear its
MW expansion coefficients (while keeping the grid refinement) with the
``GridCleaner`` (unlike the other ``TreeBuilders``, the ``GridCleaner`` will not return a
``FunctionTree`` pointer, as it would always be the same as the argument)

.. code-block:: cpp

    GridCleaner<D> C(MRA);
    C(f_tree);

This action will leave the ``FunctionTree`` in the same state as the ``GridGenerator``
(uninitialized function), and its coefficients can now be re-computed.

In certain situations it might be desireable to separate the actions of the
``TreeCalculator`` and the ``TreeAdaptor``. For this we can combine the
``GridCleaner`` with a ``WaveletAdaptor``, which will adaptively refine the
grid (based on the wavelet norm and the given precision with ``max_iter=1``)
`before` it is cleared

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
    GridCleaner<D> C(MRA, prec); // The precision parameter is passed as
    MWProjector<D> Q(MRA);       // argument to the cleaner, not the projector

    int n_nodes = 1;
    while (n_nodes > 0) {
        Q(f_tree, f_func);
        n_nodes = C(f_tree);
    }
    Q(f_tree, f_func);

