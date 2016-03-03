===========
 MRCPP API
===========

--------------
 Introduction
--------------

The main features of MRCPP are the numerical multiwavelet (MW) representations 
of functions and operators. Two integral convolution operators are implemented 
(the Poisson and Helmholtz operators), as well as the 
partial derivative and arithmetic operators. In addition
to the numerical representations there are a limited number of analytic
functions that are usually used as starting point for the numerical
computations.

------------------
 AnalyticFunction
------------------

Note that all these analytic functions needs to be projected onto the MW basis 
before they can be used in numerical computations (see below).

GaussFunc
=========

A very important analytic function is the Cartesian Gaussian

math::

    f(r) = \alpha (x-x_0)^a (y-y_0)^b (z-z_0)^c e^{-\beta |r-r_0|^2}

which is initialized in the following way (in 3D)::

    double alpha, beta;
    double pos[3] = {x_0, y_0, z_0};
    int pow[3] = {a, b, c};
    GaussFunc<3> f_func(beta, alpha, pos, pow);

GaussPoly
=========

GaussPoly is a generalization of the GaussFunc, where there is an arbitrary
polynomial in front of the exponential

math::
    f(r) = \alpha P(r-r_0) e^{-\beta |r-r_0|^2}

GaussExp
========

A GaussExp is a collection of Gaussian in the form 

math::

    G(r) = \sum_i c_i g_i(r)

where :math:`g_i` can be either GaussFunc or GaussPoly

math::
    g_i(r) =  \alpha_i P_i(r-r_i)e^{-\beta_i|r-r_i|^2}

PositionFunction
================

There is also a very simple analytic function for the Cartesian coordinate in D
dimensions relative to an origin, e.g. the z coordinate in 3D

math::
    f(x, y, z) = z - z_0

This PositionFunction is available as::

    double r_0[3] = {0.0, 0.0, 0.0};
    PositionFunction<3> r_func(r_0);
    r_func.setDir(2);

--------------------------------
 MultiResolution Analysis (MRA)
--------------------------------

In order to combine different functions and operators in mathematical operations
they must be compatible, that is, they have to be
defined on the same computational domain and constructed using the same 
polynomial basis (order and type). This information constitutes an MRA, which 
needs to be defined and passed as argument to all function and operator 
constructors, and only functions and operators with compatible MRAs can be 
combined in subsequent calculations. An MRA is defined in two steps, first the
computational domain is given by a BoundingBox (D is the dimension)::

    int n;                          // Root scale
    int l[D];                       // Translation of first box
    int nb[D];                      // Number of boxes
    NodeIndex<D> idx(n, l);         // Index defining the first box
    BoundingBox<D> world(idx, nb);

which is combined with a ScalingBasis to give an MRA::

    int k;                          // Polynomial order
    int type;                       // Polynomial type
    ScalingBasis basis(k, type);
    MultiResolutionAnalysis<D> MRA(world, basis);

Two polynomial types are supported (Legendre and Interpol), and they are 
both available at orders :math:`k=1,2,\dots,40` (note that some operators are 
constructed using intermediates of order $2k$, so in that case the maximum 
order available is :math:`k=20`).

--------------------------
 Function representations
--------------------------

MW representations of functions are called FunctionTrees, and are in principle 
available in any dimension using the template parameter D (in practice D=1,2,3).
Constructing a full grown FunctionTree involves a number of steps, including
setting up a memory allocator, constructing root nodes, building a tree
structure, distributing memory among MPI hosts and computing MW
coefficients. For this reason the FunctionTree constructor is made protected,
and all construction is done by TreeBuilder objects.

There are tree different ways of constructing MW functions (computing the 
expansion coefficients in the MW basis)

* Projection of analytic function
* Arithmetic operations
* Application of MW operator

and these will be described in the following sections. Integrals are
computed very efficiently in the orthonormal MW basis, and among the important
methods of the FunctionTree are estimating the error in the representation
(based on the wavelet norm), obtaining the squared :math:`L^2`-norm of the
function, as well as its integral and dot product with another FunctionTree
(both over the full computational domain)::

    double error = f_tree.estimateError();
    double sq_norm = f_tree.getSquareNorm();
    double integral = f_tree.integrate();
    double dot_prod = f_tree.dot(g_tree);

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

-------------
 TreeBuilder
-------------

This is the class that is responsible for the construction of 
FunctionTrees, which involves allocating memory, growing a tree structure and 
calculating MW coefficients. The TreeBuilder has two important members: a
TreeCalculator that defines how the MW coefficients are computed, and a
TreeAdaptor that defines how the tree structure is grown. There are four 
different ways of computing MW coefficients (projection, addition,
multiplication and operator application), and we have the corresponding
TreeBuilders (the MW prefix indicates that they compute MW coefficients) 

* MWProjector
* MWAdder
* MWMultiplier
* MWOperator

Each of these is a specialization of the TreeBuilder class that differs in the
type of TreeCalculator, and can be combined with any TreeAdaptor. The interface
of these classes is mainly the ``operator()``, which is overloaded with the 
proper input argument(s) and will return a pointer to the newly constructed 
FunctionTree. All TreeBuilders take an MRA as the first argument in their
constructor, and all FunctionTrees produced by this TreeBuilder will get the
same MRA.

----------------
 TreeCalculator
----------------

This class operates on the node level, computing MW coefficients based on the
proper input data (analytic functions in the case of projection,
FunctionTrees in the case of operators). The TreeCalculator is hidden within the
TreeBuilder, and is not part of its interface. There is one calculator for each 
of the MW-types of TreeBuilder:

* ProjectionCalculator
* AdditionCalculator
* MultiplicationCalculator
* OperationCalculator

-------------
 TreeAdaptor
-------------

Like the TreeCalculator, this class operates on the node level, but instead of
computing coefficients, it decides whether each node needs to be split into
:math:`2^D` children nodes. There can be different reasons for splitting nodes, 
the most important being to reduce the wavelet norm of the representation. 
There are three different TreeAdaptors: 

* WaveletAdaptor
* AnalyticAdaptor
* CopyAdaptor

where the WaveletAdaptor tests the wavelet norm, the
AnalyticAdaptor use some known information of an analytic function, and the
CopyAdaptor will copy the node structure of another tree. 

 MWProjector
=============

Given an analytic D-dimensional function f\_func, we can obtain its 
numerical MW representation by projecting it onto the MW basis. For this we 
have the MWProjector::

    MWProjector<D> Q(MRA);
    FunctionTree<D> *f_tree = Q(f_func);

The default projector will simply project the function onto the grid that is
defined by its root nodes, with no regard on grid adaptivity. If you want to 
control the accuracy of the approximation you need to add an adaptor to 
the projector::

    double prec;
    WaveletAdaptor<D> w_adaptor(prec);
    MWProjector<D> Q(MRA, w_adaptor);
    FunctionTree<D> *f_tree = Q(f_func);

The WaveletAdaptor will automatically construct the necessary grid needed to 
represent the function to the given precision, based on the wavelet norm of 
the representation. 

 Arithmetic operations
=======================

Given two functions :math:`f` and :math:`g` in MW representation 
(FunctionTrees), we can compute the sum (e.g. :math:`h = f - 2g`) or 
product (e.g. :math:`h = f\times 2g`)::

    MWAdder<D> add(MRA);
    FunctionTree<D> *h_tree = add(1.0, f_tree, -2.0, g_tree)

    double prec;
    WaveletAdaptor<D> w_adaptor(prec);
    MWMultiplier<D> mult(MRA, w_adaptor);
    FunctionTree<D> *h_tree = mult(2.0, f_tree, g_tree)

where the addition is using a CopyAdaptor that will perform the addition on the 
union grid of the input functions. The WaveletAdaptor will build the product tree
adaptively based on the wavelet norm for the multiplication. Note that any
adaptor can in principle be used for any TreeBuilder, and the CopyAdaptor is the
default for the arithmetic operations, but there are situations where the
WaveletAdaptor is appropriate (e.g. when the addition corresponds to a unitary
transformation among a set of functions).

When more than two functions are involved in the arithmetics it might
be beneficial to combine them into a single operation using the STL vector::

    vector<double> coefs;
    vector<FunctionTree<D> *> trees;

    FunctionTree<D> *h_tree = add(coefs, trees);
    FunctionTree<D> *h_tree = mult(coefs, trees);

A number of in-place operations are also available::

    f_tree *= 2.0;
    f_tree *= g_tree;
    f_tree += g_tree;
    f_tree -= g_tree;
    f_tree.square();
    f_tree.pow(3.0/2.0);
    f_tree.normalize();
    f_tree.orthogonalize(g_tree);

-------------------------
 Advanced initialization
-------------------------

The TreeBuilders, as presented above, have a clear and limited interface, but 
there is one important drawback: every operation require the construction
of a new FunctionTree from scratch (including extensive memory allocation). 
In many practical applications however (e.g. iterative algorithms), we are 
recalculating the same functions over and over, where the requirements on the
numerical grids change only little between each iteration. In such situations it 
will be beneficial to be able to reuse the existing grids without reallocating
the memory. For this purpose we have the following additional TreeBuilder 
sub-classes (the Grid prefix indicates that they do not compute MW 
coefficients):

* GridGenerator
* GridCleaner

where the former constructs empty grids from scratch and the latter clears the
MW coefficients on an existing FunctionTree. The end result is in both cases an
empty tree skeleton with no MW coefficients (undefined function).

 GridGenerator
===============

Sometimes it is useful to construct an empty grid based on some available 
information of the function that is about to be represented. This can be e.g.
that you want to copy the grid of an existing FunctionTree or that an analytic
function has more or less known grid requirements (like Gaussians). Sometimes it
is even necessary to force the grid refinement beyond the coasest scales in 
order for the WaveletAdaptor to detect a wavelet "signal" that allows it to do
its job properly (this happens for narrow Gaussians where non of the initial
quadrature points hits a function value significantly different from zero).
In such cases we use a GridGenerator to build the initial tree structure.

A special case of the GridGenerator (with no argument) corresponds to the 
default constructor of the FunctionTree::

    GridGenerator<D> G(MRA);
    FunctionTree<D> *f_tree = G();

which will construct a new FunctionTree with empty nodes (undefined
function with no MW coefs), containing only the root nodes of the given MRA.
Passing an analytic function as argument to the generator will build a grid 
based on some predefined knowledge of the function (if there are any, otherwise
it is identical to the default constructor)::

    FunctionTree<D> *f_tree = G(f_func);

while passing a FunctionTree to the generator will copy its grid::

    FunctionTree<D> *f_tree = G(g_tree);

Both of these will produce a skeleton FunctionTree with empty nodes. In order 
to define a function in the new tree it is passed as the first argument to the 
regular TreeBuilders presented above, e.g for projection::

    GridGenerator<D> G(MRA);
    MWProjector Q(MRA);
    FunctionTree<D> *f_tree = G(f_func);
    Q(*f_tree, f_func);

This will first produce an empty grid suited for representing the analytic
function f\_func and then perform the projection on the given numerical grid.
Similar notation applies for all TreeBuilders, if an undefined FunctionTree is 
given as first argument, it will not construct a new tree but perform the 
operation on the one given (the given tree is used as starting point for the 
TreeBuilder, and further grid refinements can occur if a TreeAdaptor is
present), e.g. the grid copy can be done in two steps as::

    FunctionTree<D> *f_tree = G();
    G(*f_tree, g_tree);

Actually, the effect of the GridGenerator is to `extend` the existing grid 
with any missing nodes relative to the input. This means that we can build the
union of two grids by successive application of the generator::

    G(f_tree, g_tree);
    G(f_tree, h_tree);

and one can make the grids of two functions equal to their union::

    G(f_tree, g_tree);
    G(g_tree, f_tree);


 GridCleaner
=============

Given a FunctionTree that is a valid function representation we can clear its 
MW expansion coefficients (while keeping the grid refinement) with the 
GridCleaner (unlike the other TreeBuilders, the GridCleaner will not return a 
FunctionTree pointer, as it would always be the same as the argument)::

    GridCleaner<D> C(MRA);
    C(f_tree);

This action will leave the FunctionTree in the same state as the GridGenerator
(uninitialized function), and its coefficients can now be re-computed. 

In certain situations might be desireable to separate the actions of the 
projector and the wavelet adaptor. For this we can combine the GridCleaner 
with an adaptor, which will adaptively refine the grid \emph{before} it is 
cleared::

    double prec;
    WaveletAdaptor<D> w_adaptor(prec);
    GridCleaner<D> C(MRA, w_adaptor);
    C(f_tree);

One example where this might be
useful is in iterative algorithms where you want to fix the grid size for 
all calculations within one cycle and then relax the grid in the end in 
preparation for the next iteration. The following is equivalent to the adaptive
projection above (the cleaner returns the number of new nodes that were created
in the process)::

    double prec;
    WaveletAdaptor<D> w_adaptor(prec);
    GridCleaner<D> C(MRA, w_adaptor);     // The adaptor is passed as argument
    MWProjector<D> Q(MRA);                // to the cleaner, not the projector

    int n_nodes = 1;
    while (n_nodes > 0) {
        Q(f_tree, f_func);
        n_nodes = C(f_tree);
    }
    Q(f_tree, f_func);

