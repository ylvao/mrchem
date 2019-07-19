============
Installation
============


-------------------
Build prerequisites
-------------------

The supplied ``setup`` script should be able to configure things
correctly, provided all the necessary tools are available.

Python
------

MRChem depends on some Python tools for its compilation and execution.
More specifically, users will need Python:

- for configuration using the script ``setup``.
- to run the integration tests, using `runtest <https://runtest.readthedocs.io/en/latest/>`_.
- for input parsing, using `parselglossy  <https://github.com/dev-cafe/parselglossy>`_.
- for running calculations, using the ``mrchem`` script.

Developers will additionally need Python to build the documentation locally with
Sphinx.

We **strongly** suggest not to install these Python dependencies globally, but
rather to use a local virtual environment. We provide both a ``Pipfile`` and
a ``requirements.txt`` specifying the Python dependencies.
We recommend using `Pipenv <https://pipenv.readthedocs.io/en/latest/>`_, since
it manages virtual environment and package installation seamlessly.
After installing it with your package manager, run::

    $ pipenv install

to create a virtual environment with all packages installed.

.. note::
   Developers will want to run::

      $ pipenv install --dev

   to also install development packages, such as those needed to generate
   documentation with Sphinx.

The environment can be activated with::

    $ pipenv shell

Alternatively, any Python command can be run within the virtual environment by
doing::

    $ pipenv run python -c "print('Hello, world')"

Stallo
------

Using the Intel tool chain on Stallo::

    $ module load intel/2018b
    $ module load Python/3.6.6-intel-2018b
    $ module load Eigen/3.3.5
    $ module load CMake/3.12.1-GCCcore-7.3.0
    $ module load Git/2.8.1

Fram
----

Using the Intel tool chain on Fram::

    $ module load intel/2017a
    $ module load Python/3.6.1-intel-2017a
    $ module load CMake/3.12.1
    $ module load git/2.14.1-GCCcore-6.4.0

Eigen is not available through the module system on Fram, but it will be
download during the CMake configuration step. It can also be installed manually
by the user, see the `Eigen3 <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_
home page or the ``.ci/eigen.sh`` source file for details.

-------------------------------
Obtaining and building the code
-------------------------------

The public version of MRChem is available on GitHub::

    $ git clone git@github.com:MRChemSoft/mrchem.git

MRChem also depends on `XCFun <https://github.com/dftlibs/xcfun>`_,
`Getkw <https://github.com/dev-cafe/libgetkw>`_ and `MRCPP <https://github.com/MRChemSoft/mrcpp>`_.
We have structured the build system such that these dependencies will be **fetched**
when configuring the project if they are not already installed.

To build the code with only shared memory OpenMP parallelization::

    $ cd mrchem
    $ pipenv install
    $ ./setup --prefix=<install-dir> --omp <build-dir>
    $ cd <build-dir>
    $ make
    $ make install

With the Intel tool chain on Stallo or Fram you need to specify the compilers
in the setup::

    $ ./setup --prefix=<install-dir> --omp --cxx=icpc <build-dir>

To build the code with hybrid MPI/OpenMP parallelization::

    $ ./setup --prefix=<install-dir> --omp --mpi --cxx=mpiicpc <build-dir>

The MRChem executables will be installed in ``<install-dir>/bin``.

If the dependencies are already installed, say in the folder ``<dependencies-prefix>``,
you can pass the following flag to the ``setup`` script to use them::

    $ ./setup --prefix=<install-dir> --cmake-options="-DXCFun_DIR="<dependencies-prefix>/share/cmake/XCFun" -Dgetkw_DIR="<dependencies-prefix>/share/cmake/getkw" -DMRCPP_DIR="<dependencies-prefix>/share/cmake/MRCPP""

When fetching dependencies at configuration, a shallow Git clone is performed into ``<build-dir>``.
For the example of MRCPP, sources are saved into the folders ``<build-dir>/mrcpp_sources-src`` and built into ``<build-dir>/mrcpp_source-build``.
If you want only to recompile one of the external libraries (for example MRCPP), without rebuilding from scratch, try::

   $ cd <build-dir>/mrcpp_sources-build
   $ make

Note that this will leave your build in an undefined state, since it will not try to update the MRChem parts.

-------
Testing
-------

A test suite is provided to make sure that everything compiled properly. To run
a collection of small unit tests::

    $ cd <build-dir>
    $ ctest -L unit

To run a couple of more involved integration tests::

    $ cd <build-dir>
    $ pipenv run ctest -L integration

Note how we used Pipenv to run the integration tests. This ensures that the
Python dependencies (``parselglossy`` and ``runtest``) are satisfied in a
virtual environment and available to ``ctest``.

-------------------
Running the program
-------------------

A Python input parser will be provided along with the mrchem executable::

    $ build/bin/mrchem              // Python input parser
    $ build/bin/mrchem.x            // MRChem executable

The input parser takes a single file argument. It will process the `user_input`
(usually ``mrchem.inp``) and generate a `program_input` file with a ``.json``
extension (usually ``mrchem.json``), and pass it to the main MRChem executable.
Output is written to ``stdout`` but can be redirected to an output file::

    $ pipenv run ./mrchem mrchem.inp > mrchem.out &

Note again that the ``pipenv`` environment must be activated, either `once` by
``pipenv shell`` or on `every` command using ``pipenv run``, see Python section
above. To run the program in OpenMP parallel use the environment variable
``OMP_NUM_THREADS`` (``unset OMP_NUM_THREADS`` will give you all threads
available, otherwise use ``export OMP_NUM_THREADS N``)::

    $ export OMP_NUM_THREADS 16
    $ pipenv run mrchem mrchem.inp

When you run the program in hybrid MPI/OpenMP parallel, you must run the input
parser manually first with the dryrun ``-D`` option, before launching the main
executable with ``mpirun`` (or equivalent). For 20 threads each on 5 MPI
processes::

    $ pipenv run mrchem -D mrchem.inp
    $ OMP_NUM_THREADS=20  mpirun -np 5 mrchem.x mrchem.json >mrchem.out &
