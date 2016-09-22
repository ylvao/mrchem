============
Installation
============


-------------------
Build prerequisites
-------------------

On Stallo the supplied setup script should be able to configure things
correctly, provided all the necessary modules have been loaded::

    $ module load eigen cmake python


-------------------------------
Obtaining and building the code
-------------------------------

An experimental version of MRChem is available on GitHub::

    $ git clone git@github.com:MRChemSoft/mrchem.git
    $ cd mrchem
    $ ./setup [--cc=icc --cxx=icpc]
    $ cd build
    $ make

The official MRChem main program is located in ``src/mrchem/mrchem.cpp``, whose
executable will be built in ``build/bin/mrchem.x``. Please do not change this
file unless you know what you are doing. To try out your own ideas you can
instead write a separate main program in a file
called ``mrchem.cpp`` in the ``pilot`` directory. You will find a sample code
called ``mrchem.cpp.sample`` in this directory where some of the functionality
of MRCPP is demonstrated. To activate it, rename it ``mrchem.cpp`` *before* you
run the setup script::

    $ git clone git@github.com:MRChemSoft/mrchem.git
    $ cd mrchem/pilot
    $ cp mrchem.cpp.sample mrchem.cpp
    $ cd ..
    $ ./setup [--cc=icc --cxx=icpc]
    $ cd build
    $ make

The pilot executable will now be built in ``build/pilot/mrchem-pilot.x``.
Feel free to do whatever you like with your pilot code, but it is your own
personal playground so don't add this file to git.

A test suite is provided to make sure that everything compiled properly::

    $ cd build/tests
    $ ./unit-tests


-------------------
Running the program
-------------------

A Python input parser will be provided along with the mrchem
executables both for the official program and the pilot code::

    $ build/bin/mrchem              // Python input parser
    $ build/bin/mrchem.x            // MRChem executable

    $ build/pilot/mrchem            // Python input parser
    $ build/pilot/mrchem-pilot.x    // Pilot executable

The input parser takes a single file argument (default ``mrchem.inp``),
processes the input and calls the main executable. Output is written to stdout
but can be redirected to an output file::

    $ ./mrchem mrchem.inp > mrchem.out &

A sample input file is provided for the pilot code. For the official program,
please refer to the MRChem manual or the ``example`` directory.

By following the instructions above the code will be compiled in OpenMP
parallel. To run the program in parallel use the environment variable
``OMP_NUM_THREADS`` (``unset OMP_NUM_THREADS`` will give you all threads
available, otherwise use ``export OMP_NUM_THREADS N``)::

    $ export OMP_NUM_THREADS 16
    $ ./mrchem
