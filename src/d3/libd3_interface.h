#ifndef DFTD3_INTERFACE_H
#define DFTD3_INTERFACE_H

// Use extern "C" to prevent C++ name mangling
// the actual Fortran library builds a C binding called "wrapper" that
// initializes dftd3 and computes dispersion/gradients for given coordinates.
// signature mirrors the subroutine in lib/wrapper.f90
extern "C" {
    // Fortran declaration uses `value` for the scalar integers, therefore the
    // C binding receives them by value, not by pointer.  Failing to match this
    // was causing the D3 routine to see gibberish atom counts and versions,
    // which in turn produced zero energies/gradients or crashes.
    void wrapper(int natoms,
                 double *coords,   // array natoms-by-3 in row-major order
                 int *itype,       // atomic numbers
                 const char *fun,  // functional name string (null terminated)
                 int version,      // functional version
                 int tz,           // zero/one for tz flag
                 double *edisp,    // output energy
                 double *grads);   // output gradients (3 x natoms)
}

#endif