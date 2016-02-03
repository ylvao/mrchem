/** Class to parse and process basis sets in MOLECULE or INTGRL format.
 */
#ifndef INTGRL_H
#define INTGRL_H

#include <vector>
#include <string>
#include <iostream>
#include "Atom.h"
#include "AOBasis.h"

class Intgrl {
public:
    Intgrl(const std::string &file);
    virtual ~Intgrl();

    int getNAtoms() const { return this->atoms.size(); }

    const Atom &getAtom(int i) const { return *this->atoms[i]; }
    const AOBasis &getAOBasis(int i) const { return *this->basis[i]; }

    Atom &getAtom(int i) { return *this->atoms[i]; }
    AOBasis &getAOBasis(int i) { return *this->basis[i]; }

    GaussExp<3> getAtomBasis(int i, bool norm = true) const;
    GaussExp<3> getMolBasis(bool norm = true) const;

protected:
    std::vector<Atom *> atoms;
    std::vector<AOBasis *> basis;

    void readIntgrlFile(std::iostream &ifs);
    void readContractionBlock(std::iostream &ifs, AOBasis &basis, int l);
    void readAtomBlock(std::iostream &ifs);
    void readAtomData(std::iostream &ifs, int n_atoms, double z);
};

#endif // INTGRL_H
