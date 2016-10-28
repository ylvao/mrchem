#ifndef ORBITALVECTOR_H
#define ORBITALVECTOR_H

#include <Eigen/Core>
#include <string>
#include <vector>

//#include "NonlinearMaximizer.h"
#include "Orbital.h"

class OrbitalVector {
public:
    OrbitalVector(int n_orbs);
    OrbitalVector(int n_alpha, int n_beta);
    OrbitalVector(int ne, int mult, bool rest);
    OrbitalVector(const OrbitalVector &orb_set);
    OrbitalVector& operator=(const OrbitalVector &orb_set);
    virtual ~OrbitalVector();

    void push_back(int n_orbs, int occ, int spin);
    void clear(bool free = true);

//    void writeOrbitals(const std::string &of);
//    void readOrbitals(const std::string &of);

    void normalize();
    void orthogonalize();
    void orthogonalize(OrbitalVector &phi);

    Eigen::MatrixXcd calcOverlapMatrix();
    Eigen::MatrixXcd calcOverlapMatrix(OrbitalVector &ket);
    Eigen::MatrixXcd calcOverlapMatrix_P(OrbitalVector &ket);
    Eigen::MatrixXcd calcOverlapMatrix_P_H(OrbitalVector &ket);

    int size() const { return this->orbitals.size(); }
    int getNOccupied() const;
    int getNEmpty() const;
    int getNSingly() const;
    int getNDoubly() const;
    int getNPaired() const;
    int getNAlpha() const;
    int getNBeta() const;
    int getNElectrons(int spin = Paired) const;
    int getMultiplicity() const;

    void setSpins(const Eigen::VectorXi &spins);
    void setErrors(const Eigen::VectorXd &errors);
    void setOccupancies(const Eigen::VectorXi &occ);

    Eigen::VectorXi getSpins() const;
    Eigen::VectorXd getNorms() const;
    Eigen::VectorXd getErrors() const;
    Eigen::VectorXd getSquareNorms() const;
    Eigen::VectorXi getOccupancies() const;

    bool isConverged(double prec) const;
    double calcTotalError() const;

    const Orbital *getOrbitalPtr(int i) const { return this->orbitals[i]; }
    Orbital *getOrbitalPtr(int i) { return this->orbitals[i]; }

    const Orbital &getOrbital(int i) const;
    Orbital &getOrbital(int i);

//    void replaceOrbitals(OrbitalVector &new_orbs);
    void replaceOrbital(int i, Orbital **orb);
    void replaceTrees(int i, Orbital *orb);

    int printTreeSizes() const;

    friend std::ostream& operator<<(std::ostream &o, OrbitalVector &orb_set) {
        int oldPrec = TelePrompter::setPrecision(15);
        o << "*OrbitalVector: ";
        o << std::setw(4) << orb_set.size()          << " orbitals  ";
        o << std::setw(4) << orb_set.getNOccupied()  << " occupied  ";
        o << std::setw(4) << orb_set.getNElectrons() << " electrons " << std::endl;
        o << "------------------------------";
        o << "------------------------------\n";
        o << "   n    sqNorm               Occ Spin  Error\n";
        o << "------------------------------";
        o << "------------------------------\n";
        for (int i = 0; i < orb_set.size(); i++) {
            Orbital *orb = orb_set.getOrbitalPtr(i);
            o << std::setw(4) << i;
            if (orb != 0) {
                o << *orb;
            } else {
                o << std::endl;
            }
        }
        o << "------------------------------";
        o << "------------------------------\n";
        TelePrompter::setPrecision(oldPrec);
        return o;
    }
protected:
//    Eigen::MatrixXd getSpinMatrix() const;
//    void spinCleanMatrix(Eigen::MatrixXd &M) const;

//    void separateSpinMatrix(Eigen::MatrixXd &F,
//                            Eigen::MatrixXd &P,
//                            Eigen::MatrixXd &A,
//                            Eigen::MatrixXd &B) const;
//    void collectSpinMatrix(Eigen::MatrixXd &F,
//                           Eigen::MatrixXd &P,
//                           Eigen::MatrixXd &A,
//                           Eigen::MatrixXd &B) const;
//    void separateSpinOrbitals(OrbitalVector &phi_p,
//                              OrbitalVector &phi_a,
//                              OrbitalVector &phi_b);
//    void collectSpinOrbitals(OrbitalVector &phi_p,
//                             OrbitalVector &phi_a,
//                             OrbitalVector &phi_b);
private:
    //Data
    std::vector<Orbital *> orbitals;
};

#endif // ORBITALVECTOR_H
