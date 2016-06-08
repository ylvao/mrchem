#ifndef ORBITALSET_H
#define ORBITALSET_H

#include <Eigen/Core>
#include <string>
#include <vector>

//#include "NonlinearMaximizer.h"
#include "Orbital.h"

//class OrbExp;

class OrbitalSet {
public:
    OrbitalSet(int n_orbs);
    OrbitalSet(int n_alpha, int n_beta);
    OrbitalSet(int ne, int mult, bool rest);
    OrbitalSet(const OrbitalSet &orb_set);
    virtual ~OrbitalSet();

    void push_back(int n_orbs, int occ, int spin);
    void clear();

//    OrbitalSet& operator=(const OrbitalSet &orb_set);

//    void writeOrbitals(const std::string &of);
//    void readOrbitals(const OrbitalSet &orbs);
//    void readOrbitals(const std::string &of);
//    void readOrbitals(const std::string &bf,
//                      const std::string &mo);
//    void readOrbitals(const std::string &bf,
//                      const std::string &mo_a,
//                      const std::string &mo_b);
//    void readVirtuals(const std::string &bf,
//                      const std::string &mo,
//                      int n_occ);

//    void normalize();
//    void orthogonalize(double prec);
//    void orthogonalize(double prec, OrbitalSet &orbitals);
//    Eigen::MatrixXd orthonormalize(double prec, Eigen::MatrixXd *F = 0);
//    Eigen::MatrixXd localize(double prec, Eigen::MatrixXd *F = 0);
//    Eigen::MatrixXd diagonalize(double prec, Eigen::MatrixXd *F);

//    Eigen::MatrixXd calcOverlapMatrix();
//    Eigen::MatrixXd calcOverlapMatrix(OrbitalSet &orbitals);
//    Eigen::MatrixXd calcLocalizationMatrix();
//    Eigen::MatrixXd calcOrthonormalizationMatrix();
//    Eigen::MatrixXd calcDiagonalizationMatrix(Eigen::MatrixXd F);

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

//    void rotate(double prec, const Eigen::MatrixXd &U);
//    void add(double a, OrbitalSet &set_1, double b, OrbitalSet &set_b);
//    void addInPlace(OrbitalSet &orb_set, double c = 1.0);
//    void crop(double prec = -1.0);

    const Orbital *getOrbitalPtr(int i) const { return this->orbitals[i]; }
    Orbital *getOrbitalPtr(int i) { return this->orbitals[i]; }

    const Orbital &getOrbital(int i) const;
    Orbital &getOrbital(int i);

//    void replaceOrbitals(OrbitalSet &new_orbs);
//    void replaceOrbital(int i, Orbital **orb);

//    int printTreeSizes() const;

    friend std::ostream& operator<<(std::ostream &o, OrbitalSet &orb_set) {
        int oldPrec = TelePrompter::setPrecision(15);
        o << "*OrbitalSet:  ";
        o << std::setw(4) << orb_set.size()          << " orbitals   ";
        o << std::setw(4) << orb_set.getNOccupied()  << " occupied   ";
        o << std::setw(4) << orb_set.getNElectrons() << " electrons" << std::endl;
        o << "------------------------------";
        o << "------------------------------\n";
        o << "   n    Norm                 Occ Spin  Error\n";
        o << "------------------------------";
        o << "------------------------------\n";
        for (int i = orb_set.size(); i > 0; i--) {
            Orbital *orb = orb_set.getOrbitalPtr(i-1);
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
//    OrbExp* readOrbitalExpansion(const std::string &bf, const std::string &mo);

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
//    void separateSpinOrbitals(OrbitalSet &phi_p,
//                              OrbitalSet &phi_a,
//                              OrbitalSet &phi_b);
//    void collectSpinOrbitals(OrbitalSet &phi_p,
//                             OrbitalSet &phi_a,
//                             OrbitalSet &phi_b);
private:
    //Data
    std::vector<Orbital *> orbitals;
};

/** subclass which defines the particular Gradient and Hessian
 * and other specific functions for a maximization of
 * f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 * The resulting transformation includes the orthonormalization of the orbitals.
 * For details see the tex documentation in doc directory
 *
 */
//class RR : public NonlinearMaximizer {
//public:
//    RR(OrbitalSet &orbitals);//make the matrices <i|R_x|j>,<i|R_y|j>,<i|R_z|j>
//    const Eigen::MatrixXd &getTotalU() const { return this->total_U; }
//protected:
//    int N;//number of orbitals
//    Eigen::MatrixXd r_i_orig;//<i|R_x|j>,<i|R_y|j>,<i|R_z|j>
//    Eigen::MatrixXd r_i ;// rotated  r_i_orig
//    Eigen::MatrixXd total_U;// the rotation matrix of the orbitals

//    //NB:total_U is not Unitary if the basis set is not orthonormal
//    double functional(); //the functional to maximize
//    double make_gradient();
//    double make_hessian();
//    void do_step(Eigen::VectorXd step);
//};

#endif // ORBITALSET_H
