#ifndef DENSITY_H
#define DENSITY_H

#include "FunctionTree.h"

class DensExp;

class Density : public FunctionTree<3> {
public:
    Density() { }
    Density(DensExp &orbs);

    virtual ~Density() { }

    Density(const Density &dens);
    Density &operator =(const Density &dens);

    void setNumberDensity(double numDens);

    bool saveTree(const std::string &file);
    bool loadTree(const std::string &file);
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar.register_type(static_cast<Density *>(NULL));
        ar & boost::serialization::base_object<FunctionTree<3> >(*this);
    }
};

#endif // DENSITY_H
