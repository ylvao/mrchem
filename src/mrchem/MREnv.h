#ifndef MRENV_H
#define MRENV_H

#include <Eigen/Core>

#include "parallel.h"
#include "constants.h"
#include "Timer.h"

class MREnv {
public:
    static void initializeMRCPP(int argc, char **argv);
    static void finalizeMRCPP(const Timer t);
    static void initializeMRA();
};

#endif // MRENV_H
