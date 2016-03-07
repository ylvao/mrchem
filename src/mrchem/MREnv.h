#ifndef MRENV_H
#define MRENV_H

#include <Eigen/Core>

#include "parallel.h"
#include "constants.h"

#include "TelePrompter.h"

class MREnv {
public:
    static void initializeMRCPP();
    static void finalizeMRCPP(double t);
};

#endif // MRENV_H
