#ifndef MRENV_H
#define MRENV_H

class MREnv {
public:
    static void initializeMRCPP(int k, double prec);
    static void finalizeMRCPP(double t);

private:
    static void initializeTrees(int k, int depth, double prec, int type);

};
#endif // MRENV_H
