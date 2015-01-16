#ifndef MRENV_H
#define MRENV_H

class MREnv {
public:
    static void initializeMRCPP();
    static void initializeMRCPP(int argc, char **argv, const char *fname = 0);

private:
    static void initializeTrees(int k, int depth, double prec, int type);

};
#endif // MRENV_H
