#ifndef MFG_CONSTS_H
#define MFG_CONSTS_H

#if defined _WIN32 || defined __CYGWIN__ || defined WIN32
    extern long timerStart;
#else
    extern struct timeval timerStart;
#endif

extern double SIGMA;
extern double SIGMA_SQ;
extern double TAU;
extern double A;
extern double B;
extern unsigned int NX;
extern unsigned int N_1;
extern int TIME_STEP_CNT;
extern double H;
extern double H_SQ;
extern double ALPHA_COEF;

#endif //MFG_CONSTS_H
