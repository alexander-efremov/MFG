#ifndef FEM_CIRCLE_CONSTS_H
#define FEM_CIRCLE_CONSTS_H


#if defined _WIN32 || defined __CYGWIN__ || defined WIN32
    extern long timerStart;
#else
    extern struct timeval timerStart;
#endif

extern double SIGMA;
extern double TAU;
extern double A;
extern double B;
extern unsigned int NX;
extern unsigned int NX_1;
extern int TIME_STEP_CNT;
extern double HX;
extern double U;

#endif //FEM_CIRCLE_CONSTS_H
