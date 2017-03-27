#include "consts.h"

#if defined _WIN32 || defined __CYGWIN__ || defined WIN32
long timerStart = 0;
#else

#endif

double SIGMA = 0.;
double SIGMA_SQ = 0.;
double TAU = 0.;
double A = 0.;
double B = 0.;
unsigned int NX = 0;
unsigned int N_1 = 0;
int TIME_STEP_CNT = 0;
double H = 0.;
double H_SQ = 0.;
double ALPHA_COEF = 0.;

