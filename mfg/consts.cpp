#include "consts.h"

#if defined _WIN32 || defined __CYGWIN__ || defined WIN32
long timerStart = 0;
#else

#endif

double SIGMA = 0.;
double TAU = 0.;
double A = 0.;
double B = 0.;
unsigned int NX = 0;
unsigned int NX_1 = 0;
int TIME_STEP_CNT = 0;
double HX = 0.;
double U = 0.;

