#include "consts.h"

#if defined _WIN32 || defined __CYGWIN__ || defined WIN32
long timerStart = 0;
#else
//struct timeval timerStart;
#endif

