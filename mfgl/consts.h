#ifndef FEM_CIRCLE_CONSTS_H
#define FEM_CIRCLE_CONSTS_H


#if defined _WIN32 || defined __CYGWIN__ || defined WIN32
    extern long timerStart;
#else
    extern struct timeval timerStart;
#endif


#endif //FEM_CIRCLE_CONSTS_H
