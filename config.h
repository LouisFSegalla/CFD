#ifndef CONFIG_H
#define CONFIG_H
//Macros para as constantes do problema, todas estão em unidades do sistema internacional

static const unsigned int NUM_X_CELLS = 100;
static const unsigned int NUM_GHOST_CELLS = 3;
static const unsigned int NUM_RK_STEPS = 1;
static const double STOPPING_TIME = 2.0;
static const unsigned int MAX_ITERATIONS = 1000000;
static const double MAX_X_SIDE = 1.0;
static const double MIN_X_SIDE = -1.0;
static const double ADVECTION_VEL = 1.0;
static const double COURANT_NUM = 1.0;
#endif
