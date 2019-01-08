#ifndef CONFIG_H
#define CONFIG_H
//Macros para as constantes do problema, todas estão em unidades do sistema internacional

static const unsigned int NUM_X_CELLS = 60;
static const unsigned int NUM_Y_CELLS = 80;
static const unsigned int NUM_GHOST_CELLS = 4;
static const unsigned int NUM_RK_STEPS = 3;
static const double STOPPING_TIME = 2.0;
static const unsigned int MAX_ITERATIONS = 100000;
static const double MAX_X_SIDE = 1.0;
static const double MIN_X_SIDE = -1.0;
static const double MAX_Y_SIDE = 1.0;
static const double MIN_Y_SIDE = -1.0;
static const double ADVECTION_VEL_X = 1.0;
static const double ADVECTION_VEL_Y = 1.0;
static const double COURANT_NUM = 0.5;
#endif
