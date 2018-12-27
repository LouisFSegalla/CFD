#ifndef CLASS_CELL_H
#define CLASS_CELL_H

#include "config.h"

class Cell
{
public:
    double u[NUM_RK_STEPS + 1];//talvez não dẽ
    double uEast;
    double uWest;
    double cx;
    double dx;
    double TotalFlux;
    
};

#endif
