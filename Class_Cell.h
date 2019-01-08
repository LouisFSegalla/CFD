#ifndef CLASS_CELL_H
#define CLASS_CELL_H

#include "config.h"

class Cell
{
public:
    double u[NUM_RK_STEPS + 1];//talvez não dẽ
    double uEast;
    double uWest;
    double uNorth;
    double uSouth;
    double cx;
    double cy;
    double dx;
    double dy;
    double TotalFlux;
};

#endif
