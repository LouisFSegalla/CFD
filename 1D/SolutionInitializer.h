#ifndef SOLUTIONINITIALIZER_H
#define SOLUTIONINITIALIZER_H

#include <cmath>

#include "config.h"
#include "Class_Cell.h"

void SolutionInitializerSquareWave(Cell *cells)
{
    for(int i = 0; i < NUM_X_CELLS + 2*NUM_GHOST_CELLS; i++)
    {
        if(cells[i].cx > -0.25 && cells[i].cx < 0.25)
        {
            cells[i].u[0] = 1.0;
        }
        else
        {
            cells[i].u[0] = 0.0;
        }
    }
}

void SolutionInitializerSineWave(Cell *cells)
{
    for(int i = 0; i < NUM_X_CELLS + 2*NUM_GHOST_CELLS; i++)
    {
        cells[i].u[0] = sin( M_PI * cells[i].cx);
    }
}
#endif
