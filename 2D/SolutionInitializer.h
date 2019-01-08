#ifndef SOLUTIONINITIALIZER_H
#define SOLUTIONINITIALIZER_H

#include <cmath>
#include "boost/multi_array.hpp"

#include "config.h"
#include "Class_Cell.h"

void SolutionInitializerSquareWave(boost::multi_array<Cell,2> &cells)
{
    for(int i = 0; i < cells.shape()[0]; i++)
    {
        for(int j = 0; j < cells.shape()[1]; j++)
        {
            if(cells[i][j].cx > -0.25 && cells[i][j].cx < 0.25 && cells[i][j].cy > -0.25 && cells[i][j].cy < 0.25)
            {
                cells[i][j].u[0] = 1.0;
            }
            else
            {
                cells[i][j].u[0] = 0.0;
            }
        }
    }
}

void SolutionInitializerSineWave(boost::multi_array<Cell,2> &cells)
{
    for(int i = 0; i < cells.shape()[0]; i++)
    {
        for(int j = 0; j < cells.shape()[1]; j++)
        {
            cells[i][j].u[0]= sin((cells[i][j].cx + cells[i][j].cy)*M_PI) / double(2);
        }
    }
}

#endif
