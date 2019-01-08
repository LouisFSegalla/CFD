#ifndef SOLUTIONWRITER_H
#define SOLUTIONWRITER_H

#include<iostream>
#include<fstream>
#include<string>
#include "boost/multi_array.hpp"


#include "Class_Cell.h"
#include "config.h"

void WriteSolution(boost::multi_array<Cell,2> &cells, double &time, std::string name)
{
    std::ofstream out( name , std::ofstream::out);
        
    out << "Num X cells: " << NUM_X_CELLS << std::endl;
    out << "Num Y cells: " << NUM_Y_CELLS << std::endl;
    out << "Num ghost cells: " << NUM_GHOST_CELLS << std::endl;
    out << "Time: " << time << std::endl;
        
    for(int i = 0; i < cells.shape()[0]; i++)
    {
        for(int j = 0; j < cells.shape()[1]; j++)
        {
            out << cells[i][j].cx << " " << cells[i][j].cy << " " <<cells[i][j].u[0] << std::endl;
        }
    }
    
}


#endif
