#ifndef SOLUTIONWRITER_H
#define SOLUTIONWRITER_H

#include<iostream>
#include<fstream>
#include<string>

#include "Class_Cell.h"
#include "config.h"

void WriteSolution(Cell *cells, double &time, std::string name)
{
    std::ofstream out( name , std::ofstream::out);
    
    out << "Num X cells: " << NUM_X_CELLS << std::endl;
    out << "Num ghost cells: " << NUM_GHOST_CELLS << std::endl;
    out << "Time: " << time << std::endl;
        
    for(int i = 0; i < NUM_X_CELLS + 2*NUM_GHOST_CELLS; i++)
    {
        out << cells[i].cx <<  " " << cells[i].u[0] << std::endl;
    }
    
}


#endif
