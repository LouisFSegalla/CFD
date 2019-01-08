#include<iostream>
#include<cstdlib>
#include "boost/multi_array.hpp"


using namespace std;
using namespace boost;


#include "Class_Cell.h"
#include "config.h"

#include "SolutionWriter.h"
#include "CFD.h"


int main()
{
    boost::multi_array<Cell,2> cells;
    cells.resize(extents[NUM_X_CELLS + 2*NUM_GHOST_CELLS][NUM_Y_CELLS + 2*NUM_GHOST_CELLS]);//Array de objetos da classe Cell
    
    double dx = (MAX_X_SIDE - MIN_X_SIDE) / NUM_X_CELLS;
    double dy = (MAX_Y_SIDE - MIN_Y_SIDE) / NUM_Y_CELLS;
    double time = 0.0;
      
    for(int i = 0; i < cells.shape()[0]; i++)
    {
        for(int j = 0; j < cells.shape()[1]; j++)
        {
            cells[i][j].dx = dx;
            cells[i][j].dy = dy;
            cells[i][j].cx = MIN_X_SIDE + dx*( i + 0.5 - NUM_GHOST_CELLS );
            cells[i][j].cy = MIN_Y_SIDE + dy*( j + 0.5 - NUM_GHOST_CELLS );
        }
    }
    
    SolutionInitializerSquareWave(cells);
        
    std::string name = "Solução no tempo " + std::to_string(time) + ".txt";
    WriteSolution(cells,time,name);
    
    bool LastTimeStep = false;
    
    
    //Time step calculation
    double dtx = (cells[0][0].dx / abs(ADVECTION_VEL_X));
    double dty = (cells[0][0].dy / abs(ADVECTION_VEL_Y));
    double dt = (1.0/float((1.0/float(dtx)) + (1.0/float(dty))))*COURANT_NUM;
       
    for(int RKsteps = 0; RKsteps < NUM_RK_STEPS; RKsteps++)
    {
        GhostCellsUpdater(cells,RKsteps);
        ReconstructVariablesLimitedLW(cells,RKsteps);
        CalculateFlux(cells,RKsteps);
        UpdateCellAveragesRK3(cells,RKsteps,dt);
    }
    CopyVariables(cells);
    time += dt;
    
    for(int timeIteration = 1; timeIteration < MAX_ITERATIONS; timeIteration++)
    {
        if(time + dt > STOPPING_TIME)
        {
            dt = STOPPING_TIME - time;
            LastTimeStep = true;
        }
        
        for(int RKsteps = 0; RKsteps < NUM_RK_STEPS; RKsteps++)
        {
            GhostCellsUpdater(cells,RKsteps);
            ReconstructVariablesWENO5thOrder(cells,RKsteps);
            CalculateFlux(cells,RKsteps);
            UpdateCellAveragesRK3(cells,RKsteps,dt);
        }
        CopyVariables(cells);
        time += dt;
        
        if(LastTimeStep){break;}

        if( (timeIteration) % 25 == 0.0)
        {
            name = "Solução no tempo " + std::to_string(time) + ".txt";
            WriteSolution(cells,time,name);
        }
        
    }
    name = "Solução no tempo " + std::to_string(time) + ".txt";
    WriteSolution(cells,time,name);
    ErrorCalculation(cells);
    return 0;
}
