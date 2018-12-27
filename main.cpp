#include<iostream>
#include<cstdlib>

using namespace std;

#include "Class_Cell.h"
#include "config.h"
#include "SolutionInitializer.h"
#include "SolutionWriter.h"
#include "CFD.h"


int main()
{
    Cell *cells = new Cell[NUM_X_CELLS + 2*NUM_GHOST_CELLS];//Array de objetos da classe Cell
    
    int sizeArray = NUM_X_CELLS + 2*NUM_GHOST_CELLS;
    double dx = (MAX_X_SIDE - MIN_X_SIDE) / NUM_X_CELLS;
    double time = 0.0;
    
    for(int i = 0; i < sizeArray; i++)
    {
        cells[i].dx = dx;
        cells[i].cx = MIN_X_SIDE + dx*( i + 0.5 - NUM_GHOST_CELLS );
    }
    
    SolutionInitializerSineWave(cells);
        
    std::string name = "Solução no tempo " + std::to_string(time) + ".txt";
    WriteSolution(cells,time,name);
    
    bool LastTimeStep = false;
    
    double dt = (cells[0].dx / abs(ADVECTION_VEL)) * COURANT_NUM;
    
    for(int timeIteration = 0; timeIteration < MAX_ITERATIONS; timeIteration++)
    {
        if(time + dt > STOPPING_TIME)
        {
            dt = STOPPING_TIME - time;
            LastTimeStep = true;
        }
        
        for(int RKsteps = 0; RKsteps < NUM_RK_STEPS; RKsteps++)
        {
            GhostCellsUpdater(cells,RKsteps);
            ReconstructVariables(cells,RKsteps);
            CalculateFlux(cells,RKsteps);
            UpdateCellAverages(cells,RKsteps,dt);
        }
        
        CopyVariables(cells);
        
        time += dt;
        
        if(LastTimeStep){break;}

        if( (timeIteration) % 5 == 0.0)
        {
            name = "Solução no tempo " + std::to_string(time) + ".txt";
            WriteSolution(cells,time,name);
        }
        
    }
    name = "Solução no tempo " + std::to_string(time) + ".txt";
    WriteSolution(cells,time,name);
    
    return 0;
}
