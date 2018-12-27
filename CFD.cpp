#include<cstdlib>

#include "CFD.h"

void GhostCellsUpdater(Cell *cells, int &RKstep)
{
    for(int GhostCells = 0; GhostCells < NUM_GHOST_CELLS; GhostCells++)
    {
        cells[GhostCells].u[RKstep] = cells[NUM_X_CELLS + GhostCells].u[RKstep];
        
        cells[NUM_X_CELLS + NUM_GHOST_CELLS + GhostCells].u[RKstep] = cells[NUM_GHOST_CELLS + GhostCells].u[RKstep];
    }
    
}

void ReconstructVariables(Cell *cells, int &RKstep)
{
    for(int i = NUM_GHOST_CELLS - 1; i < NUM_X_CELLS + NUM_GHOST_CELLS + 1; i++)
    {
        cells[i].uEast = cells[i].u[RKstep];
        cells[i].uWest = cells[i].u[RKstep];
    }
}

void CalculateFlux(Cell *cells, int &RKstep)
{
    for(int i = 0; i < NUM_X_CELLS + 2*NUM_GHOST_CELLS; i++)
    {
        cells[i].TotalFlux = 0.0;
    }
    
    double LeftValue = 0;
    double RightValue = 0;
    double Flux = 0;
    for(int interface = NUM_GHOST_CELLS; interface < NUM_GHOST_CELLS + NUM_X_CELLS + 1; interface++)
    {
        LeftValue = cells[interface-1].uEast;
        RightValue = cells[interface].uWest;
        
        //Lax-Friedrich
        Flux = 0.5*ADVECTION_VEL*(LeftValue + RightValue) - 0.5*abs(ADVECTION_VEL)*(RightValue - LeftValue);
        
        cells[interface - 1].TotalFlux -= Flux;
        cells[interface].TotalFlux += Flux;
        
    }
}

void UpdateCellAverages(Cell *cells, int &RKstep, double &dt)
{
    for(int i = NUM_GHOST_CELLS; i < NUM_GHOST_CELLS + NUM_X_CELLS; i++)
    {
        cells[i].u[RKstep+1] = cells[i].u[RKstep] + dt/cells[i].dx * cells[i].TotalFlux; 
    }
}

void CopyVariables(Cell *cells)
{
    for(int i = NUM_GHOST_CELLS; i < NUM_X_CELLS + NUM_GHOST_CELLS + 1; i++)
    {
        cells[i].u[0] = cells[i].u[NUM_RK_STEPS];
    }
}
