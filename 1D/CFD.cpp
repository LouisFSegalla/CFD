#include <iostream>
#include <cstdlib>
#include <algorithm>

using namespace std;

#include "CFD.h"

void GhostCellsUpdater(Cell *cells, int &RKstep)
{
    for(int GhostCells = 0; GhostCells < NUM_GHOST_CELLS; GhostCells++)
    {
        cells[GhostCells].u[RKstep] = cells[NUM_X_CELLS + GhostCells].u[RKstep];
        
        cells[NUM_X_CELLS + NUM_GHOST_CELLS + GhostCells].u[RKstep] = cells[NUM_GHOST_CELLS + GhostCells].u[RKstep];
    }
    
}

void ReconstructVariablesFirstOrder(Cell *cells, int &RKstep)
{
    for(int i = NUM_GHOST_CELLS - 1; i < NUM_X_CELLS + NUM_GHOST_CELLS + 1; i++)
    {
        cells[i].uEast = cells[i].u[RKstep];
        cells[i].uWest = cells[i].u[RKstep];
    }
}

void ReconstructVariablesBeamWarming(Cell *cells, int &RKstep)
{
    double dudx = 0;
    for(int i = NUM_GHOST_CELLS - 1; i < NUM_X_CELLS + NUM_GHOST_CELLS + 1; i++)
    {
        dudx = (cells[i].u[RKstep] - cells[i-1].u[RKstep]) / double(cells[i].dx);
        cells[i].uEast = cells[i].u[RKstep] + dudx*cells[i].dx /double(2);
        cells[i].uWest = cells[i].u[RKstep] - dudx*cells[i].dx /double(2);
    }
}

void ReconstructVariablesLaxWendroff(Cell *cells, int &RKstep)
{
    double dudx = 0;
    for(int i = NUM_GHOST_CELLS - 1; i < NUM_X_CELLS + NUM_GHOST_CELLS + 1; i++)
    {
        dudx = (cells[i+1].u[RKstep] - cells[i].u[RKstep]) / double(cells[i].dx);
        cells[i].uEast = cells[i].u[RKstep] + dudx*cells[i].dx /double(2);
        cells[i].uWest = cells[i].u[RKstep] - dudx*cells[i].dx /double(2);
    }
}

void ReconstructVariablesLimitedLW(Cell *cells, int &RKstep)
{
    double dudx = 0;
    double r = 0;
    double phi = 0;
    for(int i = NUM_GHOST_CELLS - 1; i < NUM_X_CELLS + NUM_GHOST_CELLS + 1; i++)
    {
        r = (cells[i].u[RKstep] - cells[i-1].u[RKstep]) / double(cells[i+1].u[RKstep] - cells[i].u[RKstep] + 1e-6);
        r = std::max(0.0,r);
        phi = (r*r + r) / double(r*r + 1.0);
        dudx = (cells[i+1].u[RKstep] - cells[i].u[RKstep]) / double(cells[i].dx);
        cells[i].uEast = cells[i].u[RKstep] + phi*dudx*cells[i].dx /double(2);
        cells[i].uWest = cells[i].u[RKstep] - phi*dudx*cells[i].dx /double(2);
    }
}

void ReconstructVariablesFromm(Cell *cells, int &RKstep)
{
    double dudx = 0;
    for(int i = NUM_GHOST_CELLS - 1; i < NUM_X_CELLS + NUM_GHOST_CELLS + 1; i++)
    {
        dudx = (cells[i+1].u[RKstep] - cells[i-1].u[RKstep]) / double(2*cells[i].dx);
        cells[i].uEast = cells[i].u[RKstep] + dudx*cells[i].dx /double(2);
        cells[i].uWest = cells[i].u[RKstep] - dudx*cells[i].dx /double(2);
    }
}

void ReconstructVariablesWENO5thOrder(Cell *cells, int &RKstep, double &dt)
{
    double pointMinus2 = 0, pointMinus1 = 0, point = 0, pointPlus1 = 0, pointPlus2 = 0;
//     double a = 0, u = 0, v = 0;
    double p0 = 0, p1 = 0, p2 = 0;
    double w0 = 0, w1 = 0, w2 = 0;
    double alpha0 = 0, alpha1 = 0, alpha2 = 0, alphaSum = 0;
    double IS0 = 0, IS1 = 0, IS2 = 0;
    double C0 = 0.1, C1 = 0.6, C2 = 0.3;
    double e = 1e-6;
    
//     double dudx = 0;
//     double r = 0;
//     double phi = 0;
    
    double dflux[5];
    
    for(int i = NUM_GHOST_CELLS-1; i < NUM_X_CELLS + NUM_GHOST_CELLS + 1; i++)
    {
     
        pointMinus2 = cells[i-2].u[RKstep];
        pointMinus1 = cells[i-1].u[RKstep];
        point= cells[i].u[RKstep];
        pointPlus1 = cells[i+1].u[RKstep];
        pointPlus2 = cells[i+2].u[RKstep];
        
        //Polinômios para o lado uEast
        p0 = ((2*pointPlus2) - (7*pointPlus1) + (11*point))/float(6);
        p1 = ((-1*pointPlus1) + (5*point) + (2*pointMinus1))/float(6);
        p2 = ((2*point) + (5*pointMinus1) - (pointMinus2))/float(6);
        
        //Smoothness indicator
        IS0 = (13/12)*(pointMinus2 - (2*pointMinus1) + point)*(pointMinus2 - (2*pointMinus1) + point) 
              + 0.25*(pointMinus2 - (4*pointMinus1) + 3*point)*(pointMinus2 - (4*pointMinus1) + 3*point);
        
        IS1 = (13/12)*(pointMinus1 - (2*point) + pointPlus1)*(pointMinus1 - (2*point) + pointPlus1)
              + 0.25*(pointMinus1-pointPlus1)*(pointMinus1-pointPlus1);
              
        IS2 = (13/12)*(point - (2*pointPlus1) + pointPlus2)*(point - (2*pointPlus1) + pointPlus2)
              + 0.25*(3*point - (4*pointPlus1) + pointPlus2)*(3*point - (4*pointPlus1) + pointPlus2);
              
        alpha0 = C0/double((e + IS0)*(e + IS0));
        alpha1 = C1/double((e + IS1)*(e + IS1));
        alpha2 = C2/double((e + IS2)*(e + IS2));
        alphaSum = alpha0 + alpha1 + alpha2;
        
        w0 = alpha0 / double(alphaSum);
        w1 = alpha1 / double(alphaSum);
        w2 = alpha2 / double(alphaSum);
        
        cells[i].uEast = (w0*p0 + w1*p1 + w2*p2);
                
        //Polinômios para o lado uWest
        p0 = ((-1*pointPlus2) + (5*pointPlus1) + (2*point))/float(6);
        p1 = ((2*pointPlus1) + (5*point) - pointMinus1)/float(6);
        p2 = ((11*point) - (7*pointMinus1) + (2*pointMinus2))/float(6);
        
        //Smoothness indicator
        IS0 = (13/12)*(pointMinus2 - (2*pointMinus1) + point)*(pointMinus2 - (2*pointMinus1) + point) 
              + 0.25*(pointMinus2 - (4*pointMinus1) + 3*point)*(pointMinus2 - (4*pointMinus1) + 3*point);
        
        IS1 = (13/12)*(pointMinus1 - (2*point) + pointPlus1)*(pointMinus1 - (2*point) + pointPlus1)
              + 0.25*(pointMinus1-pointPlus1)*(pointMinus1-pointPlus1);
              
        IS2 = (13/12)*(point - (2*pointPlus1) + pointPlus2)*(point - (2*pointPlus1) + pointPlus2)
              + 0.25*(3*point - (4*pointPlus1) + pointPlus2)*(3*point - (4*pointPlus1) + pointPlus2);
        
        alpha0 = C2/double((e + IS0)*(e + IS0));
        alpha1 = C1/double((e + IS1)*(e + IS1));
        alpha2 = C0/double((e + IS2)*(e + IS2));
        alphaSum = alpha0 + alpha1 + alpha2;
        
        w0 = alpha0 / double(alphaSum);
        w1 = alpha1 / double(alphaSum);
        w2 = alpha2 / double(alphaSum);
        
        cells[i].uWest = (w0*p0 + w1*p1 + w2*p2);
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

void UpdateCellAveragesEuler(Cell *cells, int &RKstep, double &dt)
{
    for(int i = NUM_GHOST_CELLS; i < NUM_GHOST_CELLS + NUM_X_CELLS; i++)
    {
        cells[i].u[RKstep+1] = cells[i].u[RKstep] + dt/cells[i].dx * cells[i].TotalFlux; 
    }
}

void UpdateCellAveragesRK2(Cell *cells, int &RKstep, double &dt)
{
    for(int i = NUM_GHOST_CELLS; i < NUM_GHOST_CELLS + NUM_X_CELLS; i++)
    {
        if(RKstep == 0)
        {
            cells[i].u[RKstep+1] = cells[i].u[RKstep] + dt/cells[i].dx * cells[i].TotalFlux; 
        }
        if(RKstep == 1)
        {
            cells[i].u[RKstep+1] = 0.5*(cells[i].u[RKstep-1] + cells[i].u[RKstep] + dt/cells[i].dx * cells[i].TotalFlux);
        }
    }
}

void UpdateCellAveragesRK3(Cell *cells, int &RKstep, double &dt)
{
    for(int i = NUM_GHOST_CELLS; i < NUM_GHOST_CELLS + NUM_X_CELLS; i++)
    {
        if(RKstep == 0)
        {
            cells[i].u[RKstep+1] = cells[i].u[RKstep] + dt/cells[i].dx * cells[i].TotalFlux; 
        }
        if(RKstep == 1)
        {
            cells[i].u[RKstep+1] = 0.25*(3*cells[i].u[RKstep-1] + cells[i].u[RKstep] + dt/cells[i].dx * cells[i].TotalFlux);
        }
        if(RKstep == 2)
        {
            cells[i].u[RKstep+1] = (cells[i].u[RKstep-2] + 2*cells[i].u[RKstep] + 2*dt/cells[i].dx * cells[i].TotalFlux)/float(3);
        }
    }
}

void CopyVariables(Cell *cells)
{
    for(int i = NUM_GHOST_CELLS; i < NUM_X_CELLS + NUM_GHOST_CELLS + 1; i++)
    {
        cells[i].u[0] = cells[i].u[NUM_RK_STEPS];
    }
}
