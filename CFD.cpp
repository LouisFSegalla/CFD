#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include "boost/multi_array.hpp"

using namespace std;
using namespace boost;

#include "CFD.h"

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

//deixando dessa maneira parece que diminui o erro, mas não sei se é correto normalizar assim 
void SolutionInitializerCircle(boost::multi_array<Cell,2> &cells, const int &raio)
{
    double centroX = cells.shape()[0] / double(2);
    double centroY = cells.shape()[1] / double(2);
    for(int i = 0; i < cells.shape()[0]; i++)
    {
        for(int j = 0; j < cells.shape()[1]; j++)
        {
            if( (i - centroX)*(i - centroX) + (j - centroY)*(j - centroY) - raio*raio  > 0)
            {
                cells[i][j].u[0] = 1;
            }
            else if((i - centroX)*(i - centroX) + (j - centroY)*(j - centroY) - raio*raio  == 0)
            {
                cells[i][j].u[0] = 0;
            }
            else
            {
                cells[i][j].u[0] = -1;
            }
        }
    }
}

void GhostCellsUpdater(boost::multi_array<Cell,2> &cells, int &RKstep)
{
    //left and right side
    for(int j = 0; j < cells.shape()[1]; j++)
    {
        for(int GhostCells = 0; GhostCells < NUM_GHOST_CELLS; GhostCells++)
        {
            cells[GhostCells][j].u[RKstep] = cells[NUM_X_CELLS + GhostCells][j].u[RKstep];
            
            cells[NUM_X_CELLS + NUM_GHOST_CELLS + GhostCells][j].u[RKstep] = cells[NUM_GHOST_CELLS + GhostCells][j].u[RKstep];
        }
    }
    
    //bottom and top
    for(int i = 0; i < cells.shape()[0]; i++)
    {
        for(int GhostCells = 0; GhostCells < NUM_GHOST_CELLS; GhostCells++)
        {
            cells[i][GhostCells].u[RKstep] = cells[i][NUM_Y_CELLS + GhostCells].u[RKstep];
            
            cells[i][NUM_Y_CELLS + NUM_GHOST_CELLS + GhostCells].u[RKstep] = cells[i][NUM_GHOST_CELLS + GhostCells].u[RKstep];
        }
    }
}

void ReconstructVariablesFirstOrder(boost::multi_array<Cell,2> &cells, int &RKstep)
{
    for(int i = NUM_GHOST_CELLS - 1; i < NUM_X_CELLS + NUM_GHOST_CELLS + 1; i++)
    {
        for(int j = NUM_GHOST_CELLS - 1; j < NUM_Y_CELLS + NUM_GHOST_CELLS + 1; j++)
        {
            cells[i][j].uEast = cells[i][j].u[RKstep];
            cells[i][j].uWest = cells[i][j].u[RKstep];
            cells[i][j].uNorth = cells[i][j].u[RKstep];
            cells[i][j].uSouth = cells[i][j].u[RKstep];
        }
    }
}

void ReconstructVariablesBeamWarming(boost::multi_array<Cell,2> &cells, int &RKstep)
{
    double dudx = 0;
    double dudy = 0;
    for(int i = NUM_GHOST_CELLS - 1; i < NUM_X_CELLS + NUM_GHOST_CELLS + 1; i++)
    {
        for(int j = NUM_GHOST_CELLS - 1; j < NUM_Y_CELLS + NUM_GHOST_CELLS + 1; j++)
        {
            dudx = (cells[i][j].u[RKstep] - cells[i-1][j].u[RKstep]) / double(cells[i][j].dx);
            dudy = (cells[i][j].u[RKstep] - cells[i][j-1].u[RKstep]) / double(cells[i][j].dy);
            cells[i][j].uEast = cells[i][j].u[RKstep] + dudx*cells[i][j].dx /double(2);
            cells[i][j].uWest = cells[i][j].u[RKstep] - dudx*cells[i][j].dx /double(2);
            cells[i][j].uNorth = cells[i][j].u[RKstep] + dudy*cells[i][j].dy /double(2);
            cells[i][j].uSouth = cells[i][j].u[RKstep] - dudy*cells[i][j].dy /double(2);
        }
    }
}

void ReconstructVariablesLaxWendroff(boost::multi_array<Cell,2> &cells, int &RKstep)
{
    double dudx = 0;
    double dudy = 0;
    for(int i = NUM_GHOST_CELLS - 1; i < NUM_X_CELLS + NUM_GHOST_CELLS + 1; i++)
    {
        for(int j = NUM_GHOST_CELLS - 1; j < NUM_Y_CELLS + NUM_GHOST_CELLS + 1; j++)
        {
            dudx = (cells[i+1][j].u[RKstep] - cells[i][j].u[RKstep]) / double(cells[i][j].dx);
            dudy = (cells[i][j+1].u[RKstep] - cells[i][j].u[RKstep]) / double(cells[i][j].dy);
            cells[i][j].uEast = cells[i][j].u[RKstep] + dudx*cells[i][j].dx /double(2);
            cells[i][j].uWest = cells[i][j].u[RKstep] - dudx*cells[i][j].dx /double(2);
            cells[i][j].uNorth = cells[i][j].u[RKstep] + dudy*cells[i][j].dy /double(2);
            cells[i][j].uSouth = cells[i][j].u[RKstep] - dudy*cells[i][j].dy /double(2);
        }
    }
}

void ReconstructVariablesLimitedLW(boost::multi_array<Cell,2> &cells, int &RKstep)
{
    double dudx = 0;
    double dudy = 0;
    double rx = 0;
    double ry = 0;
    double phix = 0;
    double phiy = 0;
    for(int i = NUM_GHOST_CELLS - 1; i < NUM_X_CELLS + NUM_GHOST_CELLS + 1; i++)
    {
        for(int j = NUM_GHOST_CELLS - 1; j < NUM_Y_CELLS + NUM_GHOST_CELLS + 1; j++)
        {
            rx = (cells[i][j].u[RKstep] - cells[i-1][j].u[RKstep]) / double(cells[i+1][j].u[RKstep] - cells[i][j].u[RKstep] + 1e-6);
            ry = (cells[i][j].u[RKstep] - cells[i][j-1].u[RKstep]) / double(cells[i][j+1].u[RKstep] - cells[i][j].u[RKstep] + 1e-6);
            rx = std::max(0.0,rx);
            ry = std::max(0.0,ry);
            phix = (rx*rx + rx) / double(rx*rx + 1.0);
            phiy = (ry*ry + ry) / double(ry*ry + 1.0);
            
            dudx = (cells[i+1][j].u[RKstep] - cells[i][j].u[RKstep]) / double(cells[i][j].dx);
            dudy = (cells[i][j+1].u[RKstep] - cells[i][j].u[RKstep]) / double(cells[i][j].dy);
            cells[i][j].uEast = cells[i][j].u[RKstep] + phix*dudx*cells[i][j].dx /double(2);
            cells[i][j].uWest = cells[i][j].u[RKstep] - phix*dudx*cells[i][j].dx /double(2);
            cells[i][j].uNorth = cells[i][j].u[RKstep] + phiy*dudy*cells[i][j].dy /double(2);
            cells[i][j].uSouth = cells[i][j].u[RKstep] - phiy*dudy*cells[i][j].dy /double(2);
        }
    }
}

void ReconstructVariablesFromm(boost::multi_array<Cell,2> &cells, int &RKstep)
{
    double dudx = 0;
    double dudy = 0;
    for(int i = NUM_GHOST_CELLS - 1; i < NUM_X_CELLS + NUM_GHOST_CELLS + 1; i++)
    {
        for(int j = NUM_GHOST_CELLS - 1; j < NUM_Y_CELLS + NUM_GHOST_CELLS + 1; j++)
        {
            dudx = (cells[i+1][j].u[RKstep] - cells[i-1][j].u[RKstep]) / double(2*cells[i][j].dx);
            dudy = (cells[i][j+1].u[RKstep] - cells[i][j-1].u[RKstep]) / double(2*cells[i][j].dy);
            cells[i][j].uEast = cells[i][j].u[RKstep] + dudx*cells[i][j].dx /double(2);
            cells[i][j].uWest = cells[i][j].u[RKstep] - dudx*cells[i][j].dx /double(2);
            cells[i][j].uNorth = cells[i][j].u[RKstep] + dudy*cells[i][j].dy /double(2);
            cells[i][j].uSouth = cells[i][j].u[RKstep] - dudy*cells[i][j].dy /double(2);
        }
    }
}


//Errado preciso arrumar a parte inicial
void ReconstructVariablesWENO5thOrder(boost::multi_array<Cell,2> &cells, int &RKstep)
{
    double pointMinus2 = 0, pointMinus1 = 0, point = 0, pointPlus1 = 0, pointPlus2 = 0;
    double a = 0, u = 0, v = 0;
    double p0 = 0, p1 = 0, p2 = 0;
    double w0 = 0, w1 = 0, w2 = 0;
    double alpha0 = 0, alpha1 = 0, alpha2 = 0, alphaSum = 0;
    double IS0 = 0, IS1 = 0, IS2 = 0;
    double C0 = 0.1, C1 = 0.6, C2 = 0.3;
    double e = 1e-6;
        
    double dflux[5];
    
    for(int i = NUM_GHOST_CELLS-1; i < NUM_X_CELLS + NUM_GHOST_CELLS + 1; i++)
    {
        for(int j = NUM_GHOST_CELLS - 1; j < NUM_GHOST_CELLS + NUM_Y_CELLS + 1; j++)
        {
            dflux[0] = (cells[i-1][j].TotalFlux - cells[i-2][j].TotalFlux);
            dflux[1] = (cells[i][j].TotalFlux - cells[i-2][j].TotalFlux)/double(2*cells[i][j].dx);
            dflux[2] = (cells[i+1][j].TotalFlux - cells[i-1][j].TotalFlux)/double(2*cells[i][j].dx);
            dflux[3] = (cells[i+2][j].TotalFlux - cells[i][j].TotalFlux)/double(2*cells[i][j].dx);
            dflux[4] = (cells[i+2][j].TotalFlux - cells[i+1][j].TotalFlux);
            
//             dflux[0] = cells[i-2][j].TotalFlux;
//             dflux[1] = cells[i-1][j].TotalFlux;
//             dflux[2] = cells[i][j].TotalFlux;
//             dflux[3] = cells[i+1][j].TotalFlux;
//             dflux[4] = cells[i+2][j].TotalFlux;
            
            a = abs(dflux[0]);
            for(int k = 1; k < 5; k++){if(dflux[k] > abs(a)){a=abs(dflux[k]);}}
            if(a == 0){a=1;}        
        
            pointMinus2 = 0.5*(cells[i-2][j].TotalFlux + a*cells[i-2][j].u[RKstep]);
            pointMinus1 = 0.5*(cells[i-1][j].TotalFlux + a*cells[i-1][j].u[RKstep]);
            point= 0.5*(cells[i][j].TotalFlux + a*cells[i][j].u[RKstep]);
            pointPlus1 = 0.5*(cells[i+1][j].TotalFlux + a*cells[i+1][j].u[RKstep]);
            pointPlus2 = 0.5*(cells[i+2][j].TotalFlux + a*cells[i+2][j].u[RKstep]);
        
                   
//             pointMinus2 = cells[i-2][j].u[RKstep];
//             pointMinus1 = cells[i-1][j].u[RKstep];
//             point= cells[i][j].u[RKstep];
//             pointPlus1 = cells[i+1][j].u[RKstep];
//             pointPlus2 = cells[i+2][j].u[RKstep];
            
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
            
            cells[i][j].uEast = (w0*p0 + w1*p1 + w2*p2);
            
            //Calculando flux west
            pointMinus2 = 0.5*(cells[i-2][j].TotalFlux - a*cells[i-2][j].u[RKstep]);
            pointMinus1 = 0.5*(cells[i-1][j].TotalFlux - a*cells[i-1][j].u[RKstep]);
            point= 0.5*(cells[i][j].TotalFlux - a*cells[i][j].u[RKstep]);
            pointPlus1 = 0.5*(cells[i+1][j].TotalFlux - a*cells[i+1][j].u[RKstep]);
            pointPlus2 = 0.5*(cells[i+2][j].TotalFlux - a*cells[i+2][j].u[RKstep]);
            
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
            
            cells[i][j].uWest = (w0*p0 + w1*p1 + w2*p2);
            
            ////////////////////////////////////////////////////////////////////////////////////////////////
            
            dflux[0] = (cells[i][j-1].TotalFlux - cells[i][j-2].TotalFlux);
            dflux[1] = (cells[i][j].TotalFlux - cells[i][j-2].TotalFlux)/double(2*cells[i][j].dy);
            dflux[2] = (cells[i][j+1].TotalFlux - cells[i][j-1].TotalFlux)/double(2*cells[i][j].dy);
            dflux[3] = (cells[i][j+2].TotalFlux - cells[i][j].TotalFlux)/double(2*cells[i][j].dy);
            dflux[4] = (cells[i][j+2].TotalFlux - cells[i][j+1].TotalFlux);
            
//             dflux[0] = cells[i][j-2].TotalFlux;
//             dflux[1] = cells[i][j-1].TotalFlux;
//             dflux[2] = cells[i][j].TotalFlux;
//             dflux[3] = cells[i][j+1].TotalFlux;
//             dflux[4] = cells[i][j+2].TotalFlux;
            
            a = abs(dflux[0]);
            for(int k = 1; k < 5; k++){if(dflux[k] > abs(a)){a=abs(dflux[k]);}}
//             if(a == 0){a=1;}        
        
            pointMinus2 = 0.5*(cells[i][j-2].TotalFlux + a*cells[i][j-2].u[RKstep]);
            pointMinus1 = 0.5*(cells[i][j-1].TotalFlux + a*cells[i][j-1].u[RKstep]);
            point= 0.5*(cells[i][j].TotalFlux + a*cells[i][j].u[RKstep]);
            pointPlus1 = 0.5*(cells[i][j+1].TotalFlux + a*cells[i][j+1].u[RKstep]);
            pointPlus2 = 0.5*(cells[i][j+2].TotalFlux + a*cells[i][j+2].u[RKstep]);
            
//             pointMinus2 = cells[i][j-2].u[RKstep];
//             pointMinus1 = cells[i][j-1].u[RKstep];
//             point= cells[i][j].u[RKstep];
//             pointPlus1 = cells[i][j+1].u[RKstep];
//             pointPlus2 = cells[i][j+2].u[RKstep];
            
            //Polinômios para o lado uSouth
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
            
            cells[i][j].uSouth = (w0*p0 + w1*p1 + w2*p2);
            
            //Calculando flux North
            pointMinus2 = 0.5*(cells[i][j-2].TotalFlux - a*cells[i][j-2].u[RKstep]);
            pointMinus1 = 0.5*(cells[i][j-1].TotalFlux - a*cells[i][j-1].u[RKstep]);
            point= 0.5*(cells[i][j].TotalFlux - a*cells[i][j].u[RKstep]);
            pointPlus1 = 0.5*(cells[i][j+1].TotalFlux - a*cells[i][j+1].u[RKstep]);
            pointPlus2 = 0.5*(cells[i][j+2].TotalFlux - a*cells[i][j+2].u[RKstep]);
            
            //Polinômios para o lado uNorth
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
            
            cells[i][j].uNorth = (w0*p0 + w1*p1 + w2*p2);
            
        }
    }
}

void CalculateFlux(boost::multi_array<Cell,2> &cells, int &RKstep)
{
    for(int i = 0; i < cells.shape()[0]; i++)
    {
        for(int j = 0; j < cells.shape()[1]; j++)
        {
            cells[i][j].TotalFlux = 0.0;
        }
    }
    
    double LeftValue = 0;
    double RightValue = 0;
    double UpperValue = 0;
    double BottomValue = 0;
    double Flux = 0;
    //Calculation for Vertical faces
    for(int j = NUM_GHOST_CELLS; j < NUM_GHOST_CELLS + NUM_Y_CELLS; j++)
    {
        for(int interface = NUM_GHOST_CELLS; interface < NUM_GHOST_CELLS + NUM_X_CELLS + 1; interface++)
        {
            LeftValue = cells[interface-1][j].uEast;
            RightValue = cells[interface][j].uWest;
            
            //Lax-Friedrich
            Flux = 0.5*ADVECTION_VEL_X*(LeftValue + RightValue) - 0.5*abs(ADVECTION_VEL_X)*(RightValue - LeftValue);
            
            Flux *= cells[0][0].dy;
            cells[interface - 1][j].TotalFlux -= Flux;
            cells[interface][j].TotalFlux += Flux;
            
        }
    }
    
    //Calculation for Horizontal faces
    for(int i = NUM_GHOST_CELLS; i < NUM_GHOST_CELLS + NUM_X_CELLS; i++)
    {
        for(int interface = NUM_GHOST_CELLS; interface < NUM_GHOST_CELLS + NUM_Y_CELLS + 1; interface++)
        {
            BottomValue = cells[i][interface-1].uNorth;
            UpperValue = cells[i][interface].uSouth;
            
            //Lax-Friedrich
            Flux = 0.5*ADVECTION_VEL_Y*(BottomValue + UpperValue) - 0.5*abs(ADVECTION_VEL_Y)*(UpperValue - BottomValue);
            
            Flux *= cells[0][0].dx;
            cells[i][interface - 1].TotalFlux -= Flux;
            cells[i][interface].TotalFlux += Flux;
            
        }
    }
}

void UpdateCellAveragesEuler(boost::multi_array<Cell,2> &cells, int &RKstep, double &dt)
{
    double area = cells[0][0].dx * cells[0][0].dy;
    for(int i = NUM_GHOST_CELLS; i < NUM_GHOST_CELLS + NUM_X_CELLS; i++)
    {
        for(int j = NUM_GHOST_CELLS; j < NUM_GHOST_CELLS + NUM_Y_CELLS; j++)
        {
            cells[i][j].u[RKstep+1] = cells[i][j].u[RKstep] + dt/(area) * cells[i][j].TotalFlux; 
        }
    }
}

void UpdateCellAveragesRK2(boost::multi_array<Cell,2> &cells, int &RKstep, double &dt)
{
    double area = cells[0][0].dx * cells[0][0].dy;
    for(int i = NUM_GHOST_CELLS; i < NUM_GHOST_CELLS + NUM_X_CELLS; i++)
    {
        for(int j = NUM_GHOST_CELLS; j < NUM_GHOST_CELLS + NUM_Y_CELLS; j++)
        {
            if(RKstep == 0)
            {
                cells[i][j].u[RKstep+1] = cells[i][j].u[RKstep] + dt/area * cells[i][j].TotalFlux; 
            }
            if(RKstep == 1)
            {
                cells[i][j].u[RKstep+1] = 0.5*(cells[i][j].u[RKstep-1] + cells[i][j].u[RKstep] + dt/area * cells[i][j].TotalFlux);
            }
        }
    }
}

void UpdateCellAveragesRK3(boost::multi_array<Cell,2> &cells, int &RKstep, double &dt)
{
    double area = cells[0][0].dx * cells[0][0].dy;
    for(int i = NUM_GHOST_CELLS; i < NUM_GHOST_CELLS + NUM_X_CELLS; i++)
    {
        for(int j = NUM_GHOST_CELLS; j < NUM_GHOST_CELLS + NUM_Y_CELLS; j++)
        {
            if(RKstep == 0)
            {
                cells[i][j].u[RKstep+1] = cells[i][j].u[RKstep] + dt/area * cells[i][j].TotalFlux; 
            }
            if(RKstep == 1)
            {
                cells[i][j].u[RKstep+1] = 0.25*(3*cells[i][j].u[RKstep-1] + cells[i][j].u[RKstep] + dt/area*cells[i][j].TotalFlux);
            }
            if(RKstep == 2)
            {
                cells[i][j].u[RKstep+1]=(cells[i][j].u[RKstep-2] + 2*cells[i][j].u[RKstep] + 2*dt/area*cells[i][j].TotalFlux)/float(3);
            }
        }
    }
}

void CopyVariables(boost::multi_array<Cell,2> &cells)
{
    for(int i = NUM_GHOST_CELLS; i < NUM_X_CELLS + NUM_GHOST_CELLS; i++)
    {
        for(int j = NUM_GHOST_CELLS; j < NUM_Y_CELLS + NUM_GHOST_CELLS; j++)
        {
            cells[i][j].u[0] = cells[i][j].u[NUM_RK_STEPS];//apenas um teste com o +=
        }
    }
}

void ErrorCalculation(boost::multi_array<Cell,2> &cell)
{
    boost::multi_array<Cell,2> Initialcell;
    Initialcell.resize(extents[cell.shape()[0]][cell.shape()[1]]);
    
    for(int i = 0; i < Initialcell.shape()[0]; i++)
    {
        for(int j = 0; j < Initialcell.shape()[1]; j++)
        {
            Initialcell[i][j].dx = cell[i][j].dx;
            Initialcell[i][j].dy = cell[i][j].dy;
            Initialcell[i][j].cx = cell[i][j].cx;
            Initialcell[i][j].cy = cell[i][j].cy;
        }
    }
    
    SolutionInitializerSquareWave(Initialcell);
    
    double ErrorL1 = 0.0;
    double ErrorL2 = 0.0;
    double ErrorLinf = 0.0;
    double Error = 0.0;
    for(int i = NUM_GHOST_CELLS; i < NUM_GHOST_CELLS + NUM_X_CELLS; i++)
    {
        for(int j = NUM_GHOST_CELLS; j < NUM_GHOST_CELLS + NUM_Y_CELLS; j++)
        {
            Error = abs(Initialcell[i][j].u[0] - cell[i][j].u[0]);
            ErrorL1 += Error;
            ErrorL2 += Error*Error;
            ErrorLinf = std::max(Error, ErrorLinf);
        }
    }
    
    ErrorL1 /= (NUM_X_CELLS * NUM_Y_CELLS);
    ErrorL2 = sqrt(ErrorL2/double(NUM_X_CELLS * NUM_Y_CELLS));
    
    std::ofstream out( "Erro.txt" , std::ofstream::out);
    out << "ErrorL1 = " << ErrorL1 << "\nErrorL2 = " << ErrorL2 << "\nErrorLinf = " << ErrorLinf;
    out.close();

}

void GradLenght(boost::multi_array<Cell,2> &cell, boost::multi_array<double,2> &grad, int &RKstep)
{
    double dx = 0, dy = 0;
    for(int i = NUM_GHOST_CELLS; i < NUM_GHOST_CELLS + NUM_X_CELLS; i++)
    {
        for(int j = NUM_GHOST_CELLS; j < NUM_GHOST_CELLS + NUM_Y_CELLS; j++)
        {
            dx = (cell[i+1][j].u[RKstep] - cell[i-1][j].u[RKstep]) / double(2);
            dy = (cell[i][j+1].u[RKstep] - cell[i][j-1].u[RKstep]) / double(2);
    
            grad[i][j] = sqrt(dx*dx + dy*dy);
        }
    }
}



