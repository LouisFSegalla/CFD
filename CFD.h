#ifndef CFD_H
#define CFD_H

#include "boost/multi_array.hpp"


#include "Class_Cell.h"
#include "config.h"

void SolutionInitializerSquareWave(boost::multi_array<Cell,2> &);

void SolutionInitializerSineWave(boost::multi_array<Cell,2> &);

void SolutionInitializerCircle(boost::multi_array<Cell,2> & , const int &);

void GhostCellsUpdater(boost::multi_array<Cell,2> &, int &);

void ReconstructVariablesFirstOrder(boost::multi_array<Cell,2> &, int &);

void ReconstructVariablesBeamWarming(boost::multi_array<Cell,2> &, int &);

void ReconstructVariablesLaxWendroff(boost::multi_array<Cell,2> &, int &);

void ReconstructVariablesLimitedLW(boost::multi_array<Cell,2> &, int &);

void ReconstructVariablesFromm(boost::multi_array<Cell,2> &, int &);

void ReconstructVariablesWENO5thOrder(boost::multi_array<Cell,2> &, int &);

void CalculateFlux(boost::multi_array<Cell,2> &, int &);

void UpdateCellAveragesEuler(boost::multi_array<Cell,2> &, int &, double &);

void UpdateCellAveragesRK2(boost::multi_array<Cell,2> &, int &, double &);

void UpdateCellAveragesRK3(boost::multi_array<Cell,2> &, int &, double &);

void CopyVariables(boost::multi_array<Cell,2> &);

void ErrorCalculation(boost::multi_array<Cell,2> &);

void UpdateCellAveragesRungeKutta2LevelSet(boost::multi_array<Cell,2> &, int &, double &);

#endif
