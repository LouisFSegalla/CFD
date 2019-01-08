#ifndef CFD_H
#define CFD_H

#include "Class_Cell.h"
#include "config.h"

void GhostCellsUpdater(Cell *, int &);

void ReconstructVariablesFirstOrder(Cell *, int &);

void ReconstructVariablesBeamWarming(Cell *, int &);

void ReconstructVariablesLaxWendroff(Cell *, int &);

void ReconstructVariablesLimitedLW(Cell *, int &);

void ReconstructVariablesFromm(Cell *, int &);

void ReconstructVariablesWENO5thOrder(Cell *, int &, double &);

void CalculateFlux(Cell *, int &);

void UpdateCellAveragesEuler(Cell *, int &, double &);

void UpdateCellAveragesRK2(Cell *, int &, double &);

void UpdateCellAveragesRK3(Cell *, int &, double &);

void CopyVariables(Cell *);

#endif
