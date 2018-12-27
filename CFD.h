#ifndef CFD_H
#define CFD_H

#include "Class_Cell.h"
#include "config.h"

void GhostCellsUpdater(Cell *, int &);

void ReconstructVariables(Cell *, int &);

void CalculateFlux(Cell *, int &);

void UpdateCellAverages(Cell *, int &, double &);

void CopyVariables(Cell *);

#endif
