#ifndef __FFRCore_h
#define __FFRCore_h

#include "vtkPolyData.h"
#include "vtkDoubleArray.h"

void ComputeFFRCore(vtkPolyData *clModel, vtkPolyData *clPoly, vtkDoubleArray *ffrArray, double leftOstium[3], double rightOstium[3]);

#endif