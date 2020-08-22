#pragma once

#include <vector>
#include <io.h>

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkMath.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkMatrix4x4.h"
#include "vtkImageReslice.h"
#include "vtkImageData.h"
#include "vtkPointLocator.h"
#include "vtkIdList.h"
#include "vtkImageInterpolator.h"

#include "vtkImageActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

namespace file {
	void GetFiles(std::string path, std::vector<std::string>& files, bool label = true);

}
void GetMeanPoint(vtkSmartPointer<vtkPoints> points, double mPoint[3]);