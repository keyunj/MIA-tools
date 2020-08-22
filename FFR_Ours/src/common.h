#pragma once

#include "Learning.h"

#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkImageData.h"
#include "vtkImageInterpolator.h"

#define MEDIAL_SURFACE_SUBDIVISIONS 2

typedef unsigned short LabelType;

namespace SmartCoronary
{
	enum LVCorLandMarkList {
		LEFT_VENTRICLE_APEX = 0,
		BASE_LEFT_END,
		BASE_RIGHT_END,
		BASE_ANTERIOR_END,
		BASE_POSTERIOR_END,
		LEFT_CORONARY_OSTIUM,
		RIGHT_CORONARY_OSTIUM,
		//////////////////////////////////////////////////////////////////////////
		NUMBER_OF_LVCOR_LANDMARKS
	};

	struct LVCorLandMark
	{
		LVCorLandMarkList mark;
		int vertexid;
		char *name;
	};

	extern const LVCorLandMark LVCorLandmarkTable[];
}

bool DetectLandmarks(vtkImageData *imageData, Learning& learn, double landmarks[][3], vtkImageInterpolator *interpolator, std::string& workPath);
bool DetectCenterline(vtkImageData *labelImage, vtkImageData* thinImage, vtkPolyData *centerlineModel, double leftOstium[3], double rightOstium[3]);
void SimplifyCenterline(vtkPolyData* clModel);
void AxisCenterline(vtkPolyData* clModel, double planenormal[3] = NULL);
void LumenWallCenterline(vtkPolyData* clModel, int centerline_components = 10);
bool DetectCenterlineLumenWall(vtkImageData* imageImage, vtkPolyData* clModel);