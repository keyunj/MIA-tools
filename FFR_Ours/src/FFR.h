#pragma once

#include <string>

#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkPolyData.h"
#include "vtkImageInterpolator.h"
#include "vtkIdFilter.h"
#include "vtkIdList.h"
#include "vtkPoints.h"

#include "vtkCellPicker.h"
#include "vtkPolyDataMapper.h"
#include "vtkImageActor.h"
#include "vtkImagePlaneWidget.h"
#include "vtkProperty.h"
#include "vtkImageMapToColors.h"
#include "vtkMapper.h"
#include "vtkImageMapper3D.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"

#include "vtkDecimatePro.h"
#include "vtkSmoothPolyDataFilter.h"

#include "Learning.h"
#include "ExtendTubeFilter.h"

#include "common.h"
#include "ExtendTubeFilter.h"
#include "ExtendIdFilter.h"
#include "Learning.h"
#include "FFRCore.h"


class FFR {
public:
	FFR();
	~FFR();

public:
	void initialize(vtkSmartPointer<vtkImageData> imageImgData, vtkSmartPointer<vtkImageData> labelImgData, std::string& workPath, char* imageBuffer, char* labelBuffer);

	bool detectLandmarks();
	bool detectCenterline();
	bool detectTube();

	void writePolyData();

	void rendering(vtkImageData* renderImage);

public:
	std::string workPath;
	double reductFactor;

	std::string imagePath;
	std::string labelPath;

	vtkSmartPointer<vtkImageData> imageImgData;
	vtkSmartPointer<vtkImageData> labelImgData;
	vtkSmartPointer<vtkImageData> thinImgData;

	vtkSmartPointer<vtkPolyData> centerlineModel;
	vtkSmartPointer<ExtendTubeFilter> centerlineTube;

	vtkSmartPointer<vtkImageInterpolator> interpolator;
	vtkSmartPointer<vtkIdFilter> centerlineId;
	
	Learning learn;

	double ostiumLandmark[SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS][3];
	//std::vector< std::pair< vtkSmartPointer<vtkSphereWidget>, vtkSmartPointer<vtkCaptionActor2D> > > landmarkWidget;
};