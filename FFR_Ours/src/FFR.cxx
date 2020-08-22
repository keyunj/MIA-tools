#include "FFR.h"

#include "vtkMath.h"
#include "vtkImageResample.h"
#include "vtkImageFlip.h"
#include "vtkLookupTable.h"
#include "vtkPointData.h"
#include "vtkXMLPolyDataWriter.h"

#include "common.h"
#include "ExtendTubeFilter.h"
#include "ExtendIdFilter.h"
#include "Learning.h"
#include "FFRCore.h"

using namespace std;

FFR::FFR()
{
}

FFR::~FFR()
{
}

void FFR::initialize(vtkSmartPointer<vtkImageData> imageImgData, vtkSmartPointer<vtkImageData> labelImgData, string& workPath, char* imageBuffer, char* labelBuffer)
{
	clock_t timeBegin, timeEnd;
	cout << endl << "Initializing ... ";
	timeBegin = clock();

	this->workPath = workPath;
	this->imageImgData = imageImgData;
	this->labelImgData = labelImgData;
	imagePath.assign(imageBuffer);
	labelPath.assign(labelBuffer);

	reductFactor = 0.5;

	interpolator = vtkSmartPointer<vtkImageInterpolator>::New();
	interpolator->SetInterpolationModeToLinear();
	interpolator->SetOutValue(-3024.0);
	interpolator->Initialize(this->imageImgData);

	timeEnd = clock();
	cout << "Done Time: " << (timeEnd - timeBegin) / 1000 << "." << (timeEnd - timeBegin) % 1000 << "s" << endl;
}

bool FFR::detectLandmarks()
{
	// ¼ì²âlandmarkµã
	bool flag = DetectLandmarks(imageImgData, learn, ostiumLandmark, interpolator, workPath);

	if (!flag) return false;

	return true;
}

bool FFR::detectCenterline()
{
	centerlineModel = vtkSmartPointer<vtkPolyData>::New();

	// ½µ²ÉÑù£¬0.5±¶
	vtkSmartPointer<vtkImageResample> resample = vtkSmartPointer<vtkImageResample>::New();
	resample->SetInputData(labelImgData);
	for (int i = 0; i < 3; i++) resample->SetAxisMagnificationFactor(i, reductFactor);
	try
	{
		resample->Update();
	}
	catch (...)
	{
		cerr << "Error occurs when resampling " << endl;
	}
	vtkSmartPointer<vtkImageData> resampleImg = resample->GetOutput();

	thinImgData = vtkSmartPointer<vtkImageData>::New();
	bool flag = DetectCenterline(resampleImg, thinImgData, centerlineModel, ostiumLandmark[SmartCoronary::LEFT_CORONARY_OSTIUM], ostiumLandmark[SmartCoronary::RIGHT_CORONARY_OSTIUM]);

 	return flag;
}

bool FFR::detectTube()
{
	cout << "Detecting coronary tube ..." << endl;

	bool flag = DetectCenterlineLumenWall(labelImgData, centerlineModel);

	if (!flag) {
		cerr << "Error occurs when detecting lumen wall " << endl;
		return false;
	}

	centerlineId = vtkSmartPointer<vtkIdFilter>::New();
	centerlineId->SetInputData(centerlineModel);
	centerlineId->PointIdsOff();
	centerlineId->CellIdsOn();
	centerlineId->FieldDataOn();
	centerlineId->SetIdsArrayName("SegmentId");

	centerlineTube = vtkSmartPointer<ExtendTubeFilter>::New();
	centerlineTube->SetInputConnection(centerlineId->GetOutputPort());
	try
	{
		centerlineTube->Update();
	}
	catch (...)
	{
		cerr << "Error occurs when detecting centerline tube structure " << endl;
	}

	vtkPolyData* centerlinePoly = centerlineTube->GetOutput(0);
	vtkSmartPointer<vtkDoubleArray> ffrArray = vtkSmartPointer<vtkDoubleArray>::New();
	ComputeFFRCore(centerlineModel, centerlinePoly, ffrArray, ostiumLandmark[SmartCoronary::LEFT_CORONARY_OSTIUM], ostiumLandmark[SmartCoronary::RIGHT_CORONARY_OSTIUM]);
	centerlinePoly->GetPointData()->SetScalars(ffrArray);
	centerlineTube->SetUpdateSegment(-1);
	centerlineTube->UpdateOutput0Off();
	try
	{
		centerlineTube->Update();
	}
	catch (...)
	{
		cerr << "Error occurs when computing FFR " << endl;
	}

	return true;
}

void FFR::writePolyData()
{
	if (centerlineTube) {
		vtkSmartPointer<vtkXMLPolyDataWriter> write = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		write->SetInputConnection(centerlineTube->GetOutputPort(2));

		std::string tubePath = imagePath.substr(0, imagePath.length() - 4);
		tubePath.append("_ffrtube.vtp");
		write->SetFileName(tubePath.c_str());

		try
		{
			write->Write();
		}
		catch (...)
		{
			cerr << "Error occurs when writing polydata " << endl;
		}
	}
}

void FFR::rendering(vtkImageData* renderImage)
{
	// Actor
	vtkSmartPointer<vtkImageActor> modelImgActor1 = vtkSmartPointer<vtkImageActor>::New();
	vtkSmartPointer<vtkImageActor> modelImgActor2 = vtkSmartPointer<vtkImageActor>::New();
	vtkSmartPointer<vtkImageActor> modelImgActor3 = vtkSmartPointer<vtkImageActor>::New();

	// Render
	vtkSmartPointer<vtkRenderer> win0Render = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderer> win1Render = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderer> win2Render = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderer> win3Render = vtkSmartPointer<vtkRenderer>::New();

	// Render Window
	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(win0Render);

	// Render Window Interaction
	vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renWinInteractor->SetRenderWindow(renWin);

	// Render Window Interactor Style
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> renWinInterSty = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	renWinInteractor->SetInteractorStyle(renWinInterSty);

	// äÖÈ¾XYZÇÐÃæ
	///. Cell Picker
	vtkSmartPointer<vtkCellPicker> modelCellPicker = vtkSmartPointer<vtkCellPicker>::New();
	modelCellPicker->SetTolerance(0.005);

	///. Plane widgets
	vtkSmartPointer<vtkImagePlaneWidget> modelImgPlaneWidX = vtkSmartPointer<vtkImagePlaneWidget>::New();
	modelImgPlaneWidX->SetInteractor(renWinInteractor);
	modelImgPlaneWidX->SetKeyPressActivationValue('x');
	modelImgPlaneWidX->SetPicker(modelCellPicker);
	modelImgPlaneWidX->RestrictPlaneToVolumeOn();
	modelImgPlaneWidX->GetPlaneProperty()->SetColor(1.0, 0.0, 0.0);
	modelImgPlaneWidX->DisplayTextOn();
	modelImgPlaneWidX->TextureInterpolateOff();
	modelImgPlaneWidX->SetResliceInterpolateToLinear();
	modelImgPlaneWidX->SetInputData(renderImage);
	modelImgPlaneWidX->SetPlaneOrientationToXAxes();
	modelImgPlaneWidX->SetSliceIndex(255);
	modelImgPlaneWidX->GetTexturePlaneProperty()->SetOpacity(1);
	modelImgPlaneWidX->On();

	vtkSmartPointer<vtkImagePlaneWidget> modelImgPlaneWidY = vtkSmartPointer<vtkImagePlaneWidget>::New();
	modelImgPlaneWidY->SetInteractor(renWinInteractor);
	modelImgPlaneWidY->SetKeyPressActivationValue('y');
	modelImgPlaneWidY->SetPicker(modelCellPicker);
	modelImgPlaneWidY->RestrictPlaneToVolumeOn();
	modelImgPlaneWidY->GetPlaneProperty()->SetColor(0.0, 1.0, 0.0);
	modelImgPlaneWidY->DisplayTextOn();
	modelImgPlaneWidY->TextureInterpolateOff();
	modelImgPlaneWidY->SetResliceInterpolateToLinear();
	modelImgPlaneWidY->SetInputData(renderImage);
	modelImgPlaneWidY->SetPlaneOrientationToYAxes();
	modelImgPlaneWidY->SetSliceIndex(255);
	modelImgPlaneWidY->On();

	vtkSmartPointer<vtkImagePlaneWidget> modelImgPlaneWidZ = vtkSmartPointer<vtkImagePlaneWidget>::New();
	modelImgPlaneWidZ->SetInteractor(renWinInteractor);
	modelImgPlaneWidZ->SetKeyPressActivationValue('z');
	modelImgPlaneWidZ->SetPicker(modelCellPicker);
	modelImgPlaneWidZ->RestrictPlaneToVolumeOn();
	modelImgPlaneWidZ->GetPlaneProperty()->SetColor(0.0, 0.0, 1.0);
	modelImgPlaneWidZ->DisplayTextOn();
	modelImgPlaneWidZ->TextureInterpolateOff();
	modelImgPlaneWidZ->SetResliceInterpolateToLinear();
	modelImgPlaneWidZ->SetInputData(renderImage);
	modelImgPlaneWidZ->SetPlaneOrientationToZAxes();
	modelImgPlaneWidZ->SetSliceIndex(150);
	modelImgPlaneWidZ->On();

	renWin->Render();

	///. Sclies
	vtkSmartPointer<vtkImageMapToColors> modelMapToColor1 = vtkSmartPointer<vtkImageMapToColors>::New();
	modelMapToColor1->PassAlphaToOutputOff();
	modelMapToColor1->SetActiveComponent(0);
	modelMapToColor1->SetOutputFormatToLuminance();
	modelMapToColor1->SetInputData(modelImgPlaneWidX->GetResliceOutput());
	modelMapToColor1->SetLookupTable((vtkScalarsToColors*)modelImgPlaneWidX->GetLookupTable());

	vtkSmartPointer<vtkImageMapToColors> modelMapToColor2 = vtkSmartPointer<vtkImageMapToColors>::New();
	modelMapToColor2->PassAlphaToOutputOff();
	modelMapToColor2->SetActiveComponent(0);
	modelMapToColor2->SetOutputFormatToLuminance();
	modelMapToColor2->SetInputData(modelImgPlaneWidY->GetResliceOutput());
	modelMapToColor2->SetLookupTable((vtkScalarsToColors*)modelImgPlaneWidX->GetLookupTable());

	vtkSmartPointer<vtkImageMapToColors> modelMapToColor3 = vtkSmartPointer<vtkImageMapToColors>::New();
	modelMapToColor3->PassAlphaToOutputOff();
	modelMapToColor3->SetActiveComponent(0);
	modelMapToColor3->SetOutputFormatToLuminance();
	modelMapToColor3->SetInputData(modelImgPlaneWidZ->GetResliceOutput());
	modelMapToColor3->SetLookupTable((vtkScalarsToColors*)modelImgPlaneWidX->GetLookupTable());

	// other polydata
	// landmark
	vtkSmartPointer<vtkPolyData> landmarkPoly = vtkSmartPointer<vtkPolyData>::New();
	landmarkPoly->Allocate();
	vtkSmartPointer<vtkIdList> lmIdList = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkPoints> lmPoints = vtkSmartPointer<vtkPoints>::New();
	vtkIdType id = lmPoints->InsertNextPoint(ostiumLandmark[SmartCoronary::LEFT_CORONARY_OSTIUM]);
	lmIdList->InsertNextId(id);
	id = lmPoints->InsertNextPoint(ostiumLandmark[SmartCoronary::RIGHT_CORONARY_OSTIUM]);
	lmIdList->InsertNextId(id);
	//for (vtkIdType ii = 0; ii < SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS; ii++) {
	//	//for (int jj = 0; jj < 3; jj++) ostiumLandmark[ii][jj] /= 3;
	//	vtkIdType id = lmPoints->InsertNextPoint(ostiumLandmark[ii]);
	//	lmIdList->InsertNextId(id);
	//}
	landmarkPoly->InsertNextCell(VTK_POLY_VERTEX, lmIdList);
	landmarkPoly->SetPoints(lmPoints);

	vtkSmartPointer<vtkPolyDataMapper> landmarkMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	landmarkMapper->SetInputData(landmarkPoly);
	vtkSmartPointer<vtkActor> landmarkActor = vtkSmartPointer<vtkActor>::New();
	landmarkActor->SetMapper(landmarkMapper);
	landmarkActor->GetProperty()->SetPointSize(3);
	landmarkActor->GetProperty()->SetColor(0.0, 1.0, 0.0);

	// centerline
	if (true) {
		vtkSmartPointer<vtkPolyDataMapper> centerlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		centerlineMapper->SetInputData(centerlineModel);
		vtkSmartPointer<vtkActor> centerlineActor = vtkSmartPointer<vtkActor>::New();
		centerlineActor->SetMapper(centerlineMapper);
		centerlineActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
		centerlineActor->GetProperty()->SetPointSize(5.0);

		win0Render->AddActor(centerlineActor);
	}

	// centerlinetube
	if (centerlineTube) {
		vtkSmartPointer<vtkLookupTable> clLookupTable = vtkSmartPointer<vtkLookupTable>::New();
		clLookupTable->SetTableRange(0.7, 1.0);
		clLookupTable->SetValueRange(1.0, 1.0);
		clLookupTable->SetHueRange(0.0, 0.6667);
		clLookupTable->SetSaturationRange(1.0, 1.0);
		clLookupTable->SetRampToLinear();
		clLookupTable->Build();

		//vtkSmartPointer<vtkDecimatePro> deci = vtkSmartPointer<vtkDecimatePro>::New();
		//deci->SetInputConnection(centerlineTube->GetOutputPort(2));
		//deci->SetTargetReduction(0.9);
		//deci->PreserveTopologyOn();
		//vtkSmartPointer<vtkSmoothPolyDataFilter> smoother = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
		//smoother->SetInputConnection(deci->GetOutputPort());
		//smoother->SetNumberOfIterations(20);

		vtkSmartPointer<vtkPolyDataMapper> centerlineTubeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		centerlineTubeMapper->SetInputConnection(centerlineTube->GetOutputPort(2));
		centerlineTubeMapper->SetColorModeToMapScalars();
		centerlineTubeMapper->SetScalarModeToUsePointData();
		centerlineTubeMapper->UseLookupTableScalarRangeOn();
		centerlineTubeMapper->SetLookupTable(clLookupTable);
		vtkSmartPointer<vtkActor> centerlineTubeActor = vtkSmartPointer<vtkActor>::New();
		centerlineTubeActor->SetMapper(centerlineTubeMapper);

		win0Render->AddActor(centerlineTubeActor);
	}

	// Set Properties
	modelImgActor1->GetMapper()->SetInputConnection(modelMapToColor1->GetOutputPort());
	modelImgActor2->GetMapper()->SetInputConnection(modelMapToColor2->GetOutputPort());
	modelImgActor3->GetMapper()->SetInputConnection(modelMapToColor3->GetOutputPort());

	win0Render->AddActor(landmarkActor);
	win0Render->SetBackground(0.1, 0.1, 0.1);
	win0Render->SetViewport(0.00, 0.00, 0.75, 1.00);

	win1Render->AddActor(modelImgActor1);
	win1Render->SetBackground(0.1, 0.0, 0.0);
	win1Render->SetViewport(0.75, 0.67, 1.00, 1.00);

	win2Render->AddActor(modelImgActor2);
	win2Render->SetBackground(0.0, 0.1, 0.0);
	win2Render->SetViewport(0.75, 0.33, 1.00, 0.67);

	win3Render->AddActor(modelImgActor3);
	win3Render->SetBackground(0.0, 0.0, 0.1);
	win3Render->SetViewport(0.75, 0.00, 1.00, 0.33);

	renWin->SetSize(1200, 900);
	renWin->AddRenderer(win1Render);
	renWin->AddRenderer(win2Render);
	renWin->AddRenderer(win3Render);

	// Rendering
	renWinInteractor->Initialize();
	renWin->Render();
	renWinInteractor->Start();
}
