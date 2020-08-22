// System
#include <sstream>
#include <fstream>
#include <vector>
#include <iostream>
#include <Windows.h>
#include <commdlg.h>
#include <direct.h>

// VTK
#include "vtkSmartPointer.h"
#include "vtkSmartVolumeMapper.h"
#include "vtkColorTransferFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkVolume.h"
#include "vtkPiecewiseFunction.h"
#include "vtkXMLPolyDataReader.h"

#include "vtkImageData.h"
#include "vtkMetaImageReader.h"
#include "vtkMetaImageWriter.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkImageProperty.h"
#include "vtkProperty.h"
#include "vtkCommand.h"

#include "vtkCellPicker.h"
#include "vtkImagePlaneWidget.h"
#include "vtkPolyDataMapper.h"
#include "vtkImageSliceMapper.h"
#include "vtkImageSlice.h"
#include "vtkLookupTable.h"
#include "vtkImageMapToColors.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"

#include "vtkRuledSurfaceFilter.h"

// Custom
#include "MyProcessing.h"

class keyCallback : public vtkCommand
{
public:
	static keyCallback *New() {
		return new keyCallback;
	}
	void Delete() {
		delete this;
	}

	virtual void Execute(vtkObject *caller, unsigned long, void*) {
		std::string key = iren->GetKeySym();
		if (key == "w"/* && ii < files.size()*/) {
			for (int ii = 0; ii < files.size(); ii++) {
				std::string img_path = "../data/" + files[ii];
				vtkSmartPointer<vtkMetaImageReader> reader = vtkSmartPointer<vtkMetaImageReader>::New();
				reader->SetFileName(img_path.c_str());
				reader->Update();
				mapper->SetInputData(reader->GetOutput());

				if (iren) iren->GetRenderWindow()->Render();
				//ii++;
				Sleep(200);
			}
		}
	}

public:
	vtkRenderWindowInteractor *iren;
	vtkSmartVolumeMapper *mapper;
	std::vector<std::string> files;
	int ii = 0;
};

int main() {
	std::string initImage = "../1000_animation.mha";

	// row image
	vtkSmartPointer<vtkMetaImageReader> modelImgReader = vtkSmartPointer<vtkMetaImageReader>::New();
	modelImgReader->SetFileName(initImage.c_str());
	modelImgReader->Update();

	vtkSmartPointer<vtkSmartVolumeMapper> volMapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();
	volMapper->SetInputData(modelImgReader->GetOutput());
	volMapper->SetBlendModeToComposite();
	volMapper->SetRequestedRenderModeToDefault();

	vtkSmartPointer<vtkPiecewiseFunction> opac = vtkSmartPointer<vtkPiecewiseFunction>::New();
	opac->AddPoint(50.0, 0.0);
	opac->AddPoint(50.1, 1.0);
	//opac->AddPoint(255.0, 1.0);

	vtkSmartPointer<vtkColorTransferFunction> comp = vtkSmartPointer<vtkColorTransferFunction>::New();
	comp->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
	comp->AddRGBPoint(100.0, 0.5, 0.5, 0.5);
	comp->AddRGBPoint(254.0, 1.0, 1.0, 1.0);
	comp->AddRGBPoint(255.0, 1.0, 0.0, 0.0);

	vtkSmartPointer<vtkVolumeProperty> volProp = vtkSmartPointer<vtkVolumeProperty>::New();
	volProp->ShadeOff();
	//volProp->SetInterpolationType(VTK_LINEAR_INTERPOLATION);
	volProp->SetColor(comp);
	volProp->SetScalarOpacity(opac);

	vtkSmartPointer<vtkVolume> vol = vtkSmartPointer<vtkVolume>::New();
	vol->SetMapper(volMapper);
	vol->SetProperty(volProp);

	//std::string coro_path = "../mesh.vtp";
	//vtkSmartPointer<vtkXMLPolyDataReader> coroReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	//coroReader->SetFileName(coro_path.c_str());
	//coroReader->Update();

	//vtkSmartPointer<vtkPolyDataMapper> coroMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	//coroMapper->SetInputData(coroReader->GetOutput());
	//
	//vtkSmartPointer<vtkActor> coroActor = vtkSmartPointer<vtkActor>::New();
	//coroActor->SetMapper(coroMapper);
	//coroActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
	
	// Render
	vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
	//ren->AddActor(coroActor);
	ren->AddVolume(vol);
	ren->SetBackground(0.1, 0.1, 0.1);
	//ren->SetViewport(0.00, 0.00, 0.50, 1.00);

	// Render Window
	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(ren);
	
	// Render Window Interaction
	vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(renWin);
	
	// Render Window Interactor Style
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> renWinInterSty = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	iren->SetInteractorStyle(renWinInterSty);
	
	// property
	vtkSmartPointer<keyCallback> keycb = vtkSmartPointer<keyCallback>::New();
	keycb->iren = iren;
	keycb->mapper = volMapper;
	file::GetFiles("../data/", keycb->files, false);
	
	iren->AddObserver(vtkCommand::CharEvent, keycb);

	iren->Initialize();
	renWin->Render();
	iren->Start();

	return 0;
}

