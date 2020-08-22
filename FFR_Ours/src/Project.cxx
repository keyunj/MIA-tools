#include <conio.h>
#include <string>
#include <iostream>
#include <Windows.h>
#include <commdlg.h>
#include <time.h>

#include "FFR.h"

#include "vtkSmartPointer.h"
#include "vtkMetaImageReader.h"
#include "vtkNIFTIImageReader.h"
#include "vtkImageData.h"
#include "vtkPolyData.h"
#include "vtkImageReslice.h"

using namespace std;

int main() {
	char orderChar;
	clock_t timeBegin, timeEnd;

	char currentPath[1000];
	GetCurrentDirectory(1000, currentPath);
	string workPath(currentPath);

	while (true) {
		cout << "Press e(or E) to exit..." << endl;
		orderChar = _getch();
		system("cls");
		fflush(stdin);

		if (orderChar == 'e' || orderChar == 'E')
			break;

		OPENFILENAME openFileName;
		ZeroMemory(&openFileName, sizeof(OPENFILENAME));

		cout << "Select Image File ... ";
		char imageBuffer[MAX_PATH] = { 0 };
		openFileName.lpstrTitle = "Select Image File";
		openFileName.lStructSize = sizeof(OPENFILENAME);
		openFileName.lpstrFilter = "MHA File(*.mha)\0*.mha\0MHD File(*.mhd)\0*.mhd\0\0";
		openFileName.nFilterIndex = 1;
		openFileName.lpstrFile = imageBuffer;
		//openFileName.lpstrFile[0] = '\0';
		openFileName.nMaxFile = sizeof(imageBuffer);
		openFileName.Flags = OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST;
		if (!GetOpenFileName(&openFileName)) {
			cout << "Can not open the file! Please retry..." << endl;
			continue;
		}
		cout << "Done" << endl;

		cout << "Select Label File ... ";
		char labelBuffer[MAX_PATH] = { 0 };
		openFileName.lpstrTitle = "Select Label File";
		openFileName.lpstrFilter = "MHA File(*.mha)\0*.mha\0MHD File(*.mhd)\0*.mhd\0\0";
		openFileName.lpstrFile = labelBuffer;
		if (!GetOpenFileName(&openFileName)) {
			cout << "Can not open the file! Please retry..." << endl;
			continue;
		}
		cout << "Done" << endl;

		vtkSmartPointer<vtkMetaImageReader> imageReader = vtkSmartPointer<vtkMetaImageReader>::New();
		vtkSmartPointer<vtkMetaImageReader> labelReader = vtkSmartPointer<vtkMetaImageReader>::New();
		vtkSmartPointer<vtkImageData> imageImgData = vtkSmartPointer<vtkImageData>::New();
		vtkSmartPointer<vtkImageData> labelImgData = vtkSmartPointer<vtkImageData>::New();

		// load data
		imageReader->SetFileName(imageBuffer);
		imageReader->Update();
		imageImgData = imageReader->GetOutput();
		labelReader->SetFileName(labelBuffer);
		labelReader->Update();
		labelImgData = labelReader->GetOutput();

		FFR FFRProc;
		FFRProc.initialize(imageImgData, labelImgData, workPath, imageBuffer, labelBuffer);

		bool flag = FFRProc.detectLandmarks();
		if (!flag) continue;

		flag = FFRProc.detectCenterline();
		if (!flag) continue;

		flag = FFRProc.detectTube();
		if (!flag) continue;
		FFRProc.writePolyData();

		FFRProc.rendering(FFRProc.imageImgData);
		//TODO: 绑定原图像与建模之后的FFR显示之间的camera参数，使得两者的视角能够保持同步
	}

	system("exit");
	return 0;
}