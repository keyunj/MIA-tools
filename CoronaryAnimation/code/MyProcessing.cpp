#include "MyProcessing.h"

namespace file {
	void GetFiles(std::string path, std::vector<std::string>& files, bool label) {
		struct _finddata_t fileInfo;
		std::string p;
		intptr_t hFile = _findfirst(p.assign(path).append("*").c_str(), &fileInfo);
		if (hFile != -1L) {
			do {
				if (fileInfo.attrib & _A_SUBDIR) {
					if (strcmp(fileInfo.name, ".") != 0 && strcmp(fileInfo.name, "..") != 0) {
						if (label) files.push_back(fileInfo.name);
					}
				}
				else {
					if (fileInfo.name != "." && fileInfo.name != "..") {
						if (!label) files.push_back(fileInfo.name);
					}
				}
			} while (_findnext(hFile, &fileInfo) == 0);
			_findclose(hFile);
		}
	}
}

void GetMeanPoint(vtkSmartPointer<vtkPoints> points, double mPoint[3]) {
	mPoint[0] = mPoint[1] = mPoint[2] = 0;
	for (vtkIdType ii = 0; ii < points->GetNumberOfPoints(); ii++) {
		mPoint[0] += points->GetPoint(ii)[0];
		mPoint[1] += points->GetPoint(ii)[1];
		mPoint[2] += points->GetPoint(ii)[2];
	}
	mPoint[0] /= points->GetNumberOfPoints();
	mPoint[1] /= points->GetNumberOfPoints();
	mPoint[2] /= points->GetNumberOfPoints();
}