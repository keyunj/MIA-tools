#include <omp.h>

#include <set>
#include <list>
#include <time.h>

#include "psimpl.h"

#include "vtkSmartPointer.h"
#include "vtkMath.h"
#include "vtkIdList.h"
#include "vtkIntArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkCellArray.h"
#include "vtkVariantArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkTransform.h"
#include "vtkWindowedSincPolyDataFilter.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkCleanPolyData.h"
#include "vtkAppendPolyData.h"
#include "vtkLandmarkTransform.h"
#include "vtkThinPlateSplineTransform.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTable.h"
#include "vtkPCAStatistics.h"
#include "vtkImageResample.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkMetaImageWriter.h"
#include "vtkTriangle.h"
#include "vtkPointLocator.h"
#include "vtkCellLocator.h"

#include "itkVTKImageToImageFilter.h"
#include "itkImageToVTKImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryFillholeImageFilter.h"
#include "itkOrImageFilter.h"
#include "itkBinaryThinningImageFilter3D.h"
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"
#include "itkResampleImageFilter.h"

#include "common.h"
#include "Learning.h"
#include "LearningImpl.h"
#include "ExtendSplineFilter.h"

namespace SmartCoronary
{
	const LVCorLandMark LVCorLandmarkTable[] = { { LEFT_VENTRICLE_APEX	,  0, "LEFT_VENTRICLE_APEX" },
	{ BASE_LEFT_END			, 85, "BASE_LEFT_END" },
	{ BASE_RIGHT_END			, 72, "BASE_RIGHT_END" },
	{ BASE_ANTERIOR_END		, 84, "BASE_ANTERIOR_END" },
	{ BASE_POSTERIOR_END		, 76, "BASE_POSTERIOR_END" },
	{ LEFT_CORONARY_OSTIUM  , -1, "LEFT_CORONARY_OSTIUM" },
	{ RIGHT_CORONARY_OSTIUM , -1, "RIGHT_CORONARY_OSTIUM" },
	};
}

void GetRotationMatrix(double axis[3], double angle, double rot[3][3])
{
	double c = cos(angle);
	double s = sin(angle);
	rot[0][0] = c + axis[0] * axis[0] * (1.0 - c);
	rot[0][1] = axis[0] * axis[1] * (1.0 - c) - s*axis[2];
	rot[0][2] = axis[0] * axis[2] * (1.0 - c) + s*axis[1];
	rot[1][0] = axis[0] * axis[1] * (1.0 - c) + s*axis[2];
	rot[1][1] = c + axis[1] * axis[1] * (1.0 - c);
	rot[1][2] = axis[1] * axis[2] * (1.0 - c) - s*axis[0];
	rot[2][0] = axis[0] * axis[2] * (1.0 - c) - s*axis[1];
	rot[2][1] = axis[1] * axis[2] * (1.0 - c) + s*axis[0];
	rot[2][2] = c + axis[2] * axis[2] * (1.0 - c);
}

void ImageGradient(vtkImageInterpolator* interpolator, const double point[3], double gradient[3], double scale = 1.0)
{
	double coord[3];
	double grad1, grad2;
	for (int j = 0; j<3; j++)
	{
		double axis[3] = { 0, 0, 0 };
		axis[j] = scale / 2.0;
		vtkMath::Add(point, axis, coord);
		interpolator->Interpolate(coord, &grad1);
		vtkMath::Subtract(point, axis, coord);
		interpolator->Interpolate(coord, &grad2);
		gradient[j] = grad1 - grad2;
	}
}

void ImageFeatures(vtkImageInterpolator* interpolator, const double point[3], const double normal[3], double thickness, cv::Mat& features)
{
	double coord[3];
	double intensity, gradient[3], gradmag;
	for (int i = 0; i<5; i++)
	{
		for (int j = 0; j<3; j++) coord[j] = thickness*normal[j] * (i - 2);
		vtkMath::Add(coord, point, coord);
		interpolator->Interpolate(coord, &intensity);
		ImageGradient(interpolator, coord, gradient);
		gradmag = vtkMath::Norm(gradient);
		features.at<float>(12 * i) = (float)intensity;
		features.at<float>(12 * i + 1) = (float)sqrt(abs(intensity));
		features.at<float>(12 * i + 2) = (float)intensity*intensity;
		features.at<float>(12 * i + 3) = (float)intensity*intensity*intensity;
		features.at<float>(12 * i + 4) = (float)gradmag;
		features.at<float>(12 * i + 5) = (float)sqrt(gradmag);
		features.at<float>(12 * i + 6) = (float)gradmag*gradmag;
		features.at<float>(12 * i + 7) = (float)gradmag*gradmag*gradmag;
		features.at<float>(12 * i + 8) = (float)gradient[0];
		features.at<float>(12 * i + 9) = (float)gradient[1];
		features.at<float>(12 * i + 10) = (float)gradient[2];
		features.at<float>(12 * i + 11) = (float)vtkMath::Dot(gradient, normal);
	}
}

void FillIntegralImage(vtkImageData* intergalImage, vtkImageData *imageData, vtkImageInterpolator* interpolator)
{
	int	 imageDims[3];
	double imageOrigins[3];
	double imageSpacings[3];

	imageData->GetDimensions(imageDims);
	imageData->GetOrigin(imageOrigins);
	imageData->GetSpacing(imageSpacings);

	intergalImage->SetExtent(0, (int)(imageDims[0] * imageSpacings[0]) + 72, 0, (int)(imageDims[1] * imageSpacings[1]) + 72, 0, (int)(imageDims[2] * imageSpacings[2]) + 72);
	intergalImage->SetOrigin(imageOrigins[0] - 36.0, imageOrigins[1] - 36.0, imageOrigins[2] - 36.0);
	intergalImage->SetSpacing(1.0, 1.0, 1.0);
	intergalImage->AllocateScalars(VTK_DOUBLE, 6);

	int* dims = intergalImage->GetDimensions();
	double* origin = intergalImage->GetOrigin();
	double* spacing = intergalImage->GetSpacing();
	double coord[3];
	int dim[3];
	int pim[7][3];
	double pixel;
	int contri;
	int sign[][4] = { { -1, 0, 0, 1 },{ 0, -1, 0, 1 },{ -1, -1, 0, -1 },{ 0, 0, -1, 1 },{ -1, 0, -1, -1 },{ 0, -1, -1, -1 },{ -1, -1, -1, 1 } };

	double *integralimage = static_cast<double *>(intergalImage->GetScalarPointer());
	int dims01 = dims[0] * dims[1] * 6;
	int dims0 = dims[0] * 6;

	for (dim[2] = 0; dim[2] < dims[2]; dim[2]++)
	{
		for (dim[1] = 0; dim[1] < dims[1]; dim[1]++)
		{
			for (dim[0] = 0; dim[0] < dims[0]; dim[0]++)
			{
				for (int k = 0; k<3; k++) coord[k] = origin[k] + spacing[k] * dim[k];
				interpolator->Interpolate(coord, &pixel);
				if (pixel<-512.0) contri = 0;
				else if (pixel<   0.0) contri = 1;
				else if (pixel< 256.0) contri = 2;
				else if (pixel< 512.0) contri = 3;
				else if (pixel< 768.0) contri = 4;
				else				  contri = 5;
				double* currhist = &integralimage[dim[2] * dims01 + dim[1] * dims0 + dim[0] * 6];
				memset(currhist, 0, 6 * sizeof(double));
				currhist[contri] = 1.0;

				for (int j = 0; j<7; j++)
				{
					bool neg = false;
					for (int k = 0; k<3; k++)
					{
						pim[j][k] = dim[k] + sign[j][k];
						neg = pim[j][k]<0;
						if (neg) break;
					}
					if (neg) continue;

					double* prevhist = &integralimage[pim[j][2] * dims01 + pim[j][1] * dims0 + pim[j][0] * 6];
					for (int i = 0; i<6; i++)
					{
						currhist[i] += sign[j][3] * prevhist[i];
					}
				}
			}
		}
	}
}

void IntegralImageHist(vtkImageData* intergalImage, int corner1[3], int corner2[3], double hist[6])
{
	int dim[3];
	double* inthist;
	for (int i = 0; i<6; i++) hist[i] = 0.0;
	int sign[][4] = { { 0, 0, 0, 1 },{ -1, 0, 0, -1 },{ 0, -1, 0, -1 },{ -1, -1, 0, 1 },{ 0, 0, -1, -1 },{ -1, 0, -1, 1 },{ 0, -1, -1, 1 },{ -1, -1, -1, -1 } };
	for (int d = 0; d<8; d++)
	{
		dim[0] = sign[d][0] == 0 ? corner2[0] : corner1[0];
		dim[1] = sign[d][1] == 0 ? corner2[1] : corner1[1];
		dim[2] = sign[d][2] == 0 ? corner2[2] : corner1[2];
		inthist = static_cast<double*>(intergalImage->GetScalarPointer(dim));
		for (int i = 0; i<6; i++) hist[i] += sign[d][3] * inthist[i];
	}
}

void HistNormalize(double hist[6])
{
	double mag = 0.0;
	for (int i = 0; i<6; i++) mag += hist[i] * hist[i];
	mag = std::max(sqrt(mag), 0.001);
	for (int i = 0; i<6; i++) hist[i] /= mag;
}

void HistSubtract(const double hist1[6], const double hist2[6], double hist3[6])
{
	for (int i = 0; i<6; i++) hist3[i] = hist1[i] - hist2[i];
}

void InsertFeature(cv::Mat& featureRow, double hist[6], int& count)
{
	HistNormalize(hist);
	for (int i = 0; i<6; i++)
	{
		featureRow.at<float>(count) = hist[i];
		count++;
	}
}

void ImageFeatures(vtkImageData* intergalImage, double coord[3], cv::Mat& featureRow)
{
	double* origin = intergalImage->GetOrigin();
	double* spacing = intergalImage->GetSpacing();

	int dim[3];
	for (int i = 0; i<3; i++) dim[i] = (int)((coord[i] - origin[i]) / spacing[i]);
	int corner1[3], corner2[3];

	double hist1[6], hist2[6];
	int count = 0;
	for (int i = 1; i <= 3; i++)
	{
		int cellsize = i * 8;

		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist1);
		InsertFeature(featureRow, hist1, count);

		corner1[0] = dim[0] - cellsize; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0]; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist1);
		corner1[0] = dim[0]; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist2, hist1, hist1);
		InsertFeature(featureRow, hist1, count);

		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1]; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist1);
		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1]; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist2, hist1, hist1);
		InsertFeature(featureRow, hist1, count);

		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2];
		IntegralImageHist(intergalImage, corner1, corner2, hist1);
		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2];
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist2, hist1, hist1);
		InsertFeature(featureRow, hist1, count);

		corner1[0] = dim[0] - cellsize * 3 / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] - cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist1);
		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist2, hist1, hist1);
		corner1[0] = dim[0] + cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize * 3 / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist1, hist2, hist1);
		InsertFeature(featureRow, hist1, count);

		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize * 3 / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] - cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist1);
		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist2, hist1, hist1);
		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] + cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize * 3 / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist1, hist2, hist1);
		InsertFeature(featureRow, hist1, count);

		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize * 3 / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] - cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist1);
		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist2, hist1, hist1);
		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] + cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize * 3 / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist1, hist2, hist1);
		InsertFeature(featureRow, hist1, count);
	}
}

void SmoothMesh1(vtkPolyData *poly, int iterations = 2)
{
	if (!poly) return;

	vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
	smoother->SetNumberOfIterations(iterations);
	smoother->BoundarySmoothingOn();
	smoother->FeatureEdgeSmoothingOff();
	smoother->NonManifoldSmoothingOn();
	smoother->NormalizeCoordinatesOn();
	smoother->GenerateErrorScalarsOff();
	smoother->GenerateErrorVectorsOff();
	smoother->SetInputData(poly);
	smoother->Update();

	poly->DeepCopy(smoother->GetOutput());
}

void SmoothMesh2(vtkPolyData *poly, int iterations = 10)
{
	if (!poly) return;

	vtkSmartPointer<vtkSmoothPolyDataFilter> smoother = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoother->SetNumberOfIterations(iterations);
	smoother->BoundarySmoothingOn();
	smoother->FeatureEdgeSmoothingOff();
	smoother->GenerateErrorScalarsOff();
	smoother->GenerateErrorVectorsOff();
	smoother->SetInputData(poly);
	smoother->Update();

	poly->DeepCopy(smoother->GetOutput());
}

void SavePolyData(vtkPolyData *poly, const char* fileName)
{
	if (!poly) return;
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(poly);
	writer->SetFileName(fileName);
	writer->SetDataModeToBinary();
	try
	{
		writer->Write();
	}
	catch (...)
	{
		std::cerr << "Error occurs when writing " << fileName << std::endl;
		return;
	}
}

/* allocate memory for an nrow x ncol matrix */
template<class TReal>
TReal **create_matrix(long nrow, long ncol)
{
	typedef TReal* TRealPointer;
	TReal **m = new TRealPointer[nrow];

	TReal* block = (TReal*)calloc(nrow*ncol, sizeof(TReal));
	m[0] = block;
	for (int row = 1; row < nrow; ++row)
	{
		m[row] = &block[row * ncol];
	}
	return m;
}

/* free a TReal matrix allocated with matrix() */
template<class TReal>
void free_matrix(TReal **m)
{
	free(m[0]);
	delete[] m;
}

void FillSumImage(vtkImageData* sumImage, vtkImageInterpolator* interpolator)
{
	int* dims = sumImage->GetDimensions();
	double* origin = sumImage->GetOrigin();
	double* spacing = sumImage->GetSpacing();
	double coord[3];
	int dim[3];
	int pim[3];
	double pixel;
	int sign[][4] = { { -1, 0, 0, 1 },{ 0, -1, 0, 1 },{ -1, -1, 0, -1 },{ 0, 0, -1, 1 },{ -1, 0, -1, -1 },{ 0, -1, -1, -1 },{ -1, -1, -1, 1 } };

	int dims01 = dims[0] * dims[1];
	double *image = static_cast<double*>(sumImage->GetScalarPointer());

	for (dim[2] = 0; dim[2] < dims[2]; dim[2]++)
	{
		for (dim[1] = 0; dim[1] < dims[1]; dim[1]++)
		{
			for (dim[0] = 0; dim[0] < dims[0]; dim[0]++)
			{
				for (int k = 0; k<3; k++) coord[k] = origin[k] + spacing[k] * dim[k];
				interpolator->Interpolate(coord, &pixel);
				double& curr = image[dim[2] * dims01 + dim[1] * dims[0] + dim[0]];
				curr = pixel;

				for (int j = 0; j<7; j++)
				{
					bool neg = false;
					for (int k = 0; k<3; k++)
					{
						pim[k] = dim[k] + sign[j][k];
						neg = pim[k]<0;
						if (neg) break;
					}
					if (neg) continue;

					curr += sign[j][3] * image[pim[2] * dims01 + pim[1] * dims[0] + pim[0]];
				}
			}
		}
	}
}

void SumImageHist(vtkImageData* sumImage, double* sumimage, int corner1[3], int corner2[3], double& sum)
{
	int* dims = sumImage->GetDimensions();
	int dims01 = dims[0] * dims[1];

	int dim[3];
	sum = 0.0;
	int sign[][4] = { { 0, 0, 0, 1 },{ -1, 0, 0, -1 },{ 0, -1, 0, -1 },{ -1, -1, 0, 1 },{ 0, 0, -1, -1 },{ -1, 0, -1, 1 },{ 0, -1, -1, 1 },{ -1, -1, -1, -1 } };
	for (int d = 0; d<8; d++)
	{
		dim[0] = sign[d][0] == 0 ? corner2[0] : (corner1[0] - 1);
		dim[1] = sign[d][1] == 0 ? corner2[1] : (corner1[1] - 1);
		dim[2] = sign[d][2] == 0 ? corner2[2] : (corner1[2] - 1);
		sum += sign[d][3] * sumimage[dim[2] * dims01 + dim[1] * dims[0] + dim[0]];
	}
}

void Hessian(vtkImageData* sumImage, double* sumimage, double coord[3], double eigvalue[3], double eigvector[3][3], int cellsize = 3)
{
	double* origin = sumImage->GetOrigin();
	double* spacing = sumImage->GetSpacing();

	int dim[3];
	for (int i = 0; i<3; i++) dim[i] = (int)((coord[i] - origin[i]) / spacing[i] + 0.5);
	int corner1[3], corner2[3];

	int cellsize3 = cellsize*cellsize*cellsize;
	int cellsizem3d2 = cellsize * 3 / 2;
	int cellsized2 = cellsize / 2;

	double sum1, sum2;
	double **hessian = create_matrix<double>(3, 3);
	//double hessian[3][3];
	{
		//dx^2
		corner1[0] = dim[0] - cellsizem3d2; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] - cellsized2 - 1; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum1);
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 - 2.0*sum2;
		corner1[0] = dim[0] + cellsized2 + 1; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] + cellsizem3d2; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		hessian[0][0] = (sum1 + sum2) / cellsize3;

		//dy^2
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] - cellsizem3d2; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] - cellsized2 - 1; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum1);
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 - 2.0*sum2;
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] + cellsized2 + 1; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] + cellsizem3d2; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		hessian[1][1] = (sum1 + sum2) / cellsize3;

		//dz^2
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] - cellsizem3d2;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] - cellsized2 - 1;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum1);
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 - 2.0*sum2;
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] + cellsized2 + 1;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] + cellsizem3d2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		hessian[2][2] = (sum1 + sum2) / cellsize3;

		//dxdy
		corner1[0] = dim[0] - cellsize + 1; corner1[1] = dim[1] - cellsize + 1; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0]; corner2[1] = dim[1]; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum1);
		corner1[0] = dim[0]; corner1[1] = dim[1]; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] + cellsize - 1; corner2[1] = dim[1] + cellsize - 1; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 + sum2;
		corner1[0] = dim[0] - cellsize + 1; corner1[1] = dim[1]; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0]; corner2[1] = dim[1] + cellsize - 1; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 - sum2;
		corner1[0] = dim[0]; corner1[1] = dim[1] - cellsize + 1; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] + cellsize - 1; corner2[1] = dim[1]; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		hessian[0][1] = (sum1 - sum2) / cellsize3;
		hessian[1][0] = hessian[0][1];

		//dxdz
		corner1[0] = dim[0] - cellsize + 1; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] - cellsize + 1;
		corner2[0] = dim[0]; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2];
		SumImageHist(sumImage, sumimage, corner1, corner2, sum1);
		corner1[0] = dim[0]; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2];
		corner2[0] = dim[0] + cellsize - 1; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] + cellsize - 1;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 + sum2;
		corner1[0] = dim[0] - cellsize + 1; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2];
		corner2[0] = dim[0]; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] + cellsize - 1;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 - sum2;
		corner1[0] = dim[0]; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] - cellsize + 1;
		corner2[0] = dim[0] + cellsize - 1; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2];
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		hessian[0][2] = (sum1 - sum2) / cellsize3;
		hessian[2][0] = hessian[0][2];

		//dydz
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] - cellsize + 1; corner1[2] = dim[2] - cellsize + 1;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1]; corner2[2] = dim[2];
		SumImageHist(sumImage, sumimage, corner1, corner2, sum1);
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1]; corner1[2] = dim[2];
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] + cellsize - 1; corner2[2] = dim[2] + cellsize - 1;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 + sum2;
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] - cellsize + 1; corner1[2] = dim[2];
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1]; corner2[2] = dim[2] + cellsize - 1;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 - sum2;
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1]; corner1[2] = dim[2] - cellsize + 1;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] + cellsize - 1; corner2[2] = dim[2];
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		hessian[1][2] = (sum1 - sum2) / cellsize3;
		hessian[2][1] = hessian[1][2];
	}
	double **eigvec = create_matrix<double>(3, 3);
	vtkMath::Jacobi(hessian, eigvalue, eigvec);
	for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			eigvector[i][j] = eigvec[i][j];

	free_matrix(hessian);
	free_matrix(eigvec);
}

template<class BinaryImageType>
void TraceCenterlineInternal(std::queue< std::pair<typename BinaryImageType::IndexType, typename BinaryImageType::IndexType> >& clqueue,
	typename BinaryImageType::Pointer rawCenterline, vtkPoints *clPoints, vtkPolyData* clModel)
{
	std::map<typename BinaryImageType::IndexType, vtkIdType, typename BinaryImageType::IndexType::LexicographicCompare> indexMap;
	typedef itk::NeighborhoodIterator< BinaryImageType > NeighborhoodIteratorType;
	typename NeighborhoodIteratorType::RadiusType nRadius; nRadius.Fill(1);
	NeighborhoodIteratorType labelIt(nRadius, rawCenterline, rawCenterline->GetRequestedRegion());

	bool inbounds;

	while (!clqueue.empty())
	{
		std::list<typename BinaryImageType::IndexType> indexList;

		typename BinaryImageType::IndexType startid = clqueue.front().first;
		indexList.push_back(startid);
		typename BinaryImageType::IndexType nextid = clqueue.front().second;
		clqueue.pop();

		std::queue<typename BinaryImageType::IndexType> sgqueue;
		sgqueue.push(nextid);
		while (!sgqueue.empty())
		{
			nextid = sgqueue.front();
			sgqueue.pop();
			indexList.push_back(nextid);
			labelIt.SetLocation(nextid);
			for (size_t i = 0; i<labelIt.Size(); i++)
			{
				if (labelIt.GetPixel(i, inbounds) > 0 && inbounds)
				{
					labelIt.SetPixel(i, 0);
					sgqueue.push(labelIt.GetIndex(i));
				}
			}
			if (sgqueue.size() != 1)
			{
				if (indexList.size() > 10 || !sgqueue.empty())
				{
					vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
					for (auto lit = indexList.begin(); lit != indexList.end(); ++lit)
					{
						if (indexMap.find(*lit) == indexMap.end())
						{
							typename BinaryImageType::PointType point;
							rawCenterline->TransformIndexToPhysicalPoint(*lit, point);
							indexMap[*lit] = clPoints->InsertNextPoint(point[0], point[1], point[2]);
						}
						idlist->InsertNextId(indexMap[*lit]);
					}
					clModel->InsertNextCell(VTK_POLY_LINE, idlist);
				}

				while (!sgqueue.empty())
				{
					clqueue.push(std::make_pair(nextid, sgqueue.front()));
					sgqueue.pop();
				}
			}
		}
	}
}

// 从初始中心线结果中检测
template<class BinaryImageType>
void TraceCenterline(typename BinaryImageType::Pointer rawCenterline, const typename BinaryImageType::IndexType& ostiumIndex, vtkPolyData* clModel)
{
	vtkSmartPointer<vtkPoints> clPoints = vtkSmartPointer<vtkPoints>::New();

	typedef itk::NeighborhoodIterator< BinaryImageType > NeighborhoodIteratorType;
	NeighborhoodIteratorType::RadiusType nRadius; nRadius.Fill(1);
	NeighborhoodIteratorType labelIt(nRadius, rawCenterline, rawCenterline->GetRequestedRegion());

	bool inbounds;

	if (ostiumIndex[0] >= 0 && ostiumIndex[1] >= 0 && ostiumIndex[2] >= 0)		//find ostium seeded centerlines
	{
		labelIt.SetLocation(ostiumIndex);
		if (labelIt.GetCenterPixel() <= 0)	return;
		labelIt.SetCenterPixel(-1);

		std::vector<typename BinaryImageType::IndexType> initdirs;
		std::vector<int>									 initsizes;
		std::vector< std::vector<typename BinaryImageType::IndexType> > initids;
		for (size_t i = 0; i<labelIt.Size(); i++)
		{
			if (labelIt.GetPixel(i, inbounds) > 0 && inbounds)
			{
				initdirs.push_back(labelIt.GetIndex(i));
				initsizes.push_back(0);
				initids.push_back(std::vector<typename BinaryImageType::IndexType>());
			}
		}
		for (size_t i = 0; i<initdirs.size(); i++)
		{
			std::queue<typename BinaryImageType::IndexType> queue;
			queue.push(initdirs[i]);
			while (!queue.empty())
			{
				typename BinaryImageType::IndexType pid = queue.front();
				queue.pop();
				labelIt.SetLocation(pid);

				labelIt.SetCenterPixel(-1);
				initsizes[i]++;
				initids[i].push_back(pid);
				for (size_t i = 0; i<labelIt.Size(); i++)
				{
					if (labelIt.GetPixel(i, inbounds) > 0 && inbounds) queue.push(labelIt.GetIndex(i));
				}
			}
		}

		auto maxdir = std::distance(initsizes.begin(), std::max_element(initsizes.begin(), initsizes.end()));
		for (size_t i = 0; i<initdirs.size(); i++)
		{
			if (i == maxdir)
			{
				for (size_t j = 0; j<initids[i].size(); j++)
				{
					labelIt.SetLocation(initids[i][j]);
					labelIt.SetCenterPixel(1);
				}
			}
			else
			{
				for (size_t j = 0; j<initids[i].size(); j++)
				{
					labelIt.SetLocation(initids[i][j]);
					labelIt.SetCenterPixel(0);
				}
			}
		}

		labelIt.SetLocation(ostiumIndex);
		labelIt.SetCenterPixel(0);
		labelIt.SetLocation(initdirs[maxdir]);
		labelIt.SetCenterPixel(0);
		std::queue< std::pair<typename BinaryImageType::IndexType, typename BinaryImageType::IndexType> > clqueue;
		clqueue.push(std::make_pair(ostiumIndex, initdirs[maxdir]));
		TraceCenterlineInternal<BinaryImageType>(clqueue, rawCenterline, clPoints, clModel);
	}
	else
	{
		for (labelIt.GoToBegin(); !labelIt.IsAtEnd(); ++labelIt)
		{
			if (labelIt.GetCenterPixel() > 0)
			{
				for (size_t i = 0; i<labelIt.Size(); i++)
				{
					if (labelIt.GetPixel(i, inbounds) > 0 && inbounds)
					{
						labelIt.SetLocation(labelIt.GetIndex());
						labelIt.SetCenterPixel(0);
						labelIt.SetLocation(labelIt.GetIndex(i));
						labelIt.SetCenterPixel(0);
						std::queue< std::pair<typename BinaryImageType::IndexType, typename BinaryImageType::IndexType> > clqueue;
						clqueue.push(std::make_pair(labelIt.GetIndex(), labelIt.GetIndex(i)));
						TraceCenterlineInternal<BinaryImageType>(clqueue, rawCenterline, clPoints, clModel);
					}
				}
			}
		}
	}

	clModel->SetPoints(clPoints);
}

void CleanCenterline(vtkPolyData* clModel)
{
	if (!clModel || clModel->GetNumberOfCells() == 0) return;

	clModel->BuildCells();
	clModel->BuildLinks();

	unsigned short ncells;
	vtkIdType	*cells;
	vtkIdType numcells = clModel->GetNumberOfCells();
	std::vector<vtkIdType> deleteCells;
	for (vtkIdType i = 0; i<numcells; i++)
	{
		vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
		clModel->GetCellPoints(i, idlist);
		vtkSmartPointer<vtkIdList> newlist = vtkSmartPointer<vtkIdList>::New();

		bool split = false;
		for (vtkIdType j = 0; j<idlist->GetNumberOfIds(); j++)
		{
			newlist->InsertNextId(idlist->GetId(j));

			clModel->GetPointCells(idlist->GetId(j), ncells, cells);
			if (ncells>1 && j != 0 && j != idlist->GetNumberOfIds() - 1)
			{
				clModel->InsertNextCell(VTK_POLY_LINE, newlist);
				newlist->Initialize();
				newlist->InsertNextId(idlist->GetId(j));
				split = true;
			}
		}
		if (split)
		{
			clModel->InsertNextCell(VTK_POLY_LINE, newlist);
			deleteCells.push_back(i);
		}
	}
	if (deleteCells.size()>0)
	{
		for (size_t i = 0; i<deleteCells.size(); i++)	clModel->DeleteCell(deleteCells[i]);
		clModel->RemoveDeletedCells();
		clModel->BuildCells();
		clModel->BuildLinks();
	}

	std::vector<int> newIds(clModel->GetNumberOfCells(), -1);
	int newid = 0;
	for (vtkIdType i = 0; i<clModel->GetNumberOfPoints(); i++)
	{
		vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
		clModel->GetPointCells(i, idlist);
		if (idlist->GetNumberOfIds() == 2)
		{
			int& id1 = newIds[idlist->GetId(0)];
			int& id2 = newIds[idlist->GetId(1)];
			if (id1 < 0 && id2 < 0)
			{
				id1 = newid;
				id2 = newid;
				newid++;
			}
			else if (id1 >= 0 && id2 >= 0)
			{
				int oldid2 = id2, oldid1 = id1;
				if (oldid1 < oldid2)
				{

					for (int k = 0; k<newIds.size(); k++) if (newIds[k] == oldid2) newIds[k] = oldid1;
				}
				else if (oldid1 > oldid2)
				{
					for (int k = 0; k<newIds.size(); k++) if (newIds[k] == oldid1) newIds[k] = oldid2;
				}
			}
			else
			{
				if (id1 < id2) id1 = id2;
				else            id2 = id1;
			}
		}
	}
	vtkSmartPointer<vtkPolyData> newModel = vtkSmartPointer<vtkPolyData>::New();
	newModel->Allocate();
	for (vtkIdType i = 0; i<clModel->GetNumberOfCells(); i++)
	{
		if (newIds[i] == -1)
		{
			vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
			clModel->GetCellPoints(i, idlist);
			newModel->InsertNextCell(VTK_POLY_LINE, idlist);
		}
		else if (newIds[i] >= 0)
		{
			int newid = newIds[i];
			std::set<vtkIdType> ids;

			for (int k = 0; k<newIds.size(); k++)
			{
				if (newIds[k] == newid)
				{
					ids.insert(k);
				}
			}

			std::deque<vtkIdType> pidlist;
			vtkSmartPointer<vtkIdList> idl = vtkSmartPointer<vtkIdList>::New();
			clModel->GetCellPoints(*ids.begin(), idl);
			for (vtkIdType l = 0; l<idl->GetNumberOfIds(); l++) pidlist.push_back(idl->GetId(l));
			newIds[*ids.begin()] = -2;
			ids.erase(ids.begin());
			std::queue<vtkIdType> fqueue, equeue;
			fqueue.push(idl->GetId(0));
			equeue.push(idl->GetId(idl->GetNumberOfIds() - 1));
			while (!fqueue.empty())
			{
				vtkIdType fid = fqueue.front();
				fqueue.pop();
				for (auto iter = ids.begin(); iter != ids.end(); iter++)
				{
					vtkSmartPointer<vtkIdList> idl = vtkSmartPointer<vtkIdList>::New();
					clModel->GetCellPoints(*iter, idl);
					if (idl->GetId(0) == fid)
					{
						for (vtkIdType l = 1; l<idl->GetNumberOfIds(); l++) pidlist.push_front(idl->GetId(l));
						fqueue.push(idl->GetId(idl->GetNumberOfIds() - 1));
						newIds[*iter] = -2;
						ids.erase(iter);
						break;
					}
					else if (idl->GetId(idl->GetNumberOfIds() - 1) == fid)
					{
						for (vtkIdType l = idl->GetNumberOfIds() - 2; l >= 0; l--) pidlist.push_front(idl->GetId(l));
						fqueue.push(idl->GetId(0));
						newIds[*iter] = -2;
						ids.erase(iter);
						break;
					}
				}
			}
			while (!equeue.empty())
			{
				vtkIdType fid = equeue.front();
				equeue.pop();
				for (auto iter = ids.begin(); iter != ids.end(); iter++)
				{
					vtkSmartPointer<vtkIdList> idl = vtkSmartPointer<vtkIdList>::New();
					clModel->GetCellPoints(*iter, idl);
					if (idl->GetId(0) == fid)
					{
						for (vtkIdType l = 1; l<idl->GetNumberOfIds(); l++) pidlist.push_back(idl->GetId(l));
						equeue.push(idl->GetId(idl->GetNumberOfIds() - 1));
						newIds[*iter] = -2;
						ids.erase(iter);
						break;
					}
					else if (idl->GetId(idl->GetNumberOfIds() - 1) == fid)
					{
						for (vtkIdType l = idl->GetNumberOfIds() - 2; l >= 0; l--) pidlist.push_back(idl->GetId(l));
						equeue.push(idl->GetId(0));
						newIds[*iter] = -2;
						ids.erase(iter);
						break;
					}
				}
			}
			vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
			for (size_t p = 0; p<pidlist.size(); p++) idlist->InsertNextId(pidlist[p]);
			newModel->InsertNextCell(VTK_POLY_LINE, idlist);
		}
	}
	newModel->SetPoints(clModel->GetPoints());
	newModel->GetPointData()->CopyAllOn();
	newModel->GetPointData()->PassData(clModel->GetPointData());
	clModel->DeepCopy(newModel);
}

void SimplifyCenterline(vtkPolyData* clModel)
{
	if (!clModel || clModel->GetNumberOfCells() == 0) return;

	vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyData> newModel = vtkSmartPointer<vtkPolyData>::New();
	newModel->Allocate();
	std::map<vtkIdType, vtkIdType> pointMap;
	for (vtkIdType i = 0; i < clModel->GetNumberOfCells(); i++)
	{
		vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
		clModel->GetCellPoints(i, idlist);
		if (idlist->GetNumberOfIds() < 3) continue;

		std::list<double> newpolyline;
		unsigned count = 0;
		double coord[3], last[3];
		double length = 0.0;
		for (vtkIdType l = 0; l < idlist->GetNumberOfIds(); l++) {
			clModel->GetPoint(idlist->GetId(l), coord);

			if (l == 0) {
				for (int kk = 0; kk < 3; kk++) newpolyline.push_back(coord[kk]);
				count = 0;
			}
			else if (l == idlist->GetNumberOfIds() - 1) {
				for (int kk = 0; kk < 3; kk++) newpolyline.push_back(coord[kk]);
				count++;
			}
			else {
				length += sqrt(vtkMath::Distance2BetweenPoints(coord, last));

				if (unsigned(length / 3) != count) {
					for (int kk = 0; kk < 3; kk++) newpolyline.push_back(coord[kk]);
					count++;
				}
			}

			for (int k = 0; k<3; k++) last[k] = coord[k];
		}

	//	if (newpolyline.empty()) {
	//		continue;
	//	}

		count++;
	//	vtkSmartPointer<vtkIdList> newlist = vtkSmartPointer<vtkIdList>::New();
	//	while(!newpolyline.empty()) {
	//		for (int k = 0; k < 3; k++) {
	//			coord[k] = newpolyline.front();
	//			newpolyline.pop_front();
	//		}
	//		vtkIdType newid = newPoints->InsertNextPoint(coord);
	//		newlist->InsertNextId(newid);
	//	}

	//	newModel->InsertNextCell(VTK_POLY_LINE, newlist);
	//}
	//newModel->SetPoints(newPoints);

	// smartcoronary
	//std::map<vtkIdType, vtkIdType> pointMap;
	//for (vtkIdType i = 0; i<clModel->GetNumberOfCells(); i++)
	//{
	//	vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
	//	clModel->GetCellPoints(i, idlist);
	//	if (idlist->GetNumberOfIds() < 3) continue;

	//	std::list<double> polyline;
	//	double coord[3], last[3];
	//	double length = 0.0;
	//	for (vtkIdType l = 0; l<idlist->GetNumberOfIds(); l++)
	//	{
	//		clModel->GetPoint(idlist->GetId(l), coord);
	//		if (l>0) length += sqrt(vtkMath::Distance2BetweenPoints(coord, last));
	//		for (int k = 0; k<3; k++)
	//		{
	//			polyline.push_back(coord[k]);
	//			last[k] = coord[k];
	//		}
	//	}
	//	unsigned count = unsigned(length / 5.0);
	//	if (count<2) count = 2;
	//	else if (count>idlist->GetNumberOfIds()) count = idlist->GetNumberOfIds();
	//	std::vector<double> newpolyline;
	//	psimpl::simplify_douglas_peucker_n<3>(polyline.begin(), polyline.end(), count, std::back_inserter(newpolyline));
	//	if (newpolyline.size() != 3 * count)
	//	{
	//		std::cerr << "Centerline simplification results wrong number of points" << std::endl;
	//		return;
	//	}

		//double begin[3], end[3], newbegin[3], newend[3];
		//clModel->GetPoint(idlist->GetId(0), begin);
		//clModel->GetPoint(idlist->GetId(idlist->GetNumberOfIds() - 1), end);
		//std::copy(newpolyline.begin(), newpolyline.begin() + 3, newbegin);
		//std::copy(newpolyline.end() - 3, newpolyline.end(), newend);
		//if (vtkMath::Distance2BetweenPoints(newbegin, begin) > 1e-6 || vtkMath::Distance2BetweenPoints(newend, end) > 1e-6)
		//{
		//	std::cerr << "Centerline simplification results wrong begin and end points" << std::endl;
		//	return;
		//}

		vtkSmartPointer<vtkIdList> newlist = vtkSmartPointer<vtkIdList>::New();
		for (size_t j = 0; j<count; j++)
		{
			//for (int k = 0; k<3; k++) coord[k] = newpolyline[3 * j + k];
			for (int k = 0; k < 3; k++) {
				coord[k] = newpolyline.front();
				newpolyline.pop_front();
			}
			if (j == 0)
			{
				if (pointMap.find(idlist->GetId(0)) != pointMap.end())
				{
					newlist->InsertNextId(pointMap[idlist->GetId(0)]);
				}
				else
				{
					vtkIdType newid = newPoints->InsertNextPoint(coord);
					newlist->InsertNextId(newid);
					pointMap[idlist->GetId(0)] = newid;
				}
			}
			else if (j == count - 1)
			{
				if (pointMap.find(idlist->GetId(idlist->GetNumberOfIds() - 1)) != pointMap.end())
				{
					newlist->InsertNextId(pointMap[idlist->GetId(idlist->GetNumberOfIds() - 1)]);
				}
				else
				{
					vtkIdType newid = newPoints->InsertNextPoint(coord);
					newlist->InsertNextId(newid);
					pointMap[idlist->GetId(idlist->GetNumberOfIds() - 1)] = newid;
				}
			}
			else
			{
				vtkIdType newid = newPoints->InsertNextPoint(coord);
				newlist->InsertNextId(newid);
			}
		}
		newModel->InsertNextCell(VTK_POLY_LINE, newlist);
	}
	newModel->SetPoints(newPoints);

	vtkSmartPointer<ExtendSplineFilter>  clSpline = vtkSmartPointer<ExtendSplineFilter>::New();
	clSpline->SetSubdivideToLength();
	clSpline->SetLength(2.0);
	clSpline->SetGenerateTCoordsToOff();
	clSpline->SetInputData(newModel);
	clSpline->Update();

	clModel->DeepCopy(clSpline->GetOutput());
	//clModel->DeepCopy(newModel);
}

void RadiusCenterline(vtkPolyData* clModel)
{
	vtkSmartPointer<vtkDoubleArray> clRadius = vtkSmartPointer<vtkDoubleArray>::New();
	clRadius->SetName("Radius");
	clRadius->SetNumberOfValues(clModel->GetNumberOfPoints());
	for (vtkIdType id = 0; id<clModel->GetNumberOfPoints(); id++)
	{
		clRadius->SetValue(id, 1.0);
	}
	clModel->GetPointData()->SetScalars(clRadius);

	vtkSmartPointer<vtkIntArray> clRadialExtent = vtkSmartPointer<vtkIntArray>::New();
	clRadialExtent->SetName("RadialExtent");
	clRadialExtent->SetNumberOfValues(clModel->GetNumberOfCells());
	for (vtkIdType id = 0; id<clModel->GetNumberOfCells(); id++)
	{
		clRadialExtent->SetValue(id, 20);
	}
	clModel->GetCellData()->AddArray(clRadialExtent);
}

void LumenWallCenterline(vtkPolyData* clModel, int centerline_components) //centerline_components must be even number
{
	vtkDoubleArray* clRadius = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("Radius"));
	if (!clRadius) return;

	vtkSmartPointer<vtkDoubleArray> clLumenRadius = vtkSmartPointer<vtkDoubleArray>::New();
	clLumenRadius->SetName("LumenRadius");
	clLumenRadius->SetNumberOfComponents(centerline_components);
	clLumenRadius->SetNumberOfTuples(clModel->GetNumberOfPoints());

	vtkSmartPointer<vtkDoubleArray> clWallThickness = vtkSmartPointer<vtkDoubleArray>::New();
	clWallThickness->SetName("WallThickness");
	clWallThickness->SetNumberOfComponents(centerline_components);
	clWallThickness->SetNumberOfTuples(clModel->GetNumberOfPoints());

	vtkSmartPointer<vtkDoubleArray> clLongitudinalAngle = vtkSmartPointer<vtkDoubleArray>::New();
	clLongitudinalAngle->SetName("LongitudinalAngle");
	clLongitudinalAngle->SetNumberOfComponents(centerline_components);
	clLongitudinalAngle->SetNumberOfTuples(clModel->GetNumberOfPoints());

	double *radii = new double[clLumenRadius->GetNumberOfComponents()];
	for (vtkIdType id = 0; id<clModel->GetNumberOfPoints(); id++)
	{
		for (int j = 0; j<clLumenRadius->GetNumberOfComponents(); j++) radii[j] = clRadius->GetValue(id);
		clLumenRadius->SetTuple(id, radii);
		for (int j = 0; j<clWallThickness->GetNumberOfComponents(); j++) radii[j] = 0.2;
		clWallThickness->SetTuple(id, radii);
		for (int j = 0; j<clLongitudinalAngle->GetNumberOfComponents(); j++) radii[j] = 0.0;
		clLongitudinalAngle->SetTuple(id, radii);
	}
	delete[] radii;
	clModel->GetPointData()->AddArray(clLumenRadius);
	clModel->GetPointData()->AddArray(clWallThickness);
	clModel->GetPointData()->AddArray(clLongitudinalAngle);
}

bool DetectCenterlineLumenWall(vtkImageData* labelImage, vtkPolyData* clModel)
{
	vtkDoubleArray *clDir = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("Dir"));
	vtkDoubleArray *clAxis1 = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("Axis1"));
	vtkDoubleArray *clAxis2 = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("Axis2"));
	vtkDoubleArray *clRadius = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("Radius"));
	vtkDoubleArray *clLumenRadius = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("LumenRadius"));
	vtkDoubleArray *clWallThickness = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("WallThickness"));
	vtkDoubleArray *clLongitudinalAngle = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("LongitudinalAngle"));
	 
	if (!clDir || !clAxis1 || !clAxis2 || !clLumenRadius || !clWallThickness || !clLongitudinalAngle) {
		cerr << "Something missing before detecting centerline lumen wall " << endl;
		return false;
	}

	clModel->BuildCells();

	double center[3], direct[3], axis1[3], axis2[3], ray[3];
	double cirstep = 2.0*M_PI / clLumenRadius->GetNumberOfComponents();
	for (vtkIdType ii = 0; ii < clModel->GetNumberOfCells(); ii++) {
		vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
		clModel->GetCellPoints(ii, idList);
		for (vtkIdType jj = 0; jj < idList->GetNumberOfIds(); jj++) {
			vtkIdType pid = idList->GetId(jj);
			clModel->GetPoint(pid, center);
			clDir->GetTuple(pid, direct);
			clAxis1->GetTuple(pid, axis1);
			clAxis2->GetTuple(pid, axis2);

			LabelType* imgPointer = (LabelType*)labelImage->GetScalarPointer();
			double avgRadius = 0;
			for (vtkIdType kk = 0; kk < clLumenRadius->GetNumberOfComponents(); kk++) {
				for (int ll = 0; ll < 3; ll++) ray[ll] = axis1[ll] * cos(kk*cirstep) + axis2[ll] * sin(kk*cirstep);
				for (double length = 0.1; length < 3; length += 0.1) {
					double backPoint[3];
					for (int ll = 0; ll < 3; ll++) backPoint[ll] = length * ray[ll] + center[ll];
					vtkIdType id = labelImage->FindPoint(backPoint);
					if (!imgPointer[id]) {
						clLumenRadius->SetComponent(pid, kk, length + 0.1);
						break;
					}
				}
				avgRadius += clLumenRadius->GetComponent(pid, kk);
			}
			avgRadius /= clLumenRadius->GetNumberOfComponents();
			clRadius->SetValue(pid, avgRadius);
		}
	}

	return true;
}

void AxisCenterline(vtkPolyData* clModel, double planenormal[3])
{
	if (planenormal) vtkMath::Normalize(planenormal);
	vtkSmartPointer<vtkDoubleArray>	clDir = vtkSmartPointer<vtkDoubleArray>::New();
	clDir->SetName("Dir");
	clDir->SetNumberOfComponents(3);
	clDir->SetNumberOfTuples(clModel->GetNumberOfPoints());
	vtkSmartPointer<vtkDoubleArray>	clAxis1 = vtkSmartPointer<vtkDoubleArray>::New();
	clAxis1->SetName("Axis1");
	clAxis1->SetNumberOfComponents(3);
	clAxis1->SetNumberOfTuples(clModel->GetNumberOfPoints());
	vtkSmartPointer<vtkDoubleArray>	clAxis2 = vtkSmartPointer<vtkDoubleArray>::New();
	clAxis2->SetName("Axis2");
	clAxis2->SetNumberOfComponents(3);
	clAxis2->SetNumberOfTuples(clModel->GetNumberOfPoints());
	double coord[3], dir[3], axis1[3], axis2[3], olddir[3], oldaxis1[3];
	double rot[3][3];
	for (vtkIdType id = 0; id<clModel->GetNumberOfCells(); id++)
	{
		vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
		clModel->GetCellPoints(id, idlist);
		for (vtkIdType j = 0; j<idlist->GetNumberOfIds(); j++)
		{
			if (j == 0)
			{
				clModel->GetPoint(idlist->GetId(j), coord);
				clModel->GetPoint(idlist->GetId(j + 1), dir);
				vtkMath::Subtract(dir, coord, dir);
				vtkMath::Normalize(dir);
				if (planenormal)
				{
					for (int k = 0; k<3; k++) axis1[k] = planenormal[k];
					vtkMath::Cross(dir, axis1, axis2);
					vtkMath::Normalize(axis2);
				}
				else
				{
					vtkMath::Perpendiculars(dir, axis1, axis2, 0.0);
				}
			}
			else if (j == idlist->GetNumberOfIds() - 1)
			{
				clModel->GetPoint(idlist->GetId(j - 1), coord);
				clModel->GetPoint(idlist->GetId(j), dir);
				vtkMath::Subtract(dir, coord, dir);
				vtkMath::Normalize(dir);
			}
			else
			{
				clModel->GetPoint(idlist->GetId(j - 1), coord);
				clModel->GetPoint(idlist->GetId(j + 1), dir);
				vtkMath::Subtract(dir, coord, dir);
				vtkMath::Normalize(dir);
			}
			if (j>0)
			{
				vtkMath::Cross(olddir, dir, axis2);
				if (vtkMath::Norm(axis2) == 0.0)
				{
					for (int k = 0; k<3; k++) axis1[k] = oldaxis1[k];
					vtkMath::Cross(dir, axis1, axis2);
					vtkMath::Normalize(axis2);
				}
				else
				{
					vtkMath::Normalize(axis2);
					double angle = acos(vtkMath::Dot(olddir, dir));
					GetRotationMatrix(axis2, angle, rot);
					vtkMath::Multiply3x3(rot, oldaxis1, axis1);
					vtkMath::Normalize(axis1);
					vtkMath::Cross(dir, axis1, axis2);
					vtkMath::Normalize(axis2);
				}
			}
			clDir->SetTuple(idlist->GetId(j), dir);
			clAxis1->SetTuple(idlist->GetId(j), axis1);
			clAxis2->SetTuple(idlist->GetId(j), axis2);
			for (int k = 0; k<3; k++)
			{
				olddir[k] = dir[k];
				oldaxis1[k] = axis1[k];

			}
		}
	}
	clModel->GetPointData()->AddArray(clDir);
	clModel->GetPointData()->AddArray(clAxis1);
	clModel->GetPointData()->AddArray(clAxis2);
}

bool DetectLandmarks(vtkImageData *imageData, Learning& learn, double landmarks[][3], vtkImageInterpolator *interpolator, std::string& workPath)
{
	clock_t timeBegin, timeEnd;

	LearningImpl *learnimpl = learn.limpl;
	cout << endl << "Loading landmark classifier ... ";
	timeBegin = clock();
	learnimpl->LoadLandmarkClassifiers(SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS, workPath);
	timeEnd = clock();
	cout << "Done Time: " << (timeEnd - timeBegin) / 1000 << "." << (timeEnd - timeBegin) % 1000 << "s" << endl;

	cout << endl << "Detecting landmarks ... ";
	timeBegin = clock();

	vtkSmartPointer<vtkImageData> integralImage = vtkSmartPointer<vtkImageData>::New();
	FillIntegralImage(integralImage, imageData, interpolator);

	int	 imageDims[3];
	double imageOrigins[3];
	double imageSpacings[3];

	imageData->GetDimensions(imageDims);
	imageData->GetOrigin(imageOrigins);
	imageData->GetSpacing(imageSpacings);

	double coord[3];
	int dim[3];
	double maxpred[SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS];
	for (int k = 0; k<SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS; k++) maxpred[k] = std::numeric_limits<double>::lowest();

	int imageDims01 = imageDims[0] * imageDims[1];
	short* imagedata = static_cast<short*>(imageData->GetScalarPointer());

	for (dim[2] = 1; dim[2] < imageDims[2] - 1; dim[2]++)
	{
		for (dim[1] = 1; dim[1] < imageDims[1] - 1; dim[1]++)
		{
			for (dim[0] = 1; dim[0] < imageDims[0] - 1; dim[0]++)
			{
				if ((dim[0] % 5 != 0 || dim[1] % 5 != 0 || dim[2] % 5 != 0))
					continue;

				short pixel = imagedata[dim[2] * imageDims01 + dim[1] * imageDims[0] + dim[0]];
				if (pixel < 0)
					continue;

				for (int k = 0; k<3; k++)
					coord[k] = imageOrigins[k] + imageSpacings[k] * dim[k];

				cv::Mat featureRow(1, 126, CV_32F);
				ImageFeatures(integralImage, coord, featureRow);
				for (int id = 0; id<SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS; id++)
				{
					float pred = learnimpl->lmBoost[id].predict(featureRow, cv::Mat(), cv::Range::all(), false, true);
					if (pred>maxpred[id])
					{
						maxpred[id] = pred;
						for (int k = 0; k<3; k++) landmarks[id][k] = coord[k];
					}
				}

			}
		}
	}

	timeEnd = clock();
	cout << "Done Time: " << (timeEnd - timeBegin) / 1000 << "." << (timeEnd - timeBegin) % 1000 << "s" << endl;
	return true;
}

// 检测中心线
bool DetectCenterline(vtkImageData *labelImage, vtkImageData* thinImage, vtkPolyData *centerlineModel, double leftOstium[3], double rightOstium[3])
{
	clock_t timeBegin, timeEnd;

	typedef itk::Image<LabelType, 3> InputImageType;
	typedef itk::Image<short, 3> BinaryImageType;

	typedef itk::VTKImageToImageFilter<InputImageType> LabelToItkFilter;
	LabelToItkFilter::Pointer labelToItkFilter = LabelToItkFilter::New();
	labelToItkFilter->SetInput(labelImage);
	try
	{
		labelToItkFilter->Update();
	}
	catch (...)
	{
		cerr << "Error occurs when trasforming from vtk to itk " << endl;
	}

	//typedef itk::CastImageFilter<ProcessImageType, BinaryImageType> InputToBinaryFilter;
	//InputToBinaryFilter::Pointer inputToBinaryFilter = InputToBinaryFilter::New();
	//inputToBinaryFilter->SetInput(labelToItkFilter->GetOutput());
	//inputToBinaryFilter->Update();
	InputImageType::Pointer itkLabelImage = labelToItkFilter->GetOutput();

	InputImageType::IndexType leftOstiumIndex, rightOstiumIndex;
	InputImageType::PointType leftOstiumPoint, rightOstiumPoint;
	for (int k = 0; k<3; k++) leftOstiumPoint[k] = leftOstium[k];
	for (int k = 0; k<3; k++) rightOstiumPoint[k] = rightOstium[k];
	bool res1 = itkLabelImage->TransformPhysicalPointToIndex(leftOstiumPoint, leftOstiumIndex);
	bool res2 = itkLabelImage->TransformPhysicalPointToIndex(rightOstiumPoint, rightOstiumIndex);
	if (!res1 || !res2) {
		cerr << "left ostium point or right ostium point detected is out of image " << endl;
		return false;
	}
	cout << "leftOstiumIndex:  " << leftOstiumIndex << endl;
	cout << "rightOstiumIndex: " << rightOstiumIndex << endl;

	// 搜索邻域大小为20
	typedef itk::ConstNeighborhoodIterator<InputImageType> NeighborhoodIteratorType;
	NeighborhoodIteratorType::RadiusType nRadius; nRadius.Fill(10);
	NeighborhoodIteratorType labelIt(nRadius, itkLabelImage, itkLabelImage->GetRequestedRegion());
	
	// 左侧
	labelIt.SetLocation(leftOstiumIndex);
	itk::OffsetValueType leftDistanceVessness = std::numeric_limits<itk::OffsetValueType>::max();
	for (size_t i = 0; i<labelIt.Size(); i++)
	{
		short pixel = labelIt.GetPixel(i);
		if (pixel)
		{
			NeighborhoodIteratorType::OffsetType offset = labelIt.GetOffset(i);
			itk::OffsetValueType dist = offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2];
			if (dist<leftDistanceVessness)
			{
				leftOstiumIndex = labelIt.GetIndex(i);
				leftDistanceVessness = dist;
			}
		}
	}
	// 右侧
	labelIt.SetLocation(rightOstiumIndex);
	itk::OffsetValueType rightDistanceVessness = std::numeric_limits<itk::OffsetValueType>::max();
	for (size_t i = 0; i<labelIt.Size(); i++)
	{
		short pixel = labelIt.GetPixel(i);
		if (pixel)
		{
			NeighborhoodIteratorType::OffsetType offset = labelIt.GetOffset(i);
			itk::OffsetValueType dist = offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2];
			if (dist<rightDistanceVessness)
			{
				rightOstiumIndex = labelIt.GetIndex(i);
				rightDistanceVessness = dist;
			}
		}
	}
	if (leftDistanceVessness > std::numeric_limits<itk::OffsetValueType>::max() / 2) {
		cerr << "Cannot find the starting point of the left coronary artery" << endl;
		return false;
	}
	if (rightDistanceVessness > std::numeric_limits<itk::OffsetValueType>::max() / 2) {
		cerr << "Cannot find the starting point of the right coronary artery" << endl;
		return false;
	}
	std::cout << "leftOstiumIndex: " << leftOstiumIndex << std::endl;
	std::cout << "rightOstiumIndex: " << rightOstiumIndex << std::endl;

	// 检测连通域
	typedef itk::ConnectedThresholdImageFilter<InputImageType, BinaryImageType> ConnectedFilter;
	ConnectedFilter::Pointer connectedFilter = ConnectedFilter::New();
	connectedFilter->SetInput(itkLabelImage);
	connectedFilter->SetConnectivity(ConnectedFilter::FullConnectivity);
	connectedFilter->AddSeed(leftOstiumIndex);
	connectedFilter->AddSeed(rightOstiumIndex);
	connectedFilter->SetUpper(1);
	connectedFilter->SetLower(1);
	connectedFilter->SetReplaceValue(1);
	try
	{
		connectedFilter->Update();
	}
	catch (...)
	{
		cerr << "Error occurs when detecting connected domains " << endl;
	}

	// 细化
	cout << endl << "Thinning label image ... ";
	timeBegin = clock();

	typedef itk::BinaryThinningImageFilter3D<BinaryImageType, BinaryImageType> ThinningFilter;
	ThinningFilter::Pointer thinningFilter = ThinningFilter::New();
	thinningFilter->SetInput(connectedFilter->GetOutput());
	thinningFilter->Update();
	BinaryImageType::Pointer itkThinImage = thinningFilter->GetOutput();

	timeEnd = clock();
	cout << "Done Time: " << (timeEnd - timeBegin) / 1000 << "." << (timeEnd - timeBegin) % 1000 << "s" << endl;

	// 细化后的图像，搜索邻域大小为5
	typedef itk::NeighborhoodIterator<BinaryImageType>	BNeighborhoodIteratorType;
	BNeighborhoodIteratorType thinIt(nRadius, itkThinImage, itkThinImage->GetRequestedRegion());

	// 左侧
	thinIt.SetLocation(leftOstiumIndex);
	leftDistanceVessness = std::numeric_limits<itk::OffsetValueType>::max();
	for (size_t i = 0; i < thinIt.Size(); i++)
	{
		short pixel = thinIt.GetPixel(i);
		if (pixel)
		{
			NeighborhoodIteratorType::OffsetType offset = thinIt.GetOffset(i);
			itk::OffsetValueType dist = offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2];
			if (dist<leftDistanceVessness)
			{
				leftOstiumIndex = thinIt.GetIndex(i);
				leftDistanceVessness = dist;
			}
		}
	}
	// 右侧
	thinIt.SetLocation(rightOstiumIndex);
	rightDistanceVessness = std::numeric_limits<itk::OffsetValueType>::max();
	for (size_t i = 0; i < thinIt.Size(); i++)
	{
		short pixel = thinIt.GetPixel(i);
		if (pixel)
		{
			NeighborhoodIteratorType::OffsetType offset = thinIt.GetOffset(i);
			itk::OffsetValueType dist = offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2];
			if (dist<rightDistanceVessness)
			{
				rightOstiumIndex = thinIt.GetIndex(i);
				rightDistanceVessness = dist;
			}
		}
	}
	if (leftDistanceVessness > std::numeric_limits<itk::OffsetValueType>::max() / 2) {
		cerr << "Cannot find the starting point of the left coronary artery after thinning" << endl;
		return false;
	}
	if (rightDistanceVessness > std::numeric_limits<itk::OffsetValueType>::max() / 2) {
		cerr << "Cannot find the starting point of the right coronary artery after thinning" << endl;
		return false;
	}
	itkLabelImage->TransformIndexToPhysicalPoint(leftOstiumIndex, leftOstiumPoint);
	itkLabelImage->TransformIndexToPhysicalPoint(rightOstiumIndex, rightOstiumPoint);
	for (int k = 0; k < 3; k++) leftOstium[k] = leftOstiumPoint[k];
	for (int k = 0; k < 3; k++) rightOstium[k] = rightOstiumPoint[k];
	cout << "leftOstiumIndex: " << leftOstiumIndex << endl;
	cout << "rightOstiumIndex: " << rightOstiumIndex << endl;

	// 构建中心线
	cout << endl << "Constructing centerline ... ";
	timeBegin = clock();

	typedef itk::ImageDuplicator< BinaryImageType > DuplicatorType;
	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(itkThinImage);
	duplicator->Update();
	BinaryImageType::Pointer rawCenterline = duplicator->GetModifiableOutput();

	BinaryImageType::IndexType ostiumIndex[2] = { leftOstiumIndex, rightOstiumIndex };
	vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();
	for (int i = 0; i < 2; i++)
	{
		vtkSmartPointer<vtkPolyData> clModel = vtkSmartPointer<vtkPolyData>::New();
		clModel->Allocate();
		TraceCenterline<BinaryImageType>(rawCenterline, ostiumIndex[i], clModel);
		CleanCenterline(clModel);
		SimplifyCenterline(clModel);
		RadiusCenterline(clModel);
		LumenWallCenterline(clModel);
		AxisCenterline(clModel);
		append->AddInputData(clModel);
	}
	try
	{
		append->Update();
	}
	catch (...)
	{
		cerr << "Error occurs when Constructing centerline " << endl;
	}
	centerlineModel->DeepCopy(append->GetOutput());

	timeEnd = clock();
	cout << "Done Time: " << (timeEnd - timeBegin) / 1000 << "." << (timeEnd - timeBegin) % 1000 << "s" << endl;

	///////////////////////////////////////////////////////////////////////
	typedef itk::ImageToVTKImageFilter<BinaryImageType> LabelToVtkFilter;
	LabelToVtkFilter::Pointer labelToVtkFilter = LabelToVtkFilter::New();
	labelToVtkFilter->SetInput(thinningFilter->GetOutput());
	try
	{
		labelToVtkFilter->Update();
	}
	catch (...)
	{
		cerr << "Error occurs when transforming from itk to vtk " << endl;
		return false;
	}
	thinImage->DeepCopy(labelToVtkFilter->GetOutput());
	//////////////////////////////////////////////////////////////////////

	return true;
}