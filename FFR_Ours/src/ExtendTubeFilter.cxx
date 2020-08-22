#include "ExtendTubeFilter.h"

#include <map>

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyLine.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkTupleInterpolator.h"
#include "vtkCardinalSpline.h"

vtkStandardNewMacro(ExtendTubeFilter);

void ExtendTubeFilter::GetRotationMatrix(double axis[3], double angle, double rot[3][3])
{
	double c = cos(angle);
	double s = sin(angle);
	rot[0][0] = c + axis[0]*axis[0]*(1.0-c);
	rot[0][1] = axis[0]*axis[1]*(1.0-c) - s*axis[2];
	rot[0][2] = axis[0]*axis[2]*(1.0-c) + s*axis[1];
	rot[1][0] = axis[0]*axis[1]*(1.0-c) + s*axis[2];
	rot[1][1] = c + axis[1]*axis[1]*(1.0-c);
	rot[1][2] = axis[1]*axis[2]*(1.0-c) - s*axis[0];
	rot[2][0] = axis[0]*axis[2]*(1.0-c) - s*axis[1];
	rot[2][1] = axis[1]*axis[2]*(1.0-c) + s*axis[0];
	rot[2][2] = c + axis[2]*axis[2]*(1.0-c);
}

ExtendTubeFilter::ExtendTubeFilter()
{
	this->Capping = 0;
	this->LongitudinalRefineSteps = 3;
	this->CircumferentialRefineSteps = 3;
	this->RadiusScale = 1.0;

	this->UpdateSegment = -1;
	this->UpdateOutput0 = 1;

	this->SetNumberOfInputPorts( 1 );
	this->SetNumberOfOutputPorts( 4 );

	this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
		"Radius");

	this->SetInputArrayToProcess(1,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
		"LumenRadius");

	this->SetInputArrayToProcess(2,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
		"WallThickness");

	this->SetInputArrayToProcess(3,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
		"LongitudinalAngle");

	this->SetInputArrayToProcess(4,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
		"Dir");

	this->SetInputArrayToProcess(5,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
		"Axis1");

	this->SetInputArrayToProcess(6,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
		"Axis2");

	out0cache = vtkPolyData::New();

}

ExtendTubeFilter::~ExtendTubeFilter()
{
	out0cache->Delete();
}

void ExtendTubeFilter::SetUpdateSegment(vtkIdType update)
{
	if(this->UpdateSegment != update) this->UpdateSegment = update;
}

int ExtendTubeFilter::ProcessRequest(vtkInformation* request,
									 vtkInformationVector** inputVector,
									 vtkInformationVector* outputVector)
{
	if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA_NOT_GENERATED()))
	{
		if( this->UpdateOutput0 == 0 )
		{
			vtkInformation* outInfo = outputVector->GetInformationObject(0);
			outInfo->Set(vtkDemandDrivenPipeline::DATA_NOT_GENERATED(), 1);
		}
		if(this->UpdateSegment>=0)
		{
			for(int i=1; i<4; i++)
			{
				vtkInformation* outInfo = outputVector->GetInformationObject(i);
				outInfo->Set(vtkDemandDrivenPipeline::DATA_NOT_GENERATED(), 1);
			}
		}
	}

	// generate the data
	if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
	{
		return this->RequestData(request, inputVector, outputVector);
	}

	if(request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
	{
		return this->RequestUpdateExtent(request, inputVector, outputVector);
	}

	// execute information
	if(request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
	{
		return this->RequestInformation(request, inputVector, outputVector);
	}

	return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

int ExtendTubeFilter::RequestData(
	vtkInformation *vtkNotUsed(request),
	vtkInformationVector **inputVector,
	vtkInformationVector *outputVector)
{
	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *out0Info = outputVector->GetInformationObject(0);
	vtkInformation *out1Info = outputVector->GetInformationObject(1);
	vtkInformation *out2Info = outputVector->GetInformationObject(2);
	vtkInformation *out3Info = outputVector->GetInformationObject(3);

	// get the input and output
	vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *output0 = vtkPolyData::SafeDownCast(out0Info->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *output1 = vtkPolyData::SafeDownCast(out1Info->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *output2 = vtkPolyData::SafeDownCast(out2Info->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *output3 = vtkPolyData::SafeDownCast(out3Info->Get(vtkDataObject::DATA_OBJECT()));

	vtkDoubleArray *clRadius=vtkDoubleArray::SafeDownCast(this->GetInputArrayToProcess(0,inputVector));
	vtkDoubleArray *clLumenRadius=vtkDoubleArray::SafeDownCast(this->GetInputArrayToProcess(1,inputVector));
	vtkDoubleArray *clWallThickness=vtkDoubleArray::SafeDownCast(this->GetInputArrayToProcess(2,inputVector));
	vtkDoubleArray *clLongitudinalAngle=vtkDoubleArray::SafeDownCast(this->GetInputArrayToProcess(3,inputVector));
	vtkDoubleArray *clDir=vtkDoubleArray::SafeDownCast(this->GetInputArrayToProcess(4,inputVector));
	vtkDoubleArray *clAxis1=vtkDoubleArray::SafeDownCast(this->GetInputArrayToProcess(5,inputVector));
	vtkDoubleArray *clAxis2=vtkDoubleArray::SafeDownCast(this->GetInputArrayToProcess(6,inputVector));
	if(!clRadius)
	{
		vtkDebugMacro("Could not find clRadius.");
		return 1;
	}
	if(!clLumenRadius)
	{
		vtkDebugMacro("Could not find clLumenRadius.");
		return 1;
	}
	if(!clWallThickness)
	{
		vtkDebugMacro("Could not find clWallThickness.");
		return 1;
	}
	if( clLumenRadius->GetNumberOfComponents() != clWallThickness->GetNumberOfComponents() || clLumenRadius->GetNumberOfComponents()<3 )
	{
		vtkDebugMacro("The number of components in clLumenRadius or clWallThickness is incorrect. ");
		return 1;
	}
	if(!clLongitudinalAngle)
	{
		vtkDebugMacro("Could not find clLongitudinalAngle.");
		return 1;
	}
	if(!clDir)
	{
		vtkDebugMacro("Could not find clDir.");
		return 1;
	}
	if(!clAxis1)
	{
		vtkDebugMacro("Could not find clAxis1.");
		return 1;
	}
	if(!clAxis2)
	{
		vtkDebugMacro("Could not find clAxis2.");
		return 1;
	}

	vtkIdType inCellId;
	vtkIdType npts=0, *pts=NULL;

	vtkPoints *out0Points = vtkPoints::New();
	vtkCellArray *out0Lines = vtkCellArray::New();
	vtkDoubleArray *out0Radius = vtkDoubleArray::New();
	out0Radius->SetName("Radius");
	out0Radius->SetNumberOfComponents(1);
	vtkDoubleArray *out0LumenRadius = vtkDoubleArray::New();
	out0LumenRadius->SetName("LumenRadius");
	out0LumenRadius->SetNumberOfComponents(clLumenRadius->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1));
	vtkDoubleArray *out0WallThickness = vtkDoubleArray::New();
	out0WallThickness->SetName("WallThickness");
	out0WallThickness->SetNumberOfComponents(clWallThickness->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1));
	vtkDoubleArray *out0LongitudinalAngle = vtkDoubleArray::New();
	out0LongitudinalAngle->SetName("LongitudinalAngle");
	out0LongitudinalAngle->SetNumberOfComponents(clLongitudinalAngle->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1));
	vtkDoubleArray *out0Dir = vtkDoubleArray::New();
	out0Dir->SetName("Dir");
	out0Dir->SetNumberOfComponents(3);
	vtkDoubleArray *out0Axis1 = vtkDoubleArray::New();
	out0Axis1->SetName("Axis1");
	out0Axis1->SetNumberOfComponents(3);
	vtkDoubleArray *out0Axis2 = vtkDoubleArray::New();
	out0Axis2->SetName("Axis2");
	out0Axis2->SetNumberOfComponents(3);
	vtkDoubleArray *out0LongiParam = vtkDoubleArray::New();
	out0LongiParam->SetName("LongiParam");
	out0LongiParam->SetNumberOfComponents(1);
	vtkDoubleArray *out0CircumParam = vtkDoubleArray::New();
	out0CircumParam->SetName("CircumParam");
	out0CircumParam->SetNumberOfComponents(out0LumenRadius->GetNumberOfComponents());

	vtkCardinalSpline *spline = vtkCardinalSpline::New();

	output0->GetCellData()->CopyAllocate(input->GetCellData());
	vtkPoints*		inPoints = input->GetPoints();
	vtkCellArray* inLines = input->GetLines();

	double *radii = new double[clLumenRadius->GetNumberOfComponents()];
	double *thickness = new double[clWallThickness->GetNumberOfComponents()];
	double *longiangle = new double[clLongitudinalAngle->GetNumberOfComponents()];
	double *refineradii = new double[clLumenRadius->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1)];
	double *refinethickness = new double[clWallThickness->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1)];
	double *refinelongiangle = new double[clLongitudinalAngle->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1)];
	double radius, coord[3], center[3];
	double ffr;
	double sdir[3], edir[3], saxis1[3], eaxis1[3], saxis2[3], eaxis2[3];
	double dir[3], axis1[3], axis2[3];
	double longiparam;
	double *circumparam = new double[clLumenRadius->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1)];
	double rot[3][3], rotaxis[3], rotangle;
	double cstep = 1.0/out0LumenRadius->GetNumberOfComponents();
	double cirstep  = 2.0*M_PI*cstep;
	double circumstep = clLumenRadius->GetNumberOfComponents()*cstep;

	for(inCellId=0, inLines->InitTraversal(); inLines->GetNextCell(npts,pts); inCellId++)
	{
		if( this->UpdateSegment >= 0 && inCellId != this->UpdateSegment ) continue;

		vtkIdList *idlist = vtkIdList::New();

		vtkTupleInterpolator* pointInterpolator = vtkTupleInterpolator::New();
		pointInterpolator->SetNumberOfComponents(3);

		vtkTupleInterpolator* radiusInterpolator = vtkTupleInterpolator::New();
		radiusInterpolator->SetNumberOfComponents(1);
		vtkTupleInterpolator* lumenRadiusInterpolator = vtkTupleInterpolator::New();
		lumenRadiusInterpolator->SetNumberOfComponents(clLumenRadius->GetNumberOfComponents());
		vtkTupleInterpolator* wallThicknessInterpolator = vtkTupleInterpolator::New();
		wallThicknessInterpolator->SetNumberOfComponents(clWallThickness->GetNumberOfComponents());
		vtkTupleInterpolator* longitudinalAngleInterpolator = vtkTupleInterpolator::New();
		longitudinalAngleInterpolator->SetNumberOfComponents(clLongitudinalAngle->GetNumberOfComponents());
		for(vtkIdType j=0; j<npts; j++)
		{
			inPoints->GetPoint(pts[j], coord);
			pointInterpolator->AddTuple(j, coord);
			radius = clRadius->GetValue(pts[j]);
			radiusInterpolator->AddTuple(j, &radius);
			clLumenRadius->GetTuple(pts[j], radii);
			lumenRadiusInterpolator->AddTuple(j, radii);
			clWallThickness->GetTuple(pts[j], thickness);
			wallThicknessInterpolator->AddTuple(j, thickness);
			clLongitudinalAngle->GetTuple(pts[j], longiangle);
			longitudinalAngleInterpolator->AddTuple(j, longiangle);
		}

		for(vtkIdType j=0; j<npts; j++)
		{
			//out0Points, out0Radius, out0LumenRadius, out0WallThickness, out0LongitudinalAngle, out0LongiParam
			pointInterpolator->InterpolateTuple(j, coord);
			idlist->InsertNextId(out0Points->InsertNextPoint(coord));
			radiusInterpolator->InterpolateTuple(j, &radius);
			out0Radius->InsertNextValue(radius);
			lumenRadiusInterpolator->InterpolateTuple(j, radii);
			InterpolateRefine(spline, radii, clLumenRadius->GetNumberOfComponents(), refineradii, this->CircumferentialRefineSteps);
			out0LumenRadius->InsertNextTuple(refineradii);
			wallThicknessInterpolator->InterpolateTuple(j, thickness);
			InterpolateRefine(spline, thickness, clWallThickness->GetNumberOfComponents(), refinethickness, this->CircumferentialRefineSteps);
			out0WallThickness->InsertNextTuple(refinethickness);
			longitudinalAngleInterpolator->InterpolateTuple(j, longiangle);
			InterpolateRefine(spline, longiangle, clLongitudinalAngle->GetNumberOfComponents(), refinelongiangle, this->CircumferentialRefineSteps);
			out0LongitudinalAngle->InsertNextTuple(refinelongiangle);
			out0LongiParam->InsertNextValue(j);

			//out0Dir, out0Axis1, out0Axis2
			if(j==0)
			{
				if(npts==2)
				{
					inPoints->GetPoint(pts[j], coord);
					inPoints->GetPoint(pts[j+1], sdir);
					vtkMath::Subtract(sdir, coord, sdir);
					vtkMath::Normalize(sdir);
					vtkMath::Perpendiculars(sdir, saxis1, saxis2, 0.0);
				}
				else
				{
					clDir->GetTuple(pts[j+1], sdir);
					clAxis1->GetTuple(pts[j+1], saxis1);
					clAxis2->GetTuple(pts[j+1], saxis2);
				}
				std::copy(sdir, sdir+3, edir);
				std::copy(saxis1, saxis1+3, eaxis1);
				std::copy(saxis2, saxis2+3, eaxis2);
			}
			else
			{
				std::copy(edir, edir+3, sdir);
				std::copy(eaxis1, eaxis1+3, saxis1);
				std::copy(eaxis2, eaxis2+3, saxis2);
				if(j<npts-2)
				{
					clDir->GetTuple(pts[j+1], edir);
					clAxis1->GetTuple(pts[j+1], eaxis1);
					clAxis2->GetTuple(pts[j+1], eaxis2);
				}
				else
				{
					std::copy(sdir, sdir+3, edir);
					std::copy(saxis1, saxis1+3, eaxis1);
					std::copy(saxis2, saxis2+3, eaxis2);
				}
			}
			out0Dir->InsertNextTuple(sdir);
			out0Axis1->InsertNextTuple(saxis1);
			out0Axis2->InsertNextTuple(saxis2);

			if(j==npts-1) break;

			if(this->LongitudinalRefineSteps>0)
			{
				if( edir[0] == sdir[0] && edir[1] == sdir[1] && edir[2] == sdir[2] )
				{
					rotaxis[0] = 0.0; rotaxis[1] = 0.0; rotaxis[2] = 0.0;
					rotangle = 0.0;
				}
				else
				{
					vtkMath::Cross(sdir, edir, rotaxis);
					vtkMath::Normalize(rotaxis);
					rotangle = acos(vtkMath::Dot(sdir, edir));
				}

				double lrs = 1.0/(this->LongitudinalRefineSteps+1);
				for(int k=1; k<=this->LongitudinalRefineSteps; k++)
				{
					double s = j+k*lrs;
					pointInterpolator->InterpolateTuple(s, coord);
					idlist->InsertNextId(out0Points->InsertNextPoint(coord));
					radiusInterpolator->InterpolateTuple(s, &radius);
					out0Radius->InsertNextValue(radius);
					lumenRadiusInterpolator->InterpolateTuple(s, radii);
					InterpolateRefine(spline, radii, clLumenRadius->GetNumberOfComponents(), refineradii, this->CircumferentialRefineSteps);
					out0LumenRadius->InsertNextTuple(refineradii);
					wallThicknessInterpolator->InterpolateTuple(s, thickness);
					InterpolateRefine(spline, thickness, clWallThickness->GetNumberOfComponents(), refinethickness, this->CircumferentialRefineSteps);
					out0WallThickness->InsertNextTuple(refinethickness);
					longitudinalAngleInterpolator->InterpolateTuple(s, longiangle);
					InterpolateRefine(spline, longiangle, clLongitudinalAngle->GetNumberOfComponents(), refinelongiangle, this->CircumferentialRefineSteps);
					out0LongitudinalAngle->InsertNextTuple(refinelongiangle);
					out0LongiParam->InsertNextValue(s);

					double angle = k*lrs*rotangle;
					GetRotationMatrix(rotaxis, angle, rot);
					vtkMath::Multiply3x3(rot, saxis1, axis1);
					vtkMath::Multiply3x3(rot, saxis2, axis2);
					vtkMath::Multiply3x3(rot, sdir, dir);
					vtkMath::Normalize(axis1);
					vtkMath::Normalize(axis2);
					vtkMath::Normalize(dir);
					out0Dir->InsertNextTuple(dir);
					out0Axis1->InsertNextTuple(axis1);
					out0Axis2->InsertNextTuple(axis2);
				}
			}
		}
		pointInterpolator->Delete();
		radiusInterpolator->Delete();
		lumenRadiusInterpolator->Delete();
		wallThicknessInterpolator->Delete();
		longitudinalAngleInterpolator->Delete();

		for(int k=0; k<out0CircumParam->GetNumberOfComponents(); k++)
		{
			circumparam[k] = k*circumstep;
		}
		out0CircumParam->InsertNextTuple(circumparam);

		if( this->UpdateSegment < 0 )
		{
			vtkIdType outcellId = out0Lines->InsertNextCell(idlist);
			output0->GetCellData()->CopyData(input->GetCellData(),inCellId,outcellId);
		}

		idlist->Delete();
	}
	delete[] radii;
	delete[] thickness;
	spline->Delete();

	if( input->GetMTime() > output0->GetMTime() ) this->UpdateOutput0 = 1;

	if( this->UpdateOutput0 )
	{
		if( this->UpdateSegment >= 0 )
		{
			vtkPoints* cachePoints = out0cache->GetPoints();
			vtkDoubleArray *cacheRadius = vtkDoubleArray::SafeDownCast(out0cache->GetPointData()->GetArray("Radius"));
			vtkDoubleArray *cacheLumenRadius = vtkDoubleArray::SafeDownCast(out0cache->GetPointData()->GetArray("LumenRadius"));
			vtkDoubleArray *cacheWallThickness = vtkDoubleArray::SafeDownCast(out0cache->GetPointData()->GetArray("WallThickness"));
			vtkDoubleArray *cacheLongitudinalAngle = vtkDoubleArray::SafeDownCast(out0cache->GetPointData()->GetArray("LongitudinalAngle"));
			vtkDoubleArray *cacheDir = vtkDoubleArray::SafeDownCast(out0cache->GetPointData()->GetArray("Dir"));
			vtkDoubleArray *cacheAxis1 = vtkDoubleArray::SafeDownCast(out0cache->GetPointData()->GetArray("Axis1"));
			vtkDoubleArray *cacheAxis2 = vtkDoubleArray::SafeDownCast(out0cache->GetPointData()->GetArray("Axis2"));
			vtkDoubleArray *cacheLongiParam = vtkDoubleArray::SafeDownCast(out0cache->GetPointData()->GetArray("LongiParam"));
			vtkDoubleArray *cacheCircumParam = vtkDoubleArray::SafeDownCast(out0cache->GetCellData()->GetArray("CircumParam"));
			inLines = out0cache->GetLines();
			for(inCellId=0, inLines->InitTraversal(); inLines->GetNextCell(npts,pts); inCellId++)
			{
				if(inCellId == this->UpdateSegment)
				{
					for(vtkIdType j=0; j<npts; j++)
					{
						cachePoints->SetPoint(pts[j], out0Points->GetPoint(j));
						cacheRadius->SetValue(pts[j], out0Radius->GetValue(j));
						cacheLumenRadius->SetTuple(pts[j], out0LumenRadius->GetTuple(j));
						cacheWallThickness->SetTuple(pts[j], out0WallThickness->GetTuple(j));
						cacheLongitudinalAngle->SetTuple(pts[j], out0LongitudinalAngle->GetTuple(j));
						cacheDir->SetTuple(pts[j], out0Dir->GetTuple(j));
						cacheAxis1->SetTuple(pts[j], out0Axis1->GetTuple(j));
						cacheAxis2->SetTuple(pts[j], out0Axis2->GetTuple(j));
						cacheLongiParam->SetValue(pts[j], out0LongiParam->GetValue(j));
					}
					cacheCircumParam->SetTuple(inCellId, out0CircumParam->GetTuple(0));
					break;
				}
			}
			output0->DeepCopy(out0cache);
		}
		else
		{
			output0->SetPoints(out0Points);
			output0->SetLines(out0Lines);
			output0->GetPointData()->AddArray(out0Radius);
			output0->GetPointData()->AddArray(out0LumenRadius);
			output0->GetPointData()->AddArray(out0WallThickness);
			output0->GetPointData()->AddArray(out0LongitudinalAngle);
			output0->GetPointData()->AddArray(out0Dir);
			output0->GetPointData()->AddArray(out0Axis1);
			output0->GetPointData()->AddArray(out0Axis2);
			output0->GetPointData()->AddArray(out0LongiParam);
			output0->GetCellData()->AddArray(out0CircumParam);
			out0cache->DeepCopy(output0);
		}
	}

	if( this->UpdateSegment < 0 )
	{
		int hasFFR = output0->GetPointData()->HasArray("FFR");
		vtkDoubleArray *out0FFR = hasFFR ? vtkDoubleArray::SafeDownCast(output0->GetPointData()->GetArray("FFR")) : NULL;
		vtkPoints *out1Points = vtkPoints::New();
		vtkCellArray *out1Strips = vtkCellArray::New();
		vtkDoubleArray *out1Param = vtkDoubleArray::New();
		out1Param->SetName("Param");
		out1Param->SetNumberOfComponents(2);
		vtkDoubleArray *out1Radius = vtkDoubleArray::New();
		out1Radius->SetName("Radius");
		out1Radius->SetNumberOfComponents(1);
		vtkDoubleArray *out1FFR = NULL;
		if(hasFFR)
		{
			out1FFR = vtkDoubleArray::New();
			out1FFR->SetName("FFR");
			out1FFR->SetNumberOfComponents(1);
		}

		vtkPoints *out2Points = vtkPoints::New();
		vtkCellArray *out2Strips = vtkCellArray::New();
		vtkDoubleArray *out2Param = vtkDoubleArray::New();
		out2Param->SetName("Param");
		out2Param->SetNumberOfComponents(2);
		vtkDoubleArray *out2Radius = vtkDoubleArray::New();
		out2Radius->SetName("Radius");
		out2Radius->SetNumberOfComponents(1);
		vtkDoubleArray *out2FFR = NULL;
		if(hasFFR) 
		{
			out2FFR = vtkDoubleArray::New();
			out2FFR->SetName("FFR");
			out2FFR->SetNumberOfComponents(1);
		}

		vtkPoints *out3Points = vtkPoints::New();
		vtkCellArray *out3Strips = vtkCellArray::New();
		vtkDoubleArray *out3Param = vtkDoubleArray::New();
		out3Param->SetName("Param");
		out3Param->SetNumberOfComponents(2);
		vtkDoubleArray *out3Radius = vtkDoubleArray::New();
		out3Radius->SetName("Radius");
		out3Radius->SetNumberOfComponents(1);
		vtkDoubleArray *out3FFR = NULL;
		if(hasFFR) 
		{
			out3FFR = vtkDoubleArray::New();
			out3FFR->SetName("FFR");
			out3FFR->SetNumberOfComponents(1);
		}

		output1->GetCellData()->CopyFieldOff("CircumParam");
		output2->GetCellData()->CopyFieldOff("CircumParam");
		output3->GetCellData()->CopyFieldOff("CircumParam");
		output1->GetCellData()->CopyAllocate(output0->GetCellData());
		output2->GetCellData()->CopyAllocate(output0->GetCellData());
		output3->GetCellData()->CopyAllocate(output0->GetCellData());

		std::vector< std::vector<vtkIdType> > endIndex;
		std::vector<std::vector<double>> endPoints1;
		std::vector<std::vector<double>> endPoints2;
		std::vector<std::vector<double>> endPoints3;
		FindJunctionSeam(input, output0, endIndex, endPoints1, endPoints2, endPoints3);
		for(int k=0; k<out0CircumParam->GetNumberOfComponents(); k++) circumparam[k] = k*circumstep;
		for(size_t i=0; i<endPoints1.size(); i++) 
		{
			out1Points->InsertNextPoint(endPoints1[i][0], endPoints1[i][1], endPoints1[i][2]);
			out1Param->InsertNextTuple2(0.0, 0.0);
			out1Radius->InsertNextValue(10.0);
			if(hasFFR) out1FFR->InsertNextValue(1.0);
		}
		for(size_t i=0; i<endPoints2.size(); i++)
		{
			out2Points->InsertNextPoint(endPoints2[i][0], endPoints2[i][1], endPoints2[i][2]);
			out2Param->InsertNextTuple2(0.0, 0.0);
			out2Radius->InsertNextValue(10.0);
			if(hasFFR) out2FFR->InsertNextValue(1.0);
		}
		for(size_t i=0; i<endPoints3.size(); i++)
		{
			out3Points->InsertNextPoint(endPoints3[i][0], endPoints3[i][1], endPoints3[i][2]);
			out3Param->InsertNextTuple2(0.0, 0.0);
			out3Radius->InsertNextValue(10.0);
			if(hasFFR) out3FFR->InsertNextValue(1.0);
		}

		for(inCellId=0, out0Lines->InitTraversal(); out0Lines->GetNextCell(npts,pts); inCellId++)
		{
			for(int j=0; j<2; j++)
			{
				int ind = 2*inCellId+j;
				if( endIndex[ind].size()==0 ) continue;
				radius = out0Radius->GetValue(pts[j==0?0:npts-1]);
				if(hasFFR) ffr = out0FFR->GetValue(pts[j==0?0:npts-1]);
				for(int k=0; k<out0LumenRadius->GetNumberOfComponents(); k++)
				{
					out1Radius->SetValue(endIndex[ind][k], radius);
					out2Radius->SetValue(endIndex[ind][k], radius);
					out3Radius->SetValue(endIndex[ind][k], radius);
					if(hasFFR)
					{
						out1FFR->SetValue(endIndex[ind][k], ffr);
						out2FFR->SetValue(endIndex[ind][k], ffr);
						out3FFR->SetValue(endIndex[ind][k], ffr);
					}
				}
			}
		}

		for(inCellId=0, out0Lines->InitTraversal(); out0Lines->GetNextCell(npts,pts); inCellId++)
		{
			std::vector<std::vector<std::vector<double>>> displacement1(2); 
			std::vector<std::vector<std::vector<double>>> displacement2(2); 
			std::vector<std::vector<std::vector<double>>> displacement3(2);
			for(vtkIdType j=1; j<npts-1; j=j+npts-3)
			{
				out0Points->GetPoint(pts[j], center);
				out0Axis1->GetTuple(pts[j], axis1);
				out0Axis2->GetTuple(pts[j], axis2);
				radius = out0Radius->GetValue(pts[j]);
				out0LumenRadius->GetTuple(pts[j], refineradii);
				out0WallThickness->GetTuple(pts[j], refinethickness);

				int end = j==1?0:1;
				int ind = 2*inCellId+end;
				if( endIndex[ind].size()>0 )
				{
					displacement1[end].resize(out0LumenRadius->GetNumberOfComponents());
					displacement2[end].resize(out0LumenRadius->GetNumberOfComponents());
					displacement3[end].resize(out0LumenRadius->GetNumberOfComponents());
					for(int k=0; k<out0LumenRadius->GetNumberOfComponents(); k++)
					{
						displacement1[end][k].resize(3);
						for(int l=0; l<3; l++) displacement1[end][k][l] = center[l] + radius * this->RadiusScale * ( cos(k*cirstep)*axis1[l] + sin(k*cirstep)*axis2[l] ) - endPoints1[endIndex[ind][k]][l];

						displacement2[end][k].resize(3);
						for(int l=0; l<3; l++) displacement2[end][k][l] = center[l] + refineradii[k] * this->RadiusScale * ( cos(k*cirstep)*axis1[l] + sin(k*cirstep)*axis2[l] ) - endPoints2[endIndex[ind][k]][l];

						displacement3[end][k].resize(3);
						for(int l=0; l<3; l++) displacement3[end][k][l]  = center[l] + (refineradii[k]+refinethickness[k]) * this->RadiusScale * ( cos(k*cirstep)*axis1[l] + sin(k*cirstep)*axis2[l] ) - endPoints3[endIndex[ind][k]][l];
					}
				}
			}

			vtkIdList *idlist1prev = vtkIdList::New();
			vtkIdList *idlist1curr = vtkIdList::New();
			vtkIdList *idlist2prev = vtkIdList::New();
			vtkIdList *idlist2curr = vtkIdList::New();
			vtkIdList *idlist3prev = vtkIdList::New();
			vtkIdList *idlist3curr = vtkIdList::New();
			out0CircumParam->GetTuple(inCellId, circumparam);
			for(vtkIdType j=0; j<npts; j++)
			{
				if( j==1 || j==npts-2 ) continue;

				out0Points->GetPoint(pts[j], center);
				out0Dir->GetTuple(pts[j], dir);
				out0Axis1->GetTuple(pts[j], axis1);
				out0Axis2->GetTuple(pts[j], axis2);
				radius = out0Radius->GetValue(pts[j]);
				if(hasFFR) ffr = out0FFR->GetValue(pts[j]);

				out0LumenRadius->GetTuple(pts[j], refineradii);
				out0WallThickness->GetTuple(pts[j], refinethickness);
				longiparam = out0LongiParam->GetValue(pts[j]);

				if( (j==0 || j==npts-1) && endIndex[2*inCellId+(j==0?0:1)].size()>0 )
				{
					int ind = 2*inCellId+(j==0?0:1);
					for(int k=0; k<out0LumenRadius->GetNumberOfComponents(); k++)
					{
						idlist1curr->InsertNextId(endIndex[ind][k]);
						idlist2curr->InsertNextId(endIndex[ind][k]);
						idlist3curr->InsertNextId(endIndex[ind][k]);
					}
				}
				else
				{
					for(int k=0; k<out0LumenRadius->GetNumberOfComponents(); k++)
					{
						for(int l=0; l<3; l++) coord[l] = center[l] + radius * this->RadiusScale * ( cos(k*cirstep)*axis1[l] + sin(k*cirstep)*axis2[l] );

						if( j<4 && j<npts/2 && endIndex[2*inCellId].size()>0 )
						{
							for(int l=0; l<3; l++) coord[l] += 1.0/j * displacement1[0][k][l];
						}
						if( j>npts-5 && j>=npts/2 && endIndex[2*inCellId+1].size()>0 )
						{
							for(int l=0; l<3; l++) coord[l] += 1.0/(npts-1-j) * displacement1[1][k][l];
						}
						idlist1curr->InsertNextId(out1Points->InsertNextPoint(coord));
						out1Param->InsertNextTuple2(longiparam, circumparam[k]);
						out1Radius->InsertNextValue(radius);
						if(hasFFR) out1FFR->InsertNextValue(ffr);

						for(int l=0; l<3; l++) coord[l] = center[l] + refineradii[k] * this->RadiusScale * ( cos(k*cirstep)*axis1[l] + sin(k*cirstep)*axis2[l] );
						if( j<4 && j<npts/2 && endIndex[2*inCellId].size()>0 )
						{
							for(int l=0; l<3; l++) coord[l] += 1.0/j * displacement2[0][k][l];
						}
						if( j>npts-5 && j>=npts/2 && endIndex[2*inCellId+1].size()>0 )
						{
							for(int l=0; l<3; l++) coord[l] += 1.0/(npts-1-j) * displacement2[1][k][l];
						}
						idlist2curr->InsertNextId(out2Points->InsertNextPoint(coord));
						out2Param->InsertNextTuple2(longiparam, circumparam[k]);
						out2Radius->InsertNextValue(refineradii[k]);
						if(hasFFR) out2FFR->InsertNextValue(ffr);

						for(int l=0; l<3; l++) coord[l] = center[l] + (refineradii[k]+refinethickness[k]) * this->RadiusScale * ( cos(k*cirstep)*axis1[l] + sin(k*cirstep)*axis2[l] );
						if( j<4 && j<npts/2 && endIndex[2*inCellId].size()>0 )
						{
							for(int l=0; l<3; l++) coord[l] += 1.0/j * displacement3[0][k][l];
						}
						if( j>npts-5 && j>=npts/2 && endIndex[2*inCellId+1].size()>0 )
						{
							for(int l=0; l<3; l++) coord[l] += 1.0/(npts-1-j) * displacement3[1][k][l];
						}
						idlist3curr->InsertNextId(out3Points->InsertNextPoint(coord));
						out3Param->InsertNextTuple2(longiparam, circumparam[k]);
						out3Radius->InsertNextValue(refineradii[k]+refinethickness[k]);
						if(hasFFR) out3FFR->InsertNextValue(ffr);
					}
				}

				if(j!=0)
				{
					vtkIdType outcellId;
					outcellId = out1Strips->InsertNextCell(2*(idlist1curr->GetNumberOfIds()+1));
					output1->GetCellData()->CopyData(output0->GetCellData(),inCellId,outcellId);
					for(int k=0; k<idlist1curr->GetNumberOfIds(); k++)
					{
						out1Strips->InsertCellPoint(idlist1curr->GetId(k));
						out1Strips->InsertCellPoint(idlist1prev->GetId(k));
					}
					out1Strips->InsertCellPoint(idlist1curr->GetId(0));
					out1Strips->InsertCellPoint(idlist1prev->GetId(0));

					outcellId = out2Strips->InsertNextCell(2*(idlist2curr->GetNumberOfIds()+1));
					output2->GetCellData()->CopyData(output0->GetCellData(),inCellId,outcellId);
					for(int k=0; k<idlist2curr->GetNumberOfIds(); k++)
					{
						out2Strips->InsertCellPoint(idlist2curr->GetId(k));
						out2Strips->InsertCellPoint(idlist2prev->GetId(k));
					}
					out2Strips->InsertCellPoint(idlist2curr->GetId(0));
					out2Strips->InsertCellPoint(idlist2prev->GetId(0));

					outcellId = out3Strips->InsertNextCell(2*(idlist3curr->GetNumberOfIds()+1));
					output3->GetCellData()->CopyData(output0->GetCellData(),inCellId,outcellId);
					for(int k=0; k<idlist3curr->GetNumberOfIds(); k++)
					{
						out3Strips->InsertCellPoint(idlist3curr->GetId(k));
						out3Strips->InsertCellPoint(idlist3prev->GetId(k));
					}
					out3Strips->InsertCellPoint(idlist3curr->GetId(0));
					out3Strips->InsertCellPoint(idlist3prev->GetId(0));
				}
				idlist1prev->DeepCopy(idlist1curr);
				idlist2prev->DeepCopy(idlist2curr);
				idlist3prev->DeepCopy(idlist3curr);
				idlist1curr->Reset();
				idlist2curr->Reset();
				idlist3curr->Reset();
			}
			idlist1prev->Delete();
			idlist1curr->Delete();
			idlist2prev->Delete();
			idlist2curr->Delete();
			idlist3prev->Delete();
			idlist3curr->Delete();
		}
		output1->SetPoints(out1Points); out1Points->Delete();
		output2->SetPoints(out2Points); out2Points->Delete();
		output3->SetPoints(out3Points); out3Points->Delete();
		output1->GetPointData()->AddArray(out1Param); out1Param->Delete();
		output2->GetPointData()->AddArray(out2Param); out2Param->Delete();
		output3->GetPointData()->AddArray(out3Param); out3Param->Delete();
		output1->GetPointData()->SetScalars(out1Radius); out1Radius->Delete();
		output2->GetPointData()->SetScalars(out2Radius); out2Radius->Delete();
		output3->GetPointData()->SetScalars(out3Radius); out3Radius->Delete();
		if(hasFFR)
		{
			output1->GetPointData()->SetScalars(out1FFR); out1FFR->Delete();
			output2->GetPointData()->SetScalars(out2FFR); out2FFR->Delete();
			output3->GetPointData()->SetScalars(out3FFR); out3FFR->Delete();
		}
		output1->SetStrips(out1Strips); out1Strips->Delete();
		output2->SetStrips(out2Strips); out2Strips->Delete();
		output3->SetStrips(out3Strips); out3Strips->Delete();
	}

	delete[] refineradii;
	delete[] refinethickness;
	delete[] circumparam;
	out0Points->Delete();
	out0Lines->Delete();
	out0Radius->Delete();
	out0LumenRadius->Delete();
	out0WallThickness->Delete();
	out0Dir->Delete();
	out0Axis1->Delete();
	out0Axis2->Delete();
	out0LongiParam->Delete();
	out0CircumParam->Delete();

	return 1;
}

void ExtendTubeFilter::InterpolateRefine(vtkCardinalSpline *spline, double *in, int insize, double *out, int refinesteps)
{
	spline->RemoveAllPoints();
	for(int i=0; i<insize; i++)
	{
		spline->AddPoint(i, in[i]);
	}
	spline->ClosedOn();
	spline->Compute();
	int ct=0;
	for(int i=0; i<insize; i++)
	{
		out[ct++] = spline->Evaluate(i);
		double rs = 1.0/(refinesteps+1);
		for(int j=1; j<=refinesteps; j++)
			out[ct++] = spline->Evaluate(i+j*rs);
	}
}

void ExtendTubeFilter::FindJunctionSeam(vtkPolyData *input, vtkPolyData *output0, std::vector< std::vector<vtkIdType> >& endIndex, std::vector<std::vector<double>>& endPoints1, std::vector<std::vector<double>>& endPoints2, std::vector<std::vector<double>>& endPoints3)
{
	endIndex.clear();
	endIndex.resize(output0->GetNumberOfLines()*2);
	endPoints1.clear();
	endPoints2.clear();
	endPoints3.clear();

	output0->BuildLinks();
	vtkDoubleArray *out0LumenRadius = vtkDoubleArray::SafeDownCast(output0->GetPointData()->GetArray("LumenRadius"));
	vtkDoubleArray *out0WallThickness = vtkDoubleArray::SafeDownCast(output0->GetPointData()->GetArray("WallThickness"));

	double center[3];
	vtkIdType cIds[4];
	int ends[4]; 
	double dir[4][3], axis1[4][3], axis2[4][3];
	double edir[4][3];
	double radius[4], *lumenradius[4], *wallthickness[4];
	std::vector<std::vector<int>> leftindexs(4);
	std::vector<std::vector<int>> rightindexs(4);
	std::vector<std::vector<vtkIdType>> seamindex(4);
	for(int k=0; k<4; k++)
	{
		lumenradius[k] = new double[out0LumenRadius->GetNumberOfComponents()];
		wallthickness[k] = new double[out0WallThickness->GetNumberOfComponents()];
		leftindexs[k].resize(out0LumenRadius->GetNumberOfComponents()/2 + 1);
		rightindexs[k].resize(out0LumenRadius->GetNumberOfComponents()/2 + 1);
		seamindex[k].resize(out0LumenRadius->GetNumberOfComponents()/2 + 1);
	}

	for(vtkIdType id=0; id<input->GetNumberOfPoints(); id++)
	{
		vtkIdList *cellIds = vtkIdList::New();
		input->GetPointCells(id, cellIds);
		if( cellIds->GetNumberOfIds() > 2 )
		{
			input->GetPoint(id, center);
			if( cellIds->GetNumberOfIds() == 3 )
			{
				for(int k=0; k<3; k++) 
				{
					cIds[k] = cellIds->GetId(k);
					ends[k] = PointInCellEnd(input, id, cIds[k]);
					GetPointData(output0, cIds[k], ends[k], radius[k], lumenradius[k], wallthickness[k], dir[k], axis1[k], axis2[k]);
					if(ends[k]) for(int d=0; d<3; d++) edir[k][d] = -dir[k][d];
					else        for(int d=0; d<3; d++) edir[k][d] =  dir[k][d];
				}
				double on[3], o12[3], o23[3], o31[3];
				FindBifurcationAxes(edir[0], edir[1], edir[2], on, o12, o23, o31);
				for(int k=0; k<3; k++) 
				{
					FindSeamIndex(on, axis1[k], axis2[k], out0LumenRadius->GetNumberOfComponents(), ends[k], leftindexs[k], rightindexs[k]);
				}
				FindSeamPoint(seamindex, endPoints1, endPoints2, endPoints3, leftindexs, rightindexs, out0LumenRadius->GetNumberOfComponents(), center, on, o12, o23, o31, radius, lumenradius, wallthickness);
				for(int k=0; k<3; k++)  endIndex[2*cIds[k]+ends[k]].resize(out0LumenRadius->GetNumberOfComponents());
				for(int c=0; c<out0LumenRadius->GetNumberOfComponents()/2; c++)
				{
					endIndex[2*cIds[0]+ends[0]][leftindexs[0][c]] = seamindex[2][c];
					endIndex[2*cIds[1]+ends[1]][leftindexs[1][c]] = seamindex[0][c];
					endIndex[2*cIds[2]+ends[2]][leftindexs[2][c]] = seamindex[1][c];
				}
				for(int c=1; c<=out0LumenRadius->GetNumberOfComponents()/2; c++)
				{
					endIndex[2*cIds[0]+ends[0]][rightindexs[0][c]] = seamindex[0][c];
					endIndex[2*cIds[1]+ends[1]][rightindexs[1][c]] = seamindex[1][c];
					endIndex[2*cIds[2]+ends[2]][rightindexs[2][c]] = seamindex[2][c];
				}

			}
		}
		cellIds->Delete();
	}

	for(int k=0; k<4; k++)
	{
		delete[] lumenradius[k];
		delete[] wallthickness[k];
	}
}

inline void FindBisectAxis(double on[3], double v1[3], double v2[3], double o12[3])
{
	double t[3];
	double angle=vtkMath::AngleBetweenVectors(v1, v2)/2.0;
	double costheta = cos(0.5*angle);
	double sintheta = sin(0.5*angle);
	double quat[4];
	quat[0] = costheta;
	quat[1] = on[0]*sintheta;
	quat[2] = on[1]*sintheta;
	quat[3] = on[2]*sintheta;

	double mat[3][3];
	vtkMath::QuaternionToMatrix3x3(quat, mat);
	vtkMath::Multiply3x3(mat, v1, o12);
	vtkMath::Cross(v1, o12, t);
	if( vtkMath::Dot(on, t) < 0.0 )
	{
		vtkMath::MultiplyScalar(o12, -1.0);
	}
	vtkMath::Normalize(o12);
}

void ExtendTubeFilter::FindBifurcationAxes(double s1[3], double s2[3], double s3[3], double on[3], double o12[3], double o23[3], double o31[3])
{
	double v1[3], v2[3], v3[3];
	vtkMath::Subtract(s2, s1, v1);
	vtkMath::Subtract(s3, s1, v2);
	vtkMath::Cross(v1,v2, on);
	vtkMath::Normalize(on);
	vtkMath::ProjectVector(s1, on, v1); vtkMath::Subtract(s1, v1, v1); vtkMath::Normalize(v1);
	vtkMath::ProjectVector(s2, on, v2); vtkMath::Subtract(s2, v2, v2); vtkMath::Normalize(v2);
	vtkMath::ProjectVector(s3, on, v3); vtkMath::Subtract(s3, v3, v3); vtkMath::Normalize(v3);
	FindBisectAxis(on, v1, v2, o12);
	FindBisectAxis(on, v2, v3, o23);
	FindBisectAxis(on, v3, v1, o31);
}

void ExtendTubeFilter::FindSeamIndex(double on[3], double axis1[3], double axis2[3], int components, int end, std::vector<int>& leftindex, std::vector<int>& rightindex)
{
	double angle = atan2(vtkMath::Dot(on, axis2), vtkMath::Dot(on, axis1));
	if(angle < 0.0) angle = angle + 2.0*M_PI;
	int onk = (int( angle * components/(2.0*M_PI) + 0.5) + components ) % components;
	int halfcomponents = components/2 + 1;
	if(end)
	{
		for(int i=0; i<halfcomponents; i++) 
		{
			leftindex[i] = (onk-i+components)%components;
			rightindex[i] = (onk+i)%components;
		}
	}
	else
	{
		for(int i=0; i<halfcomponents; i++) 
		{
			leftindex[i] = (onk+i)%components;
			rightindex[i] = (onk-i+components)%components;
		}
	}
}

void ExtendTubeFilter::FindSeamPoint(std::vector<std::vector<vtkIdType>>& seamindex,
	std::vector<std::vector<double>>& endPoints1,  
	std::vector<std::vector<double>>& endPoints2,  
	std::vector<std::vector<double>>& endPoints3,
	std::vector<std::vector<int>>& leftindexs,
	std::vector<std::vector<int>>& rightindexs,
	int components,
	double center[3], double on[3], double o12[3], double o23[3], double o31[3],
	double radius[4], double *lumenradius[4], double *wallthickness[4])
{
	std::vector<double> coord(3);
	std::vector<vtkIdType> eindex(2);
	double sradius = (radius[0] + radius[1] + radius[2])/3.0;
	for(int k=0; k<3; k++) coord[k] = center[k] + on[k]*sradius;
	eindex[0] = endPoints1.size();
	endPoints1.push_back(coord);
	for(int k=0; k<3; k++) coord[k] = center[k] - on[k]*sradius;
	eindex[1] = endPoints1.size();
	endPoints1.push_back(coord);
	sradius = (lumenradius[0][leftindexs[0][0]] + lumenradius[1][leftindexs[1][0]] + lumenradius[2][leftindexs[2][0]] )/3.0;
	for(int k=0; k<3; k++) coord[k] = center[k] + on[k]*sradius;
	endPoints2.push_back(coord);
	sradius = (lumenradius[0][leftindexs[0][leftindexs.size()-1]] + lumenradius[1][leftindexs[1][leftindexs.size()-1]] + lumenradius[2][leftindexs[2][leftindexs.size()-1]] )/3.0;
	for(int k=0; k<3; k++) coord[k] = center[k] - on[k]*sradius;
	endPoints2.push_back(coord);
	sradius = (lumenradius[0][leftindexs[0][0]] + wallthickness[0][leftindexs[0][0]] + lumenradius[1][leftindexs[1][0]] + wallthickness[1][leftindexs[1][0]] + lumenradius[2][leftindexs[2][0]] + wallthickness[2][leftindexs[2][0]])/3.0;
	for(int k=0; k<3; k++) coord[k] = center[k] + on[k]*sradius;
	endPoints3.push_back(coord);
	sradius = (lumenradius[0][leftindexs[0][leftindexs.size()-1]] + wallthickness[0][leftindexs[0][leftindexs.size()-1]] + lumenradius[1][leftindexs[1][leftindexs.size()-1]] + wallthickness[1][leftindexs[1][leftindexs.size()-1]] + lumenradius[2][leftindexs[2][leftindexs.size()-1]] + wallthickness[2][leftindexs[2][leftindexs.size()-1]])/3.0;
	for(int k=0; k<3; k++) coord[k] = center[k] - on[k]*sradius;
	endPoints3.push_back(coord);
	
	double circumstep = 2.0*M_PI/components;
	for(int i=0; i<seamindex[0].size(); i++)
	{
		if( i==0 )							for(int k=0; k<3; k++) seamindex[k][i] = eindex[0];
		else if( i==seamindex[0].size()-1 ) for(int k=0; k<3; k++) seamindex[k][i] = eindex[1];
		else
		{
			sradius = (radius[0] + radius[1])/2.0;
			for(int k=0; k<3; k++) coord[k] = center[k] + sradius * ( on[k] * cos(i*circumstep) + o12[k] * sin(i*circumstep) );
			seamindex[0][i] = endPoints1.size();
			endPoints1.push_back(coord);
			sradius = (radius[1] + radius[2])/2.0;
			for(int k=0; k<3; k++) coord[k] = center[k] + sradius * ( on[k] * cos(i*circumstep) + o23[k] * sin(i*circumstep) );
			seamindex[1][i] = endPoints1.size();
			endPoints1.push_back(coord);
			sradius = (radius[2] + radius[0])/2.0;
			for(int k=0; k<3; k++) coord[k] = center[k] + sradius * ( on[k] * cos(i*circumstep) + o31[k] * sin(i*circumstep) );
			seamindex[2][i] = endPoints1.size();
			endPoints1.push_back(coord);

			sradius = (lumenradius[0][rightindexs[0][i]] + lumenradius[1][leftindexs[1][i]])/2.0;
			for(int k=0; k<3; k++) coord[k] = center[k] + sradius * ( on[k] * cos(i*circumstep) + o12[k] * sin(i*circumstep) );
			endPoints2.push_back(coord);
			sradius = (lumenradius[1][rightindexs[1][i]] + lumenradius[2][leftindexs[2][i]])/2.0;
			for(int k=0; k<3; k++) coord[k] = center[k] + sradius * ( on[k] * cos(i*circumstep) + o23[k] * sin(i*circumstep) );
			endPoints2.push_back(coord);
			sradius = (lumenradius[2][rightindexs[2][i]] + lumenradius[0][leftindexs[0][i]])/2.0;
			for(int k=0; k<3; k++) coord[k] = center[k] + sradius * ( on[k] * cos(i*circumstep) + o31[k] * sin(i*circumstep) );
			endPoints2.push_back(coord);

			sradius = (lumenradius[0][rightindexs[0][i]] + wallthickness[0][rightindexs[0][i]] + lumenradius[1][leftindexs[1][i]] + wallthickness[1][leftindexs[1][i]])/2.0;
			for(int k=0; k<3; k++) coord[k] = center[k] + sradius * ( on[k] * cos(i*circumstep) + o12[k] * sin(i*circumstep) );
			endPoints3.push_back(coord);
			sradius = (lumenradius[1][rightindexs[1][i]] + wallthickness[1][rightindexs[1][i]] + lumenradius[2][leftindexs[2][i]] + wallthickness[2][leftindexs[2][i]])/2.0;
			for(int k=0; k<3; k++) coord[k] = center[k] + sradius * ( on[k] * cos(i*circumstep) + o23[k] * sin(i*circumstep) );
			endPoints3.push_back(coord);
			sradius = (lumenradius[2][rightindexs[2][i]] + wallthickness[2][rightindexs[2][i]] + lumenradius[0][leftindexs[0][i]] + wallthickness[0][leftindexs[0][i]])/2.0;
			for(int k=0; k<3; k++) coord[k] = center[k] + sradius * ( on[k] * cos(i*circumstep) + o31[k] * sin(i*circumstep) );
			endPoints3.push_back(coord);
		}
	}
}

int ExtendTubeFilter::PointInCellEnd(vtkPolyData *input, vtkIdType pId, vtkIdType cId)
{
	vtkIdType npts, *pts;
	input->GetCellPoints(cId, npts, pts);
	if( pId == pts[0] ) 	return 0;
	else if( pId == pts[npts-1] ) return 1;
	else 
	{
		std::cerr << "PointInCellEnd error" << std::endl;
		return -1;
	}
}

void ExtendTubeFilter::GetPointData(vtkPolyData *output0, vtkIdType cId, int end, double &radius, double *lumenradius, double *wallthickness, double dir[3], double axis1[3], double axis2[3] )
{
	vtkDoubleArray *out0Radius = vtkDoubleArray::SafeDownCast(output0->GetPointData()->GetArray("Radius"));
	vtkDoubleArray *out0LumenRadius = vtkDoubleArray::SafeDownCast(output0->GetPointData()->GetArray("LumenRadius"));
	vtkDoubleArray *out0WallThickness = vtkDoubleArray::SafeDownCast(output0->GetPointData()->GetArray("WallThickness"));
	vtkDoubleArray *out0Dir = vtkDoubleArray::SafeDownCast(output0->GetPointData()->GetArray("Dir"));
	vtkDoubleArray *out0Axis1 = vtkDoubleArray::SafeDownCast(output0->GetPointData()->GetArray("Axis1"));
	vtkDoubleArray *out0Axis2 = vtkDoubleArray::SafeDownCast(output0->GetPointData()->GetArray("Axis2"));

	vtkIdType npts, *pts;
	output0->GetCellPoints(cId, npts, pts);
	if( end == 0 ) 
	{
		vtkIdType nextId = pts[1];
		radius = out0Radius->GetValue(nextId);
		out0LumenRadius->GetTuple(nextId, lumenradius);
		out0WallThickness->GetTuple(nextId, wallthickness);
		out0Dir->GetTuple(nextId, dir);
		out0Axis1->GetTuple(nextId, axis1);
		out0Axis2->GetTuple(nextId, axis2);
	}
	else if( end == 1 )
	{
		vtkIdType prevId = pts[npts-2];
		radius = out0Radius->GetValue(prevId);
		out0LumenRadius->GetTuple(prevId, lumenradius);
		out0WallThickness->GetTuple(prevId, wallthickness);
		out0Dir->GetTuple(prevId, dir);
		out0Axis1->GetTuple(prevId, axis1);
		out0Axis2->GetTuple(prevId, axis2);
	}
	else 
	{
		std::cerr << "GetPointData error" << std::endl;
	}
}

