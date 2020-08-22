#ifndef __ExtendTubeFilter_h
#define __ExtendTubeFilter_h

#include "vtkFiltersCoreModule.h"
#include "vtkPolyDataAlgorithm.h"

#define M_PI	3.14159265358979323846

class vtkCellArray;
class vtkCellData;
class vtkDataArray;
class vtkFloatArray;
class vtkPointData;
class vtkPoints;
class vtkCardinalSpline;

class /*VTKFILTERSCORE_EXPORT*/ ExtendTubeFilter : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(ExtendTubeFilter,vtkPolyDataAlgorithm);

  static ExtendTubeFilter *New();

  vtkSetClampMacro(LongitudinalRefineSteps, int, 0, 10);
  vtkGetMacro(LongitudinalRefineSteps, int);

  vtkSetClampMacro(CircumferentialRefineSteps, int, 0, 10);
  vtkGetMacro(CircumferentialRefineSteps, int);

  vtkSetMacro(RadiusScale, double);
  vtkGetMacro(RadiusScale, double);

  virtual void SetUpdateSegment(vtkIdType update);

  // Description:
  // Turn on/off whether to cap the ends with polygons. Initial value is off.
  vtkSetMacro(Capping,int);
  vtkGetMacro(Capping,int);
  vtkBooleanMacro(Capping,int);

  vtkSetMacro(UpdateOutput0,int);
  vtkGetMacro(UpdateOutput0,int);
  vtkBooleanMacro(UpdateOutput0,int);

protected:
  ExtendTubeFilter();
  ~ExtendTubeFilter();

  virtual int ProcessRequest(vtkInformation*,
                             vtkInformationVector**,
                             vtkInformationVector*);

  // Usual data generation method
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  int Capping; //control whether tubes are capped
  int LongitudinalRefineSteps;
  int CircumferentialRefineSteps;
  double RadiusScale;

  vtkIdType UpdateSegment;
  int UpdateOutput0;

  // Helper methods
  void InterpolateRefine(vtkCardinalSpline *spline, double *in, int insize, double *out, int refinesteps);
  void FindJunctionSeam(vtkPolyData *input, vtkPolyData *output0, std::vector< std::vector<vtkIdType> >& endIndex, std::vector<std::vector<double>>& endPoints1, std::vector<std::vector<double>>& endPoints2, std::vector<std::vector<double>>& endPoints3);
  void FindBifurcationAxes(double s1[3], double s2[3], double s3[3], double on[3], double o12[3], double o23[3], double o31[3]);
  void FindSeamIndex(double on[3], double axis1[3], double axis2[3], int components, int end, std::vector<int>& leftindex, std::vector<int>& rightindex);
  void ExtendTubeFilter::FindSeamPoint(std::vector<std::vector<vtkIdType>>& seamindex,
									   std::vector<std::vector<double>>& endPoints1,  
									   std::vector<std::vector<double>>& endPoints2,  
									   std::vector<std::vector<double>>& endPoints3,
									   std::vector<std::vector<int>>& leftindexs,
									   std::vector<std::vector<int>>& rightindexs,
									   int components,
									   double center[3], double on[3], double o12[3], double o23[3], double o31[3],
									   double radius[4], double *lumenradius[4], double *wallthickness[4]);
  int PointInCellEnd(vtkPolyData *input, vtkIdType pId, vtkIdType cId);
  void GetPointData(vtkPolyData *output0, vtkIdType cId, int end, double &radius, double *lumenradius, double *wallthickness, double dir[3], double axis1[3], double axis2[3] );

private:
  ExtendTubeFilter(const ExtendTubeFilter&);  // Not implemented.
  void operator=(const ExtendTubeFilter&);  // Not implemented.
  
  void GetRotationMatrix(double axis[3], double angle, double rot[3][3]);

  vtkPolyData *out0cache;
};

#endif
