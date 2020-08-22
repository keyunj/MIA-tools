#ifndef __ExtendIdFilter_h
#define __ExtendIdFilter_h

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkDataSetAlgorithm.h"
#include "vtkSetGet.h"

class /*VTKFILTERSCORE_EXPORT*/ ExtendIdFilter : public vtkDataSetAlgorithm
{
public:
  vtkTypeMacro(ExtendIdFilter,vtkDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct object with PointIds and CellIds on; and ids being generated
  // as scalars.
  static ExtendIdFilter *New();

  // Description:
  // Enable/disable the generation of point ids. Default is on.
  vtkSetMacro(PointIds,int);
  vtkGetMacro(PointIds,int);
  vtkBooleanMacro(PointIds,int);

  // Description:
  // Enable/disable the generation of point ids. Default is on.
  vtkSetMacro(CellIds,int);
  vtkGetMacro(CellIds,int);
  vtkBooleanMacro(CellIds,int);

  // Description:
  // Set/Get the flag which controls whether to generate scalar data
  // or field data. If this flag is off, scalar data is generated.
  // Otherwise, field data is generated. Default is off.
  vtkSetMacro(FieldData,int);
  vtkGetMacro(FieldData,int);
  vtkBooleanMacro(FieldData,int);

  // Description:
  // Set/Get the name of the Ids array if generated. By default the Ids
  // are named "ExtendIdFilter_Ids", but this can be changed with this function.
  vtkSetStringMacro(IdsArrayName);
  vtkGetStringMacro(IdsArrayName);

protected:
  ExtendIdFilter();
  ~ExtendIdFilter();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  int PointIds;
  int CellIds;
  int FieldData;
  char *IdsArrayName;

private:
  ExtendIdFilter(const ExtendIdFilter&);  // Not implemented.
  void operator=(const ExtendIdFilter&);  // Not implemented.
};

#endif


