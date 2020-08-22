#include "FFRCore.h"

#include <queue>
#include <set>
#include <tuple>
#include "vtkMath.h"
#include "vtkSmartPointer.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkTree.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"
#include "vtkIdList.h"

void LengthCenterline(vtkTree * tree, vtkDoubleArray *length)
{
	for(vtkIdType id=0; id<tree->GetNumberOfVertices(); id++)
	{
		if( id == tree->GetRoot() ) length->SetValue(id, 0.0);
		else 
		{
			double dist = 0.0;
			if( tree->GetParent(id) != tree->GetRoot() ) 
			{
				double coord1[3], coord2[3];
				tree->GetPoints()->GetPoint(id, coord1);
				tree->GetPoints()->GetPoint(tree->GetParent(id), coord2);
				dist = vtkMath::Distance2BetweenPoints(coord1, coord2);
			}
			length->SetValue(id, dist);
		}
	}
}

void SurfaceAreaCenterline(vtkTree * tree, vtkDoubleArray *surfacearea, vtkDoubleArray *radius)
{
	for(vtkIdType id=0; id<tree->GetNumberOfVertices(); id++)
	{
		if( id == tree->GetRoot() ) surfacearea->SetValue(id, 0.0);
		else 
		{
			if( tree->GetNumberOfChildren(id) != 0 ) 
			{
				surfacearea->SetValue(id, 0.0);
			}
			else
			{
				surfacearea->SetValue(id, 3.1416*radius->GetValue(id)*radius->GetValue(id));
			}
		}
	}
}

void InitRPArray(vtkTree *tgraph)
{
	if( tgraph->GetVertexData()->GetArray("R") == NULL )
	{
		vtkSmartPointer<vtkDoubleArray>  Rarray = vtkSmartPointer<vtkDoubleArray>::New();
		Rarray->SetName("R");
		Rarray->SetNumberOfValues(tgraph->GetNumberOfVertices());
		tgraph->GetVertexData()->AddArray(Rarray);
	}

	if( tgraph->GetVertexData()->GetArray("P") == NULL )
	{
		vtkSmartPointer<vtkDoubleArray>  Parray = vtkSmartPointer<vtkDoubleArray>::New();
		Parray->SetName("P");
		Parray->SetNumberOfValues(tgraph->GetNumberOfVertices());
		tgraph->GetVertexData()->AddArray(Parray);
	}

	vtkDoubleArray *Rarray = vtkDoubleArray::SafeDownCast(tgraph->GetVertexData()->GetArray("R"));
	vtkDoubleArray *Parray = vtkDoubleArray::SafeDownCast(tgraph->GetVertexData()->GetArray("P"));
	for(vtkIdType i=0; i<tgraph->GetNumberOfVertices(); i++)
	{
		Rarray->SetValue(i, 0.0);
		Parray->SetValue(i, 0.0);
	}
}

void UpdateRArray(vtkTree *tgraph, vtkIdType root, double alpha, double beta, vtkDoubleArray *radius, vtkDoubleArray *length, vtkDoubleArray *surfacearea, vtkDoubleArray *rarray)
{
	if( tgraph->GetNumberOfChildren(root)==0 ) 
	{
		double sa = surfacearea->GetValue(root);
		if(sa<0.0001) sa = 3.1416*radius->GetValue(root)*radius->GetValue(root);
		rarray->SetValue(root, beta/sa );
		return;
	}

	for(int i=0; i<tgraph->GetNumberOfChildren(root); i++)
	{
		vtkIdType cid = tgraph->GetChild(root, i);
		UpdateRArray(tgraph, cid, alpha, beta, radius, length, surfacearea, rarray);
	}

	if( tgraph->GetNumberOfChildren(root) == 1 )
	{
		vtkIdType cid = tgraph->GetChild(root, 0);
		double rR = rarray->GetValue(cid) + alpha * length->GetValue(cid) / pow(radius->GetValue(cid), 7.0);
		rarray->SetValue(root, rR);
	}
	else
	{
		double rinv = 0.0;
		for(int i=0; i<tgraph->GetNumberOfChildren(root); i++)
		{
			vtkIdType cid = tgraph->GetChild(root, i);
			rinv += 1.0/rarray->GetValue(cid);
		}
		rarray->SetValue(root, 1.0/rinv);
	}
}

void UpdatePArray(vtkTree *tgraph, vtkIdType root, vtkDoubleArray *rarray, vtkDoubleArray *parray)
{
	if( root == tgraph->GetRoot() ) 
	{
		parray->SetValue(root, 1.0 );
	}
	else
	{
		vtkIdType pid = tgraph->GetParent(root);
		if( tgraph->GetNumberOfChildren(pid) == 1 )
		{
			parray->SetValue( root, parray->GetValue(pid) * rarray->GetValue(root) / rarray->GetValue(pid) );
		}
		else
		{
			parray->SetValue( root, parray->GetValue(pid) );
		}
	}

	for(int i=0; i<tgraph->GetNumberOfChildren(root); i++)
	{
		UpdatePArray(tgraph, tgraph->GetChild(root, i), rarray, parray);
	}
}

void UpdateRPArray(vtkTree *tgraph, double alpha, double beta)
{
	vtkDoubleArray *radius = vtkDoubleArray::SafeDownCast(tgraph->GetVertexData()->GetArray("Radius"));
	vtkDoubleArray *length = vtkDoubleArray::SafeDownCast(tgraph->GetVertexData()->GetArray("Length"));
	vtkDoubleArray *surfacearea = vtkDoubleArray::SafeDownCast(tgraph->GetVertexData()->GetArray("SurfaceArea"));
	vtkDoubleArray *rarray = vtkDoubleArray::SafeDownCast(tgraph->GetVertexData()->GetArray("R"));
	vtkDoubleArray *parray = vtkDoubleArray::SafeDownCast(tgraph->GetVertexData()->GetArray("P"));
	UpdateRArray(tgraph, tgraph->GetRoot(), alpha, beta, radius, length, surfacearea, rarray);
	UpdatePArray(tgraph, tgraph->GetRoot(), rarray, parray);
}

//bool CheckLicense()
//{
//	char *product = "SmartCoronary";
//	char *version = "1.0";
//	int  stat, days;
//	stat = rlmez_checkout(product, version, &days);
//	if (stat)
//	{
//		char errstring[RLM_ERRSTRING_MAX];
//		printf("checkout of %s failed, status: %d: %s\n", product, stat,
//			rlmez_errstring(stat, errstring));
//		return false;
//	}
//	return true;
//}

void ComputeFFRCore(vtkPolyData *clModel, vtkPolyData *clPoly, vtkDoubleArray *ffrArray, double leftOstium[3], double rightOstium[3])
{
	//if( !CheckLicense() ) return;

	printf("Start ComputeFFRCore\n");

	// Find points closest to left and right ostiums
	clModel->BuildCells();
	clModel->BuildLinks();
	clPoly->BuildCells();
	clPoly->BuildLinks();
	double coord[3];
	double leftOstiumDist = std::numeric_limits<double>::max();
	double rightOstiumDist = std::numeric_limits<double>::max();
	vtkIdType leftOstiumPointIndex = -1;
	vtkIdType rightOstiumPointIndex = -1;
	vtkIdType leftOstiumGraphIndex = -1;
	vtkIdType rightOstiumGraphIndex = -1;
	vtkIdType leftOstiumCellIndex = -1;
	vtkIdType rightOstiumCellIndex = -1;
	for(vtkIdType id=0; id<clModel->GetNumberOfCells(); id++)
	{
		vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
		clModel->GetCellPoints(id, idlist);
		if(idlist->GetNumberOfIds() < 2) continue;
		for(int j=0; j<idlist->GetNumberOfIds(); j=j+idlist->GetNumberOfIds()-1)
		{
			clModel->GetPoint(idlist->GetId(j), coord);
			vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
			clPoly->GetPointCells(idlist->GetId(j), cellIds);
			if( cellIds->GetNumberOfIds() != 1 ) continue;
			double dist = vtkMath::Distance2BetweenPoints(leftOstium, coord);
			if( dist < leftOstiumDist )
			{
				leftOstiumDist = dist;
				leftOstiumCellIndex = id;
				vtkSmartPointer<vtkIdList> il = vtkSmartPointer<vtkIdList>::New();
				clPoly->GetCellPoints(id, il);
				leftOstiumPointIndex = il->GetId(j);
			}
			dist = vtkMath::Distance2BetweenPoints(rightOstium, coord);
			if( dist < rightOstiumDist )
			{
				rightOstiumDist = dist;
				rightOstiumCellIndex = id;
				vtkSmartPointer<vtkIdList> il = vtkSmartPointer<vtkIdList>::New();
				clPoly->GetCellPoints(id, il);
				rightOstiumPointIndex = il->GetId(j);
			}
		}
	}		

	// only use radius
	vtkDoubleArray *clRadius = vtkDoubleArray::SafeDownCast(clPoly->GetPointData()->GetArray("Radius"));
	if(!clRadius)
	{
		std::cerr << "Cannot find radius array" << std::endl;
		return;
	}
	vtkSmartPointer<vtkIdTypeArray> polyPointIds = vtkSmartPointer<vtkIdTypeArray>::New();
	polyPointIds->SetName("PointIds");
	polyPointIds->SetNumberOfValues(clPoly->GetNumberOfPoints());
	for(vtkIdType id=0; id<clPoly->GetNumberOfPoints(); id++)
	{
		polyPointIds->SetValue(id, -1);
	}

	vtkSmartPointer<vtkMutableDirectedGraph> dgraph = vtkSmartPointer<vtkMutableDirectedGraph>::New();
	vtkSmartPointer<vtkPoints> graphPoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkDoubleArray> graphRadius = vtkSmartPointer<vtkDoubleArray>::New();
	graphRadius->SetName("Radius");
	graphRadius->SetNumberOfComponents(1);
	std::queue< std::tuple<vtkIdType, vtkIdType, vtkIdType> > transverseQueue;
	if( sqrt(leftOstiumDist) < 5.0 )
	{
		clPoly->GetPoint(leftOstiumPointIndex, coord);
		vtkIdType newId = dgraph->AddVertex();
		leftOstiumGraphIndex = newId;
		graphPoints->InsertNextPoint(coord);
		polyPointIds->SetValue(leftOstiumPointIndex, newId);
		graphRadius->InsertNextValue(clRadius->GetValue(leftOstiumPointIndex));
		transverseQueue.push(std::make_tuple(newId, leftOstiumPointIndex, leftOstiumCellIndex));
	}
	if( sqrt(rightOstiumDist) < 5.0 )
	{
		clPoly->GetPoint(rightOstiumPointIndex, coord);
		vtkIdType newId = dgraph->AddVertex();
		rightOstiumGraphIndex = newId;
		graphPoints->InsertNextPoint(coord);
		polyPointIds->SetValue(rightOstiumPointIndex, newId);
		graphRadius->InsertNextValue(clRadius->GetValue(rightOstiumPointIndex));
		transverseQueue.push(std::make_tuple(newId, rightOstiumPointIndex, rightOstiumCellIndex));
	}
	if( leftOstiumPointIndex >=0 && rightOstiumPointIndex >=0 )
	{
		double rootcoord[3] = {0.0, 0.0, 0.0};
		clPoly->GetPoint(leftOstiumPointIndex, coord);
		vtkMath::Add(rootcoord, coord, rootcoord);
		clPoly->GetPoint(rightOstiumPointIndex, coord);
		vtkMath::Add(rootcoord, coord, rootcoord);
		vtkMath::MultiplyScalar(rootcoord, 0.5);
		vtkIdType rootId = dgraph->AddVertex();
		graphPoints->InsertNextPoint(rootcoord);
		graphRadius->InsertNextValue( (clRadius->GetValue(leftOstiumPointIndex)+clRadius->GetValue(rightOstiumPointIndex))/2.0 );

		dgraph->AddEdge(rootId, leftOstiumGraphIndex);
		dgraph->AddEdge(rootId, rightOstiumGraphIndex);
	}
	std::set<vtkIdType> visitedCells;    //avoid revisiting a cell (e.g. there is a loop in the centerline).
	while( !transverseQueue.empty() )
	{
		std::tuple<vtkIdType, vtkIdType, vtkIdType> qitem = transverseQueue.front();
		transverseQueue.pop();
		vtkIdType attachId = std::get<0>(qitem);
		vtkIdType startId = std::get<1>(qitem);
		vtkIdType cellId = std::get<2>(qitem);
		if( visitedCells.find(cellId) != visitedCells.end() ) 
		{
			std::cerr << "CellId: " << cellId << " already visited, which may be caused by a loop." << std::endl;
			continue;
		}
		visitedCells.insert(cellId);

		vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();
		clPoly->GetCellPoints(cellId, pointIds);
		if( startId == pointIds->GetId(0) )
		{
			for(int i=1; i<pointIds->GetNumberOfIds(); i++)
			{
				clPoly->GetPoint(pointIds->GetId(i), coord);
				vtkIdType newId = dgraph->AddVertex();
				graphPoints->InsertNextPoint(coord);
				polyPointIds->SetValue(pointIds->GetId(i), newId);
				graphRadius->InsertNextValue(clRadius->GetValue(pointIds->GetId(i)));
				dgraph->AddEdge(attachId, newId);
				attachId = newId;
			}
			vtkSmartPointer<vtkIdList> pIds = vtkSmartPointer<vtkIdList>::New();
			clModel->GetCellPoints(cellId, pIds);
			vtkIdType sId = pIds->GetId(pIds->GetNumberOfIds()-1);
			vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
			clModel->GetPointCells(sId, cellIds);
			for(int i=0; i<cellIds->GetNumberOfIds(); i++)
			{
				vtkIdType cid = cellIds->GetId(i);
				if( cid != cellId )
				{
					vtkSmartPointer<vtkIdList> cIds = vtkSmartPointer<vtkIdList>::New();
					clModel->GetCellPoints(cid, cIds);
					if( cIds->GetId(0) == sId )
					{
						clPoly->GetCellPoints(cid, cIds);
						startId = cIds->GetId(0);
					}
					else if( cIds->GetId(cIds->GetNumberOfIds()-1) == sId )
					{
						clPoly->GetCellPoints(cid, cIds);
						startId = cIds->GetId(cIds->GetNumberOfIds()-1);
					}
					else
					{
						std::cerr << "clModel and clPoly topological mismatch detected" << std::endl;
						continue;
					}
					polyPointIds->SetValue(startId, polyPointIds->GetValue(pointIds->GetId(pointIds->GetNumberOfIds()-1)));
					transverseQueue.push(std::make_tuple(attachId, startId, cid));
				}
			}
		}
		else if( startId == pointIds->GetId(pointIds->GetNumberOfIds()-1) )
		{
			for(int i=pointIds->GetNumberOfIds()-2; i>=0; i--)
			{
				clPoly->GetPoint(pointIds->GetId(i), coord);
				vtkIdType newId = dgraph->AddVertex();
				graphPoints->InsertNextPoint(coord);
				polyPointIds->SetValue(pointIds->GetId(i), newId);
				graphRadius->InsertNextValue(clRadius->GetValue(pointIds->GetId(i)));
				dgraph->AddEdge(attachId, newId);
				attachId = newId;
			}
			vtkSmartPointer<vtkIdList> pIds = vtkSmartPointer<vtkIdList>::New();
			clModel->GetCellPoints(cellId, pIds);
			vtkIdType sId = pIds->GetId(0);
			vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
			clModel->GetPointCells(sId, cellIds);
			for(int i=0; i<cellIds->GetNumberOfIds(); i++)
			{
				vtkIdType cid = cellIds->GetId(i);
				if( cid != cellId )
				{
					vtkSmartPointer<vtkIdList> cIds = vtkSmartPointer<vtkIdList>::New();
					clModel->GetCellPoints(cid, cIds);
					if( cIds->GetId(0) == sId )
					{
						clPoly->GetCellPoints(cid, cIds);
						startId = cIds->GetId(0);
					}
					else if( cIds->GetId(cIds->GetNumberOfIds()-1) == sId )
					{
						clPoly->GetCellPoints(cid, cIds);
						startId = cIds->GetId(cIds->GetNumberOfIds()-1);
					}
					else
					{
						std::cerr << "clModel and clPoly topological mismatch detected" << std::endl;
						continue;
					}
					polyPointIds->SetValue(startId, polyPointIds->GetValue(pointIds->GetId(0)));
					transverseQueue.push(std::make_tuple(attachId, startId, cid));
				}
			}

		}
		else
		{
			std::cerr << "startId is not equal to either beginning or end of the cell" << std::endl;
		}
	}
	dgraph->SetPoints(graphPoints);
	dgraph->GetVertexData()->AddArray(graphRadius);

	vtkSmartPointer<vtkTree> tree = vtkSmartPointer<vtkTree>::New();
	if(!tree->CheckedShallowCopy(dgraph))
	{
		std::cerr << "Cannot convert to a tree" << std::endl;
		return;
	}

	vtkSmartPointer<vtkDoubleArray> length = vtkSmartPointer<vtkDoubleArray>::New();
	length->SetName("Length");
	length->SetNumberOfValues(tree->GetNumberOfVertices());
	LengthCenterline(tree, length);
	tree->GetVertexData()->AddArray(length);

	vtkDoubleArray *radius = vtkDoubleArray::SafeDownCast(tree->GetVertexData()->GetArray("Radius"));
	vtkSmartPointer<vtkDoubleArray> surfaceArea = vtkSmartPointer<vtkDoubleArray>::New();
	surfaceArea->SetName("SurfaceArea");
	surfaceArea->SetNumberOfValues(tree->GetNumberOfVertices());
	SurfaceAreaCenterline(tree, surfaceArea, radius);
	tree->GetVertexData()->AddArray(surfaceArea);

	InitRPArray(tree);
	UpdateRPArray(tree, 10.0, 100000.0);

	ffrArray->SetName("FFR");
	ffrArray->SetNumberOfValues(clPoly->GetNumberOfPoints());
	vtkDoubleArray *parray = vtkDoubleArray::SafeDownCast(tree->GetVertexData()->GetArray("P"));
	for(vtkIdType id=0; id<clPoly->GetNumberOfPoints(); id++)
	{
		vtkIdType pid = polyPointIds->GetValue(id);
		if( pid >=0 )
		{
			ffrArray->SetValue(id, parray->GetValue(pid));
		}
		else
		{
			ffrArray->SetValue(id, 1.0);
		}
	}
}

