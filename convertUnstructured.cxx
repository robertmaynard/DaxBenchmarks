#include <stdlib.h>
#include <assert.h>

#include <iostream>
#include <vector>

#include <vtkNew.h>
#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkImageResample.h>
#include <vtkNrrdReader.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkHexahedron.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkCell.h>
#include <vtkIdList.h>

#include "zcurve.h"

using namespace std;

static vtkSmartPointer<vtkImageData>
readData(std::string file,  double resampleRatio=1.0) // 1.0
{
  std::cout << "reading file: " << file << std::endl;
  std::cout << "Resample rate: " << resampleRatio << std::endl;
  vtkNew<vtkNrrdReader> reader;
  reader->SetFileName(file.c_str());
  reader->Update();

  //re-sample the dataset
  vtkNew<vtkImageResample> resample;
  resample->SetInputConnection(reader->GetOutputPort());
  resample->SetAxisMagnificationFactor(0,resampleRatio);
  resample->SetAxisMagnificationFactor(1,resampleRatio);
  resample->SetAxisMagnificationFactor(2,resampleRatio);

  resample->Update();

  vtkImageData *newImageData = vtkImageData::SafeDownCast(resample->GetOutputDataObject(0));
#if 0
  //take ref
  vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
  vtkImageData *newImageData = vtkImageData::SafeDownCast(resample->GetOutputDataObject(0));
  image.TakeReference( newImageData );
  image->Register(NULL);

  //now set the buffer
  vtkDataArray *newData = image->GetPointData()->GetScalars();
  dax::Scalar* rawBuffer = reinterpret_cast<dax::Scalar*>( newData->GetVoidPointer(0) );
  buffer.resize( newData->GetNumberOfTuples() );
  std::copy(rawBuffer, rawBuffer + newData->GetNumberOfTuples(), buffer.begin() );
#endif
  return newImageData;
}

vtkSmartPointer<vtkUnstructuredGrid>
convertSimple(vtkImageData *imageData)
{
	// get dim
	int dim[3];
	imageData->GetDimensions(dim);
	cout << "resampled data width: " << dim[0] << "," << dim[1] << "," << dim[2] << endl;
	int size = dim[0] * dim[1] * dim[2];

	//
	// assign all points
	//
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(size);
	int x,y,z;
	int i=0;
	for (z=0; z<dim[2]; z++)
		for (y=0; y<dim[1]; y++)
			for (x=0; x<dim[0]; x++)
			{
				double p[3];
				p[0] = x; p[1] = y; p[2] = z;
				points->SetPoint(i++, p);
			}

	vtkSmartPointer<vtkCellArray> hexs =
	    vtkSmartPointer<vtkCellArray>::New();
	//hexs->SetNumberOfCells((vtkIdType)size); <- size non-equal cause problem
	for (z=0; z<dim[2]-1; z++)
		for (y=0; y<dim[1]-1; y++)
			for (x=0; x<dim[0]-1; x++)
			{
				// Setup the coordinates of eight points
				// (the two faces must be in counter clockwise order as viewd from the outside)
				vtkSmartPointer<vtkHexahedron> hex =
				    vtkSmartPointer<vtkHexahedron>::New();
#define GET_ID(x,y,z) ( (x) + dim[0]*( (y) + dim[1]*(z) ))
				hex->GetPointIds()->SetId(0,GET_ID(x  , y  , z));
				hex->GetPointIds()->SetId(1,GET_ID(x+1, y  , z));
				hex->GetPointIds()->SetId(2,GET_ID(x+1, y+1, z));
				hex->GetPointIds()->SetId(3,GET_ID(x  , y+1, z));
				hex->GetPointIds()->SetId(4,GET_ID(x  , y  , z+1));
				hex->GetPointIds()->SetId(5,GET_ID(x+1, y  , z+1));
				hex->GetPointIds()->SetId(6,GET_ID(x+1, y+1, z+1));
				hex->GetPointIds()->SetId(7,GET_ID(x  , y+1, z+1));
#undef GET_ID
				hexs->InsertNextCell(hex);
			}

	vtkSmartPointer<vtkUnstructuredGrid> uGrid =
	    vtkSmartPointer<vtkUnstructuredGrid>::New();
	uGrid->SetPoints(points);
	uGrid->SetCells(VTK_HEXAHEDRON, hexs);
	uGrid->GetPointData()->AddArray( imageData->GetPointData()->GetScalars() );

	return uGrid;
}


vtkSmartPointer<vtkUnstructuredGrid>
reorder_zcurve( vtkSmartPointer<vtkUnstructuredGrid> inGrid, vtkSmartPointer<vtkImageData> image )
{
	vtkSmartPointer<vtkUnstructuredGrid> outGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	int dim[4];
	dim[3] = 1;
	image->GetDimensions(dim);

	cout << "Generating order" << endl;
	vector<int> mappingID;
	gen_zcurve(dim, mappingID);
	int size = mappingID.size();
	printf("%d %d %d %d , size=%d\n", dim[0], dim[1], dim[2], dim[3], size);

	int dim1[4];
	dim1[0] = dim[0] -1;
	dim1[1] = dim[1] -1;
	dim1[2] = dim[2] -1;
	dim1[3] = 1;
	vector<int> mappingID1;
	gen_zcurve(dim1, mappingID1);
	int size1 = mappingID1.size();


	//
	// start reordering points
	//
	cout << "reordering points" << endl;

	vtkSmartPointer<vtkPoints> inPoints = inGrid->GetPoints();
	vtkSmartPointer<vtkPoints> outPoints = vtkSmartPointer<vtkPoints>::New();
	int i;
	assert(inPoints->GetNumberOfPoints() == size);
	outPoints->SetNumberOfPoints(size);
	for (i=0; i<size; i++)
	{
		double *p = inPoints->GetPoint(mappingID[i]);
		outPoints->SetPoint(i, p);
		assert(mappingID[i] < size);
		//printf("Point: %lf %lf %lf\n", p[0], p[1], p[2]);
	}
	outGrid->SetPoints(outPoints);

	cout << "reordering field" << endl;
	vtkFloatArray *inField = vtkFloatArray::SafeDownCast( inGrid->GetPointData()->GetArray(0) );
	vtkSmartPointer<vtkFloatArray> outField = vtkSmartPointer<vtkFloatArray>::New();
	outField->SetName( inField->GetName() );
	//outField->SetNumberOfValues(size);
	assert(inField->GetNumberOfValues() == size);
	for (i=0; i<size; i++)
	{
		outField->InsertNextValue(inField->GetValue(mappingID[i]));
		assert(inField->GetValue(mappingID[i]) < size);
		//printf("outField: inField(%d)=%d\n", mappingID[i], inField->GetValue(mappingID[i]));
	}
	//inField->Delete();
	outGrid->GetPointData()->AddArray( outField );


	cout << "set and reassign cells" << endl;
	vector<int> mappingPos(size);
	for (i=0; i<size; i++)
		mappingPos[ mappingID[i] ] = i;
	vtkSmartPointer<vtkCellArray> outHexs =
			vtkSmartPointer<vtkCellArray>::New();
	assert(inHexs->GetNumberOfCells() == size1);
	//outHexs->SetNumberOfCells(size1);
	for (i=0; i<size1; i++)
	{
		//vtkSmartPointer<vtkHexahedron> outHex = vtkHexahedron::New();
		vtkSmartPointer<vtkHexahedron> hex = vtkHexahedron::SafeDownCast( inGrid->GetCell( mappingID1[i] ) );
		for (int j=0; j<8; j++)
		{
			hex->GetPointIds()->SetId( j, mappingPos[ hex->GetPointId(j) ] );
			//printf("SetId[%d]=mapping[inHex(%d)]\n", j, inHex->GetPointId(j));
			assert (inHex->GetPointId(j) < size);
		}
		outHexs->InsertNextCell(hex);
		//hex->Delete();
	}
	outGrid->SetCells(VTK_HEXAHEDRON, outHexs);

	return outGrid;
}

void writeData(vtkSmartPointer<vtkUnstructuredGrid> uGrid, string str)
{
	cout << "Writing file " << str << endl;
	// Write file
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
		vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(str.c_str());
	writer->SetInputData(uGrid);
	writer->Write();
}

int main(int argc, const char **argv)
{
	if (argc==1) {
		cout << "Usage: convertUnstructured filename resample_ratio zcurve=1" << endl;
		exit(0);
	}
	const char *filename = argv[1];
	float ratio = 1.f;
	if (argc>2) ratio = atof(argv[2]);
	bool use_zcurve = false;
	if (argc>3) use_zcurve = atoi(argv[3]);

	vtkSmartPointer<vtkImageData> imageData = readData(filename, ratio);

	vtkSmartPointer<vtkUnstructuredGrid> uGrid = convertSimple(imageData);

	char s[100];
	sprintf(s, "%s_r%g.vtu", filename, ratio);
	writeData(uGrid, s);

	if (use_zcurve) {
		vtkSmartPointer<vtkUnstructuredGrid> uGridr = reorder_zcurve(uGrid, imageData);

		char s[100];
		sprintf(s, "%s_r%gz.vtu", filename, ratio );
		writeData(uGridr, s);
	}

	return 0;
}
