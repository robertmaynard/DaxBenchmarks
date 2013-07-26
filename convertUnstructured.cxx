#include <stdlib.h>

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

	//
	// assign all points
	//
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	int x,y,z;
	for (z=0; z<dim[2]; z++)
		for (y=0; y<dim[1]; y++)
			for (x=0; x<dim[0]; x++)
			{
				float p[3];
				p[0] = x; p[1] = y; p[2] = z;
				points->InsertNextPoint(p);
			}

	vtkSmartPointer<vtkCellArray> hexs =
	    vtkSmartPointer<vtkCellArray>::New();
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

	vtkSmartPointer<vtkHexahedron> hex =
	    vtkSmartPointer<vtkHexahedron>::New();

	vtkSmartPointer<vtkUnstructuredGrid> uGrid =
	    vtkSmartPointer<vtkUnstructuredGrid>::New();
	uGrid->SetPoints(points);
	uGrid->SetCells(VTK_HEXAHEDRON, hexs);
	uGrid->GetPointData()->AddArray( imageData->GetPointData()->GetScalars() );

	return uGrid;
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
	if (argc<2) {
		cout << "Usage: convertUnstructured filename resample_ratio" << endl;
		exit(0);
	}
	const char *filename = argv[1];
	float ratio = atof(argv[2]);
	vtkSmartPointer<vtkImageData> imageData = readData(filename, ratio);

	vtkSmartPointer<vtkUnstructuredGrid> uGrid = convertSimple(imageData);

	char s[100];
	sprintf(s, "%s_r%g.vtu", filename, ratio );
	writeData(uGrid, s);

	return 0;
}
