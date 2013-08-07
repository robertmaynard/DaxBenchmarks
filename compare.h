
#include "compare_mc.h"
#include "compare_thresh.h"


#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkImageResample.h>
#include <vtkNrrdReader.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkNew.h>

#include <iostream>
#include <vector>

#include "SharedStatus.h"
#include "tlog/tlog.h"

static const int NUM_TRIALS = 10;

static vtkSmartPointer<vtkImageData>
ReadData(std::vector<dax::Scalar> &buffer, std::string file,  double resampleSize=1.0)
{
  //make sure we are testing float benchmarks only
  assert(sizeof(float) == sizeof(dax::Scalar));

  std::cout << "reading file: " << file << std::endl;
  vtkNew<vtkNrrdReader> reader;
  reader->SetFileName(file.c_str());
  reader->Update();

  //re-sample the dataset
  vtkNew<vtkImageResample> resample;
  resample->SetInputConnection(reader->GetOutputPort());
  resample->SetAxisMagnificationFactor(0,resampleSize);
  resample->SetAxisMagnificationFactor(1,resampleSize);
  resample->SetAxisMagnificationFactor(2,resampleSize);

  resample->Update();

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

  return image;
}


static vtkSmartPointer<vtkUnstructuredGrid>
ReadData_uGrid(std::string file)
{
  //make sure we are testing float benchmarks only
  assert(sizeof(float) == sizeof(dax::Scalar));

  std::cout << "reading file: " << file << std::endl;
  vtkNew<vtkXMLUnstructuredGridReader> reader;
  reader->SetFileName(file.c_str());
  reader->Update();

  //take ref
  vtkSmartPointer<vtkUnstructuredGrid> uGrid = reader->GetOutput(0);
  //std::cout << uGrid->GetPointData() << std::endl;

  return uGrid;
}


void WriteData(vtkSmartPointer<vtkPolyData> uGrid, std::string str)
{
	cout << "Writing file " << str << endl;
	// Write file
	vtkNew<vtkXMLPolyDataWriter> writer;
	writer->SetFileName(str.c_str());
	writer->SetInputData(uGrid);
	writer->Write();
}


int RunComparison(std::string device, std::string file, int pipeline, double resample_ratio=1.0)
{
  // Jimmy added
  SharedStatus::init();
  tlog = new TLog();
  tlog->regThread("Main");

  vtkSmartPointer< vtkUnstructuredGrid > image = ReadData_uGrid(file);

  //get dims of image data
  int dims[3];
  //image->GetDimensions(dims);

  //pipeline 1 is equal to threshold
  if(pipeline <= 1)
  {
#if 0
    const bool WithPointResolution=true;

    std::cout << "Warming up the machine" << std::endl;
    //warm up the card, whatever algorithm we run first will get a performance
    //hit for the first 10 iterations if we don't run something first
    RunDaxThreshold(dims,buffer,device,NUM_TRIALS,!WithPointResolution,true);

    //print out header of csv
    std::cout << "Benchmarking Threshold" << std::endl;

    std::cout << "DaxNoResolution,Accelerator,Time,Trial" << std::endl;
    RunDaxThreshold(dims,buffer,device,NUM_TRIALS,!WithPointResolution);

    std::cout << "DaxResolution,Accelerator,Time,Trial" << std::endl;
    RunDaxThreshold(dims,buffer,device,NUM_TRIALS, WithPointResolution);

    std::cout << "Piston,Accelerator,Time,Trial" << std::endl;
    RunPistonThreshold(dims,buffer,device,NUM_TRIALS);

    if(device == "Serial")
      {
      std::cout << "VTK,Accelerator,Time,Trial" << std::endl;
      RunVTKThreshold(image,NUM_TRIALS);
      }
#endif
  }
  else //marching cubes
  {

    const bool WithPointResolution=true;

    std::cout << "Warming up the machine" << std::endl;
    //warm up the card, whatever algorithm we run first will get a performance
    //hit for the first 10 iterations if we don't run something first
    RunDaxMarchingCubes(image,device,NUM_TRIALS,WithPointResolution,true);

    std::cout << "Benchmarking Marching Cubes" << std::endl;

    // Jimmy added
    //std::cout << "VTKDax,Accelerator,Time,Trial" << std::endl;
    //RunVTKDaxMarchingCubes(device, image, NUM_TRIALS);

    std::cout << "DaxNoResolution,Accelerator,Time,Trial" << std::endl;
    RunDaxMarchingCubes(image,device,NUM_TRIALS,!WithPointResolution);

    std::cout << "DaxResolution,Accelerator,Time,Trial" << std::endl;
    RunDaxMarchingCubes(image,device,NUM_TRIALS, WithPointResolution);

    std::cout << "Piston,Accelerator,Time,Trial" << std::endl;
    //RunPistonMarchingCubes(image,device,NUM_TRIALS);

    if(device == "Serial")
      {
      std::cout << "Serial,Accelerator,Time,Trial" << std::endl;
      RunVTKMarchingCubes(image,NUM_TRIALS);
      }

#if 0
    // debug result
    vtkNew<vtkPolyData> polyData;
    RunDaxMarchingCubes(image, device, 1, true, true, polyData.GetPointer());
    WriteData(polyData.GetPointer(), "out.vtp");
#endif

  }

  // Jimmy
  delete tlog;
  SharedStatus::getInstance()->print();

  return 0;
}
