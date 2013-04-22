
#include "compare_mc.h"
#include "compare_thresh.h"


#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkImageResample.h>
#include <vtkNrrdReader.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <iostream>
#include <vector>

static const int NUM_TRIALS = 20;

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


int RunComparison(std::string device, std::string file, int pipeline)
{

  std::vector<dax::Scalar> buffer;
  double resample_ratio = 1.0; //full data
  vtkSmartPointer< vtkImageData > image = ReadData(buffer, file, resample_ratio);

  //get dims of image data
  int dims[3]; image->GetDimensions(dims);

  //pipeline 1 is equal to threshold
  if(pipeline <= 1)
  {
    //print out header of csv
    std::cout << "Benchmarking Threshold" << std::endl;

    bool WithPointResolution=true;
    std::cout << "DaxNoResolution,Accelerator,Time,Trial" << std::endl;
    RunDaxThreshold(dims,buffer,device,NUM_TRIALS,!WithPointResolution);

    std::cout << "DaxResolution,Accelerator,Time,Trial" << std::endl;
    RunDaxThreshold(dims,buffer,device,NUM_TRIALS, WithPointResolution);

    std::cout << "Piston,Accelerator,Time,Trial" << std::endl;
    RunPistonThreshold(image,device,NUM_TRIALS);

    if(device == "Serial")
      {
      std::cout << "VTK,Accelerator,Time,Trial" << std::endl;
      RunVTKThreshold(image,NUM_TRIALS);
      }
  }
  else //marching cubes
  {

    std::cout << "Benchmarking Marching Cubes" << std::endl;

    bool WithPointResolution=true;

    std::cout << "DaxNoResolution,Accelerator,Time,Trial" << std::endl;
    RunDaxMarchingCubes(dims,buffer,device,NUM_TRIALS,!WithPointResolution);

    std::cout << "DaxResolution,Accelerator,Time,Trial" << std::endl;
    RunDaxMarchingCubes(dims,buffer,device,NUM_TRIALS, WithPointResolution);

    std::cout << "Piston,Accelerator,Time,Trial" << std::endl;
    RunPistonMarchingCubes(image,device,NUM_TRIALS);

    if(device == "Serial")
      {
      RunVTKMarchingCubes(image,NUM_TRIALS);
      }
  }

  return 0;
}
