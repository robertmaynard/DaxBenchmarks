
#include <dax/CellTag.h>
#include <dax/CellTraits.h>

#include <dax/cont/ArrayHandle.h>
#include <dax/cont/GenerateTopology.h>
#include <dax/cont/Scheduler.h>
#include <dax/cont/Timer.h>
#include <dax/cont/UniformGrid.h>
#include <dax/cont/UnstructuredGrid.h>
#include <dax/worklet/Threshold.h>

#include <vtkContourFilter.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkThreshold.h>
#include <vtkTrivialProducer.h>

#ifdef PISTON_ENABLED
#include "pistonImage3d.h"
#include "piston-threshold_geometry.h"
#endif PISTON_ENABLED

#include <vector>

static const dax::Scalar MIN_VALUE=0.07;
static const dax::Scalar MAX_VALUE=1.0;



static void RunDaxThreshold(int dims[3], std::vector<dax::Scalar>& buffer,
                                std::string device, int MAX_NUM_TRIALS,
                                bool enablePointResolution,
                                bool silent=false)
{
  dax::cont::UniformGrid<> grid;
  grid.SetExtent(dax::make_Id3(0, 0, 0), dax::make_Id3(dims[0]-1, dims[1]-1, dims[2]-1));

  typedef dax::cont::GenerateTopology<dax::worklet::ThresholdTopology> GenerateT;
  typedef GenerateT::ClassifyResultType  ClassifyResultType;

  //construct the scheduler that will execute all the worklets
  dax::cont::Scheduler<> scheduler;
  for(int i=0; i < MAX_NUM_TRIALS; ++i)
    {
    dax::cont::Timer<> timer;

    dax::cont::ArrayHandle<dax::Scalar> field = dax::cont::make_ArrayHandle(buffer);

    //construct the two worklets that will be used to do the marching cubes
    dax::worklet::ThresholdClassify<dax::Scalar> classifyWorklet(MIN_VALUE,MAX_VALUE);
    dax::worklet::ThresholdTopology generateWorklet;

    //run the first step
    ClassifyResultType classification; //array handle for the first step classification
    scheduler.Invoke(classifyWorklet, grid, field, classification);

    //construct the topology generation worklet
    GenerateT generate(classification,generateWorklet);
    generate.SetRemoveDuplicatePoints(enablePointResolution);

    //run the second step
    dax::cont::UnstructuredGrid<dax::CellTagHexahedron> outGrid;

    //schedule threshold
    scheduler.Invoke(generate, grid, outGrid);

    double time = timer.GetElapsedTime();
    if(!silent)
      std::cout << "Dax," << device << "," << time << "," << i << std::endl;
    }
}

static void RunVTKThreshold(vtkImageData* image, int MAX_NUM_TRIALS)
{
  vtkNew<vtkTrivialProducer> producer;
  producer->SetOutput(image);
  producer->Update();

  for(int i=0; i < MAX_NUM_TRIALS; ++i)
    {

    vtkNew<vtkThreshold> threshold;
    threshold->SetInputConnection(producer->GetOutputPort());

    dax::cont::Timer<> timer;

    threshold->SetPointsDataTypeToFloat();
    threshold->AllScalarsOn();
    threshold->ThresholdBetween(MIN_VALUE, MAX_VALUE);

    threshold->Update();

    double time = timer.GetElapsedTime();
    std::cout << "VTK,Serial," << time << "," << i << std::endl;
    }
}


static void RunPistonThreshold(int dims[3], std::vector<dax::Scalar>& buffer,
                               std::string device, int MAX_NUM_TRIALS)
{
#ifdef PISTON_ENABLED
  typedef piston::threshold_geometry< piston_scalar_image3d > TG;

  for (int i=0; i < MAX_NUM_TRIALS; ++i)
    {
    //piston moves memory when constructing the marching cubes object
    dax::cont::Timer<> timer;
    piston_scalar_image3d pimage(dims[0],dims[1],dims[2],buffer);
    TG threshold(pimage,MIN_VALUE, MAX_VALUE);


    threshold();
    double time = timer.GetElapsedTime();
    std::cout << "Piston," << device << "," << time << "," << i << std::endl;
    }
#endif PISTON_ENABLED
}
