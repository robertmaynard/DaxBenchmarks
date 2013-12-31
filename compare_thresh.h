
#include <dax/CellTag.h>
#include <dax/CellTraits.h>

#include <dax/cont/ArrayHandle.h>
#include <dax/cont/DispatcherGenerateTopology.h>
#include <dax/cont/DispatcherMapCell.h>
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
#endif

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

  typedef dax::cont::DispatcherGenerateTopology<
                                  dax::worklet::ThresholdTopology> DispatchTopo;
  typedef DispatchTopo::CountHandleType  CountHandleType;

  //construct the dispatcher for the map cell
  typedef dax::cont::DispatcherMapCell<
            dax::worklet::ThresholdCount<dax::Scalar> > DispatchThresholdCount;
  dax::worklet::ThresholdCount<dax::Scalar> countWorklet(MIN_VALUE,MAX_VALUE);

  for(int i=0; i < MAX_NUM_TRIALS; ++i)
    {
    dax::cont::Timer<> timer;

    dax::cont::ArrayHandle<dax::Scalar> field = dax::cont::make_ArrayHandle(buffer);

    //run the first step
    CountHandleType count; //array handle for the first step count
    DispatchThresholdCount( countWorklet ).Invoke(grid, field, count);

    //run the second step
    dax::cont::UnstructuredGrid<dax::CellTagHexahedron> outGrid;
    DispatchTopo generate(count);
    generate.SetRemoveDuplicatePoints(enablePointResolution);
    generate.Invoke(grid, outGrid);

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
#endif
}
