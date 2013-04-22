
#include <dax/CellTag.h>
#include <dax/CellTraits.h>

#include <dax/cont/ArrayHandle.h>
#include <dax/cont/GenerateInterpolatedCells.h>
#include <dax/cont/Scheduler.h>
#include <dax/cont/Timer.h>
#include <dax/exec/CellField.h>
#include <dax/math/Trig.h>
#include <dax/math/VectorAnalysis.h>
#include <dax/worklet/MarchingCubes.h>

#include <vtkContourFilter.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkTrivialProducer.h>

#ifdef PISTON_ENABLED
#include "pistonImage3d.h"
#include <piston/marching_cube.h>
#endif

#include <vector>

static const dax::Scalar ISO_VALUE=0.07;

namespace dax { namespace worklet {

struct Normals: public dax::exec::WorkletMapCell
{
typedef void ControlSignature(Topology, Field(In,Point),
                  Field(Out));
typedef _3 ExecutionSignature(_2);

DAX_EXEC_EXPORT
dax::Vector3 operator()(const dax::exec::CellField<
          dax::Vector3,
          dax::CellTagTriangle >& coordinates) const
  {
  const dax::Vector3 a = coordinates[1] - coordinates[0];
  const dax::Vector3 b = coordinates[2] - coordinates[0];
  return dax::math::Normal( dax::math::Cross(a,b) );
  }
};


} }


static void RunDaxMarchingCubes(int dims[3], std::vector<dax::Scalar>& buffer,
                                std::string device, int MAX_NUM_TRIALS,
                                bool enablePointResolution,
                                bool silent=false)
{
  dax::cont::UniformGrid<> grid;
  grid.SetExtent(dax::make_Id3(0, 0, 0), dax::make_Id3(dims[0]-1, dims[1]-1, dims[2]-1));

  typedef dax::cont::GenerateInterpolatedCells<dax::worklet::MarchingCubesTopology> GenerateIC;
  typedef GenerateIC::ClassifyResultType  ClassifyResultType;

  //construct the scheduler that will execute all the worklets
  dax::cont::Scheduler<> scheduler;


  for(int i=0; i < MAX_NUM_TRIALS; ++i)
    {
    dax::cont::Timer<> timer;

    dax::cont::ArrayHandle<dax::Scalar> field = dax::cont::make_ArrayHandle(buffer);

    //construct the two worklets that will be used to do the marching cubes
    dax::worklet::MarchingCubesClassify classifyWorklet(ISO_VALUE);
    dax::worklet::MarchingCubesTopology generateWorklet(ISO_VALUE);
    dax::worklet::Normals normWorklet;

    //run the first step
    ClassifyResultType classification; //array handle for the first step classification
    scheduler.Invoke(classifyWorklet, grid, field, classification);

    //construct the topology generation worklet
    GenerateIC generate(classification,generateWorklet);
    generate.SetRemoveDuplicatePoints(enablePointResolution);

    //run the second step
    dax::cont::UnstructuredGrid<dax::CellTagTriangle> outGrid;

    //schedule marching cubes worklet generate step, saving
    scheduler.Invoke(generate, grid, outGrid, field);

    //compute the normals of each output triangle

    dax::cont::ArrayHandle<dax::Vector3> normals;
    scheduler.Invoke(normWorklet, outGrid, outGrid.GetPointCoordinates(), normals);

    double time = timer.GetElapsedTime();
    if(!silent)
      std::cout << "Dax," << device << "," << time << "," << i << std::endl;
    }
}

static void RunVTKMarchingCubes(vtkImageData* image, int MAX_NUM_TRIALS)
{
  vtkNew<vtkTrivialProducer> producer;
  producer->SetOutput(image);
  producer->Update();

  for(int i=0; i < MAX_NUM_TRIALS; ++i)
    {

    vtkNew<vtkContourFilter> marching;
    marching->SetInputConnection(producer->GetOutputPort());

    dax::cont::Timer<> timer;

    marching->ComputeGradientsOff();
    marching->ComputeNormalsOn();
    marching->ComputeScalarsOn();
    marching->SetNumberOfContours(1);
    marching->SetValue(0, ISO_VALUE);

    marching->Update();

    double time = timer.GetElapsedTime();
    std::cout << "VTK,Serial," << time << "," << i << std::endl;
    }
}


static void RunPistonMarchingCubes(int dims[3], std::vector<dax::Scalar>& buffer,
                                   std::string device, int MAX_NUM_TRIALS)
{
#ifdef PISTON_ENABLED
  typedef piston::marching_cube< piston_scalar_image3d,
                                 piston_scalar_image3d > MC;

  for (int i=0; i < MAX_NUM_TRIALS; ++i)
    {
    dax::cont::Timer<> timer;
    piston_scalar_image3d pimage(dims[0],dims[1],dims[2],buffer);
    //piston moves memory when constructing the marching cubes object
    MC marching(pimage,pimage,ISO_VALUE);

    marching();
    double time = timer.GetElapsedTime();
    std::cout << "Piston," << device << "," << time << "," << i << std::endl;
    }
#endif
}
