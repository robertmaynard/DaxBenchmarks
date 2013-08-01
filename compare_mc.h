
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

#include <vtkMarchingCubes.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkTrivialProducer.h>

#include <vtkToDax/CellTypeToType.h>
#include <vtkToDax/DataSetTypeToType.h>
#include <vtkToDax/DataSetConverters.h>
#include <vtkUnstructuredGrid.h>
#include <vtkHexahedron.h>
#include <vtkDataArray.h>

#ifdef PISTON_ENABLED
#include "pistonImage3d.h"
#include <piston/marching_cube.h>
#endif

// Jimmy added:
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkDataSet.h>
#include <vtkDaxMarchingCubes.h>
#include "vtkDaxPolyDataNormals.h"
#include <vtkXMLPolyDataWriter.h> // debug
#include <daxToVtk/DataSetConverters.h>
#include <daxToVtk/CellTypeToType.h>
#include <vtkToDax/DataSetConverters.h>
#include <vtkToDax/Containers.h>

#include <vector>
#include <string.h>

#include "SharedStatus.h"
#include "tlog/tlog.h"
#include "tlog/tlogDefs.h"

//#define DEBUG_NORMAL

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
#ifdef _PROFILE
	static long long pretime=0;
	if (Timer::getTimeMS()-pretime > 10) {
		tlog->logMessage("dax worklet");
		pretime = Timer::getTimeMS();
	}
#endif

  const dax::Vector3 a = coordinates[1] - coordinates[0];
  const dax::Vector3 b = coordinates[2] - coordinates[0];
  return dax::math::Normal( dax::math::Cross(a,b) );
  }

};


} }


/////////////////////////////////////////////////////////

bool GetPointNormals(vtkPolyData* polydata)
{
  //std::cout << "In GetPointNormals: " << polydata->GetNumberOfPoints() << std::endl;
  //std::cout << "Looking for point normals..." << std::endl;

  // Count points
  vtkIdType numPoints = polydata->GetNumberOfPoints();
  //std::cout << "There are " << numPoints << " points." << std::endl;

  // Count triangles
  vtkIdType numPolys = polydata->GetNumberOfPolys();
  //std::cout << "There are " << numPolys << " polys." << std::endl;

  ////////////////////////////////////////////////////////////////
  // Double normals in an array
  vtkDoubleArray* normalDataDouble =
    vtkDoubleArray::SafeDownCast(polydata->GetPointData()->GetArray("Normals"));

  if(normalDataDouble)
    {
    int nc = normalDataDouble->GetNumberOfTuples();
    //std::cout << "There are " << nc
    //        << " components in normalDataDouble" << std::endl;
    return true;
    }

  ////////////////////////////////////////////////////////////////
  // Double normals in an array
  vtkFloatArray* normalDataFloat =
    vtkFloatArray::SafeDownCast(polydata->GetPointData()->GetArray("Normals"));

  if(normalDataFloat)
    {
    int nc = normalDataFloat->GetNumberOfTuples();
    //std::cout << "There are " << nc
    //        << " components in normalDataFloat" << std::endl;
    return true;
    }

  ////////////////////////////////////////////////////////////////
  // Point normals
  vtkDoubleArray* normalsDouble =
    vtkDoubleArray::SafeDownCast(polydata->GetPointData()->GetNormals());

  if(normalsDouble)
    {
    //std::cout << "There are " << normalsDouble->GetNumberOfComponents()
    //          << " components in normalsDouble" << std::endl;
    return true;
    }

  ////////////////////////////////////////////////////////////////
  // Point normals
  vtkFloatArray* normalsFloat =
    vtkFloatArray::SafeDownCast(polydata->GetPointData()->GetNormals());

  if(normalsFloat)
    {
    //std::cout << "There are " << normalsFloat->GetNumberOfComponents()
    //          << " components in normalsFloat" << std::endl;
    return true;
    }

  /////////////////////////////////////////////////////////////////////
  // Generic type point normals
  vtkDataArray* normalsGeneric = polydata->GetPointData()->GetNormals(); //works
  if(normalsGeneric)
    {
    //std::cout << "There are " << normalsGeneric->GetNumberOfTuples()
    //          << " normals in normalsGeneric" << std::endl;

    double testDouble[3];
    normalsGeneric->GetTuple(0, testDouble);

    //std::cout << "Double: " << testDouble[0] << " "
    //          << testDouble[1] << " " << testDouble[2] << std::endl;

    // Can't do this:
    /*
    float testFloat[3];
    normalsGeneric->GetTuple(0, testFloat);

    std::cout << "Float: " << testFloat[0] << " "
              << testFloat[1] << " " << testFloat[2] << std::endl;
    */
    return true;
    }


  // If the function has not yet quit, there were none of these types of normals
  std::cout << "Normals not found!" << std::endl;
  return false;

}



//convert an unstructured grid type
namespace vtkToDax{
#if 0
//template<> struct CellTypeAndDataType<VTK_POLY_DATA,VTK_HEXAHEDRON>{enum{Valid=1};};

template<typename CellTypeToTypeDef> struct DataSetTypeToType<CellTypeToTypeDef,vtkHexahedron>
{
  typedef CellTypeToTypeDef CellTypeToType;
  typedef typename CellTypeToTypeDef::DaxCellType DaxCellType;
  enum {VTKDataSetType=VTK_UNSTRUCTURED_GRID};
  enum {Valid=(CellTypeAndDataType<VTK_UNSTRUCTURED_GRID,CellTypeToTypeDef::VTKCellType>::Valid)};
  typedef dax::cont::UnstructuredGrid<DaxCellType,
          vtkToDax::vtkTopologyContainerTag<CellTypeToType>,
          vtkToDax::vtkPointsContainerTag>
          DaxDataSetType;
};
#endif
} // namespace

///////////////////////////////////////////////////


static void RunDaxMarchingCubes(vtkUnstructuredGrid *data,
                                std::string device, int MAX_NUM_TRIALS,
                                bool enablePointResolution,
                                bool silent=false,
                                vtkPolyData *vtkPolyOut = NULL)
{
	//
	// convert data to dax
	//
	//now set the buffer
	vtkDataArray *newData = data->GetPointData()->GetArray(0);
	dax::Scalar* rawBuffer = reinterpret_cast<dax::Scalar*>( newData->GetVoidPointer(0) );
	std::vector<dax::Scalar> buffer( newData->GetNumberOfTuples() );
	std::copy(rawBuffer, rawBuffer + newData->GetNumberOfTuples(), buffer.begin() );

	typedef vtkToDax::CellTypeToType<vtkHexahedron> VTKCellTypeStruct;
	typedef vtkToDax::DataSetTypeToType<VTKCellTypeStruct, vtkUnstructuredGrid> DataSetTypeToTypeStruct;
	typedef DataSetTypeToTypeStruct::DaxDataSetType InputDataSetType;
	InputDataSetType topology =
		  vtkToDax::dataSetConverter(data, DataSetTypeToTypeStruct() );


  //
  // original codes
  //
  typedef dax::cont::GenerateInterpolatedCells<dax::worklet::MarchingCubesGenerate> GenerateIC;
  typedef GenerateIC::ClassifyResultType  ClassifyResultType;

  //construct the scheduler that will execute all the worklets
  dax::cont::Scheduler<> scheduler;

  // tlog
  int evt_id = tlog->createEventID(device.c_str(), 0, 255, 0);


  for(int i=0; i < MAX_NUM_TRIALS; ++i)
    {
		dax::cont::Timer<> timer;
		char s[100];
		sprintf(s, "Dax round %d, enablePointResolution=%d", i, enablePointResolution);
		tlog->startEvent(evt_id, s);

		dax::cont::ArrayHandle<dax::Scalar> field = dax::cont::make_ArrayHandle(buffer);

		//construct the two worklets that will be used to do the marching cubes
		dax::worklet::MarchingCubesClassify classifyWorklet(ISO_VALUE);
		dax::worklet::MarchingCubesGenerate generateWorklet(ISO_VALUE);
		dax::worklet::Normals normWorklet;

		//run the first step
		ClassifyResultType classification; //array handle for the first step classification
		scheduler.Invoke(classifyWorklet, topology, field, classification);

		//construct the topology generation worklet
		GenerateIC generate(classification,generateWorklet);
		generate.SetRemoveDuplicatePoints(enablePointResolution);

		//run the second step
		dax::cont::UnstructuredGrid<dax::CellTagTriangle,
			vtkToDax::vtkTopologyContainerTag< vtkToDax::CellTypeToType<vtkTriangle>  >,
			vtkToDax::vtkPointsContainerTag > outGrid;

		//schedule marching cubes worklet generate step, saving
		scheduler.Invoke(generate, topology, outGrid, field);

		//compute the normals of each output triangle

		dax::cont::ArrayHandle<dax::Vector3> normals;
		scheduler.Invoke(normWorklet, outGrid, outGrid.GetPointCoordinates(), normals);

		double time = timer.GetElapsedTime();
		tlog->endEvent(evt_id);

		//if(!silent)
			std::cout << "Dax," << device << "," << time << "," << i << std::endl;

	    // timing
		if (!silent) {
			if (enablePointResolution)
				SharedStatus::getInstance()->dax_mc_res_time.push_back(time);
			else
				SharedStatus::getInstance()->dax_mc_nores_time.push_back(time);
		}

		// covert the mesh into vtk format for debugging
		if (vtkPolyOut) {
			daxToVtk::dataSetConverter(outGrid, vtkPolyOut);
			std::cout << "vtkPoly Points: " << vtkPolyOut->GetNumberOfPoints() << ", polys: " << vtkPolyOut->GetNumberOfPolys() << std::endl;
		}

    }
}

static void RunVTKMarchingCubes(vtkUnstructuredGrid* image, int MAX_NUM_TRIALS)
{
#if 0
  vtkNew<vtkTrivialProducer> producer;
  producer->SetOutput(image);
  producer->Update();

  int evt_id = tlog->createEventID("VTK", 0, 255, 0);

  for(int i=0; i < MAX_NUM_TRIALS; ++i)
    {
		char s[100];
		sprintf(s, "VTK Serial round %d", i);
		tlog->startEvent(evt_id, s);

		vtkNew<vtkMarchingCubes> marching;
		marching->SetInputConnection(producer->GetOutputPort());

		dax::cont::Timer<> timer;

		marching->ComputeGradientsOff();
		marching->ComputeNormalsOn();
		marching->ComputeScalarsOn();
		marching->SetNumberOfContours(1);
		marching->SetValue(0, ISO_VALUE);

		marching->Update();

		double time = timer.GetElapsedTime();
		tlog->endEvent(evt_id);

		std::cout << "VTK,Serial," << time << "," << i << std::endl;


	    // timing
		SharedStatus::getInstance()->vtk_mc_time.push_back(time);

    }
#endif
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



static void RunVTKDaxMarchingCubes(std::string &device, vtkUnstructuredGrid* image, int MAX_NUM_TRIALS)
{
#if 0
	bool hasPointNormals;
  vtkNew<vtkTrivialProducer> producer;
  producer->SetOutput(image);
  producer->Update();

  int evt_id = tlog->createEventID(device.c_str(), 0, 128, 0);
  int evtMC = tlog->createEventID("Update MC filter", 0, 255, 0);
  int evtNorm = tlog->createEventID("Update Norm filter", 0, 255, 0);

#ifdef DEBUG_NORMAL

  // Count points
  vtkIdType numPoints = image->GetNumberOfPoints();
  std::cout << "There are " << numPoints << " points." << std::endl;

#endif

  for(int i=0; i < MAX_NUM_TRIALS; ++i)
    {
	char s[100];
	sprintf(s, "vtkDax round %d", i);
	tlog->startEvent(evt_id, s);
    dax::cont::Timer<> timer, timer_mc;

#ifdef DEBUG_NORMAL
        cout << "new marching cubes" << endl;
#endif
	vtkNew<vtkDaxMarchingCubes> marching;
    marching->SetInputConnection(producer->GetOutputPort());


    marching->ComputeGradientsOff();
    marching->ComputeNormalsOff();
    marching->ComputeScalarsOff();
    marching->SetNumberOfContours(1);
    marching->SetValue(0, ISO_VALUE);

	tlog->startEvent(evtMC);
    marching->Update(); //this calls marching::RequestData
    tlog->endEvent(evtMC);

#ifdef DEBUG_NORMAL
    // debug
    hasPointNormals = GetPointNormals(marching->GetOutput());
    std::cout << "Marching cubes result has point normals? " << hasPointNormals << std::endl;
#endif
    SharedStatus::getInstance()->vtkdax_mc_time.push_back(timer_mc.GetElapsedTime());

    //
    // Normal Computation
    //
    dax::cont::Timer<> timer_norm;
#if 1  // normal computation using Dax
    vtkNew<vtkDaxPolyDataNormals> normalGenerator;
	normalGenerator->SetInputConnection(marching->GetOutputPort());

	tlog->startEvent(evtNorm);
	normalGenerator->Update();
	tlog->endEvent(evtNorm);
#elif 1 // normal computation using VTK
    vtkNew<vtkPolyDataNormals> normalGenerator;
    normalGenerator->SetInputConnection(marching->GetOutputPort());
    normalGenerator->ComputePointNormalsOn();
    normalGenerator->ComputeCellNormalsOff();
    normalGenerator->Update();
#else // Skip normal computation
    marching->Update();
#endif
    SharedStatus::getInstance()->vtkdax_norm_time.push_back(timer_norm.GetElapsedTime());

    double time = timer.GetElapsedTime();
    tlog->endEvent(evt_id);

    std::cout << "VTKDax" << "," << device << "," << time << "," << i << std::endl;
    SharedStatus::getInstance()->vtkdax_total_time.push_back(time);

#ifdef DEBUG_NORMAL
    // debug
    vtkPolyData *polydata = normalGenerator->GetOutput();
    hasPointNormals = GetPointNormals(polydata);
    std::cout << "vtkDaxPolyDataNormals result has point normals? " << hasPointNormals << std::endl << std::endl;
#endif

    }
#endif
}

template <typename T> struct MarchingCubesOuputType
	{
  	typedef dax::CellTagTriangle type;
	};
// debugging the result
static void genVTKDaxMarchingCubes(vtkUnstructuredGrid* data)
{
#if 0
	//
	// convert vtk data to dax
	//
	//now set the buffer
	vtkDataArray *newData = data->GetPointData()->GetArray(0);
	dax::Scalar* rawBuffer = reinterpret_cast<dax::Scalar*>( newData->GetVoidPointer(0) );
	std::vector<dax::Scalar> buffer( newData->GetNumberOfTuples() );
	std::copy(rawBuffer, rawBuffer + newData->GetNumberOfTuples(), buffer.begin() );

	typedef vtkToDax::CellTypeToType<vtkHexahedron> VTKCellTypeStruct;
	typedef vtkToDax::DataSetTypeToType<VTKCellTypeStruct, vtkUnstructuredGrid> DataSetTypeToTypeStruct;
	typedef DataSetTypeToTypeStruct::DaxDataSetType InputDataSetType;
	InputDataSetType topology =
		  vtkToDax::dataSetConverter(data, DataSetTypeToTypeStruct() );


	//
	// run dax
	//
	dax::cont::Scheduler<> scheduler;
  	typedef dax::cont::GenerateInterpolatedCells<dax::worklet::MarchingCubesGenerate> GenerateIC;
  	typedef GenerateIC::ClassifyResultType  ClassifyResultType;
	dax::cont::ArrayHandle<dax::Scalar> field = dax::cont::make_ArrayHandle(buffer);

	//construct the two worklets that will be used to do the marching cubes
	dax::worklet::MarchingCubesClassify classifyWorklet(ISO_VALUE);
	dax::worklet::MarchingCubesGenerate generateWorklet(ISO_VALUE);
	dax::worklet::Normals normWorklet;

	//run the first step
	ClassifyResultType classification; //array handle for the first step classification
	scheduler.Invoke(classifyWorklet, topology, field, classification);

	//construct the topology generation worklet
	GenerateIC generate(classification,generateWorklet);
	generate.SetRemoveDuplicatePoints(true);

	//run the second step
    //dax::cont::UnstructuredGrid<dax::CellTagTriangle> outGrid;
    typedef typename MarchingCubesOuputType< typename VTKCellTypeStruct::DaxCellType >::type OutCellType;
	dax::cont::UnstructuredGrid<OutCellType,
                 vtkToDax::vtkTopologyContainerTag<VTKCellType>,
                 vtkToDax::vtkPointsContainerTag > outGrid;

	//schedule marching cubes worklet generate step, saving
	scheduler.Invoke(generate, topology, outGrid, field);

	//compute the normals of each output triangle

	dax::cont::ArrayHandle<dax::Vector3, vtkToDax::vtkArrayContainerTag<vtkFloatArray> > normals;
	scheduler.Invoke(normWorklet, outGrid, outGrid.GetPointCoordinates(), normals);


	//
	// convert to vtk
	//
	vtkNew<vtkPolyData> output ;

    daxToVtk::dataSetConverter(outGrid,output.GetPointer());

	//calling GetPortalControl brings memory back from execution to control env
	//so we should time this too.
	vtkDataArray *dataArray = normals.GetPortalControl().GetVtkData();
	dataArray->SetName("Normals");
	output->GetPointData()->AddArray(dataArray);

	vtkNew<vtkXMLPolyDataWriter> writer;
    writer->SetInputData(output.GetPointer());
    writer->SetFileName("out.vtp");
    writer->Write();

    cout << "polydata saved" << endl;
#endif
}
