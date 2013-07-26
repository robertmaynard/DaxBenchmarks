cmake_minimum_required(VERSION 2.8)

project(mccompare CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${mccompare_SOURCE_DIR})
option(ENABLE_PROFILING "Enable thread loggging for profiling" OFF)

add_definitions(-DNDEBUG)
if (ENABLE_PROFILING)
	add_definitions(-D_PROFILE)
endif ()

find_package(Dax REQUIRED)
find_package(VTK 6.1
  COMPONENTS
  vtkCommonCore
  vtkCommonDataModel
  vtkCommonExecutionModel
  vtkCommonMisc
  vtkFiltersCore
  vtkFiltersGeometry
  vtkIOLegacy
  vtkIOImage
  vtkImagingCore
  vtkAcceleratorsDax # Jimmy added 
  REQUIRED)

DaxConfigureSerial(REQUIRED)
#DaxConfigureTBB(REQUIRED)
DaxConfigureCuda()

include(${VTK_USE_FILE})

#
# tools
#
vtk_module_load (vtkIOXML)

add_executable(convertUnstructured 
	convertUnstructured.cxx)
target_link_libraries (convertUnstructured
  vtkCommonCore
  vtkCommonDataModel
  vtkCommonExecutionModel
  vtkCommonMisc
  vtkFiltersCore
  vtkFiltersGeometry
  vtkImagingCore
  vtkIOImage
  vtkIOLegacy
  vtkIOXML
)


#
# Benchmark
#


if (false) ################ SKIP ###############

set_source_files_properties(compare.h
  PROPERTIES HEADER_FILE_ONLY TRUE
  )

set(srcs
  compare.h
  compare_mc.h
  compare_thresh.h
  ArgumentsParser.cxx #Jimmy added this line and below:
  vtkDaxPolyDataNormals.h
  vtkDaxPolyDataNormalsImpl.h
  SharedStatus.cxx
  tlog/tlogDefs.cpp
  )

option(ENABLE_PISTON "Benchmark piston comparison" OFF)
if(${ENABLE_PISTON})

 find_path( PISTON_INCLUDE
    NAMES piston/piston_math.h
    DOC "Piston headers"
    )
 include_directories(${PISTON_INCLUDE})
 add_definitions("-DPISTON_ENABLED")

endif()


#-----------------------------------------------------------------------------
# setup timing libs we need to link too
set(DAX_TIMING_LIBS)
if(UNIX AND NOT APPLE)
  set(DAX_TIMING_LIBS rt)
endif()


if (false)
add_executable(BenchmarkSerial
  ${srcs}
  mainSerial.cxx
  vtkDaxPolyDataNormals.cxx
  )

set_target_properties(BenchmarkSerial PROPERTIES COMPILE_DEFINITIONS "DAX_DEVICE_ADAPTER=DAX_DEVICE_ADAPTER_SERIAL")

target_link_libraries(BenchmarkSerial
  vtkAcceleratorsDax # Jimmy added  
  vtkCommonCore
  vtkCommonDataModel
  vtkCommonExecutionModel
  vtkCommonMisc
  vtkFiltersCore
  vtkFiltersGeometry
  vtkImagingCore
  vtkIOImage
  vtkIOLegacy
  )
endif () #serial

if(Dax_TBB_FOUND)
  #Add TBB version
  add_executable(BenchmarkTBB
    ${srcs}
    mainTBB.cxx
    vtkDaxPolyDataNormals.cxx
    )

  set_target_properties(BenchmarkTBB PROPERTIES COMPILE_DEFINITIONS "DAX_DEVICE_ADAPTER=DAX_DEVICE_ADAPTER_TBB")


  target_link_libraries(BenchmarkTBB
    ${DAX_TIMING_LIBS}
    vtkCommonCore
    vtkCommonDataModel
    vtkCommonExecutionModel
    vtkCommonMisc
    vtkFiltersCore
    vtkFiltersGeometry
    vtkImagingCore
    vtkIOImage
    	vtkIOLegacy
    vtkAcceleratorsDax # Jimmy added 
    )
endif()

if(Dax_Cuda_FOUND)
  #Add CUDA version
  cuda_add_executable(BenchmarkCuda
    ${srcs}
    main.cu
    vtkDaxPolyDataNormals.cu
    )

  set_target_properties(BenchmarkCuda PROPERTIES COMPILE_DEFINITIONS "DAX_DEVICE_ADAPTER=DAX_DEVICE_ADAPTER_CUDA")


  target_link_libraries(BenchmarkCuda
    vtkCommonCore
    vtkCommonDataModel
    vtkCommonExecutionModel
    vtkCommonMisc
    vtkFiltersCore
    vtkFiltersGeometry
    vtkImagingCore
    vtkIOImage
    vtkIOLegacy
    vtkAcceleratorsDax # Jimmy added 
    )
endif()



endif () #-------------------
