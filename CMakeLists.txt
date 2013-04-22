cmake_minimum_required(VERSION 2.8)

project(mccompare CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${mccompare_SOURCE_DIR})

add_definitions(-DNDEBUG)

find_package(Dax REQUIRED)
find_package(VTK 6.0
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
  REQUIRED)

DaxConfigureSerial(REQUIRED)
DaxConfigureOpenMP(REQUIRED)
DaxConfigureTBB(REQUIRED)
DaxConfigureCuda(REQUIRED)

include(${VTK_USE_FILE})

set_source_files_properties(compare.h
  PROPERTIES HEADER_FILE_ONLY TRUE
  )

set(srcs
  compare.h
  compare_mc.h
  compare_thresh.h
  ArgumentsParser.cxx
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


add_executable(BenchmarkSerial
  ${srcs}
  mainSerial.cxx
  )

target_link_libraries(BenchmarkSerial
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

#Add TBB version
add_executable(BenchmarkTBB
  ${srcs}
  mainTBB.cxx
  )

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
  )

#Add CUDA version
cuda_add_executable(BenchmarkCuda
  ${srcs}
  main.cu
  )

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
  )
