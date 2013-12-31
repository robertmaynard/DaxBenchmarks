## Data Analysis at Extreme ##

The Dax Toolkit is an open source C++ header only library that provides a collection of data analysis and visualization algorithms that run well on multi- and many-core processors.  The Dax Toolkit also makes it easy to design other visualization algorithms that work well on such processors. We currently support CUDA, OpenMP, Intel TBB, and Serial Execution.

### How ###

The Dax Toolkit uses fine-grained concurrency for data analysis and visualization algorithms.
The basic computational unit of the Dax Toolkit is a worklet, a function that implements the algorithm’s behavior on an element of a mesh (that is, a point, edge, face, or cell) or a small local neighborhood.

The worklet is constrained to be serial and stateless; it can access only the element passed to and from the invocation. With this constraint, the serial worklet function can be concurrently scheduled on an unlimited number of threads without the complications of threads or race conditions.

Although worklets are not allowed communication, many visualization algorithms require operations such as variable array packing and coincident topology resolution that intrinsically require significant coordination among threads. Dax enables such algorithms by classifying and implementing the most common and versatile communicative operations into worklet types which are managed by the Dax dispatcher.

## Getting Dax ##


The Dax repository is located at [https://github.com/kitware/DaxToolkit](https://github.com/kitware/DaxToolkit)

Dax dependencies are:


+  [CMake 2.8.8](http://cmake.org/cmake/resources/software.html)
+  [Boost 1.49.0](http://www.boost.org) or greater
+  [Cuda Toolkit 4+](https://developer.nvidia.com/cuda-toolkit) or [Thrust 1.4 / 1.5](https://thrust.github.com)
   depending if you want Cuda and/or OpenMP

```
git clone git://github.com/Kitware/DaxToolkit.git dax
mkdir dax-build
cd dax-build
cmake-gui ../dax
```

A detailed walk-through of installing and building Dax can be found on our [Install page](http://www.daxtoolkit.org/index.php/Building_the_Dax_Toolkit)

A walk-through of integrating Dax into your project can be found on our
[Integration page] (http://www.daxtoolkit.org/index.php/How_To_Use_Dax_In_Your_Project)


## How To Use Bencmarks ##

Each program has two consitent arguments, file and pipeline. File points
to a nrrd file, and pipeline can be 1 for threshold or 2 for marching cubes

For tbb you can also set cores to specify the number of threads to use.

```
./BenchmarkTBB --file=./data.nhdr --pipeline=2 --cores=64

```




## License ##

Copyright (c) Kitware, Inc.
All rights reserved.
[See LICENSE.txt for details](LICENSE.txt).

Copyright 2011 Sandia Corporation.
Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
the U.S. Government retains certain rights in this software.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
[See the copyright file for more information](LICENSE.txt).

