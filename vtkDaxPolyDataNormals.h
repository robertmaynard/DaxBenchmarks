//=============================================================================
//
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2012 Sandia Corporation.
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//
//=============================================================================

#ifndef vtkDaxPolyDataNormal_H
#define vtkDaxPolyDataNormal_H

#include <vtkPolyDataNormals.h>
#include <dax/cont/Timer.h>
#include "vtkAcceleratorsDaxModule.h" //required for correct implementation



class VTKACCELERATORSDAX_EXPORT vtkDaxPolyDataNormals : public vtkPolyDataNormals
{
public:
  vtkTypeMacro(vtkDaxPolyDataNormals,vtkPolyDataNormals)
  //void PrintSelf(ostream& os, vtkIndent indent);
  static vtkDaxPolyDataNormals* New();

protected:
  vtkDaxPolyDataNormals();
  ~vtkDaxPolyDataNormals();


  // convenience method
  virtual int RequestInformation(vtkInformation* request,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector);

  // Description:
  // This is called by the superclass.
  // This is the method you should override.
  virtual int RequestUpdateExtent(vtkInformation*,
                                  vtkInformationVector**,
                                  vtkInformationVector*);

  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  vtkDaxPolyDataNormals(const vtkDaxPolyDataNormals&); //Not implemented
  void operator=(const vtkDaxPolyDataNormals&); // Not implemented
  dax::cont::Timer<> timer;
};

#endif // vtkDaxPolyDataNormal_H
