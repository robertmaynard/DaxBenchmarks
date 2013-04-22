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

#define DAX_DEVICE_ADAPTER DAX_DEVICE_ADAPTER_SERIAL

//I can't get piston to run serially, so make sure we don't compile it.
// The primary issue is that piston uses device space tags, and thrust
// cpp backend is labeled with a host space tag
#undef PISTON_ENABLED

#include "ArgumentsParser.h"
#include "compare.h"


int main(int argc, char* argv[])
  {
  dax::testing::ArgumentsParser parser;
  if (!parser.parseArguments(argc, argv))
    {
    return 1;
    }

  const std::string file = parser.file();
  const int pipeline = parser.pipeline();
  RunComparison("Serial", file, pipeline);
  return 0;
}
