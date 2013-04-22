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

#include "ArgumentsParser.h"

#include <dax/testing/OptionParser.h>
#include <iostream>
#include <sstream>
#include <string>

enum  optionIndex { UNKNOWN, HELP, FILEPATH, PIPELINE, CORES};
const dax::testing::option::Descriptor usage[] =
{
  {UNKNOWN,   0,"" , ""    ,      dax::testing::option::Arg::None, "USAGE: example [options]\n\n"
                                                                    "Options:" },
  {HELP,      0,"h" , "help",    dax::testing::option::Arg::None, "  --help, -h  \tPrint usage and exit." },
  {FILEPATH,      0,"", "file",      dax::testing::option::Arg::Optional, "  --file  \t nrrd file to read." },
  {PIPELINE,  0,"", "pipeline",  dax::testing::option::Arg::Optional, "  --pipeline  \t What pipeline to run ( 1 threshold, 2 marching cubes)." },
  {CORES,  0,"", "cores",        dax::testing::option::Arg::Optional, "  --cores  \t number of cores to use, only valid for TBB." },
  {UNKNOWN,   0,"",  "",         dax::testing::option::Arg::None, "\nExample:\n"
                                                                   " example --file=./test --pipeline=1\n"},
  {0,0,0,0,0,0}
};


//-----------------------------------------------------------------------------
dax::testing::ArgumentsParser::ArgumentsParser():
  File(""),
  Pipeline(THRESHOLD),
  Cores(-1)
{
}

//-----------------------------------------------------------------------------
dax::testing::ArgumentsParser::~ArgumentsParser()
{
}

//-----------------------------------------------------------------------------
bool dax::testing::ArgumentsParser::parseArguments(int argc, char* argv[])
{

  argc-=(argc>0);
  argv+=(argc>0); // skip program name argv[0] if present

  dax::testing::option::Stats  stats(usage, argc, argv);
  dax::testing::option::Option* options = new dax::testing::option::Option[stats.options_max];
  dax::testing::option::Option* buffer = new dax::testing::option::Option[stats.options_max];
  dax::testing::option::Parser parse(usage, argc, argv, options, buffer);

  if (parse.error())
    {
    delete[] options;
    delete[] buffer;
    return false;
    }

  if (options[HELP] || argc == 0)
    {
    dax::testing::option::printUsage(std::cout, usage);
    delete[] options;
    delete[] buffer;

    return false;
    }

  if ( options[FILEPATH] )
    {
    std::string sarg(options[FILEPATH].last()->arg);
    std::stringstream argstream(sarg);
    argstream >> this->File;
    }


  if ( options[CORES] )
    {
    std::string sarg(options[CORES].last()->arg);
    std::stringstream argstream(sarg);
    argstream >> this->Cores;
    }

  if ( options[PIPELINE] )
    {
    std::string sarg(options[PIPELINE].last()->arg);
    std::stringstream argstream(sarg);
    int pipelineflag = 0;
    argstream >> pipelineflag;
    if (pipelineflag == THRESHOLD)
      {
      this->Pipeline = THRESHOLD;
      }
    else if (pipelineflag == MARCHING_CUBES)
      {
      this->Pipeline = MARCHING_CUBES;
      }
    else
      {
      std::cerr << "Incorrect pipeline choice: " << pipelineflag << std::endl;
      std::cerr << "Threshold is : " << THRESHOLD  << std::endl;
      std::cerr << "Marching Cubes is : " << MARCHING_CUBES  << std::endl;
      }
    }

  delete[] options;
  delete[] buffer;
  return true;
}
