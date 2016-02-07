/*
 *  File - parsers.h 
 *  Project - Entropy_parsing
 *
 *  Declares some functions used to parse data files into data
 *  structures and report errors.
 *
 *  Copyright (c) 2001, 2002 Ilya Nemenman. 
 */

#ifndef ENTR_PARSER_H
#define ENTR_PARSER_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "EntrData.h"
#include "mpi_map.h"
#include "node.h"
#include "simple_dicts.h"
#include "simple_recorder.h"
#include "parallel_recorder.h"
#include "nsb.h"

extern ENTRDATA::VERBOSITY vlevel;

namespace ENTRDATA{

  void error(const char*, const char* ="");

  void Stream2Dict2Data(VEDP&, std::istream&, const EntrDict*, Recorder*, 
			unsigned long, unsigned long, int,
			STARTPOINT, double,
			unsigned long&, unsigned long&);
  void Data2Cts2File(const VEDP&, const std::string&, const std::string&, int, 
		     int, STARTPOINT, DOCALCS, DOPART, const std::string&, VNsbSet&);
  void ProcessingCycle(const std::string&, unsigned long, unsigned long, 
		       IMPLEDATA, int, int, STARTPOINT, DOCALCS, DOPART,
		       double, int, std::vector<std::string>*, 
		       const std::string&,EntrDict*, Recorder*, VNsbSet&);
  void ProcessFile(const std::string&, unsigned long, unsigned long, 
		   IMPLEDATA, const std::string&, const std::string&, int, 
		   int, WRITELEVEL, int,
		   STARTPOINT, DOCALCS, DOPART, double);
  void ChooseDictRec(const std::string&, const std::string&, int, WRITELEVEL,
		     EntrDict*&, Recorder*&);
  void KillDictRec(EntrDict*, Recorder*);
  void ChooseData(const IMPLEDATA, 
		  const std::vector<std::string>&, VEDP&, int);
  void KillData(VEDP&);
}

#endif
