/*
 *  File - mpi_map.h
 *  Project - Entropy_parsing
 * 
 *  Declaration of the multiple precision integer map class. This is
 *  to be used as a realization of a EntrDatas -- main data class.
 *
 *  Copyright (c) 2001,2002 Ilya Nemenman. All rights reserved.
 *
 */

#ifndef MPI_MAP_H
#define MPI_MAP_H

#include <map>
#include <vector>
#include <math.h>
//#include <gmpxx.h>

#include "../include/utility.h"
#include "EntrData.h"
#include "mpi.h"

namespace ENTRDATA{
  using namespace GMP;
  
  typedef std::map<mpi,long> DATAMAP;
  //    typedef map<mpz_class,long> DATAMAP;
  typedef DATAMAP::const_iterator DMCI;
  typedef std::vector<DATAMAP*> vDATAMAP;
  typedef vDATAMAP::const_iterator vDMCI;
  typedef vDATAMAP::iterator vDMI;


  class mpi_map: public EntrData{
  private:
    vDATAMAP	data;	        // vector of maps, each counting words of
				// different length 
    void Print(std::ostream&) const;	// printing the data structures
    void RecordList(const lint&, int) 
      throw(EntrData_range_error, EntrData_bad_alloc);	
				// recording list -- int=0 -- at the
				// length of the list; 
				// int!=0 -- at all length up to the
				// length of the list 
    counts* MakeCounts(int) const throw(EntrData_bad_alloc); 
				// recording data into "counts" class 

  public:
    mpi_map(const int sz =2, const std::string newname = "", 
	    const int length=1) throw(EntrData_bad_alloc);
				// construction specifying the
				// alphabet size, and the 
				// initial size of the map
    ~mpi_map();
    
    void RecordAll(const lint& v) throw(EntrData_range_error, EntrData_bad_alloc)	
    {RecordList(v,1);}	        // recording list in the map at all depth
    void Record(const lint& v) throw(EntrData_range_error, EntrData_bad_alloc)
    {RecordList(v,0);}		// recording list in the map only at
				// the depth of its length

  };	
}



#endif
