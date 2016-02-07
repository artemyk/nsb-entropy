/*
 *  File - mpi_map.cpp
 *  Project - Entropy_parsing
 *
 *  Defining the mpi_map class, one of possible implementations for
 *  EntrData class.
 *
 *  Copyright (c) 2001, 2002 Ilya Nemenman. 
 *
 */


#include "mpi_map.h"

/* 
   We initialize the vectors of maps to its initially	
   specified size
 */
ENTRDATA::mpi_map::mpi_map(const int sz, const std::string newname, const int length)
  throw(EntrData_bad_alloc): EntrData(sz, newname), data(length){
  try{
    for(vDMI p=data.begin(); p!=data.end(); ++p)
      *p = new DATAMAP();
  }
  catch(std::bad_alloc){
    throw EntrData_bad_alloc();
  }
}

/*
 * Destructor has to delete all of the maps stored in data
 */
ENTRDATA::mpi_map::~mpi_map(){
  for(vDMCI p=data.begin(); p!=data.end();p++)
    delete *p;			// deleting data objects
}


/*
 * Recording a list into the data structure.  We first make the key
 * from the list utilizing the fact that the alphabet size is ASize(),
 * and then increment "(data[i])[key]"
 *
 * Note that if we are recording, for example, list of length 2, we
 * keep 3 maps: for empty list, length 1, and length 2
 *
 * We keep all the maps up to the depth of the longest recorded list,
 * however not all of them will be filled.  Correspondingly, counts at
 * non-filled depths will be zero.
 */
void ENTRDATA::mpi_map::RecordList(const lint& v, int what) 
  throw(EntrData_range_error, EntrData_bad_alloc){

  // first clear all counts since they are not valid any more
  ClearCounts();

  mpi key((unsigned long)0);
  mpi t((unsigned long)1);
  //mpz_class key(0);
  //mpz_class t(0);

  if(data.size()<(v.size()+1)){	// need to resize the "data" list
    unsigned int p=data.size();
    data.resize(v.size()+1);	// resizing the vector of pointers to maps
    try{
      for(; p<data.size(); ++p){
	data[p] = new DATAMAP(); // adding maps at newly created positions
      }
    } 
    catch(std::bad_alloc){
      throw EntrData_bad_alloc();	
    }
  }
        
  // recording into "data"
  if(what || (v.size()==0)) (*data[0])[key]++; // recording the empty
					       // list if recording at all levels,
					       // or at zero's depth
  {
    unsigned int j =0; 
    for(LCI i = v.begin() ; i!=v.end(); ++i){
      // if bad range for the next list element
      if((*i<0)||(*i>=ASize())) 
	throw EntrData_range_error("Index" +num2str(*i) +
				   " exceeds alphabet size " +num2str(ASize()) + ".\n");
      //key += t* (*i);
      //t *= ASize();
      key.addmul_ui(t, *i);
      t.mul_ui(ASize());
      if(what || (v.size() == j+1)){ // record list if recording at all levels or 
				// at approptiate depth
	(*data[j+1])[key]++;	// empty list goes in 0, list of length 1 goes in map #1, etc
      }
      j++;
    }
  }
}    //    out << "ASize = " << ed.ASize() << "\n";


/*
 * printing the map data
 */
void ENTRDATA::mpi_map::Print(std::ostream& out) const {
  out<<"# Created by "<< ENTRDATA_VER(EDVERNAME)<<", " <<  ENTRDATA_VER(EDVERVER)<<"\n";
  out<<"# name: ASize\n";
  out<<"# type: scalar\n";
  out<<ASize()<<"\n";
  out<<"# name: struct_name\n";
  out<<"# type: string\n";
  out<<"# elements: 1\n";
  out<<"# length: "<<strlen(GetName())<<"\n";
  out<<GetName()<<"\n";
  out<<"# name: len\n";
  out<<"# type: matrix\n";
  out<<"# rows: 1\n";
  out<<"# columns: " << data.size()<<"\n";
  for(unsigned long u=0;u<data.size();u++)
    out<<" "<<u;
  out<<"\n";
  std::string str;
  for(unsigned long u=0; u<data.size();u++){
    unsigned int sz=(*data[u]).size();
    out<<"# name: counts"<<u<<"\n";
    out<<"# type: matrix\n";
    out<<"# rows: 1\n";
    out<<"# columns: "<<sz<<"\n";
    for(DMCI p=(*data[u]).begin(); p!=(*data[u]).end(); ++p)
      out<<" "<<p->second;
    out<<"\n";
    out<<"# name: words"<<u<<"\n";
    out<<"# type: string\n";
    out<<"# elements: "<<(*data[u]).size()<<"\n";
    for(DMCI p=(*data[u]).begin(); p!=(*data[u]).end(); ++p){
      out<<"# length: "<< u<<"\n";
      str = (p->first).to_string(ASize());
      for(unsigned int ii=0;ii<(u-str.size());ii++)
	out<<"0";
      out<<str<<"\n";
    }

    //    out<< "Data structure " << GetName() <<":\n";
    //    out<< "   Word length = " << u << ".\n";
    //for(DMCI p=(*data[u]).begin(); p!=(*data[u]).end(); ++p)
    //  out<<"\t"<<p->first<<"\t"<<p->second<<"\n";
  }
}

/*
 * The function creates the "counts" object and counts events
 * in the map_mpi into it.
 */ 
ENTRDATA::counts* ENTRDATA::mpi_map::MakeCounts(int depth) const throw(EntrData_bad_alloc){
  counts* c;
  try{
    c= new counts(pow((double) ASize(), (double) depth));
  }
  catch (std::bad_alloc){
    throw EntrData_bad_alloc();
  };

  if(depth<= ((int)data.size()-1)){ // something is recorded only at this depth
    for(DMCI p = (*data[depth]).begin(); p != (*data[depth]).end(); ++p) {
				// do the actual counting and recording
      (*c)(p->second)+=1.0;
    }
  }

  c->NEvents();			// to force recalculation of counts[0]

  mcts[depth]=c;
  return c;
}

    
