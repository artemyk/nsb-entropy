/*
 *  File - EntrData.h
 *  Project - Entropy_parsing
 *
 *  Declaration of the main data class (EntrData) interface and the
 *  dictionary class (EntrDict) interface, as well as other main data
 *  types and constants.
 *
 *  Copyright (c) 2001-2013 Ilya Nemenman. All rights reserved.
 *
 */

#ifndef ENTRDATA_H
#define ENTRDATA_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <list>
#include "EntrData_except.h"
#include "counts.h"

namespace ENTRDATA{
  enum ENTRDATA_VER_CALL {EDVERNAME, EDVERVER, EDVERDATE, EDVERAUTH, EDVERLIC, EDVERALL}; 


  inline std::string ENTRDATA_VER(ENTRDATA_VER_CALL type){
    switch (type){
    case EDVERNAME:
      return std::string("NSB entropy estimator collection");
      break;
    case EDVERVER:
      return std::string("Ver. 1.14");
      break;
    case EDVERDATE:
      return std::string("09/04/2013");
      break;
    case EDVERAUTH:
      return std::string("2002-2013 by Ilya Nemenman (ilya@menem.com)");
      break;
    case EDVERLIC:
      return std::string("This software is distributed under GPL 3.0") +
	std::string(" and comes with ABSOLUTELY NO WARRANTY.\n") +
	std::string("This is free software, and you are welcome to ") +
	std::string("redistribute it under conditions specified in the License.");
    default:			// EDVERALL
      break;
    }
    return ENTRDATA_VER(EDVERNAME) + ", " + ENTRDATA_VER(EDVERVER) + " (" 
      + ENTRDATA_VER(EDVERDATE) + ")" + ".\n" + "Copyright (C) " 
      + ENTRDATA_VER(EDVERAUTH) + ".\n" +ENTRDATA_VER(EDVERLIC);
  }

  // verbosity
  enum VERBOSITY{VERBMIN, VERBMED, VERBMAX};

  class EntrData;
  class EntrDict;

  typedef std::list<int> lint;
  typedef lint::const_iterator LCI;  
  inline std::ostream& operator<<(std::ostream& out, lint& word){
    for(LCI p=word.begin();p!=word.end();++p)
      out<<(*p)<<" ";
    return out;
  }

  // the following vector is used to represent the occurences data
  // for different possible subsamples of a full sample, depending
  // on the phse (starting point)
  typedef std::vector<EntrData*> VEDP;
  
  // different possible implementations of EntrData
  enum IMPLEDATA {EDATAMPI, EDATANODE};
  
  // different possible ways to write the data in
  enum WRITELEVEL {WRITEONE, WRITEALL};

  // for shift != 1, do we need to get counts for every possible
  // start position, or just the first one?
  // I also use values of startpoint >1 to indicate how many
  // starts are required
  typedef int STARTPOINT;
  const STARTPOINT MANYSTARTS=0;
  const STARTPOINT ONESTART=1;

  // do or not do nsb calculations
  enum DOCALCS {YESNSB, YESCOUNT, NSBNONE};
  // partiioning data into regions based on sampling before doing NSB
  // calculations
  enum DOPART {YESPART, NOPART};

  // Map connecting the depth of counting with the pointer to the
  // appropriate "counts" structure
  typedef std::map<int, counts*> mpcounts;
  typedef mpcounts::iterator mPCI;
  typedef mpcounts::const_iterator mPCCI;


  /* 
   *  Main data storage class. This is a virtual class declaring
   *  what any implementation of the EntrData children must
   *  provide. At present, there are two working implementations for
   *  the EntrData -- one done thrugh "mpi_map" class, and the other
   *  through "node" class. "Mpi_map" keeps different maps for each
   *  word-length recorded into it. This is not optimal since it
   *  does not take into the account the fact that short words not
   *  occuring means that some long words will not occur either. The
   *  "node" implementation takes this into account. However,
   *  overhead for storing nodes is far greater and outweights
   *  achieved memory economy. Thus at the moment, "node" is not
   *  developed any more as it is a lot more memory consuming.
   */    
  class EntrData{
  private:
    int         asize;		// alphabet size
    std::string name;		// name of the data structure

    virtual void Print(std::ostream&) const =0; // virtual print function

  protected:
    mutable mpcounts     mcts;	// map of counts
    void ClearCounts() const;	// removing all the counts (to be, e.g.,
				// used when new data is added)
    virtual counts* MakeCounts(int) const throw(EntrData_bad_alloc) =0; // recording data into 
				// "counts" class;
				// int -- the length of words to be traced; 
				// if there are no words of this length recorded,
				// empty counts are returned    

  public:
    EntrData(const int sz, const std::string& newname = "") 
      throw (EntrData_bad_alloc): asize(sz), name(newname) {}
    virtual ~EntrData(){ClearCounts();}
    int	ASize() const {return asize;}
    virtual void Record(const lint&) =0; // recording a list into the data structure at 
				// the depth of the list 
    virtual void RecordAll(const lint&)=0; // recording at all depth up to the 
				// length of the list
    const char* GetName() const {return name.c_str();} // returning the name of the data
				// structure
    const counts* count(int) const throw(EntrData_bad_alloc); // returning
				// counts produced by MakeCounts
    friend std::ostream& operator<<(std::ostream&, const EntrData&);	// outputting the data
  };	
    
  // Operator form for printing.
  inline std::ostream& operator<<(std::ostream& out, const EntrData& ed){
    ed.Print(out);
    return out;
  }
  

  
  
  /* 
   *  Virtual dictionary class. All other types of dictionaries have
   *  to conform to this prototype. The dictionary has to be able to
   *  initialize itself, accept input from an input stream, and
   *  return an outputa character, tell if it's ready to
   *  output a vector of translated characters, and finally output
   *  it. Note that this vector is a same-time value of some
   *  multi-dimensional signal, but NOT successive time values of
   *  the same one-dimensional signal.
   */
  
  class EntrDict{
  private:
    int asize;			// the alphabet size 
    int def;			// see Def()
    std::string name;		// name of the alphabet
  protected:			
    // children have access to the following functions
    void SetASize(int sz){asize=sz;}
    void SetName(const std::string& na){name=na;}
    void SetDef(int d){def=d;}

  public:
    EntrDict(const int as, const int d, const std::string& newname="") throw() : 
      asize(as), def(d), name(newname){} // constructor 
    virtual ~EntrDict(){}	// destructor 
    
    const char* GetName() const{ // description of the
      return name.c_str();	// dictionary
    };				// which exception decratation here?

    virtual int Translate(std::istream&) const throw (EntrData_failed_transl) =0; 
				// translating the next symbols in the stream
    void Start(std::istream& in, unsigned long n=0) 
      throw (EntrData_failed_transl){ // preparing the stream
				//and then skipping a few data from it
      Prepare(in);
      Skip(in,n);
    }
    virtual void Prepare(std::istream&) const throw(EntrData_failed_transl)=0;
				// preparing the stream

    void Skip(std::istream& in, unsigned long n=1)
      const throw (EntrData_failed_transl){ // skipping
				// symbols from the stream
      for(unsigned long i =0; i<n; i++){
	Translate(in);
      }
    }
    int ASize() const throw() {return asize;}
    int operator[](std::istream& in) const throw (EntrData_failed_transl){
				// translation as an operator
      return Translate(in);
    }

    // the default value returned by the dictionary if reading is
    // successfull, but something wrong is read
    int Def() const throw(){return def;}


  };
}

#endif
