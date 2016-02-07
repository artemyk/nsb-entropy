/*
 *  File - recorder.h
 *  Project - Entropy_parsing
 *
 *  Declares the recorder class, that records a translated word in 
 *  the appropriate data structures.
 *
 *  Copyright (c) 2002 Ilya Nemenman. All rights reserved.
 *
 */

#ifndef RECORDER_H
#define RECORDER_H

#include "EntrData.h"
#include "EntrData_except.h"
#include "../include/random.h"



namespace ENTRDATA{
  
  class Recorder{
  private:
    int depth;			// requested recording depth
    WRITELEVEL all;             // do we record at all depth, or at
				// "depth" only
    virtual unsigned int Shorten() =0; // decreasing the stored data
				// words
  protected:
    int maxrecs;	        // number of possible recordings in
				// the file being read
 
  public:
    Recorder(int d, WRITELEVEL wl=WRITEONE,int m=0) :  depth(d), all(wl),maxrecs(m) {} 
    virtual ~Recorder() {}
    virtual unsigned int WLength() const =0; // length of data words
				// curretly stored in the recorder
    virtual int Ready() const=0; // check if recorder is ready to record
    virtual int Full() const=0;	// recorder cannot accept any more data
    virtual void UpdateWord(int) =0; // adding a new symbol to the word
    virtual void Record(VEDP*,STARTPOINT,double) =0; // recording the word to the appropriate 
				// elements of VEDP 
    virtual int FinishRecord(VEDP*,STARTPOINT, double); // finishing up recording
   
    virtual std::ostream& print(std::ostream&)=0; // outputting the
				// data
    friend std::ostream& operator<<(std::ostream&, Recorder&); // outputting the data
    
    unsigned int get_depth()const{return (unsigned int) depth;}
    WRITELEVEL get_all()const{return all;}
    long to_record()const{return maxrecs;} // how many recordings should
				// we prepare to read?
  };	
    
  // Operator form for printing.
  inline std::ostream& operator<<(std::ostream& out, Recorder& rec){
    rec.print(out);
    return out;
  }


}

#endif
