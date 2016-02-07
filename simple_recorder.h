/*
 *  File - simple_recorder.h
 *  Project - Entropy_parsing
 *
 *  Declares the simple_recorder class, that records a translated word in 
 *  the appropriate data structures in a simple successive form allowing
 *  for different overlaps between successive words.
 *
 *  Copyright (c) 2002 Ilya Nemenman. All rights reserved.
 *
 */

#ifndef SIMPLE_RECORDER_H
#define SIMPLE_RECORDER_H

#include "recorder.h"
#include "../include/random.h"

namespace ENTRDATA{

  class Simple_Recorder : public Recorder{
  private:
    lint word;			// the word currently being recorded 
    unsigned int to_rec;	// number of the data structure, in which the data
				// will be recorded the next time

    unsigned int Shorten(){	// decreasing the length of stored data words
      if (WLength()>0)  word.pop_front();
      return WLength();
    }
    void RecordSingle(EntrData*);  // recording the word in one particular 
                                   // data structure

  public:
    Simple_Recorder(int d, WRITELEVEL wl=WRITEONE) : 
    Recorder(d,wl,1), word(0), to_rec(0) {} 
    ~Simple_Recorder() {}

    void Record(VEDP*, STARTPOINT, double); // recording the word to the appropriate 
                        // elements of VEDP 
    unsigned int WLength() const{ // current length of the translation word
      return word.size();
    }
    int Ready()const{return (WLength()==get_depth());} // can record if have enough data
    int Full()const{return 0;}	// this recorder is never full
    void UpdateWord(int);	// adding a new symbol to the word
    std::ostream& print(std::ostream& out){ // printing
      out << word;
      return out;
    }
  };



}

#endif
