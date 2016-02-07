/*
 *  File - simple_recorder.h
 *  Project - Entropy_parsing
 *
 *  Declares the parallel_recorder class, that records a translated
 *  words from a set of parallel stream into
 *  the appropriate data structures in a simple successive form allowing
 *  for different overlaps between successive words.
 *  Importantly, equal time slices through stream are *horizontal*.
 *
 *  Copyright (c) 2004 Ilya Nemenman. All rights reserved.
 *
 */

#ifndef PARALLEL_RECORDER_H
#define PARALLEL_RECORDER_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "recorder.h"
#include "../include/random.h"

/* 
 * For parallel recording we record ALL slices. The probability to skip
 * or the recording everyth shift's number is for recording within the slice.
 * In such cases we skip some data in a slice, but still do all slices.
 *
 */

namespace ENTRDATA{
  typedef std::vector<lint> vlint; // vector of word-lists, to be used
				// for recording parallel words

  class Parallel_Recorder : public Recorder{
  private:
    vlint words;		// the words currently being recorded 
    int next_word;		// the word into which the next translated
				// symbol will be recorded; if this is
				// zero, all words must have the same
				// length 

    unsigned int Shorten();	// decresing the length of stored data words
    void RecordSingle(EntrData*, int); // recording the word with a
				// given number in one particular data
				// structure
  public:
    Parallel_Recorder(int,const std::string&,WRITELEVEL) throw(EntrData_failed_transl);
    ~Parallel_Recorder() {}


    void Record(VEDP*, STARTPOINT, double); // recording the word to
				// the appropriate elements of VEDP
    unsigned int WLength() const{ // current length of the translation word
      return words[0].size();	// all lists are to be of the same length
    }
    int Ready() const{ 		// recorder is ready when each word is
				// of proper length, and the next
				// word to record is 0
      return ((next_word==0)&&(WLength()==get_depth()));
    }
    int Full() const{		// this recorder is full when it's
				// ready (one slice through data is done)
      return Ready();
    }
    void UpdateWord(int);	// adding a new symbol to the word
    int FinishRecord(VEDP*, STARTPOINT, double); // redefine finishing

    std::ostream& print(std::ostream& out){ // printing
      for(int i=0;(unsigned int)i<words.size();i++)
	out << words[i];
      return out;
    }
  };



}

#endif
