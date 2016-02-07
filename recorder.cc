/*
 *  File - recorder.cc
 *  Project - Entropy_parsing
 *
 *  Definitions of the functions from the recorder class.
 *
 *  Copyright (c) 2002 Ilya Nemenman. All rights reserved.
 *
 */
#include "recorder.h"
#include "EntrData_except.h"






/*
 * Finishing up recording. If WRITELEVEL all, then 
 * we need to record the stuff which is in the word at
 * shorter depths than the maximum one.
 * Return the number of times "Record" has been called
 */
int ENTRDATA::Recorder::FinishRecord(VEDP* data, STARTPOINT st, double prob){
  int i=0;
  if(all){			// only then need recording clean-up
    while (Shorten() >0){	// shorten
      Record(data, st, prob);	// record shortened version
      i++;			// make a not of # of recordings
    }
  }
  return i;
}
