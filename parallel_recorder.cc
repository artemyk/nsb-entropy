/*
 *  File - simple_recorder.cc
 *  Project - Entropy_parsing
 *
 *  Definition of functions for the simple recorder class.
 *
 *  Copyright (c) 2002 Ilya Nemenman. All rights reserved.
 *
 */

#include "parallel_recorder.h"

/*
  getting the number of columns (streams) from a data file;
  preparing recording to depth d; number of repetitions to skip is "skip"
*/
ENTRDATA::Parallel_Recorder::Parallel_Recorder(int d, const std::string& fname, WRITELEVEL wl) 
  throw(EntrData_failed_transl): Recorder(d, wl),words(0),next_word(0){
  //opening file    
  std::string infname = fname + ".txt";
  std::ifstream in(infname.c_str());
  if(!in) throw EntrData_failed_transl("Cannot open input file "+ infname);	


  std::string str;		// temorary input string
  do{
    in>>str;
  }while((str!="rows:")&&(in));
  in>>maxrecs;
  if (!in)
    throw EntrData_failed_transl("Could not read the number of rows (different slices).");
  do{
    in>>str;
  }while((str!="columns:")&&(in));
  // getting the alphabet size
  int ncols;
  in>>ncols;
  if (!in)
    throw EntrData_failed_transl("Could not read the number of columns (streams).");
  words.resize(ncols);			// resizing to correct # of
					// streams
}



/*
 * Recording the word in the appropriate elements of
 * data structures. The number of the data structures
 * is equal to the displacement between successive words 
 * that have to be recorded into the same data structure,
 * and we record the words circulating between the structures.
 * This way one pass through the file allows to record n different 
 * sets of every n-th word. Each set has a different phase
 * (starting point).
 */
void ENTRDATA::Parallel_Recorder::Record(VEDP* data, 
				       STARTPOINT st, double prob){
  // create a random number generator needed to decide whether to record 
  // a word or not
  static rnd_gen_class rnd_gen(gsl_rng_mt19937);  
  
  // We have a set of data structures recorded in "data". Further we
  // have a vector of words to be recorded. We need to check which of
  // the words goes into which of the structures. So we set up a
  // counter for which EntrData element is to be recorded in next, and
  // for each recorded word, we do:
  int to_rec =0;		// next EntrData to be recorded into
  for(int i=0; (unsigned int)i<words.size();i++){
    // if we do all possible recordings (MANYSTEPS) or if we are 
    // now at the recording #0, do the recording
    if ((st == MANYSTARTS) || (to_rec==0))
      if (rnd_gen.uniform()<=prob) // accept datum with some probability
	RecordSingle((*data)[to_rec], i );
    to_rec = ( (unsigned int)(to_rec+1) == (data->size())) ? 0 : to_rec+1;
  }
}


/*
 * Updating the translation words
 */
void ENTRDATA::Parallel_Recorder::UpdateWord(int newint){
  // add the newly translated symbol to the correct word
  words[next_word].push_back(newint);
  if(words[next_word].size() > get_depth()) // if the word reached needed
    words[next_word].pop_front();		// depth, pop old int from front
  
  // increment next_word counter, or zero it, if the end reached
  next_word = ( (unsigned int)(next_word+1) == words.size()) ? 0 : next_word+1;
}


/*
 * Recording the word in the single data structure
 */
void ENTRDATA::Parallel_Recorder::RecordSingle(EntrData* ed, int num){
  if(get_all())
    ed->RecordAll(words[num]);	//recording at all levels
  else
    ed->Record(words[num]);	//recording at one level
}


/*
  shortening all words 
*/
unsigned int ENTRDATA::Parallel_Recorder::Shorten(){
  if (WLength()>0){
    if(next_word){		// if next word is not zeroth, need
				// to even the words length
      for(;next_word>0;next_word--)
	words[next_word-1].pop_back();
    }
    for(int i=0; i <(int)words.size();i++)
      words[i].pop_front();
  }
  return WLength();
}



/*
 * Finishing up recording. If WRITELEVEL all, then 
 * we need to record the stuff which is in the word at
 * shorter depths than the maximum one.
 * Return the number of times "Record" has been called
 */
int ENTRDATA::Parallel_Recorder::FinishRecord(VEDP* data, STARTPOINT st, double prob){
  if(get_all()&&(Shorten()>0)){	// only then need recording clean-up
    Record(data, st, prob);	// record shortened version
    return 1;			// 1 recording done
  }
  return 0;			// no recordings done
}
