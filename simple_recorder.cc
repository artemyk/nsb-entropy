/*
 *  File - simple_recorder.cc
 *  Project - Entropy_parsing
 *
 *  Definition of functions for the simple recorder class.
 *
 *  Copyright (c) 2002 Ilya Nemenman. All rights reserved.
 *
 */

#include "simple_recorder.h"

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
void ENTRDATA::Simple_Recorder::Record(VEDP* data, 
				       STARTPOINT st, double prob){
  // create a random number generator needed to decide whether to record 
  // a word or not
  static rnd_gen_class rnd_gen(gsl_rng_mt19937);  

  // if prob!=1.0, then we are doing probabilistics recording. In this
  // case, get a random # for each of the available recording objects
  if(prob<1.0){ 		// this is probabilistic recording

    // with prob<1.0, manystarts is not possible; we always record the 0'th start 
    // only with the requested probability
    for(unsigned int i=0;i<data->size();i++){
      if(rnd_gen.uniform()<=prob)
	RecordSingle((*data)[i]);
    }
  } // otherwise, it's recording with SHIFT
  // if we do all possible recordings (MANYSTEPS) or if we are 
  // now at the recording # with corresponding  VEDP existing
  else if ((st == MANYSTARTS) || ((unsigned int)to_rec<data->size())){
    RecordSingle((*data)[to_rec]);
  }
  to_rec = ( (unsigned int)(to_rec+1) == (data->size())) ? 0 : to_rec+1;
}



/*
 * Updating the translation word
 */
void ENTRDATA::Simple_Recorder::UpdateWord(int newint){

  // add the newly translated symbol to the word
  word.push_back(newint);
  if(word.size() > get_depth()) // if the word reached needed
    word.pop_front();		        // depth, pop old int from front
}


/*
 * Recording the word in the single data structure
 */
void ENTRDATA::Simple_Recorder::RecordSingle(EntrData* ed){
  if(get_all()){
    ed->RecordAll(word); //recording at all levels
  }
  else{
    ed->Record(word);    //recording at one level
  }
}
