/*
  Function definition for the main EntrData class
*/


#include "EntrData.h"
/*
   Removing all the counts (to be used when adding new data
   or deconstructing data objects)
*/
void ENTRDATA::EntrData::ClearCounts() const{
  if(!mcts.empty()){
    // removing the counts 
    for(mPCCI p=mcts.begin(); p!=mcts.end();p++)
      if (p->second != NULL)
	delete (p->second);	// deleting data objects
				// does this make p->second equal to NULL?
    
    // now need to resize the map to zero size
    mcts.clear();
  }
}


/*
  Returning counts produced by MakeCounts
*/
const ENTRDATA::counts* ENTRDATA::EntrData::count(int len) 
  const throw(EntrData_bad_alloc){

  if (mcts[len] == NULL)
    MakeCounts(len);		// create counts if absent
  return mcts[len];
}

