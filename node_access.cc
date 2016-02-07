/*
 *  File - node_access.cpp
 *  Project - Entropy_parsing
 *
 *  Interface functions for the node class.
 *
 *  Created by nemenman on Mon Nov 05 2001.
 *  Copyright (c) 2001 Ilya Nemenman. All rights reserved.
 *
 */
 
#include "../include/utility.h"
#include "node.h"

/*
 * Incrementing the visitation count
 */
void ENTRDATA::node::IncVisited(const long v) throw(EntrData_range_error){
  ClearCounts();
  visited+=v;
  if (visited <=0) throw 
    EntrData_range_error("Negative visitation count " + num2str(visited) +".\n");
}
    

/*
 * Decrementing the visitation count
 */    
void ENTRDATA::node::DecVisited(const long v) throw(EntrData_range_error) {
    IncVisited(-v);
}


/*
 * Parsing the word (list) into the tree. 
 * Arguments: current position in the word, next after the end of the word.
 * 
 * The function increments the visitation count by one, and then calls itself on
 * the appropriate child with the same word shortened by one.
 */
void ENTRDATA::node::RecordList(const lint& word, LCI pos) 
throw (EntrData_range_error,EntrData_bad_alloc){
    visited++;
    ClearCounts();		// clear the counts when making recordings

    if (pos!=word.end()){		// this is not the last data point in the parsed list
        int act_child= *pos;		// active (to be worked upon) child
        AddChild(act_child);		// adding the child #word[pos]
        
        children[act_child]->RecordList(word, ++pos);
    };
    // if it was the last data point, just exit
}


/*
 * Forward printing of the tree to the output stream
 */
void ENTRDATA::node::Print(std::ostream& out) const{
    out<<GetVisited()<<" ";
    for(int i=0;i<ASize();i++){
        if(children[i])
            children[i]->Print(out);
	else
	    out << "0 ";
    };
    out<<"\n";
}


    
/*
 * Prunning the tree of the "hanging" nodes with zero visitation counter. The root of the tree 
 * is not prunned, even if its counter is zero.
 */
void ENTRDATA::node::prune(){
  ClearCounts();
  for(int i=0;i<ASize();i++)
    if(children[i])				// children[i] exists
      if(children[i]->GetVisited() ==0){ 
	delete children[i];			// delete if needed
	children[i]=NULL;
      }
      else
	children[i]->prune();			// else prune the child
}

