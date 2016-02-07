/*
 *  File - node_count.cpp
 *  Project - tree
 *
 *  Counting nodes and their frequencies.
 *
 *  Created by nemenman on Thu Nov 08 2001.
 *  Copyright (c) 2001 Ilya Nemenman. All rights reserved.
 *
 */

#include <math.h>


#include "node.h"

/*
 * The function creates the "counts" object and counts events
 * in the node into it.
 */ 
ENTRDATA::counts* ENTRDATA::node::MakeCounts(int depth) const throw(EntrData_bad_alloc){
    counts* p;
    try{
	p= new counts(pow((double) ASize(), (double) depth));
    }
    catch (std::bad_alloc){
	throw EntrData_bad_alloc();
    };

    JustCount(p, depth);			// do the actual counting and recording
    
    p->NEvents();				// to force recalculation of counts[0]
    mcts[depth]=p;
    return p;
}

/*
 * The function that does the actual counting and recording.
 */
void ENTRDATA::node::JustCount(counts* p, int depth) const{
    if(!depth){			// we are at the appropriate depth and need start writing
	(*p)(GetVisited()) += 1.0;
    }
    else				// decrease depth and start counting children
	for(int i=0; i<ASize();i++)
	    if(children[i])		// this child exists
		children[i]->JustCount(p, depth-1);
}
