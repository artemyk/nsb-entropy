/*
 *  File - node.h
 *  Project - Entropy_parsing
 * 
 *  Declaration of the node class.
 *
 *  Created by nemenman on Fri Nov 02 2001.
 *  Copyright (c) 2001, 2002 Ilya Nemenman. All rights reserved.
 *
 */

#ifndef NODE_H
#define NODE_H

#include <list>
#include <iostream>

#include "EntrData.h"

namespace ENTRDATA{
  class node: public EntrData{
  private:
    node**	children;	// pointer to the array of pointers to children
    long visited;		// the number of times this particular node was visited

    void AddChild(const int, const int v=0) 
      throw (EntrData_range_error, EntrData_bad_alloc);	
                                // creation of the child numbered n; by default 
                                // the child's been visted v=0 times
    void IncVisited(const long v =1) throw (EntrData_range_error); 
				// increments the visitation count
    void DecVisited(const long v =1) throw (EntrData_range_error); 
				// decrements the visitation count	
    void JustCount(counts*, int) const;	// the actual function to
				// record data into "counts" class
    void Print(std::ostream&) const;	// printing the tree in forward direction
    void RecordList(const lint&, LCI) throw (EntrData_range_error,EntrData_bad_alloc);	
				// parsing a list and recording it into a tree
    counts* MakeCounts(int) const throw(EntrData_bad_alloc);
				// recording data into "counts" class

  public:
    node(const int=2, const std::string= "", const int=0) 
      throw (EntrData_bad_alloc);// construction specifying the
				// maximum number of children
				// (default -- binary tree), and
				// the number of times this
				// particular node was visited
    ~node();
    
    int  GetVisited() const {return visited;} // gets the visitation count
    void prune();		// prunning the tree below the current node
    void RecordAll(const lint& v) throw (EntrData_range_error,EntrData_bad_alloc)
      {RecordList(v,v.begin());} // recording list in the tree
    void Record(const lint& v) throw (EntrData_range_error,EntrData_bad_alloc)
      {RecordList(v,v.begin());} // for trees, there's no diference in all-level and
				// single level recording

    };	
}



#endif
