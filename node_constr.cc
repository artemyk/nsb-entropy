/*
 *  File - node_constr.cpp
 *  Project - Entropy_parsing
 * 
 *  Description of the construction/destruction methods for the 
 *  class node (synonym to tree). Creation and destruction of children
 *  of this node.
 *
 *  Created by nemenman on Fri Nov 02 2001.
 *  Copyright (c) 2001, 2002 Ilya Nemenman. All rights reserved.
 *
 */

#include <string>

#include "../include/utility.h"
#include "node.h"


/*
 * Constructor gets the maximum number of children (alphabet size),
 * name of the structure, and the visitation count. Then it creates an
 * array of pointers to the children, and sets the pointers to NULL.
 * The children themselves are not created
 *
 * A bad_alloc_node exception is thrown if allocation fails.
 */ 
ENTRDATA::node::node(const int nchild,  const std::string newname, const int v) 
  throw (EntrData_bad_alloc): EntrData(nchild, newname), visited(v) {
    try{ 
      children = new node* [ASize()];
    }
    catch (std::bad_alloc& exc){
      throw EntrData_bad_alloc();
    };
    
    for(int i=0; i<ASize(); i++){
      children[i] = NULL;	// C++ does not guarantee that _new_ returns cleared 
                                // memory; we clan it
    };
}


/* 
 * The destructor calls destructor for all the children (if they are present)
 * and then clears the memory allocated for the array of pointers to children
 */
ENTRDATA::node::~node(){
  for(int i=0; i<ASize(); i++)
    if(children[i]) delete children[i];
  
  delete [] children;
}


/* 
 * Creating a child #n, and setting its visitation count to v.
 * If the child already exists, do nothing.
 *
 * Exception range_error_node is thrown if n >ASize() or n<0;
 */
void ENTRDATA::node::AddChild(const int n, const int v) 
  throw (EntrData_range_error,EntrData_bad_alloc) {
  if((n>=ASize()) || (n<0)) {
    throw EntrData_range_error("node::AddChild: child number " + num2str(n) +" out of range.\n");
  };
    
  if(!children[n]) 		// child does not exist
    children[n] = new node(ASize(), GetName(), v);
}
