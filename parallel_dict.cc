/*
 *  File - parallel_dict.cpp
 *  Project - Entropy_parsing
 *
 *  Dictionary routines.
 *
 *  Created by nemenman Sep 2003
 *  Copyright (c) 2003 by Ilya Nemenman. All rights reserved.
 *
 */
 
#include "simple_dicts.h"

/*
  initializing the dictionary
*/
ENTRDATA::Dict_numb::Dict_numb(const std::string& fname) 
  throw(EntrData_failed_transl): EntrDict(0,0, "num00"){
  //opening file    
  std::string infname = fname + ".txt";
  std::ifstream in(infname.c_str());
  std::string str;		// temorary input string
  if(!in) throw EntrData_failed_transl("Cannot open input file "+ infname);	
  
  // reading stuff at the beginning of the file
  do{
    in>>str;
  }while((str!="scalar")&&(in));
  // getting the alphabet size
  int sz;
  in>>sz;
  if (!in)
    throw EntrData_failed_transl("Could not read the alphabet size from the file.");
  if (sz<=1)
    throw EntrData_failed_transl("Invalid alphabet size.");
  SetASize(sz);
  SetName("num"+num2str(sz));
  SetDef(0);
}


/*
  getting the alphabet size and the number of repeats (streams)
  from the in-file, and returning them through referenced arguments
*/
void ENTRDATA::Dict_numb::Prepare(std::istream& in) const throw(EntrData_failed_transl){
  std::string str;

  // reading crap in the beginning
  do{
    in>>str;
  }while((str!="columns:")&&(in));
  // getting the alphabet size
  int ncols;
  in>>ncols;
  if (!in)
    throw EntrData_failed_transl("Could not read the input stream header.");
}


/* 
 * Translation of chars into ints
 */
int ENTRDATA::Dict_numb::Translate(std::istream& in) 
  const throw (EntrData_failed_transl){

  int intin;

  // if we can get the next char from the string
  if(in>>intin){
    if((intin>=ASize())||(intin<0))
      return Def();
    else 
      return intin;
  }
  else
    // if we cannot get the next integer, throw corresponding exception:
    // could not translate
    throw EntrData_failed_transl("Could not read an integer from the stream.");
}



