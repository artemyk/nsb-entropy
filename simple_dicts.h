/*
 *  File - simple_dicts.h
 *  Project - Entropy_parsing
 * 
 *  Declaration of simple many-to-one disctionaries with different
 *  alphabets.
 *
 *  Copyright (c) 2001, 2002 Ilya Nemenman. All rights reserved.
 *
 */

#ifndef SIMPLE_DICTS_H
#define SIMPLE_DICTS_H

#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include "EntrData.h"
#include "utility.h"

namespace ENTRDATA{

   
  // Map classes used for the m2o dictionary.
  typedef std::map<char,int> DictMap;
  typedef DictMap::const_iterator DictMapCI;

  //different possible m2o dictionaries
  enum DIFFM2O {BINARY, TETRA, OCTA, GENOME, ENGL3, ENGL29, ENGL38};
  const int M2OSizes[] = {2, 4, 8, 4, 3, 29, 38};
  const std::string M2ONames[] = 
	    {"2", "4", "8", "gene", "e3", "e29", "e38"};
  const std::string ParName("par");
  const std::string NumName("num");

  
  /*
   * Dictionary classes types for ChooseDict function.
   */
  enum DICTCLASS {M2O,PAR,NUM};	// different possible dictionary
				// Many-to-one, parallel, and very similar
				// number 
				// dictionaries are implemented.

  struct DictParams{		// parameters of dictionaries 
    DICTCLASS type;		// the class they belong to
    int       param;		// extra description within the class
  };


  /*
   *  Many-to-one dictionary class. This class translates input symbol
   *  by symbol. Different standard mappings are possible.
   */
  class Dict_m2o: public EntrDict{
  private:
    DictMap alphabet;		// the alphabet

  public:
    Dict_m2o(DIFFM2O = ENGL29);	// exceptions?
    ~Dict_m2o(){}
    int Translate(std::istream&) const throw (EntrData_failed_transl);	
				// translating the char according to the alphabet
    void Prepare(std::istream& in) const throw(){} // no preparation is needed for
				// this type of dictionary 
 
  };
  

  /*
    Number dictionary, where each nonnegative integer # from the input
    is just passed on without translation.
  */
  class Dict_numb : public EntrDict{
  private:
  public:
    Dict_numb(const std::string&) 
      throw(EntrData_failed_transl); // intitializing 
				// the father EntrDict with
				// some absured values, to be redone
				// from this->Prepare(..)
    ~Dict_numb(){}
    int Translate(std::istream&) const throw (EntrData_failed_transl);	
				// translating the char according to
				// the alphabet 
    void Prepare(std::istream& in) const throw(EntrData_failed_transl); 
				// preparing the file for reading

  };  

}
#endif
