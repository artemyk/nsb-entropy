/*
 *  File - Dict_m2o.cpp
 *  Project - Entropy_parsing
 *
 *  Dictionary routines.
 *
 *  Created by nemenman on Tue Nov 20 2001.
 *  Copyright (c) 2001, 2003 by Ilya Nemenman. All rights reserved.
 *
 */
 
#include "simple_dicts.h"

/*
 * Constructing the dictionary of predefined types
 * default dictionary is ENGL29, this is also a value
 * to be adopted for incorrect values of alh
 */
ENTRDATA::Dict_m2o::Dict_m2o(DIFFM2O alph):
  EntrDict(M2OSizes[alph], M2OSizes[alph]-1, M2ONames[alph]){

  switch (alph){
  case BINARY:
    alphabet['0'] = 0;
    alphabet['1'] = 1;
    break;
  case TETRA:
    alphabet['0'] = 0;
    alphabet['1'] = 1;
    alphabet['2'] = 2;
    alphabet['3'] = 3;
    break;
  case OCTA:
    alphabet['0'] = 0;
    alphabet['1'] = 1;
    alphabet['2'] = 2;
    alphabet['3'] = 3;
    alphabet['4'] = 4;
    alphabet['5'] = 5;
    alphabet['6'] = 6;
    alphabet['7'] = 7;
    break;
  case GENOME:
    alphabet['A'] = 0;
    alphabet['a'] = 0;
    alphabet['C'] = 1;
    alphabet['c'] = 1;
    alphabet['G'] = 2;
    alphabet['g'] = 2;
    alphabet['T'] = 3;
    alphabet['t'] = 3;
    break;
  case ENGL3:			// alph size: 1 wovel, 1 consonant, 1 everything else
    alphabet['a'] = 0;
    alphabet['A'] = 0; 
    alphabet['b'] = 1;
    alphabet['B'] = 1;
    alphabet['c'] = 1;
    alphabet['C'] = 1; 
    alphabet['d'] = 1;
    alphabet['D'] = 1;
    alphabet['e'] = 0;
    alphabet['E'] = 0; 
    alphabet['f'] = 1;
    alphabet['F'] = 1;
    alphabet['g'] = 1;
    alphabet['G'] = 1; 
    alphabet['h'] = 1;
    alphabet['H'] = 1;
    alphabet['i'] = 0;
    alphabet['I'] = 0; 
    alphabet['j'] = 1;
    alphabet['J'] = 1;
    alphabet['k'] = 1;
    alphabet['K'] = 1; 
    alphabet['l'] = 1;
    alphabet['L'] = 1;
    alphabet['m'] = 1;
    alphabet['M'] = 1; 
    alphabet['n'] = 1;
    alphabet['N'] = 1;
    alphabet['o'] = 0;
    alphabet['O'] = 0; 
    alphabet['p'] = 1;
    alphabet['P'] = 1;
    alphabet['q'] = 1;
    alphabet['Q'] = 1; 
    alphabet['r'] = 1;
    alphabet['R'] = 1;
    alphabet['s'] = 1;
    alphabet['S'] = 1; 
    alphabet['t'] = 1;
    alphabet['T'] = 1;
    alphabet['u'] = 0;
    alphabet['U'] = 0; 
    alphabet['v'] = 1;
    alphabet['V'] = 1;
    alphabet['w'] = 1;
    alphabet['W'] = 1; 
    alphabet['x'] = 1;
    alphabet['X'] = 1;
    alphabet['y'] = 1;
    alphabet['Y'] = 1; 
    alphabet['z'] = 1;
    alphabet['Z'] = 1;
    alphabet['0'] = 2;		// digits and all the rest are 2
    break;
  case ENGL38:			// alphabet size -- 26 roman letters (case-insensitive),
				// + ten digits, plus space, plus all the rest 
    alphabet['a'] = 0;
    alphabet['A'] = 0; 
    alphabet['b'] = 1;
    alphabet['B'] = 1;
    alphabet['c'] = 2;
    alphabet['C'] = 2; 
    alphabet['d'] = 3;
    alphabet['D'] = 3;
    alphabet['e'] = 4;
    alphabet['E'] = 4; 
    alphabet['f'] = 5;
    alphabet['F'] = 5;
    alphabet['g'] = 6;
    alphabet['G'] = 6; 
    alphabet['h'] = 7;
    alphabet['H'] = 7;
    alphabet['i'] = 8;
    alphabet['I'] = 8; 
    alphabet['j'] = 9;
    alphabet['J'] = 9;
    alphabet['k'] = 10;
    alphabet['K'] = 10; 
    alphabet['l'] = 11;
    alphabet['L'] = 11;
    alphabet['m'] = 12;
    alphabet['M'] = 12; 
    alphabet['n'] = 13;
    alphabet['N'] = 13;
    alphabet['o'] = 14;
    alphabet['O'] = 14; 
    alphabet['p'] = 15;
    alphabet['P'] = 15;
    alphabet['q'] = 16;
    alphabet['Q'] = 16; 
    alphabet['r'] = 17;
    alphabet['R'] = 17;
    alphabet['s'] = 18;
    alphabet['S'] = 18; 
    alphabet['t'] = 19;
    alphabet['T'] = 19;
    alphabet['u'] = 20;
    alphabet['U'] = 20; 
    alphabet['v'] = 21;
    alphabet['V'] = 21;
    alphabet['w'] = 22;
    alphabet['W'] = 22; 
    alphabet['x'] = 23;
    alphabet['X'] = 23;
    alphabet['y'] = 24;
    alphabet['Y'] = 24; 
    alphabet['z'] = 25;
    alphabet['Z'] = 25;
    alphabet['0'] = 26;
    alphabet['1'] = 27; 
    alphabet['2'] = 28;
    alphabet['3'] = 29;
    alphabet['4'] = 30;
    alphabet['5'] = 31; 
    alphabet['6'] = 32;
    alphabet['7'] = 33;
    alphabet['8'] = 34;
    alphabet['9'] = 35; 
    alphabet[' '] = 36;
    alphabet['.'] = 37;
    break;
  default:			//case ENGL29
				// alphabet size -- 26 roman letters (case-insensitive),
				// + digits, plus space, plus all the rest 
    alphabet['a'] = 0;
    alphabet['A'] = 0; 
    alphabet['b'] = 1;
    alphabet['B'] = 1;
    alphabet['c'] = 2;
    alphabet['C'] = 2; 
    alphabet['d'] = 3;
    alphabet['D'] = 3;
    alphabet['e'] = 4;
    alphabet['E'] = 4; 
    alphabet['f'] = 5;
    alphabet['F'] = 5;
    alphabet['g'] = 6;
    alphabet['G'] = 6; 
    alphabet['h'] = 7;
    alphabet['H'] = 7;
    alphabet['i'] = 8;
    alphabet['I'] = 8; 
    alphabet['j'] = 9;
    alphabet['J'] = 9;
    alphabet['k'] = 10;
    alphabet['K'] = 10; 
    alphabet['l'] = 11;
    alphabet['L'] = 11;
    alphabet['m'] = 12;
    alphabet['M'] = 12; 
    alphabet['n'] = 13;
    alphabet['N'] = 13;
    alphabet['o'] = 14;
    alphabet['O'] = 14; 
    alphabet['p'] = 15;
    alphabet['P'] = 15;
    alphabet['q'] = 16;
    alphabet['Q'] = 16; 
    alphabet['r'] = 17;
    alphabet['R'] = 17;
    alphabet['s'] = 18;
    alphabet['S'] = 18; 
    alphabet['t'] = 19;
    alphabet['T'] = 19;
    alphabet['u'] = 20;
    alphabet['U'] = 20; 
    alphabet['v'] = 21;
    alphabet['V'] = 21;
    alphabet['w'] = 22;
    alphabet['W'] = 22; 
    alphabet['x'] = 23;
    alphabet['X'] = 23;
    alphabet['y'] = 24;
    alphabet['Y'] = 24; 
    alphabet['z'] = 25;
    alphabet['Z'] = 25;
    alphabet['0'] = 26;
    alphabet['1'] = 26; 
    alphabet['2'] = 26;
    alphabet['3'] = 26;
    alphabet['4'] = 26;
    alphabet['5'] = 26; 
    alphabet['6'] = 26;
    alphabet['7'] = 26;
    alphabet['8'] = 26;
    alphabet['9'] = 26; 
    alphabet[' '] = 27;
    alphabet['.'] = 28;
    break;
  }
}

/* 
 * Translation of chars into ints
 */
int ENTRDATA::Dict_m2o::Translate(std::istream& in) const throw (EntrData_failed_transl){
  char a;


  // if we can get the next char from the string
  if(in.get(a)){
    DictMapCI res=alphabet.find(a); // look for the char key in the map
    if (res == alphabet.end())	// if not found,
      return Def();		// return the default value
    else			// if found,
      return res->second;	// return what's found
  }
  else
    // if we cannot get the next char, throw corresponding exception:
    // could not translate
    throw EntrData_failed_transl("Could not read from the stream.");
}



