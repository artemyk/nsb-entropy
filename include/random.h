#ifndef RANDOM_NUMBERS
#define RANDOM_NUMBERS
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include "specfunctions/specfunctions.h"



////////////////////////////////////
//random number generator class/////
////////////////////////////////////
/*
  Random number generator may be of two types
  (1) standard gsl_*** generators supplied by GSL (with all of 
      the corresponding subtypes)
  (2) generator which just reads [presumably U(0,1)] numbers from some 
      text file 
*/

class rnd_gen_class
{
private:
  int which_type;		// 0 -- gsl generator; 
				// 1 -- readthrough generator
  int state;			// 0 -- ok to generate random #'s
				// 1 -- some error
  unsigned long produced;	// # of randoms already generated
  gsl_rng* rnd_gen;		// random number generator
  std::ifstream* input;		// pointer to the input file stream for a
				// reathrough generator
  std::string in_name;		// input file name

public:

  // constructor for standard GSL generators
  rnd_gen_class(const gsl_rng_type* =gsl_rng_mt19937);

  // constructor for the read-from-file generator
  rnd_gen_class(const std::string&, const unsigned int=0);

  //destructor
  ~rnd_gen_class();

  // checking the generator status
  int get_state() const {return state;}; 
  unsigned long how_many() const {return produced;};

  //getting a uniform [0,1) deviate from the generator
  double uniform();

  //getting a normal deviate around zero with std. dev. sigma
  double normal(const double);

};

#endif
