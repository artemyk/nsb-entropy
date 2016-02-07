/* 
   File  -- random.cc
   Definition of random number function for the 
   read-through random number generator
*/

#include "random.h"

  /* constructor for standard GSL random generators;
     gets the desired generator type (T) on the input; 
     allocates and seeds the genrator from system time;
     default generator type is gsl_rng_mt19937
  */
rnd_gen_class::rnd_gen_class(const gsl_rng_type * T ): 
  which_type(0), state(0),produced(0), input(NULL),in_name("") {
  rnd_gen = gsl_rng_alloc (T);	// generator allocation
  if (rnd_gen)			// non-zero pointer, OK allocation
    gsl_rng_set (rnd_gen, time(NULL)); // initialization with time
  else			// zero pointer, bad allocation
    state = 1;
}

/* constructor for the read-thorugh generator
   open the input file; check for errors; initialize (default first=0)*/
rnd_gen_class::rnd_gen_class(const std::string& fname, const unsigned int first): 
  which_type(1), state(0),produced(0), rnd_gen(NULL), in_name(fname){
  input = new std::ifstream(in_name.c_str());
  if (!(*input)){		// bad initialization
    state = 1;
  }
  else				// correct initialization
    for(unsigned int i=0; i<first; i++)	// now do seeding
      if (!state)		// by means of reading the "first" 
	uniform();		// number of symbols
}


/* 
   destructor
   need to check if the generator is file or GSL and destruct 
   differently
*/
rnd_gen_class::~rnd_gen_class(){
  if (!which_type)		// GSL generator 
    gsl_rng_free (rnd_gen); 
  else				// file generator
    delete input;
}

/* getting a uniform deviate */
double rnd_gen_class::uniform() {
  double rand_number;

  if (!state)			// no error in the status
    if (!which_type){		// gsl generator
      rand_number = gsl_rng_uniform (rnd_gen);
      produced++;
    }
    else			// read-from-file generator
      if((*input)>>rand_number){ // OK reading from file
	produced ++;
      }
      else{			// reading was bad -- probably either 
				// EOF or misformatting;
	if (input->eof()){	// for EOF, reopen file
	  delete input;
	  input = new std::ifstream(in_name.c_str());
	  state = 0;
	  if (!(*input)){	// bad initializatin
	    state = 1;
	    rand_number =0.0;	// will return 0
	  }
	  else			// ok initialization, get a new random #
	    rand_number = uniform();
	}
	else{			// for misformatting, change 
	  state = 1;		// the state to error
	  rand_number = 0.0;
	}
      }
  else				// error in the status
    rand_number= 0.0;		// return 0

  return rand_number;
} 

/* returning a normal deviate */
double rnd_gen_class::normal(const double sigma) { 
    if (!state)			// ok status
      if (!which_type){		// gsl generator
	produced++;
	return gsl_ran_gaussian (rnd_gen, sigma);
      }
      else{			// file generator
	double x(uniform()-0.5); // getting a uniform deviate
	// and transforming it to normal
	return sigma*M_SQRT2*dierfc(2.0*fabs(x))* ((x>0)? 1.0 : -1.0);
      }
    else			// not ok status
      return 0.0;
  } 
