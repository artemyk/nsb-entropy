/*
 *  File - mpi.h
 *  Project - Entropy_parsing
 * 
 *  Declaration of the mpi -- multiple precision integer class.
 *
 *  Copyright (c) 2001, 2002 Ilya Nemenman. All rights reserved.
 *
 */

#ifndef KEY_MPI_H
#define KEY_MPI_H

#include <string>
#include <iostream>
#include <gmp.h>

namespace GMP{
    class mpi	{
    private:
	mpz_t data;		// the real multiple precision data
    
    public:  
 
	// data access functions
	void get_data(mpz_t out) const //returning data in
				       //pre-initialized mpz
	    {mpz_set(out, data); return;} 
	void get_data_init(mpz_t out) const //initializing mpz and
					    //returning data in it
	    {mpz_init_set(out, data); return;}
	mpz_t* get_mpz()  {return &data;} // returning the same mpz as
					  // in "data", not its copy.
	
	// constructors and destructors
	mpi(unsigned long t=0){mpz_init_set_ui(data,t);} // initializing 
				// and setting from unsigned long
	mpi(const mpi& src){src.get_data_init(data);} //copy constructor; 
				//initializes and sets from another key_mpi
	mpi(const mpz_t src){mpz_init_set(data,src);} //copying and constructing 
				//from an mpz 
	~mpi(){mpz_clear(data);} //destructor -- we must clear all mpz's

	void copy(const mpz_t i){mpz_set(data,i);} //copying from an mpz
	
	// output functions
	std::string to_string(int base) const { // outputting mpi to a string
	    std::string s = mpz_get_str(NULL, base, data);	    
	    return s;
	};
    
	// math functions
	void addmul_ui(mpi m, const unsigned long u)
	  {mpz_addmul_ui(data,*(m.get_mpz()),u);} // data = data +m*u
	void mul_ui(const unsigned long u)
	  {mpz_mul_ui(data,data,u);} // data = data*u
    
    	// output operators
	friend std::ostream& operator<<(std::ostream&, const mpi&);
	
	//math operators
	void operator+=(unsigned long& x){
	    mpz_add_ui(data, data, x);
	}    
	    
	friend bool operator==(const mpi&, const mpi&);
	friend bool operator<(const mpi&, const mpi&);
    };

    inline std::ostream& operator<<(std::ostream& os, const mpi& x) {
      std::string tmp= x.to_string(10);
	return os << tmp;
    }

    //is lhs equal ro rhs?
    inline bool operator==(const mpi& lhs, const mpi& rhs){
	return mpz_cmp(lhs.data,rhs.data);
    }

    //is lhs less than rhs?
    inline bool operator<(const mpi& lhs, const mpi& rhs){
	return (mpz_cmp(lhs.data,rhs.data)<0);
    }
}

#endif


