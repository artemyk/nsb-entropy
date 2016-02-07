/*
 *  File - counts.cpp
 *  Project - Entropy_parsing
 *
 *  Definitions for the counts class.
 *
 *  Copyright (c) 2001, 2002 Ilya Nemenman. All rights reserved.
 *
 */

#include "counts.h"


/* 
 * Simple subscripting. If accessing data[0] (events that did not
 * occur), need to recalculate it first, if necessary (this field
 * usually contains excessively large numbers, and they are better
 * calculated as (maximum # of events) - (# of events that occured),
 * instaed of ++'ing for each not occured event.
 */
double ENTRDATA::counts::operator[](const long& key) const throw (EntrData_bad_alloc){
  try{
    if (key)	//key !=0
	return data[key];
    else{
	if(!isvalid) calculate();
	return data[0];
    }
  }
  catch(std::bad_alloc){
    throw EntrData_bad_alloc();
  }
}

/*
 * Subscripting operator, which allows to change the data.  Note that
 * changing counts[key] makes counts[0] as well as "ne" invalid (see
 * comments for operator[]).  Thus every access to this operator must
 * be followed by a reset of the "isvalid" flag.
 */
double& ENTRDATA::counts::operator()(const long& key) throw (EntrData_bad_alloc){
    if(!key)		 	// if trying to access key == 0
	if(!isvalid) 		// and if invalid
	    calculate();	// then need to recalculate

    isvalid = 0;
    try{
	return data[key];
    }
    catch(std::bad_alloc){
	throw EntrData_bad_alloc();
    };
}


/*
 * Number of events in the "counts"
 */
double ENTRDATA::counts::NEvents()const{
    if(isvalid)
	return ne;
    else{
	calculate();
	return ne;
    };
}


/*
 * Calculating the number of events in "counts" and the value of counts[0]
 */
void ENTRDATA::counts::calculate()const{
    if (isvalid) return; 	// no recalculation is necessary
    
    double occured=0.0;		// total number of events that occured
				// at least once. To set data[0] we
				// need this (see comments for
				// operator []).
    ne = 0.0;
    // starting from events that occured at least once
    // for this we make sure that the entry for events 
    // that did not occure exist
    // if it exists, don't change it so far; else if 
    // it doesn't -- insert it
    data[0];
    for(CIMCO p=(++data.begin()); p!=data.end();++p){
	ne += p->first * p->second;
	occured += (double)p->second;
    }
    
    data[0] = mne - occured;
    
    // don't do this check; always keep data[0] element
    //if(!floor(data[0]))		// if 0 events with zero count, erase this entry
    //	data.erase(data.find(0));
	
    isvalid=1;
}
    
/*
 * Printing the counts.
 */
std::ostream& ENTRDATA::operator<<(std::ostream& out, const counts& cts){
  //by default, newprec=14
  int oldprec = out.precision();
  out<<std::setprecision(cts.precision());
  for(CIMCO p=cts.data.begin(); p!=cts.data.end();p++)
    out<<p->first<<"\t"<<p->second<<"\n";
  out<<std::setprecision(oldprec);
  return out;
}
