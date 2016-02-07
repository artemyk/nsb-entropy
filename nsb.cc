/*
 *  File - nsb.cc
 *  Project - Entropy_parsing
 *
 *  Estimating entropies from counts, creation of appropriate class and 
 *  other auxilliary work.
 *
 *  Copyright (c) 2003 Ilya Nemenman. All rights reserved.
 *
 */

#include "nsb.h"
// number of iterations
const int ENTRDATA::NsbCalc::maxcounter=40;

/*
  transforming new counts into NsbCalc and adding them to NsbSet
*/
void ENTRDATA::NsbSet::add(const counts& cts, const Prior* pp, DOPART pa, double prec) 
  throw(EntrData_bad_alloc, EntrData_badnum){
 
  if(cts.NEvents()>0){
    try{
      NsbCalc* nsb = new NsbCalc(cts, *this, *interp, pp, pa, prec);
      
      if(nsb->copied_from()==NULL)
	data.push_back(nsb);	// do not need to delete nsb here,
				// just record the pointer to it
      else{
	data.push_back(nsb->copied_from()); // record the pointer to
				// precalculated nsb
	nsb->copied_from()->addref();
	delete nsb;		// delete current nsb
      }
    }
    catch (std::bad_alloc){
      throw EntrData_bad_alloc();
    }
  }
}



/*
  Comparing if two nx, kx sets are the same
  return 1  -- the same; not 1 -- not the same
 */
int ENTRDATA::operator==(const NsbCalc& a, const NsbCalc& b){
  int sz=a.size();

  if (sz != b.size()) return 0; // sizes are unequal
  if (sz==0) return 1;		// if both are of zero length, they 
				// are equal
  if (a.get_K()!=b.get_K()) return 0; // alphabet sizes unequal
  if (a.get_N()!=b.get_N()) return 0; // number of counts different
  if (a.get_part()!=b.get_part()) return 0; // partitioning strategy is different

  if (gsl_fcmp(a.GetErr(),b.GetErr(),a.GetErr()))
    return 0;			// requested error tolerances must
				// also be the same; gsl_fcmp returns
				// 0 if arguments are the same

  if (a.pr->GetName() != b.pr->GetName()) 
    return 0;			// comparing if prior names are same
  if (gsl_fcmp(a.pr->GetMin(), b.pr->GetMin(), a.GetErr()))
    return 0;			// now comapring sameness of min and max
  if (gsl_fcmp(a.pr->GetMax(), b.pr->GetMax(), a.GetErr()))
    return 0;

  int eq=1;			// assume they are equal
  // notice that we compare only rounded nx's and kx's
  // fractional counts and occupancy will be rounded!!!
  for(int i=0; ((i<sz) && eq); i++){
    eq &= (round(a.nx[i])==round(b.nx[i]));
    eq &= (round(a.kx[i])==round(b.kx[i]));
  } // this will come out 1 only if both arrays are term by term
				// equal 
  return eq; 
}



/*
  copying all calculated data from another NsbCalc structure
  this may be a non-neeeded function as I later erase duplicate
  objects, leaving pointers to the first copy instead
*/
void ENTRDATA::NsbCalc::CopyCalcs(NsbCalc& old){
  Snsb = old.get_Snsb();
  dSnsb= old.get_dSnsb();
  Sml  = old.get_Sml();
  Scl  = old.get_Scl();
  Bcl  = old.get_Bcl();
  xicl = old.get_xicl();
  dxicl = old.get_dxicl();
  Sas  = old.get_Sas();
  dSas = old.get_dSas();


  N  = old.get_N();
  K1 = old.get_K1();
  K  = old.get_K();
  K2 = old.get_K2();
  part = old.get_part();
}

/*
  Writing results for the case of K<=1
*/
void ENTRDATA::NsbCalc::one_bin(){
  Snsb = 0.0;
  dSnsb= 0.0;
  Sml  = 0.0;
  Scl  = 0.0;
  Bcl  = 0.0;
  xicl = 0.0;
  dxicl = 0.0;
  Sas  = 0.0;
  dSas = 0.0;
}



/*
  Constructor of the NsbCalc from counts.
  We do the following things:
  1) Initialize simple constants;
  2) Construct nx and kx arrays (checking for allocation problems)
  3) Check if same nx, kx have already appeared and, if yes,
     repopulate difficult calculational results from there
  4) Do calculations, if the data has not appeared before
  Later another routine will/won't add this new objects to NsbSet
 */
ENTRDATA::NsbCalc::NsbCalc(const counts& cts, const NsbSet& comp, BXIset& bxiall,
			   const Prior* pp, DOPART pa, double newerr)
  throw(EntrData_bad_alloc, EntrData_badnum) : 
  refs(1), kxng1(NULL), nxng1(NULL), sz(cts.size()-1), szng1(0), 
  K(cts.Cardinality()), maxent(log(K)),  N(cts.NEvents()),
  interp(NULL), pr(pp), part(pa), same_as(NULL), 
  err(newerr), Snsb(0.0), dSnsb(0.0), Sml(0.0), Scl(0.0), 
  Sas(0.0), dSas(0.0), Bcl(0.0), xicl(0.0), dxicl(0.0), 
  mlog(0.0), warncode(NSB_OK) {	
  // note that sz is set to cts.size()-1 so that we do not record
  // nx=0 element

  // resize arrays to the appropriate size 
  try{
    nx = new double[sz];
    kx = new double[sz];
  }
  catch(std::bad_alloc){
    throw EntrData_bad_alloc();
  } 

  // copy data from counts
  {
    int i =0;
    for(CIMCO p= cts.begin(); p!=cts.end(); p++){
      if((p->first) != 0){	// don't copy data with nx==0
	nx[i] = (double) p->first;
	kx[i] = p->second;
	i++;
      }
    }
  }

  // check if this data has been analysed before and, if yes, copy
  // data from there, and reset
  for(ciLNC p = comp.begin(); ((p!=comp.end()) && (same_as==NULL)); p++){
    if ((*(*p))==(*this)){
      same_as = *p;		// these counts are not new
      CopyCalcs(*(*p));		// copy results
    }
  }
  
  CommonConstructor(bxiall);
}



/*
  Constructor of the NsbCalc from two vectors, not associated with the counts object.
  We do the following things:
  1) Initialize simple constants;
  2) Construct nx and kx arrays (checking for allocation problems)
  3) No check for if this data has been analyzed before is done.
  4) Do calculations.
  Later another routine will/won't add this new objects to NsbSet

  This routine has a large overap with the other constructor. At some
  point they need to be merged.
 */
ENTRDATA::NsbCalc::NsbCalc(const vecdoub& _kx, const vecdoub& _nx, double _K, 
			   BXIset& bxiall, const Prior* pp, DOPART pa, double newerr)
  throw(EntrData_bad_alloc, EntrData_badnum) : 
  refs(1), kxng1(NULL), nxng1(NULL), sz(_kx.size()), szng1(0), 
  K(_K), maxent(log(K)), N(0),  interp(NULL), pr(pp), part(pa), same_as(NULL),   
  err(newerr), Snsb(0.0), dSnsb(0.0), Sml(0.0), Scl(0.0), 
  Sas(0.0), dSas(0.0),
  Bcl(0.0), xicl(0.0), dxicl(0.0), mlog(0.0), warncode(NSB_OK) {	
  // note that sz is set to cts.size()-1 so that we do not record
  // nx=0 element

  // resize arrays to the appropriate size 
  try{
    nx = new double[sz];
    kx = new double[sz];
  }
  catch(std::bad_alloc){
    throw EntrData_bad_alloc();
  } 

  // counting number of samples in this data set
  // and copying data
  for(unsigned long i=0;i<sz;i++){
    N+= _kx[i]*_nx[i];
    nx[i] = _nx[i];
    kx[i] = _kx[i];
  }
  CommonConstructor(bxiall);
}



/**********************
  Common part of both constructors
*/
void ENTRDATA::NsbCalc::CommonConstructor(BXIset& bxiall){
  // check for K<=1, and for N==0
  if(K<=1.0){
    one_bin();			// set calculated values to zero and return
    return;
  }
  if(N<=1.0){
    one_bin();
    warning(NSB_NOEVENTS);	// issue a warning
    return;
  }
  if(same_as==NULL){		// if the obj has not appeared before
    if(part==NOPART){
      interp = bxiall.get_spline(K, err); // setting spline pointers.
      DoCalcs();
    }else{			// YESPART
      DoCalcsPart(bxiall);
    }
  } 
}


/*
  The destructor for the NsbCalc object. Need to deallocate all
  arrays.
 */
ENTRDATA::NsbCalc::~NsbCalc(){
  delete[] nx;
  delete[] kx;
  delete[] nxng1;		// can delete even if they are
  delete[] kxng1;		// unset
  delete pr;			// remove the prior
}



/*
  Printing all NsbSet's. Fname has the actual file name + extra tags. 
  "names" is the names of different phases. Start depth is the
  smallest depth in VNsbSet
*/
void ENTRDATA::VNsbSet::Print(const std::string& fname, 
			      const std::vector<std::string>& names)const{
  int nlength=size()/stride;	// # of different length's recorded
  for(int i=0;i<nlength;i++){
    int depth=i+start_depth;	// actual depth being analyzed
    for(int j=0;j<stride;j++){
      //opening file    
      int ind=i*stride+j;	// index of current record
      std::string outfname = fname+names[j] + "_" + depth + "_entr.txt";
      std::ofstream out(outfname.c_str());
      if(!out) throw EntrData_bad_alloc("Error working with file " +outfname);
      out<< "# Created by "<< ENTRDATA_VER(EDVERNAME)<<", " <<  ENTRDATA_VER(EDVERVER)<<"\n";
      out<<std::setprecision(14); // 14 digits precision
	
      if((*this)[ind].size()>0){
	out<< "# name: N\n";
	out<< "# type: matrix\n";
	out<< "# rows: 1\n";
	out<< "# columns: "<< (*this)[ind].size() <<"\n";
	for(ciLNC p=(*this)[ind].begin(); p!=(*this)[ind].end();p++)
	  out<< (*p)->get_N() << " ";
	out<< "\n";
	
	out<< "# name: len\n";
	out<< "# type: matrix\n";
	out<< "# rows: 1\n";
	out<< "# columns: 1\n";
	out<< depth << "\n";	      // we assume that the depth
				      // is the same for the
				      // whole NsbSet
	
	out<< "# name: Snsb\n";
	out<< "# type: matrix\n";
	out<< "# rows: 1\n";
	out<< "# columns: "<< (*this)[ind].size() <<"\n";
	for(ciLNC p=(*this)[ind].begin(); p!=(*this)[ind].end();p++)
	  out<< (*p)->get_Snsb() << " ";
	out<< "\n";
	
	out<< "# name: dSnsb\n";
	out<< "# type: matrix\n";
	out<< "# rows: 1\n";
	out<< "# columns: "<< (*this)[ind].size() <<"\n";
	for(ciLNC p=(*this)[ind].begin(); p!=(*this)[ind].end();p++)
	  out<< (*p)->get_dSnsb() << " ";
	out<< "\n";
	
	out<< "# name: Sml\n";
	out<< "# type: matrix\n";
	out<< "# rows: 1\n";
	out<< "# columns: "<< (*this)[ind].size() <<"\n";
	for(ciLNC p=(*this)[ind].begin(); p!=(*this)[ind].end();p++)
	  out<< (*p)->get_Sml() << " ";
	out<< "\n";
	
	out<< "# name: Sas\n";
	out<< "# type: matrix\n";
	out<< "# rows: 1\n";
	out<< "# columns: "<< (*this)[ind].size() <<"\n";
	for(ciLNC p=(*this)[ind].begin(); p!=(*this)[ind].end();p++)
	  out<< (*p)->get_Sas() << " ";
	out<< "\n";
	
	out<< "# name: dSas\n";
	out<< "# type: matrix\n";
	out<< "# rows: 1\n";
	out<< "# columns: "<< (*this)[ind].size() <<"\n";
	for(ciLNC p=(*this)[ind].begin(); p!=(*this)[ind].end();p++)
	  out<< (*p)->get_dSas() << " ";
	out<< "\n";

	out<< "# name: Scl\n";
	out<< "# type: matrix\n";
	out<< "# rows: 1\n";
	out<< "# columns: "<< (*this)[ind].size() <<"\n";
	for(ciLNC p=(*this)[ind].begin(); p!=(*this)[ind].end();p++)
	  out<< (*p)->get_Scl() << " ";
	out<< "\n";
	
	out<< "# name: Bcl\n";
	out<< "# type: matrix\n";
	out<< "# rows: 1\n";
	out<< "# columns: "<< (*this)[ind].size() <<"\n";
	for(ciLNC p=(*this)[ind].begin(); p!=(*this)[ind].end();p++)
	  out<< (*p)->get_Bcl() << " ";
	out<< "\n";
	
	out<< "# name: xicl\n";
	out<< "# type: matrix\n";
	out<< "# rows: 1\n";
	out<< "# columns: "<< (*this)[ind].size() <<"\n";
	for(ciLNC p=(*this)[ind].begin(); p!=(*this)[ind].end();p++)
	  out<< (*p)->get_xicl() << " ";
	out<< "\n";
	
	out<< "# name: dxicl\n";
	out<< "# type: matrix\n";
	out<< "# rows: 1\n";
	out<< "# columns: "<< (*this)[ind].size() <<"\n";
	for(ciLNC p=(*this)[ind].begin(); p!=(*this)[ind].end();p++)
	  out<< (*p)->get_dxicl() << " ";
	out<< "\n";
	
      }
    }
  }
}
