/*
 *  File - parsers.cpp 
 *  Project - Entropy_parsing
 *
 *  Defines some functions used to parse data files into data
 *  structures.
 *
 *  Copyright (c) 2001-2013 Ilya Nemenman. 
 */


// note about exceptions: the first routine being called from main()
// -- ProcessFile has try-catch block for bad_alloc in it. We do not
// insert such checks in routines called from ProcessFile.

#include "parsers.h"


namespace ENTRDATA{
  /*
   * Error output function
   */
  void error(const char* p1, const char* p2){
    std::cerr << p1 << ' '<< p2 << '\n';
    exit(1);
  }
  




  /*
   *  Takes an input stream and records it into all of the data
   *  structures provided. Tries to recorde nsymb symbols in these
   *  data structures (skipping the first "skip" symbols produced by
   *  the dictionary), but if reading from the stream fails, records
   *  less, and returns the number of recorded symbols. Variable
   *  "depth" controls the length of the words to be studied, and
   *  "all" specifies if the words are recorded at all length up to
   *  "depth", or only at "depth". Return the number of symbols
   *  actually read and actually recorded.
   */
  void Stream2Dict2Data(VEDP& out, std::istream& in, const EntrDict* di, 
			Recorder* rec, unsigned long skip, 
			unsigned long nsymb, int depth, 
			STARTPOINT st, double prob,
			unsigned long& have_read, unsigned long& have_rec){

    int tr;			//translated symbol
    
    if((have_read>=nsymb)&&(nsymb!=0)){
      return;			// if already read what's needed
    }

    //recording after the first symbols are skipped
    try{			//try to translate the next symbol
      do{
      	tr = (*di)[in];
	rec->UpdateWord(tr);
	have_read++;    //we've read one more symbol
	
	// is the word appropriate length for recording, or should
	// we grow it further?
	if (rec->Ready()){ 
	  rec->Record(&out, st, prob);
	  have_rec++;
	}
      }while(((have_read<nsymb) || (!nsymb)) && (!rec->Full()) ); 
				// read until have recorded amount, 
				// or, if nsymb==0, read all the way
				// to the end; we must stop if the
				// recorder says that it is full

      //if stopped for any reason but Full(), then FinishRecord()
      if(!rec->Full()) have_rec+=rec->FinishRecord(&out,st,prob);
    }
    catch (EntrData_failed_transl){
      // if cannot read, then FinshRecord 
      have_rec+=rec->FinishRecord(&out,st,prob);
    }

    // remove!!!
    //    std::cout<<(*(out[0]))<<"\n";
  }




  /* 
   * Counting data structures and writing them to file. Each data
   * structure in the supplied vector is counted at depth between start
   * and end and recorded into corresponding files.
   */
  void Data2Cts2File(const VEDP& data, const std::string& fname, const std::string& tag,
		     int depth_start, int depth_end, STARTPOINT st,
		     DOCALCS ca, DOPART pa, const std::string& prior, VNsbSet& results){
    const counts* cts;		// won't need to delete the counts
				// since they are a part of the
				// EntrData object
    
    Prior* pr=NULL;		// pointer to a prior (NULL so that we
				// can erase it easily)
    // depth_start and depth_end show the depths which are to be
    // recorded to counts and eventually to NSB and file;
    // results.SDepth() and results.EDepth() are the start and end
    // depths saved in results.

    for(unsigned int j=0; j < (data.size());j++){ // writing out each data structure
      if((j==0)||(st!=ONESTART)){ // only needed in this case
	for(int depth=depth_start; depth<=depth_end; depth++){ 
	  // access index to "Results"
	  int ind = (depth-results.SDepth())*results.Stride()+j;
	  //writing out each required depth
	  if(ca==NSBNONE){	// no counting needed, just output frequencies of occurencies
	    std::string outname=fname + "_"+tag +"_"+
	      ((data[j])->GetName()) + "_full_" + depth + ".txt";
	    std::ofstream out(outname.c_str());
	    if(!out) error("Cannot open output file", outname.c_str());	
	    out<< *(data[j]);
	    if(!out) error("Could not output full counts with word identities.");
	  }else{         	// need to create counts
	    cts = (data[j])->count(depth);
	    if(cts->NEvents()>0.0){ // only do the rest if the data is 
				// there for this structure, 
				// not for empty structures
	      if(ca==YESCOUNT){	// don't need computations, just
				// output counts
		//creating out-file name, and opening the file
		std::string outname=fname + "_"+tag +
		  ((std::string((data[j])->GetName())).compare("") ? 
		   (std::string((data[j])->GetName())) + "_" : "") 
		   + depth + ".txt";
		std::ofstream out(outname.c_str());
		if(!out) error("Cannot open output file", outname.c_str());	
		out<< *cts;
		if(!out) error("Could not output counts.");
	      }
	      else{			// need to do computations
		// making the prior
		if(prior==UniName)
		  pr = new UniPrior(depth*log(data[j]->ASize()));
		else if (prior==BndName){
		  double a=0.0;
		  double b=depth*log(data[j]->ASize());
		  // as this NsbCalc is not yet attached to its
		  // respective NsbSet, the size of the current NsbSet is
		  // the slice number we are working on
		  int slice = results[ind].size();
		  
		  // loop over partitioning N-mer into N-n and n mers and
		  // bounding entropy above and below
		  for(int part=1; part<depth; part++){
		    // 1 stands for the long part of the word, 2 for its
		    // compliment
		    // ind are the access indices for the long part and
		    // the compliment
		    int ind1 = (depth-part-results.SDepth())*results.Stride()+j;
		    int ind2 = (part-results.SDepth())*results.Stride()+j;
		    // mean and std are for the corresponding NSB estimates
		    double mean1=0.0;
		    double mean2=0.0;
		    double std1=0.0;
		    double std2=0.0;
		    // if such index exists, and its NsbSet s of OK size
		    if(ind1>=0){
		      if(results[ind1].size()>=slice){
			mean1= (results[ind1][slice])->get_Snsb();
			std1 = (results[ind1][slice])->get_dSnsb();
		      }
		    }
		    if(ind2>=0){
		      // getting data from the complimentary size record
		      // shifted by depth-part, which is the length of
		      // the first record
		      if(results[ind2].size()>=slice+depth-part+1){
			// this works for parallel recording, when the
			// entropy of word of lenth N is < than that of
			// N-m  and m taken sequentially, and now the
			// index of ...+depth-part points at this
			// seqeuntial addition.
			mean2= (results[ind2][slice+depth-part])->get_Snsb();
			std2 = (results[ind2][slice+depth-part])->get_dSnsb();
		      }
		      else if(results[ind2].size()>slice){
			// this will work for regular recording, when the
			// size of all entropy vectors is 1 (no slices,
			// and we estimate along the file, no across). in
			// this case we always get Snsb from the same
			// spot as the current length, ind, entry in the
			// NSB natrix. in
			// addition, here we should not take sequential
			// N-m and m, but rather the only ones, the first
			// ones
			mean2= (results[ind2][slice])->get_Snsb();
			std2 = (results[ind2][slice])->get_dSnsb();
		      }
		    }
		    // we make the uniform box with 3.0 stddevs on each side
		    a = GSL_MAX(a,GSL_MAX(mean1-3.0*std1,mean2-3.0*std2)); 
		    // if the actual entries did not exist, need to fill
		    // entropy with n*ASize
		    mean1 = (mean1==0.0)?part*log(data[j]->ASize()):mean1;
		    mean2 = (mean2==0.0)?part*log(data[j]->ASize()):mean2;
		    // and now get the b from these means
		    b = GSL_MIN(b,mean1+mean2+3.0*sqrt(gsl_pow_2(std1)+gsl_pow_2(std2)));
		  }//for loop over different depths
		  
		  
		  // if we are doing serial recording, then we can use
		  // the fact that S(N) is a convex function to create
		  // possibly tighter lower and upper bounds on S
		  // we  look at depths of A and B (see notebook) and
		  // build linear bounds for entropy
		  
		  //one must be careful here when applying bounds to the
		  //cases when the estimator underestimates!!!
		  for(int A=1;A<depth;A++){
		    double mean1=0.0; // means and stds
		    double mean2=0.0;
		    double std1=0.0;
		    double std2=0.0;
		    // these bounds are very tight, so I put them up
		    // to 4 standard
		    // deviations instead of 3 above
		    int ind1 = (A-results.SDepth())*results.Stride()+j;
		    if(ind1>=0){	// this object exists
		      if(results[ind1].size()==1){ // serial, not
			// parallel recording,
			// and this bound is
			// valid only for serial
			// getting mean snd std
			mean1= (results[ind1][slice])->get_Snsb();
			std1 = (results[ind1][slice])->get_dSnsb();
			for(int B=1;B<A;B++){
			  int ind2 = (B-results.SDepth())*results.Stride()+j;
			  if(ind2>=0){ // if this object exists
			    // getting mean and std
			    mean2= (results[ind2][slice])->get_Snsb();
			    std2 = (results[ind2][slice])->get_dSnsb();
			    const double coef1=((double)(depth-B))/(A-B);
			    const double coef2=((double)(depth-A))/(B-A);
			    b= GSL_MIN(b,  mean1*coef1 + mean2*coef2 +
				       4.0*sqrt(gsl_pow_2(std1*coef1) + 
						gsl_pow_2(std2*coef2)));
			  } // if ind2 exists
			}   // for over ind 2
		      }     // if ind1 exists
		    }	// check if many slices, or just one slice
		  }		// for over ind1
		
		  pr = new Prior(a,b);
		  
		}
		else		
		  error("Undefined prior requested.");
		
		if(depth>0){
		  results[ind].add(*cts, pr, pa, 1e-6); // if we are finishing,
		  // counts for large length
		  // may not exist (NEvents==0)
		}			// then add() will ignore them
		else{
		  // prior needs to be erased only when we did not count;
		  // otherwise prior gets attached to an nsb object and gets
		  // erased with it
		  delete pr;
		  pr=NULL;
		}
	      } // end of doing nsb calculations
	    }	// end of checking if the counts are not empty
	  } // end of checking if counts need to be evaluated
	} // end of doing loop over depths
      }	// end of checking if need to work with the given phase
    } // end of doing the loop over phases
  }
  




  /*
   * The function does one pass through the data file and records the
   * data from it into output files consistently with the required
   * depth's of recording (from start to end) and required
   * implementation of EntrData structures.
   */
  void ProcessingCycle(const std::string& fname, unsigned long skip, 
		       unsigned long nsymb, IMPLEDATA impl, 
		       int start_depth, int end_depth, STARTPOINT st, 
		       DOCALCS ca, DOPART pa, double prob,
		       int VEDPsize, std::vector<std::string>& names, 
		       const std::string& prior,
		       EntrDict* di, Recorder* rec, VNsbSet& results){

    unsigned long have_read=0;	// how many symbols did we actually read
    unsigned long old_read=0;	// # read before last Stream2... call
    unsigned long have_rec=0;	// how many symbols did we actually recorded
    unsigned long old_rec=0;	// recorded at a previous call

    VEDP vedp(VEDPsize);	//vector of data structures

    if(end_depth<start_depth) error("First n-gramm is longer than the last one.");


    //opening file    
    std::string infname = fname + ".txt";
    std::ifstream in(infname.c_str());
    if(!in) error("Cannot open input file", infname.c_str());	
    if(end_depth!= (int)rec->get_depth()) 
      error("Recorder depth and requested depth mismatched");

    //skipping the first symbols
    try{
      di->Start(in, skip);
      have_read=skip;
    }	// preparing input file
    //
    // for all processsing (including parallel) we skip the first
    // "skip" symbols, but not the first "skip" repetitions
    //
    catch (EntrData_failed_transl){ //if skipping failed, return
      return;
    }
    
    int last_out=-1;		// counter used to output every
				// 10000'th recording for vlevel=VERBMIN
    // tag for saved file names
    std::string tag= ((ca==YESNSB)? prior + "_" : std::string("")) +
      (((nsymb==0) && (skip==0)) ? "" : std::string("k")+skip+"n"+nsymb +"_") +
      ((prob>=1.0) ? "" : std::string("p")+prob+"_");
    do{
      old_read=have_read;
      old_rec=have_rec;
      //creating data structures
      ChooseData(impl, names, vedp, di->ASize());
      //reading at end_depth using required WRITELEVEL
      Stream2Dict2Data(vedp, in, di, rec, skip, nsymb, end_depth, st, prob,
		       have_read, have_rec); // after this function,
					     // we have either
					     // exhausted the whole
					     // instream, or have
					     // done one row for
					     // parallel recording
      
      //counting and writing; this is not needed if this is the last
      //run through the cycle, and no new recordings were done 
      if(old_rec<have_rec) {
	if ((vlevel>VERBMIN)||((int)(have_read/10000)>last_out)){
	  last_out=have_read/10000;
	  std::cout << "Read "<< have_read << " symbols. ";
	  std::cout << "Recorded "<< have_rec << " translations.\n";
	}
	Data2Cts2File(vedp, fname, tag, start_depth, end_depth, st, ca, pa, prior, results);
      }
      
      KillData(vedp);
    }while((old_read<have_read)||(old_rec<have_rec)); // read while
				// you can still can read or record
				// something 

    if(ca==YESNSB)  results.Print(fname+"_"+tag,names);
  }



  /*
   * Gets all the control switches, calls correct dictionary and data
   * structure creation utilities, then parses data files in a correct
   * way, writes data, and cleans after itself.
   */
  void ProcessFile(const std::string& fname, unsigned long skip, unsigned long nsymb, 
		   IMPLEDATA impl, const std::string& dict, const std::string& prior,
		   int start_depth, int end_depth, 
		   WRITELEVEL all, int shift, STARTPOINT st, DOCALCS ca, DOPART pa, double prob){

    EntrDict* di;		//pointer to the dictionary
    Recorder* rec;              //pointer to the recorder


    // all necessary checks
    if ((start_depth ==0)&&(shift==0))
      error("ERROR: Don't know how to shift words of length 0 by their length.");

    // a "matrix" of results (long sequences of entropy
    // calculations). The element, which corresponds to the depth of
    // i and the startpoint j is accessed by [i*maxVEDPsize+j]
    int maxVEDPsize = (shift>0)? shift : end_depth+abs(shift); // check based on shift
    maxVEDPsize = (prob<1.0)? (int)round(1.0/prob) : maxVEDPsize; // check based on prob
    if(st>1) maxVEDPsize=st;				     // if
				// specific # of records requested, choose it
    BXIset interp_all(12);
    VNsbSet results(interp_all,maxVEDPsize*(end_depth-start_depth+1),
		    maxVEDPsize,0,start_depth,end_depth);



    try{
      for(int depth = (all==WRITEALL?end_depth:start_depth);depth<=end_depth;depth++){
	if(vlevel>=VERBMIN)
	  std::cout<< "Processing word lengths " << 
	    (all==WRITEALL?  (std::string("") + start_depth + std::string(" to ") + 
			      end_depth +".\n") : 
	     std::string("")+depth+". \n");

	ChooseDictRec(dict, fname, depth, all, di, rec);
	results.reserve_all(rec->to_record());
       
	// how many recordings will we have at the same time?
	int VEDPsize = (shift>0)? shift : depth+abs(shift); // check based on shift
	VEDPsize = (prob<1.0)? (int)round(1.0/prob) : VEDPsize;  // check based on prob
	if(st>1) VEDPsize=st;				    // if
				// specific # of records requested, choose it
	
	std::vector<std::string> names(maxVEDPsize); // names of data structures
	for (int i=0; i<VEDPsize; i++){
	  // data structures (and outputs too) will be labeled my the 
	  // name of the dictionary and their shift
	  names[i] =  std::string(di->GetName()) + "_mf" +VEDPsize+"f" + i;
	}

	ProcessingCycle(fname, skip, nsymb, impl, (all==WRITEALL?start_depth:depth),
			(all==WRITEALL?end_depth:depth), st, ca, pa, prob, 
			VEDPsize, names, prior, di, rec, results);
	KillDictRec(di, rec);
      }
    }
    catch(std::bad_alloc){
      throw EntrData_bad_alloc();
    }
  }
  

  
  /*
   *  Choosing the dictionary and the recorder. Gets dictionary type 
   *  and the working file name in, creates the appropriate dictionary
   *  and recorder, calculates needed number of data structures, and 
   *  creates their names.
   *
   *  No exceptions check done.
   */
  void ChooseDictRec(const std::string& dict, const std::string& fname, 
		     int depth, WRITELEVEL all,
		     EntrDict*& di, Recorder*& rec){

    typedef std::map<std::string, DictParams> MSDP;
    typedef MSDP::iterator MSDPI;
    
    MSDP which;			// properties of differenet dictionaries
    which[M2ONames[BINARY]].type  = M2O; // populating the description of different
    which[M2ONames[BINARY]].param = BINARY; // dictionaries with real data
    which[M2ONames[TETRA]].type  = M2O;	
    which[M2ONames[TETRA]].param = TETRA;	
    which[M2ONames[OCTA]].type  = M2O;	
    which[M2ONames[OCTA]].param = OCTA;	
    which[M2ONames[GENOME]].type  = M2O;
    which[M2ONames[GENOME]].param = GENOME;
    which[M2ONames[ENGL3]].type = M2O;
    which[M2ONames[ENGL3]].param = ENGL3;
    which[M2ONames[ENGL29]].type  = M2O;
    which[M2ONames[ENGL29]].param = ENGL29;
    which[M2ONames[ENGL38]].type  = M2O;
    which[M2ONames[ENGL38]].param = ENGL38;
    which[ParName].type  = PAR;
    which[ParName].param = 0;
    which[NumName].type  = NUM;
    which[ParName].param = 0;
    
    MSDPI current = which.find(dict); // currently requested dictionary
				// do not use subscripting, as it may
				// insert items into the map
    if(current == which.end())	// did not find correct key
      error("Incorrect dictionary type specified: ", dict.c_str());
    


    switch (current->second.type){
    case M2O:
      di = new Dict_m2o( (DIFFM2O) current->second.param); // creating a new dictionary
				// with correct alphabet
      //for M2O dictionaries, the simple recorder is needed
      rec = new Simple_Recorder(depth,all);
      break;
    case PAR:
      di = new Dict_numb(fname);
      rec = new Parallel_Recorder(depth,fname,all);
      break;
    case NUM:
      di = new Dict_numb(fname);
      rec = new Simple_Recorder(depth,all);
      break;
    default:			// we should never get to this point
      error("This error must never happen!");
      break;
    };
    

  }

  /*
   *  Removing dictionary-related allocated memory 
   */
  void KillDictRec(EntrDict* di, Recorder* rec){
    delete(di); di=NULL;
    delete(rec); rec=NULL;
  }
  
  

  /*
   *  Creating data related structures. We choose which implementation
   *  of EntrData to utilize, create the needed vector of pointers to
   *  such implementations, and finally create the implementatioins
   *  themselves with required names and alphabet sizes
   *
   * No exceptions check is done by this routine.
   */
  void ChooseData(const IMPLEDATA impl, 
		  const std::vector<std::string>& names, VEDP& vedp, int asize){
    
    switch (impl){
    case EDATANODE:
      for(int i=0; i<(int)vedp.size();i++)
	vedp[i] = new node(asize, names[i]);
      break;
      
    default:			//case EDATAMPI
      for(int i=0; i<(int)vedp.size();i++)
	vedp[i] = new mpi_map(asize, names[i]);
      break;
    }
  }

  
  /*
   *  Removing data-related allocated memory 
   */
  void KillData(VEDP& vedp){
    for(unsigned int i =0; i< (vedp.size());i++){
      delete vedp[i];
      vedp[i]=NULL;
    }
  }
}

