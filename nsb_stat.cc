/*
 *  File - nsb_stat.cc
 *  
 *  Main file that creates nsb_stat program for parsing of the data files.
 *
 *  Created by nemenman on Tue Nov 06 2001.
 *  Copyright (c) 2001-2013 Ilya Nemenman. All rights reserved.
 *
 */


#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <unistd.h>

#include "Entropy.h"

using namespace ENTRDATA;

VERBOSITY vlevel;

/*
 * String constants, command line arguments
 */

const char use_mpi  = 'm';	// maps
const char use_node = 't';	// trees
const char use_wall = 'a';	// all levels at a time writing
const char use_wone = 'o';	// one level at a time writing
const char phase    = 'u';	// record mUltiple phases, or just the zero phase (Y/N) 
const char use_prior= 'i';	// prior
const char calcs    = 'c'; 	// calculate or not 
const char dname    = 'd';      // dictionary name
const char sdepth   = 's';      // starting depth (word length)
const char edepth   = 'e';      // ending depth 
const char sk       = 'k';	// first symbols to skip
const char ns       = 'n';	// symbold to translate
const char sft      = 'f';      // shift
const char part     = 'r';	// partitioning counts
const char pr       = 'p';	// probability to accept a word
const char vrb      = 'v';	// verbosity
const char hlp      = 'h';      // help
const char* optstring = "mtaou:i:c:r:d:s:e:k:n:f:p:v:h";
const std::string YES("Y");  
const std::string NO("N");

void usage(){
  std::cout<<"usage:  nsb-entropy [-oamt] [-dDICT -iPRIOR ";
  std::cout<<"-c0/1/2 -rY/N -sSTART -eEND -kSKIP ";
  std::cout<<"-nNUMBER -fSHIFT -pPROB -uY/N/# -v0/1/2] file\n";
}

void help(){
  std::cout<<"The program counts occurences of n-gramms in a text file\n";
  std::cout<<"and estimates entropies using NSB method if needed.\n";
  usage();
  std::cout<<"\nOptions controlling the ouput:\n";
  std::cout<<"  DICT   - dictionary to use; the following values are supported now:\n";
  std::cout<<"       2   - binary (every symbol that is not 0 is treated as 1);\n";
  std::cout<<"       4   - tetrary (0, 1, 2, and everything else is 3);\n";
  std::cout<<"       8   - octary (0 through 6, and everything else is 7);\n";
  std::cout<<"       e3  - English, 3 symbols (vowels, consonants, and everything else);\n";
  std::cout<<"       e29 - (default) English, 29 symbols (26 letters, one symbol for digits,\n"; 
  std::cout<<"            one for space, and one for everything else);\n";
  std::cout<<"       e38 - English, 38 symbols (26 letters, 10 digits, space, and everything\n";
  std::cout<<"             else).\n";
  std::cout<<"       gene- genetic alphabet (ACGT, or their small versions).\n";
  std::cout<<"       num - stream consists of integers separated by spaces.\n";
  std::cout<<"       par - parallel stream processing with an integer number dictionary.\n";
  std::cout<<"  PRIOR  - encodes a priori information about the entropy:\n";
  std::cout<<"       uni - (default) uniform prior;\n";
  std::cout<<"       bnd - prior for an n-gramm is obtained by looking at all combinations of\n";
  std::cout<<"             estimates for n-i and i-grams (if these available) and limiting the\n";
  std::cout<<"             entropy to (variance-enlarged) windows allowed by these combinations.\n";
  std::cout<<"  c0/1/2 - (default 0) which calculation is to be performed: 0 - full NSB entropy\n";
  std::cout<<"           analysis, 1 - counting the frequency of occurence of sample occurences\n";
  std::cout<<"           2 - counting the frequency of samples.\n";
  std::cout<<"           If DICT=par, this will be set to 0 and any other input will be ignored.\n";
  std::cout<<"  rY/N   - (default N) whether to pArtition the bins into (a) bins that each have a\n";
  std::cout<<"           substantial fraction of all counts, (b) well sampled bins (where each\n";
  std::cout<<"           bin, though well sampled, has small fraction of all counts), and\n";
  std::cout<<"           bins with bad sampling. If this is set to Y, the prior is automatically\n";
  std::cout<<"           reset to uni.\n";
  std::cout<<"  START  - (default 0) starting n-gramm length, nonnegative integer.\n";
  std::cout<<"  END    - (default 100) ending n-gramm length, nonnegative integer.\n";
  std::cout<<"  SKIP   - (default 0) number of symbols to skip in the beginning of the file.\n";
  std::cout<<"           nonnegative integer.\n";
  std::cout<<"  NUMBER - (default 0=all) number of symbols to analyze.\n";
  std::cout<<"           nonnegative integer.\n";
  std::cout<<"  SHIFT  - (default 1) specifies the shift between successive words. For\n";
  std::cout<<"           positive SHIFT, starting position of the next word is  determined\n";
  std::cout<<"           by adding SHIFT to the starting position of the previous one. For\n";
  std::cout<<"           negative or zero SHIFT, next starting position is calculated as the\n";
  std::cout<<"           end of the previous word, plus one, plus abs(SHIFT). If PROB<1.0 (see\n";
  std::cout<<"           below), then the requested value of SHIFT is disregarded, and SHIFT=1.\n";
  std::cout<<"  PROB   - (default 1.0) a real number between 0 and 1 that specifies the \n";
  std::cout<<"           probability of acceptance of each word for recording.\n";
  std::cout<<"  uY/N/# - (default N) mUltiple or single counting. This parameter is used only\n";
  std::cout<<"           for SHIFT!=1 or PROB!=1.0. If N, then for non-trivial SHIFT, the\n";
  std::cout<<"           program record only the words with starting positions of N*SHIFT+1.\n";
  std::cout<<"           For non-trivial PROB, only one random set of samples is recorded.\n";
  std::cout<<"           If Y, then for nontrivial SHIFT, we  make SHIFT different records, and\n";
  std::cout<<"           words with their starts at N*SHIFT+i will be in those records. For \n";
  std::cout<<"           non-trivial PROB, we subsample round(1/PROB) times from the data file.\n";
  std::cout<<"           Finally, for #, we make # different records instead of SHIFT or 1/PROB.\n";
  std::cout<<"  v0/1/2 - (default 2) verbosity level, from lowest to highest.\n";
  std::cout<<"\nOptions controlling execution:\n";
  std::cout<<"  m or t - (default m) specifies if Maps or Trees are used to handle the data\n";
  std::cout<<"           internally; `t' is usually slightly faster, but consumes a lot more\n";
  std::cout<<"           memory, and thus should be used for short n-gramms or data sets only;\n";
  std::cout<<"  o or a - (default o) Specifies if n-gramms of different size are processed\n";
  std::cout<<"           One-by-one, or All together; 'a' requires just one pass through the\n";
  std::cout<<"           data file and is faster, but consumes a lot more memory, and thus\n";
  std::cout<<"           should be used for small n-gramm ranges only. All-level writing \n";
  std::cout<<"           should only be used if SHIFT<=START. In other regimes selection of\n";
  std::cout<<"           samples for counting is ambiguos, and the results are unpredictable.\n";
  std::cout<<"           Further, for PRIOR=bnd, it is essential to have this option to 'o',\n";
  std::cout<<"           otherwise it is equivalent to PRIOR=uni.\n"; 
  std::cout<<"\nNote that filename argument `file' must not have an extension. Text format is\n";
  std::cout<<"assumed for it, and the `.txt' extension will be appended internally.\n\n"; 
}


/*
 * Main
 */
int main(int argc, char* argv[]){

  // default values for the arguments
  IMPLEDATA  impl = EDATAMPI;
  WRITELEVEL all  = WRITEONE;
  STARTPOINT st   = ONESTART;
  DOCALCS    ca   = YESNSB;
  DOPART     pa   = NOPART;
  VERBOSITY  verb = VERBMED;
  std::string dict = "e29";
  std::string prior = "uni";
  int start_depth = 0;
  int end_depth = 100;
  unsigned long skip = 0;
  unsigned long nsymb = 0;
  int shift = 1;
  std::string fname="";
  double prob = 1.0;

  int curr_opt;			// current option
  long tmp_opt;			// current option argument

  // My logo output    
  std::cout << "This is " << ENTRDATA_VER(EDVERALL) <<"\n";

  while((curr_opt = getopt(argc, argv, optstring))!= -1){
    switch (curr_opt){
    case part:
      pa = ((optarg==YES)? YESPART: NOPART);
      break;
    case calcs:
      {std::istringstream in(optarg);
	int x;
	in >> x;
	ca=(DOCALCS) x;
      }
      if(ca<YESNSB) ca=YESNSB;
      if(ca>NSBNONE) ca=NSBNONE;
      break;
    case dname:
      dict = optarg;
      if(dict==ParName)
	impl = EDATAMPI;
      break;
    case vrb:
      {std::istringstream in(optarg);
	int x;
	in >> x;
	verb=(VERBOSITY) x;
      }
      if(verb<VERBMIN) verb=VERBMIN;
      if(verb>VERBMAX) verb=VERBMAX;
      vlevel=verb;
      break;
    case use_prior:
      prior = optarg;
      break;
    case use_mpi:
      impl = EDATAMPI;
      break;
    case use_node:
      impl = EDATANODE;
      break;
    case use_wall:
      all = WRITEALL;
      break;
    case use_wone:
      all = WRITEONE;
      break;
    case phase:
      if(optarg==YES)
	st=MANYSTARTS;
      else if (optarg==NO)
	st=ONESTART;
      else
	{std::istringstream in(optarg);
	  in >> st;
	}
      break;
    case sdepth: 
      {std::istringstream in(optarg);
	in >> start_depth;}
      // no negative scan depths
      start_depth = (start_depth>0)? start_depth : 0;
      break;
    case edepth:
      {std::istringstream in(optarg);
      in >> end_depth;}      
      // no negative scan depths
      end_depth = (end_depth>0)? end_depth : 0;
      break;
    case ns:
      {std::istringstream in(optarg);
      in >> tmp_opt;}
      // no negative scan depths
      nsymb = (tmp_opt>0)? (unsigned long)tmp_opt : 0;
      break;
    case sk:
      {std::istringstream in(optarg);      
      in >> tmp_opt;}
      // no negative scan depths
      skip = (tmp_opt>0)? (unsigned long)tmp_opt : 0;
      break;
    case sft:
      {std::istringstream in(optarg);
      in >> shift;}
      break;
    case pr:
      {std::istringstream in(optarg);
      in >> prob;}
      if((prob<=0.0)||(prob>1.0)) // is this a valid probability?
	prob = 1.0;
      break;
    case hlp:
      help();
      exit(0);
      break;
    default:
      usage();
      exit(0);
    }
  }
  // overwriting input values as required
  if(dict==ParName)  {ca=YESNSB;}
  if(pa==YESPART) {prior = "uni";}

  if(optind != (argc-1)){	// file name must be in the last argv; if it is not,
    usage();			// some options are bad
    exit(0);
  }

  fname = argv[optind];		// it must be the file name

  if(prob<1.0)			// nonzero value of shift for
    shift = 1;	 // rarification of data is possible only when prob is
				// not used for the same (i.e.,
				// PROB!=0)

  std::cout<< "Executing with parameters:\n";
  std::cout<< "\tData file name: " << fname <<".txt\n";
  std::cout<< "\t- " << ((skip)? (std::string("will be skipped till symbol ") + skip + " and then read"):
		    "will be read from the start") << " ";
  std::cout<<((!nsymb) ? "to the end" :
		   (std::string("for the next ") + nsymb + " symbols")) << ";\n";
  std::cout<< "\t- will be scanned at word length from " << 
    start_depth << " to " <<end_depth<<";\n";
  std::cout<< "\t- recorderd with ";
  if(st == MANYSTARTS) std::cout<<"all possible";
  else std::cout<<st;
  std::cout<< " start point(s);\n";
  std::cout<< "\t- using the dictionary '" << dict <<"';\n";
  std::cout<< "\t- using shift of " << shift <<";\n"; 
  std::cout<< "\t- with probability of accepting words " << prob << ".\n";
  std::cout<< "\tThe data will " << ((ca==YESNSB)? "be NSB analyzed" : "") << 
    ((ca==YESCOUNT)? "be counted in a compact form" : "") <<
    ((ca==NSBNONE)? "be counted in a full form" : "") << " with\n";
  std::cout<< "\t- prior set to '" << prior <<"';\n";
  std::cout<< "\t- partitioning of data will " << ((pa==YESPART)? "" : "not ") << "be performed.\n";
  std::cout<< "\tInternal settings:\n";
  std::cout<< "\t- using " << ((impl==EDATANODE)? "tree" : "map") << " implementation;\n";
  std::cout<< "\t- using " << ((all ==WRITEONE)? "one" : "multiple") << " level writing.\n\n";

  if (all && (shift!=1)){
    std::cout <<"WARNING: Should not be using multiple level writing with SHIFT=" <<shift<<".\n";
  }


  // we do no error check here for proper file name 
  // this will be checked by working routines later
  try{
    ProcessFile(fname, skip, nsymb, impl, dict, prior,
		start_depth, end_depth, all, shift, st, ca, pa, prob);
    }
    catch(std::range_error& e){
      error(e.what());
    }
    catch(EntrDataException& e){
      error(e.what());
    }    
    catch(...){
      error("Unknown error has occured.");
    }
}

