#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include "../include/random.h"

#include "../include/utility.h"
int main () {
  std::string fname("");
  unsigned long maxrnd;
  std::cout<<"Enter randoms name: ";
  std::cin>>fname;
  rnd_gen_class rgen(fname);
  std::cout<<"Enter # of randoms to generate: ";
  std::cin>>maxrnd;

  for(unsigned long i =0; i<maxrnd; i++){
    std::cout<<rgen.uniform()<<"\t"<<rgen.how_many()<<"\t"<<rgen.get_state()<<"\n";
  }
  return 0;
}
