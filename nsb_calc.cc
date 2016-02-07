/*
 *  File - nsb_calc.cc
 *  Project - Entropy_parsing
 *
 *  Estimating entropies from counts, doing actual numerics.
 *
 *  Copyright (c) 2003 Ilya Nemenman. All rights reserved.
 *
 */


#include "nsb.h"
#include "nsb_aux.h"

// nsb computational warnings names
const char* ENTRDATA::warn_names[] = 
  { "All OK. Recovered.",
    "No coincidences.", 
    "Number of events is <= 1.0. ",
    "All data coincide.",
    "Series expansion for B0 did not converge. Will probably recover.",
    "Newton-Raphson search for B0 did not converge. Will probably recover." ,
    "Newton-Raphson ERROR in B0 calculation: wrong sign. May recover.",
    "Newton-Raphson ERROR in Bcl calculation: wrong sign. Possibly serious.",
    "Newton-Raphson search for Bcl did not converge.",
    "Numerical evaluation of an integral did not converge.",
    "Over 100 coincidences. Consider using UNIform prior."
  };

// nsb computational errors names 
const char* ENTRDATA::badnum_names[] = 
  { "Improper call order; internal data uninitialized. ", 
    "GSL error. ",
    "Requested entropy larger than log(K) or smaller than 0. ",
    "Bad parameter ranges."};


/*********************************************
  calculating nx and kx for samples with nx>1;
  calculating K1 and K2
 */
void ENTRDATA::NsbCalc::make_ng1()
  throw(EntrData_badnum,EntrData_bad_alloc) {
  if (round(N)<=1.0) warning(NSB_NOEVENTS);

  try{
    int* tmp;			// variable to store nx>1.0
    tmp = new int[sz];		// finding those nx>1.0 (use 1.05 to
    for(unsigned long i=0;i<sz;i++){	// make sure that reals don't cause
      tmp[i] = (nx[i]>1.05);	// problem)
    };			
    
    szng1= vec_sum(sz,tmp);	// # of elements>0; this may be 0 for
				// # no coincidences
    // a safety check; remove all ng1 arrays;
    // remember that calling delete to NULL pointer does nothing so
    // we don't check for this (the pointers are initialized to
    // NULL by the constructor)
    delete[] nxng1;
    delete[] kxng1;
    if(szng1){			// if size>0
      nxng1 = new double[szng1];
      kxng1 = new double[szng1];
      
      unsigned long j =0;
      for(unsigned long i =0; i<sz;i++){
	if(tmp[i]){
	  nxng1[j]=nx[i];
	kxng1[j]=kx[i];
	j++;
	};
      };
    };
    delete[] tmp;		// getting rid of tmp
  } catch (std::bad_alloc){
    throw EntrData_bad_alloc();
  };

  K1= vec_sum(sz,kx);		// private variables in the NsbCalc class
  K2= vec_sum(szng1,kxng1);
}    



/*********************************************
   mean value of the posterior entropy for given counts
   and a given value of \beta (the input parameter B)
*/
double ENTRDATA::NsbCalc::meanS(const double& B) const{
  // remember that no nx==0 checking is needed, since these bins are
  // not copied in the constructor; and the other checks are done in
  // make_ng1()

  if(gsl_isinf(B)) return maxent; // infinite Beta -> maximum entropy
  const double ovrNB = 1/(N+B);
  const double BoK = B/K;
  double prod= tri_dot_f_offset2(sz, nx, nx, kx, 
				 BoK+1.0, BoK, &gsl_sf_psi);
  // sum is written in a form to avoid the lossof precision as K->Inf
  return psi(N+B+1.0) - ovrNB*prod - 
    B*ovrNB*(1.0-K1/K)*psi(BoK+1.0);
}



/*********************************************
   mean value of the posterior entropy squared for given counts
   and a given value of \beta (the input parameter B)
*/
double ENTRDATA::NsbCalc::meanS2(const double& B) const
  throw(EntrData_bad_alloc){
  // remember that no nx==0 checking is needed

  if(gsl_isinf(B)) return gsl_pow_2(maxent);
  try{
    double f =0.0;		// the variable to return later

    // temporary variables
    const double BoK = B/K;
    const double p0NB2 = psi(N+B+2.0);
    const double p1NB2 = psi1(N+B+2.0);
    const double pb1 = psi(BoK + 1.0) - p0NB2;
    // temporary arrays
    double* pnxb1;		// defined differently from Octave's code
    double* nxb;
    pnxb1 = new double[sz];
    nxb   = new double[sz];
    for (unsigned long i=0;i<sz;i++){
      nxb[i] = nx[i]+BoK;
      pnxb1[i] = psi(nxb[i] + 1.0) - p0NB2;
    };

    // i, j term, summing over all i and j (including i==j terms,
    // thus overcounting, and then correcting for it)  

    // ni*nj ~= 0 contribution
    for(unsigned long i=0;i<sz;i++)
      for(unsigned long j=0;j<sz;j++)
	f += nxb[i]*pnxb1[i]*kx[i] * nxb[j]*pnxb1[j]*kx[j] -
	  nxb[i]*kx[i] * nxb[j]*kx[j] * p1NB2; 
    //ni*b contribution 
    for(unsigned long i=0;i<sz;i++)
      f += 2.0*B*(1-K1/K)* nxb[i]*(pnxb1[i]*pb1 - p1NB2) * kx[i];
    // b*b contribution 
    f += (1-K1/K)*(1-(K1+1.0)/K)*B*B*(pb1*pb1-p1NB2);
    // correcting for overcounting
    for(unsigned long i=0;i<sz;i++)
      f -= gsl_pow_2(nxb[i]*pnxb1[i])*kx[i] - nxb[i]*nxb[i]*kx[i]*p1NB2;

    //-----------------------------------------------------
    // i term
    
    // ni contribution
    for(unsigned long i=0;i<sz;i++)
      f += (nxb[i]*(nxb[i]+1.0) * 
	    (gsl_pow_2(psi(nxb[i]+2.0) - p0NB2) +	    
	     psi1(nxb[i]+2.0) - p1NB2))*kx[i];

    // b contribution
    f += B*(1-K1/K)*(1.0+BoK) * (gsl_pow_2(psi(2.0+BoK)-p0NB2)
				   + psi1(BoK+2.0) - p1NB2);

    //-----------------------------------------------------
    // normalizing
    f /= ((N+B)*(N+B+1.0));
    
    delete[] pnxb1;
    delete[] nxb;

    return f;
  }
  catch (std::bad_alloc){
    throw EntrData_bad_alloc();
  }
}




/*********************************************
  Computes the "action" (the negative logarithms of the 
  evidence) for the integral over \xi (see NSB method for 
  calculating entropies of discrete pdf's). Does not include
  the contribution from the prior over \xi. Note that even though 
  the integration variable is \xi, the argument of this function 
  is B.
*/
double ENTRDATA::NsbCalc::mlog_evidence(const double& B) const 
  throw(EntrData_badnum, EntrData_bad_alloc){

  double f=0.0;

  if(B<=0.0) 
    return std::numeric_limits<double>::infinity(); // infinite evidence for negative B
  if(gsl_isinf(B)) 
    return std::numeric_limits<double>::infinity(); // or for infinite B

  const double BoK= B/K;
  if(szng1)			// if there are coincidences
    f +=  -dot_f_offset(szng1, nxng1, kxng1, BoK, &gsl_sf_lngamma)
      + K2*gsl_sf_lngamma(1+BoK);
    
  // Need to calcuate f += -K1*log(B) +gammaln(B+N) - gammaln(B)
  // But to avoid (big number) - (big number) = lost precision
  // problem, need to treat different regimes of N and B
  //differently
    
  // First regime, aymptotically large B and B/N: B>100 and N<0.01*B;
  // here we can expand gammaln-gammaln = psi*N + psi_1/2*N^2+...
  const int large = (B> GSL_MAX(100.0, 100.0*N));

  // polygamma(n,B) ~B^(-n); thus we expand in (N/B);
  // which of expansion parameters is the worst? how many 
  // series terms will we need? as seen below, the leading
  // term in f for N==K1 (worst case) is psi_asymp(B)/N ~ N/B.
  // we need 10^(-15) precision relative to that term. Note
  // that the series expansion has the form 
  // f = leading + psi_1/2!*N^2 +psi_2/3!*N^3 +... =
  //   = leading + (N/B +(N/B)^2+...)
  if(large){
    const int nterms  = (int)ceil(fabs((-15.0 -log10(N))/log10(N/B))) + 1;
    double* ifac;
    try{
      ifac = new double[nterms];
      fact_frac(nterms, 2.0, ifac); // populating ifac with 1/factorial
      f += series_f(nterms, ifac, N, 2, &gsl_sf_psi_n, B, 1);
      delete[] ifac;
    }
    catch(std::bad_alloc){
      throw EntrData_bad_alloc();
    };
    f += (N-K1)*log(B) + psi_asymp(B)*N;
  }
  else {			// no asymptotic expansion needed
    f += - K1*log(B) + gsl_sf_lngamma(B+N) - gsl_sf_lngamma(B);
  };

  return f;
}




/*********************************************
  The function finds the position of the minimum of the a posteriori
  evidence and the variance around it. The integration variable is 
  \xi -- the a-priori value of the entropy.
   Used data:
     data - map of "event occured data.first times" -> "data.second
            events like that" (nx -> kx)
   Changed on output:
     Bcl     - the classical value of B;
     xicl    - the classical value of xi;
     dxicl   - std. dev. near the classical value;
   Returned:
     error code; 0 - all ok;
                 1 - no coincidences; saddle point
                     evaluation invalid (wide variance);
                 2 - all data coincides; saddle point
                     evaluation invalid - Bcl close to zero
                 3 - no convergence in Newton-Raphson root
                     finding of Bcl;
*/
ENTRDATA::NSB_WARN ENTRDATA::NsbCalc::max_evidence() {
  if (round(K1)==round(N)){	// no coincidences
    Bcl =  std::numeric_limits<double>::infinity();
    xicl =  std::numeric_limits<double>::infinity();
    warning(NSB_NOCOINC);
    dxicl = GSL_NAN;
  }
  else if (round(K1) == 1){	// all data coincides
    Bcl = 0.0;
    xicl= 0.0;
    warning(NSB_ALLCOINC);
    dxicl = GSL_NAN;
  }
  else {		        // some non-trivial value of B0, Bcl has to be
				// calcuated 
    double B0= 0.0;
    // summing the series
    { 
      const int order = 10;	// we will calculate B to this order in ep
      const double N2 = gsl_pow_2(N);
      const double N3 = gsl_pow_3(N);
      const double N4 = gsl_pow_4(N);
      const double N5 = gsl_pow_5(N);
      const double N6 = gsl_pow_6(N);
      const double N7 = gsl_pow_7(N);
      const double N8 = gsl_pow_8(N);
      const double N9 = gsl_pow_9(N);
      const double N10 = gsl_pow_int(N,10);
      const double N11 = gsl_pow_int(N,11);
      
      const double ovrN = 1/N;
      const double Nm   = N-1.0;
      const double Nm2  = gsl_pow_2(Nm);
      const double Nm3  = gsl_pow_3(Nm);
      const double Nm4  = gsl_pow_4(Nm);
      const double Nm5  = gsl_pow_5(Nm);
      const double Nm6  = gsl_pow_6(Nm);
      const double Nm7  = gsl_pow_7(Nm);
      const double Nm8  = gsl_pow_8(Nm);
      const double Nm9  = gsl_pow_9(Nm);
      const double Nm10 = gsl_pow_int(Nm,10);
      
      const double b[] = {    // coefficients of the expansion of B in
				// powers of \epsilon, calculated by
				// Mathematica
	(-1.0 + N)/(2.0*N) ,	// b(-1)
	(-2.0 + ovrN)/3.0 ,	// b(0)
	(2.0 + N - N2)/(9.0*N - 9.0*N2) , // b(1)
	(2.0*(2.0 - 3.0*N - 3.0*N2 + 2.0*N3))/(135.0*Nm2*N) , // b(2)
	(4.0*(22.0 + 13.0*N - 12*N2 - 2.0*N3 + N4))/(405.0*Nm3*N) , // b(3)
	(4.0*(-40.0 + 58.0*N + 65.0*N2 - 40.0*N3 - 5.0*N4 +	
	      2.0*N5))/(1701.0*Nm4*N) , // b(4)
	(4.0*(-9496.0 - 6912.0*N + 5772.0*N2 + 2251.0*N3 - 
	      1053.0*N4 - 87.0*N5 + 29.0*N6))/(42525.0*Nm5*N) , // b(5)
	(16.0*(764.0 - 1030.0*N - 1434.0*N2 + 757.0*N3 + 295.0*N4 
	       - 111.0*N5 - 7.0*N6 + 2.0*N7))/(18225.0*Nm6*N) , // b(6)
	(16.0*(167000.0 + 142516.0*N - 108124.0*N2 - 66284.0*N3 
	       + 26921.0*N4 + 7384.0*N5 - 2326.0*N6 - 116.0*N7 
	       + 29.0*N8))/(382725*Nm7*N) , // b(7)
	(16.0*(-17886224.0 + 22513608.0*N + 37376676.0*N2 
	       - 17041380.0*N3 - 11384883.0*N4 + 3698262.0*N5 
	       + 846930.0*N6 - 229464.0*N7 - 9387.0*N8 + 2086*N9))
	/(37889775.0*Nm8*N), 	// b(8)
	(16.0*(-4166651072.0 - 3997913072.0*N + 2783482560.0*N2 +
	       2290151964.0*N3 - 803439834.0*N4 - 395614251.0*N5 +
	       108055443.0*N6 + 20215218.0*N7 - 4805712.0*N8 - 165395.0*N9
	       + 33079.0*N10))/(795685275*Nm9*N), // b(9)
	(32.0*(52543486208.0 - 62328059360.0*N - 118489458160.0*N2 +
	       47185442088.0*N3 + 44875379190.0*N4 - 12359832987.0*N5 -
	       5400540075.0*N6 + 1272974916.0*N7 + 200644800.0*N8 - 
	       42495955.0*N9 - 1255067.0*N10 +
	       228194.0*N11))/(14105329875.0*Nm10*N)}; // b(10)
      
      // calculating the value of B0 as a series exansion in powers of 
      // epsilon = (N-K1)/N
      B0 = N*series(order+2, b, (N-K1)/N, -1);
    };	    // all coefficients for the expansion removed from memory 
      
    if(B0<0) {	       // bad, but can still recover precision with NR
				// root finding 
      B0     = err;		// assigning 0+precision error to B0
      warning(NSB_SER_B0);
    };
      
    // using Newton-Raphson to polish the value of B0 and find
    // B0 to the desired precision.
    // The equation we are solving is:
    //   K1/B0 - \Sum_{I=0}^{N-1} 1/(B0+I) \Equiv F(B0) = 0
    // The derivative is:
    //   dF/dB0 = - K1/B0^2 + \Sum_{I=0}^{N-1} 1/(B0+I)^2
    {
      double dB0= 0.0;
      double F  = 0.0;
      double dF = 0.0;
      int counter = 0;
      do{
	counter++;
	F   =  K1/B0 + gsl_sf_psi(B0) - gsl_sf_psi(B0+N);
	dF  =  - K1/gsl_pow_2(B0) + gsl_sf_psi_n(1, B0) -
	  gsl_sf_psi_n(1,B0+N);
	dB0 =  - F/dF;
	B0 += dB0;
	if(B0<=0.0){
	  B0=err;
	  dB0=0.0;
	  warning(NSB_NR_B0_SIGN);
	};
      }while ((fabs(dB0) > fabs(B0*err)) && (counter<=maxcounter));
	
      if(counter>=maxcounter)
	warning(NSB_NR_B0);
    };				// clearing all "while" variables

    Bcl = B0;			// take the calculated value 
				// as the first approx for Bcl
    {
      const int order_K = 4;  // number of terms in the series, orders
      // up to 10 are in the file other_orders_K.m
      
      // temporary variables
      const double pg1B0  = gsl_sf_psi_n(1, B0);
      const double pg1NB0 = gsl_sf_psi_n(1, N+B0);
      const double denum  = K1/gsl_pow_2(B0) - pg1B0 + pg1NB0; // denumerator
      const double pg2B0  = gsl_sf_psi_n(2, B0);
      const double pg2NB0 = gsl_sf_psi_n(2, N+B0);
      const double pg21   = gsl_sf_psi_n(2,1);
      const double pg3B0  = gsl_sf_psi_n(3, B0);
      const double pg3NB0 = gsl_sf_psi_n(3, N+B0);
      const double pg4B0  = gsl_sf_psi_n(4, B0);
      const double pg4NB0 = gsl_sf_psi_n(4, N+B0);
      
      const double f0     = dot_f(szng1, nxng1, kxng1, &gsl_sf_psi);
      const double d1f0   = dot_f_int(szng1, nxng1, kxng1, &gsl_sf_psi_n, 1);
      const double d2f0   = dot_f_int(szng1, nxng1, kxng1, &gsl_sf_psi_n, 2);
      const double d3f0   = dot_f_int(szng1, nxng1, kxng1, &gsl_sf_psi_n, 3);
      
      const double B02 = gsl_pow_2(B0);
      const double B03 = gsl_pow_3(B0);
      const double B04 = gsl_pow_4(B0);
      const double B05 = gsl_pow_5(B0);

      const double PI_2 = gsl_pow_2(M_PI);
      const double PI_4 = gsl_pow_2(PI_2);
	
      // all 4 expansion orders
      double b[order_K];
      b[0] =  B02*(M_EULER*K2 + f0) / (B02*denum);
      const double b02=gsl_pow_2(b[0]);
      const double b03=gsl_pow_3(b[0]);
      const double b04=gsl_pow_4(b[0]);
      b[1] = (K2*PI_2*B0 - (6.0*K1*b02)/B03 - 3*b02*pg2B0 +
	      3.0*b02*pg2NB0 - 6.0*B0*d1f0)/(-6.0*denum);
      const double b12=gsl_pow_2(b[1]);
      b[2] = (K2*PI_2*b[0] + (6.0*K1*b03)/B04 -(12*K1*b[0]*b[1])/B03 + 
	      3.0*K2*B02*pg21 - 6.0*b[0]*b[1]*pg2B0 + 6.0*b[0]*b[1]*pg2NB0 - 
	      b03*pg3B0 + b03*pg3NB0 - 6.0*b[0]*d1f0 - 3.0*B02*d2f0)/ 
	(-6.0*denum);
      b[3] =  -(-(K2*PI_4*B03)/90.0 + (K1*b04)/B05 - (K2*PI_2*b[1])/6.0 - 
		(3.0*K1*b02*b[1])/B04 + (K1*b12)/B03 +  
		(2.0*K1*b[0]*b[2])/B03 - K2*B0*b[0]*pg21 + 
		((b12 + 2.0*b[0]*b[2])*pg2B0)/2.0 
		- ((b12 + 2*b[0]*b[2])*pg2NB0)/2.0 + 
		(b02*b[1]*pg3B0)/2.0 - (b02*b[1]*pg3NB0)/2.0 +
		(b04*pg4B0)/ 24.0 - (b04*pg4NB0)/24.0 +  b[1]*d1f0 + 
		B0*b[0]*d2f0 + (B03*d3f0)/6.0)/(-denum);
      
      Bcl += series(order_K, b, 1/K, 1); // getting the expansion
    };		// clearing all temporary variables for this expansion
      
    const double Ksq  = gsl_pow_2(K);
    // using Newton-Raphson to polish the value of Bcl and find
    // Bcl to the desired precision.
    {
      int counter = 0;
      double dBcl = 0.0;
      double F    = 0.0;
      double dF   = 0.0;
      do{
	counter++;
	F = 1/K*dot_f_offset(szng1, nxng1, kxng1, Bcl/K, &gsl_sf_psi) 
	  - K2/K*gsl_sf_psi(1.0+Bcl/K)
	  + K1/Bcl + gsl_sf_psi(Bcl) - gsl_sf_psi(Bcl+N);
	dF= 1/Ksq* dot_f_int_offset(szng1, nxng1, kxng1, Bcl/K, gsl_sf_psi_n, 1)
	  - K2/Ksq*psi1(1+Bcl/K) - K1/gsl_pow_2(Bcl) + psi1(Bcl) 
	  - psi1(Bcl +N);
	dBcl =  - F/dF;
	Bcl+=dBcl;
	if(Bcl<=0.0){
	  dBcl=0.0;
	  Bcl=err;
	  warning(NSB_NR_BCL_SIGN);
	};
      }while ((fabs(dBcl) > fabs(Bcl*err)) && (counter<=maxcounter));
      
      if(counter>=maxcounter)
	warning(NSB_NR_BCL);
      else if ((warncode == NSB_NR_B0)||
	       (warncode == NSB_NR_B0_SIGN)||
	       (warncode == NSB_SER_B0))
	
	warning(NSB_OK);	// method seems to have recovered
    };			     // clear all unnecessary vars from memory
    
    // Bcl is already set; setting the other output
    
    const double dBcl = 1/Ksq *
      dot_f_int_offset(szng1, nxng1, kxng1, Bcl/K, &gsl_sf_psi_n, 1)
      - K2/Ksq*psi1(1.0+Bcl/K) 
      - K1/gsl_pow_2(Bcl) + psi1(Bcl) - psi1(Bcl +N);
    

    // get a possible mismatch here with max_evidence() from octave code
    // check again if this turns out serious. the mismatch is in 
    // dxicl
    xicl= xi_KB(Bcl);
    dxicl = pow(-dBcl/gsl_pow_2(dxi_KB(Bcl)),-0.5);
  };
  
  return warncode;
}



/*************************************
  As noticed during the sfly034 analysis, bins with large counts
  screw up calculations of entropy for less populated counts. So we
  partition counts into three blocks -- bins that take more than
  10000 samples or simultaneously more than 5% of the total sample
  and more than 1000 samples; the second grup is the bins with >50
  counts, and the third group is all other bins. We then calculate
  entropy for each group separately and get the total entropy out of
  this. 
*/
ENTRDATA::NSB_WARN ENTRDATA::NsbCalc::DoCalcsPart(BXIset& bxiall)   
  throw(EntrData_badnum, EntrData_bad_alloc) {

  std::vector<double> _kx1(sz);
  std::vector<double> _nx1(sz);
  std::vector<double> _kx2(sz);
  std::vector<double> _nx2(sz);
  std::vector<double> _kx3(sz);
  std::vector<double> _nx3(sz);
  

  // data partitioned into 3 regions
  unsigned long j1=0;		// how many unique counts in first,
  unsigned long j2=0;		// second, and third regions
  unsigned long j3=0;
  double alphsize1=0.0;		// number of such bins (aphabet size)
  double alphsize2=0.0;
  double alphsize3=0.0;
  double events1=0.0;		// number of samples in each of the groups
  double events2=0.0;
  double events3=0.0;
  for(unsigned long i=0;i<sz;i++){
    if((nx[i]>10000)||((nx[i]>1000)&&(nx[i]>0.05*N))){
      _kx1[j1]=kx[i];
      _nx1[j1]=nx[i];
      events1 += _kx1[j1]*_nx1[j1];
      j1++;
      alphsize1+=kx[i];
    }else if(nx[i]>50){
      _kx2[j2]=kx[i];
      _nx2[j2]=nx[i];
      events2 += _kx2[j2]*_nx2[j2];
      j2++;
      alphsize2+=kx[i];
    }else{
      _kx3[j3]=kx[i];
      _nx3[j3]=nx[i];
      events3+= _kx3[j3]*_nx3[j3];
      j3++;
    }
  }
  _kx1.resize(j1);
  _nx1.resize(j1);
  _kx2.resize(j2);
  _nx2.resize(j2);
  _kx3.resize(j3);
  _nx3.resize(j3);
  alphsize3=K-alphsize1-alphsize2;
  UniPrior* pp1= new UniPrior(log(alphsize1));
  UniPrior* pp2= new UniPrior(log(alphsize2));
  UniPrior* pp3= new UniPrior(log(alphsize3));

  NsbCalc  nsb1(_kx1, _nx1, alphsize1, bxiall, pp1, NOPART, err);
  NsbCalc  nsb2(_kx2, _nx2, alphsize2, bxiall, pp2, NOPART, err);
  NsbCalc  nsb3(_kx3, _nx3, alphsize3, bxiall, pp3, NOPART, err);
  
  // we now need the entropy of choice between the three classes
  std::vector<double> _kx_choice(3);
  std::vector<double> _nx_choice(3);
  _kx_choice[0]=1;
  _kx_choice[1]=1;
  _kx_choice[2]=1;
  _nx_choice[0]=events1;
  _nx_choice[1]=events2;
  _nx_choice[2]=events3; 
  UniPrior* pp4= new UniPrior(log(3));
  NsbCalc  nsb4(_kx_choice, _nx_choice, 3, bxiall, pp4, NOPART, err);

  // from now on, events are actually fractions
  events1=events1/N;
  events2=events2/N;
  events3=events3/N;
  
  Snsb = events1*nsb1.get_Snsb() + events2*nsb2.get_Snsb()+events3*nsb3.get_Snsb() + nsb4.get_Snsb();
  //see notebook june 24, 2005 for below
  dSnsb= sqrt(gsl_pow_2(nsb1.get_Snsb()-nsb3.get_Snsb())*events1*(1-events1)/N +
	      gsl_pow_2(nsb2.get_Snsb()-nsb3.get_Snsb())*events2*(1-events2)/N +
	      gsl_pow_2(events1*nsb1.get_dSnsb()) +
	      gsl_pow_2(events2*nsb2.get_dSnsb()) +
	      gsl_pow_2(events3*nsb3.get_dSnsb()) +
	      gsl_pow_2(nsb4.get_dSnsb()));
  Sml = events1*nsb1.get_Sml() + events2*nsb2.get_Sml()+events3*nsb3.get_Sml() + nsb4.get_Sml();
  // for other quantities, we record their value from the third class
  // of bins, those that are undersampled
  Scl  = nsb3.get_Scl();
  Bcl  = nsb3.get_Bcl();
  xicl = nsb3.get_xicl();
  dxicl = nsb3.get_dxicl();
  Sas  = nsb3.get_Sas();
  dSas = nsb3.get_dSas();

  K1 = 0;
  K2 = 0;
  
  return GSL_MAX(GSL_MAX(GSL_MAX(nsb1.get_warncode(),nsb2.get_warncode()), nsb3.get_warncode()),nsb4.get_warncode());

}


/*********************************************
  Main calculational routine for the NSB method; equivalent to
  find_nsb_entropy in the Octave code.
  Calculate (change) the variables:
     Snsb   - entropy estimate by the NSB method, scalar;
     dSnsb  - the standard deviation of the estimate;
     Scl    - entropy at the saddle point;
     dScl   - standard deviation at the saddle point;
     xicl   - saddle point;
     Sml    - maximum likelihood (naive) entropy estimate;
     warncode - takes values of type NSB_WARN, see descriprion in
     nsb.h
*/
ENTRDATA::NSB_WARN ENTRDATA::NsbCalc::DoCalcs() 
  throw(EntrData_badnum, EntrData_bad_alloc) {
    
  const int numint=3;		// number of integrands
  // three different integrands and their names 
  gsl_function ints[numint];
  ints[0].function = &g_int_1;
  ints[0].params = this;
  ints[1].function = &g_int_S;
  ints[1].params = this;
  ints[2].function = &g_int_S2;
  ints[2].params = this;
  const char* msgs[numint]={"normalization","S", "S^2"};
  double vals[numint];		// values of integrals
  double abserr[numint];	// estimate of absolute errors

  // preparing for saddle point integration
  make_ng1();			// filling in n>1 arrays
  max_evidence();		// finding the saddle
  Scl=meanS(Bcl);


  // now, suppose we use "bnd" prior. The problem is then that for
  // large N (or rather large N-K1), the value ox xicl is not the
  // same as the value of Scl or <S>, so that the bound which is on
  // xi is wrong, and instead it should be on Scl or <S>. So xi must
  // be such that the corresponding S(xi,N) is within bounds. But
  // this means that the prior is determined after observations are
  // made, and the procedure is not Bayes kosher (see more details in
  // the notebook). We still do this, but only if the number of
  // coincidences is small, less than 100 or so.

  // we start with edges err away from 
  double edges[]={GSL_MAX(err,pr->GetMin()), GSL_MIN(maxent-err,pr->GetMax())}; 
  // then we build a spline that maps between xicl and mean S
  {			      
    const int npts = 40;	// number of points in spline
				// relating xi and S
    double x[npts];		// this will be mean S
    double y[npts];		// this will be xi
    for(int i=0;i<npts;i++){
      y[i]=err+(maxent-2*err)/(npts-1)*i;
      x[i]=meanS(B_xiK(y[i]));
    }
    spline xi_S(x,y,npts,gsl_interp_linear);

    // setting worst case integration limits  
    // if the limits are not the 0, logK boundaries, then we need to
    // get xi that corresponds to them
    // then we need to do a bunch of sanity cheks on the values of edges
    if (edges[0]>err){
      edges[0]= GSL_MAX(err, GSL_MIN(maxent-err,xi_S.eval(pr->GetMin())));
    }
    if (edges[1]<maxent-err){
      edges[1]= GSL_MIN(maxent-err,GSL_MAX(err,xi_S.eval(pr->GetMax()))); 
    }
    if(edges[0]>=edges[1]){
      edges[0]=err;
      edges[1]=maxent-err;
    }
  }
  if(xicl<maxent-err){		// only think of removing bound if
				// xicl is smaller than the maximum allowed
    // if bad xicl
    if((xicl<edges[0])||(xicl > edges[1])){
      edges[0]=err;
      edges[1]=maxent-err;
    } 
    // if a lot of coincidences -- one may choose to lift the edge
    // bouds, we don't do it here, and leave bounds even the number of
    // coincidences is large, so that xi and S are not the same, and
    // the bounds are not kosher
    if((N-K1)>100){
      if (pr->GetName()!=UniName){ // consider switching prior
	std::cout<<"  "<<warn_names[NSB_MANY_COINC]<<'\n';
      }
    }
  }


  double xilim[2];		// actual integration boundaries
  double delta = 0.0; // delta is the interval around the peak on which
  // gaussian approximation falls to "err/2" on each side


  if(warncode!=NSB_OK){
    // some problem with finiding saddle value; switch to full range
    if(vlevel>VERBMIN) std::cout<<"  Switching to integration over the full range due to: " << 
      warn_names[warncode] <<'\n';
    // integrating either over the whole range
    xilim[0]  = edges[0];
    xilim[1]  = edges[1];
    mlog = mlog_evidence(B_xiK(xilim[0])); // don't know the value at
					   // the saddle 
    // now recurse through the entire domain to find smaller mlog
    for(double x=xilim[0];x<=xilim[1];x+=(xilim[1]-xilim[0])/100){
      mlog=GSL_MIN(mlog,mlog_evidence(B_xiK(x)));
    }
  }
  else{
    if(vlevel>VERBMED) std::cout << "  Expect:  S = " << Scl << " +- " << dxicl << ". ";
    mlog = mlog_evidence(Bcl); // value at the saddle
    delta = dierfc(err/2.0)*M_SQRT2; // not sure about this since I
				// don't quite know how dierfc
				// is defined; but the details
				// aren't important here
    
    if(vlevel>VERBMED) std::cout<<"Integrating around the peak.\n";
  };

  workspace_AQ ws;		// creating integration workspace
  for(int i=0; i<numint;i++){
    // if integrating around the peak, need the limits of integration
    if(delta>0.0){
      // the value of the integrand at xi_cl
      const double cent = (*(ints[i].function))(xicl,ints[i].params);
      // if the integral was purely Gaussian, the value at +-delta
      // would've been of required precision 1/sqrt(2pi)*
      // exp(-(delta^2)/2). As a safety, we require the value at the
      // limits of integration to be ten times smaller.
      const double good = 0.1*cent*exp(-(delta*delta)/2)*1/sqrt(2*M_PI);
      // increase the integration window, until the integrand is small 
      // at the edges; start with xicl +- delta*dxicl and expand by +-
      // 0.5*dxicl
      // make sure that the limits don't go outside the edges

      xilim[0] = GSL_MAX(xicl-delta*dxicl, edges[0]);
      xilim[1] = GSL_MIN(xicl+delta*dxicl, edges[1]);

      double limval[2] = {(*(ints[i].function))(xilim[0],ints[i].params),
			  (*(ints[i].function))(xilim[1],ints[i].params)};

      for(int j=0; j<2;j++){	// do a loop over lower and upper
				// limits
	double window = 0.0;	// limits = (delta+window)*dxicl
	while((limval[j]>good) && (xilim[j]!=edges[j])){
	  if(window>10.0){window *=1.2;} // fast growth for large windows
	  else{window +=0.5;};	// regular growth
	  // increase the corresponding limit, but make sure we are
	  // still in range
	  if(!j){		// for left limit
	    xilim[0] = GSL_MAX(edges[0], xicl - (delta+window)*dxicl);
	  } 
	  else {		// for right limit
	    xilim[1] = GSL_MIN(edges[1], xicl + (delta+window)*dxicl);
	  };
	  limval[j] = (*(ints[i].function))(xilim[j],ints[i].params);
	};
      };
    };				// finished calculating the limits
    
    if(xilim[0]>=xilim[1]){	// a crude hack to check the problem
      xilim[0]=edges[0];
      xilim[1]=edges[1];
    }


    if(vlevel>VERBMED)std::cout<< "    Doing " << msgs[i]<< " integral. " <<
      "Limits: " <<  xilim[0] << " < xi < " << xilim[1] <<".\n";

    int ec=gsl_integration_qag(&ints[i], xilim[0], xilim[1], 0.0, err, 
			       ws.get_size(), GSL_INTEG_GAUSS21,  ws.get_ws(),
			       &vals[i], &abserr[i]);
    if(abserr[i]/vals[i]>err) warning(NSB_INT_NOCONV);
    
    if(ec!=GSL_SUCCESS)		// error handling
      throw EntrData_badnum(NSB_GSL, ec);
  };				// end of the for loop over int_1,
				// int_S, int_S2

  Snsb = vals[1]/vals[0];
  dSnsb = sqrt(vals[2]/vals[0] - gsl_pow_2(Snsb));

  if(vlevel>VERBMIN)std::cout << "  Found S = " << Snsb << " +- " << dSnsb << ".\n";

  // calculating Sml
  for(unsigned long i =0; i<sz;i++)
    Sml -= nx[i]/N*log(nx[i]/N)*kx[i];

  double D=N-K1;		// # of coincidences
  if(D){
    Sas = -psi(1.0) -log(2.0) +2.0*log(N)-psi(D);
    dSas= sqrt(psi1(D));
  }
  else{
    Sas =  std::numeric_limits<double>::infinity();
    dSas=  std::numeric_limits<double>::infinity();
  }
  return warncode;
}

//
// integrands
//


/*******************************************
  Integrand for the normalization integral; the only parameter here
  is the value of the argument; the rest is implicit through class
  structures.
  This includes prior over xi;
*/
double ENTRDATA::NsbCalc::int_1(const double& xi) const{
  return  exp(-mlog_evidence(B_xiK(xi)) + mlog)*prior_xi(xi);
}



/*******************************************
   Same for the S integral.
*/
double ENTRDATA::NsbCalc::int_S(const double& xi)const{
  double B=B_xiK(xi);
  return  exp(-mlog_evidence(B) + mlog)*prior_xi(xi)*meanS(B);
}



/*******************************************
   Same for the S^2 integral.
*/
double ENTRDATA::NsbCalc::int_S2(const double& xi)const{
  double B=B_xiK(xi);
  return  exp(-mlog_evidence(B) + mlog)*prior_xi(xi)*meanS2(B);
}



/***********************************************/  
// integration routines envelopes used to call the member
// integration routines from nsb class objects; g_ in the name
// stands for "global"
double ENTRDATA::g_int_1(double xi,  void* obj){
  NsbCalc* p = (NsbCalc*) obj;
  return p->int_1(xi);
}
double ENTRDATA::g_int_S(double xi,  void* obj){
  NsbCalc* p = (NsbCalc*) obj;
  return p->int_S(xi);
}
double ENTRDATA::g_int_S2(double xi, void* obj){
  NsbCalc* p = (NsbCalc*) obj;
  return p->int_S2(xi);
}



