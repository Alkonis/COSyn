#include "COfit.h"
#include <iostream>
#include <cstdlib>
//optimizations so far:  read in files ONCE in FitData() rather than in each iteration
//                       removal of redundant declarations
//                       memory allocation of rate_eqtn variables in access order
//                       adoption of an optimized linear algebra computation library (armadillo)
//                       some optimizations to calculation order in big loops



//To do:  Fix transposition (later)
//        benchmark rate_eqtn solving
//        Properly skip fluorcalc, etc
//        add input file
//        fix g-factors (?)
//        fix grid

using namespace std;
using namespace idlarma;

FitData* data;

constexpr double FitData::vib_einA[];
constexpr double FitData::Bv12[];
constexpr double FitData::Bv13[];

int FitData::numGuesses;

//returns an array of size numGuesses filled with random numbers on the interval (offset, modulus+offset)
//used to generate the matrix of guess data
double* FitData::FillArray(int modulus, int offset)
{
cout << "numguesses:  " << numGuesses << endl;
  cout << "Filling array..." << endl;
  double* array;
  array  = new double[numGuesses];
cout << "Made array..." << endl;
  for(int i=0; i<numGuesses;i++) {
    array[i]=rand() % modulus + offset;
  }
  cout << "Returning array..." << endl;
  return array;
}

int FitData::runTrial() {

//============================
//Divide into annulli
//=============================
  r=disk_in;  

  HD100546_luminosity = HD100546_luminosity*rel_lum;
  ranulli.push_back(disk_in);
  if (r<0.1)
  {
    while ((ranulli.back() < 1) && (ranulli.back() < disk_out))
    {
      r_index_a=r_index_a+1;
      ranulli.push_back(disk_in+0.01*r_index_a);
    }
  }

  maxra = ranulli.back();
  if ((maxra < 1) && (maxra >= .1))
  {
    while ((ranulli.back() < 1.0) || (ranulli.back() < disk_out))
    {
     r_index_b++;
     ranulli.push_back(maxra+0.1*r_index_b);
    }
  }
  maxrb=ranulli.back();
  if ((maxrb <= disk_out) && (maxrb >= 1.0))
  {
    while  (ranulli.back() < disk_out)
    {
      r_index_c++;
      ranulli.push_back(maxrb+1.0*r_index_c);
    }
  }

  steps=ranulli.size();  //check that this size function is correct
  
  vec rdisk= zeros<vec>(steps);
 
  for (int i=0; i < steps; i++ )
  {
    rdisk.at(i)=ranulli.at(i)*1.496E13;
  }
 
//benchmark here!  Make sure these results are the same from the IDL code and the C++ code!
  

//=========================================
//  CollData setup
//========================================
  CollData* d = new CollData(this->layers, this->steps);
  for (auto i=0; i<steps;i++) {
    for (auto j=0; j<layers;j++) {
      for (auto k=0; k<10;k++) {
        for (auto q=0; q<11; q++)
        {
          d->rate_eqtn(i,0)(k+11,q,j)=einA.at(k,q);//k+11, q,j,i,   k ,q  //why no zero here
        }
      }
    }
  }

  for (auto j=0; j<steps; j++) 
  {
    d->rate_eqtn.at(j,0)(span(1),span(0),span::all).fill(vib_einA[0]);   
    d->rate_eqtn.at(j,0)(span(10),span(10),span::all).fill(-vib_einA[9]); 

    for (auto i=0; i<9; i++)
    {
      d->rate_eqtn.at(j,0)(span(i+11),span(i+11),span::all).fill(-sum(einA.row(i)));
    }
   

    for (int z=0; z<9; z++)
    {
      for (int q=0; q<layers; q++) 
      {
        d->rate_eqtn.at(j,0).at(z+1,z+1,q)=-vib_einA[z];  
        d->rate_eqtn.at(j,0).at(z+2,z+1,q)=vib_einA[z+1]; 
      }
    }
   
    for (int q=0; q<layers; q++)
    {
      d->rate_eqtn.at(j,0).subcube(span(1),span(0),span::all).fill(vib_einA[0]);  //1,0,all,k
      d->rate_eqtn.at(j,0).subcube(span(9),span(9),span::all).fill(-vib_einA[8]);  //9,9,all,k
    }
   
  }

  d->rate_eqtn2=d->rate_eqtn;   //Clone colldata for the second pass without collisions
   
  cerr << "CollData generated." << endl;

  T_rot_fl = T_rot0_fl *pow((1.5E13/rdisk),T_rot_alpha_fl);
  
  ivec T_rot_index=where(T_rot_fl, [] (double datum) {return datum >= 3500;});

//cout << "T_rot_cnt:  " << T_rot_cnt << endl;
cout << "T_rot_fl:  " << T_rot_fl << endl;
cout << "T_rot_index:  " << T_rot_index << endl;
  if (T_rot_index.at(0) != -1) {
    for (auto elem : T_rot_index )
    {
      T_rot_fl[elem]=3500;
    }
  }


  T_rot_cl=T_rot0_cl * pow((1.496E13/rdisk), T_rot_alpha_cl);
  H_den=H_den0 * pow((1.496E13/rdisk),H_den_alpha);
  T_rot_index = where(T_rot_cl, [] (double datum) {return datum > 3500;});
//  T_rot_cnt=T_rot_index.size();
cout << "T_rot_index:  " << T_rot_index << endl;
  if (T_rot_index.at(0) != -1) {
    for (auto elem : T_rot_index)
    {
      T_rot_cl[elem]=3500;
    }
  }

  cout << "T_rot_cl:" << endl << T_rot_cl << endl;
  cout << "T_rot_fl:" << endl << T_rot_fl << endl;
  cout << "H_den:"  << endl << H_den << endl;

  
 //   CollData d = new CollData(this->layers, this->steps);
//  if (rel_lum <= 1e-3) goto skip_fluorcalc;
//{
//
    vec b_tot = sqrt(2*1.38e-16*6.02e23*T_rot_fl/28 + pow(v_turb,2));



    vec tau = linspace<vec>(0,19.99,2000); 
    vec F_tau=zeros<vec>(2000);
    //==============CALCULATE DFDT BY INTEGRAATION USING TRAPEZOID RULE================
    cout << "b_tot" << endl;
    cout << b_tot;
    
    cout << "tau:";
    cout << tau;
 
    cerr << "Integrating...";
    double sum;
    double n = 1000;
    double a = 0;
    double b = 5;

    double delta=(b-a)/n;

    for (int i=0; i<2000; i++) 
    {
      sum=0.5 * ( (1 - exp(-tau(i)*exp(-pow(a,2)))) + (1 - exp(-tau(i)*exp(-pow(b,2)))) );
      for (int j=0; j<n; j++)  
      {
        auto x=a+delta*j;
        sum=sum+(1-exp(-tau(i)*exp(-pow(x,2))));
      }
      sum = sum*delta;
      F_tau(i)=sum;
    }
    
    vec dFdt=deriv(tau,F_tau);

//    cout << "F_tau:  " << endl;    
//   cout << F_tau << endl;
//    cout << "dFdt  :" << endl;
//    cout << dFdt << endl;



      cube Nv_coll;
      mat tot_col_fluor;
      mat tot_col_fluor_back;                      //move, perhaps, to colldata?

      cube Nv_nocoll;
      mat tot_col_fluor_nocoll;



    cerr << "Beginning collisional computations..." << endl;

//==========================================================
//   BIG COLLISION LOOP
//=========================================================
    for (int coll_loop=0; coll_loop<2; coll_loop++) 
    {
       cerr << "in Coll_Loop :  " << coll_loop << endl;
      mat U;
      mat s;
      mat v;
       
      vec sol;
      ivec Nv_index;


      double gsum;

      for (int k=0; k<steps; k++) 
      {

        mat Fuv = HD100546_luminosity / ( 4*3.1415926535897*pow(rdisk[k],2) );

	cerr <<"Steps loop, k=" << k << endl;
       
       //=============COLLISIONS=============

       if (coll_loop != 1)
       {
	 auto k_HCO_dn = (7.57e-15*T_rot_cl[k]/(1-exp(-3084/T_rot_cl[k])))*H_den[k];
	 auto k_HCO_up=k_HCO_dn*exp(-3084/T_rot_cl[k]);

	 d->rate_eqtn.at(k,0).subcube(span(0),span(0),span::all)-=k_HCO_up;
	 d->rate_eqtn.at(k,0).subcube(span(1),span(0),span::all)+=k_HCO_dn;

	 for (int i=1; i<9; i++) 
	 {
	   d->rate_eqtn.at(k,0).subcube(span(i-1),span(i),span::all)+=k_HCO_up;
	   d->rate_eqtn.at(k,0).subcube(span(i  ),span(i),span::all) = d->rate_eqtn.at(k  ,0).subcube(span(i),span(i),span::all)-k_HCO_dn-k_HCO_up;
	   d->rate_eqtn.at(k,0).subcube(span(i+1),span(i),span::all)+=k_HCO_dn;
	 }

       }

	 //===================ULTRAVIOLET PUMPING========================

	 for (int j=0; j<layers; j++)
	 {

	   if (j>0)
	   {
	     for (int i=0; i<10; i++)
	     {
	      d->tau_0.subcube(span(i),span::all,span(j))=7.7585e12 * 0.02654 * arma::sum(d->Nv.slice(k).submat(span(0,11),span(0,j)),1) % fXA.submat(span(i),span::all).t() % lam_ang.submat(span(i),span::all).t()*1e-8/(sqrt(3.1415926535897)*b_tot[k]);
	     }
	   }
           auto max_iter=d->tau_0.slice(j).row(0).n_elem;
	   for (auto ii=0; ii<max_iter; ii++)
	   {
	     if (d->tau_0.at(0,ii,j) < 20)
	     {
	       ivec dFdt_0_index=where(tau, [&] (double elem) {return elem == round(d->tau_0.at(0,ii,j)*10)/10;});
	       auto count = dFdt_0_index.size();
//	       cout << "dFdt_0_index:"  << dFdt_0_index <<  "  count:  " << count << endl; 
	       if (dFdt_0_index.at(0) != -1)
	       {
	         d->dFdt_0.subcube(span::all,span(ii),span(j)).fill(dFdt.at(dFdt_0_index.at(0)));    //weird line!!!!!
	       }

	     }
	     else 
	     {
	       d->dFdt_0.slice(j).row(ii).fill((1/(2*d->tau_0.at(0,ii,j))*sqrt(log(d->tau_0.at(0,ii,j)))));    //(span::all,span(i),span(j))=1/(2*tau_0.at(0,ii,j))*sqrt(alog(tau_0.at(0,ii,j)));   CHECK THIS WITH SOURCE
	     }
	   }
//           cout << "dFdt_0.slice(j):  " << d->dFdt_0.slice(j) << endl;

	   d->dwdn.slice(j)=d->dFdt_0.slice(j)*.02654*2.0 % (lam_ang*1e-4) % (lam_ang*1e-8) % fXA/(sqrt(3.1415926535897)*c*1e5);

//cout << "dwdn.slice(j):  " << endl;
//cout << d->dwdn.slice(j) << endl;

	   for (int ii=0; ii<10; ii++) 
	   {
//	     cout << "ii,j,k (" << ii << "," << j << "," << k << ")" << endl;
//             cout << (d->dwdn.slice(j).row(ii) * 3.1415926535897 % Fuv.row(ii) / (hc * wavenum.row(ii))).t() << endl;
	     d->g.at(ii,0).slice(k).col(j) = (d->dwdn.slice(j).row(ii) * 3.1415926535897 % Fuv.row(ii) / (hc * wavenum.row(ii))).t();  //check this to be sure the constants are filling in right...
           } 
       
//=======================
//  Add in G-factors
//======================
	   gsum=0;

	   for (int ii=0; ii<10; ii++)
	   { 
	     gsum+=d->g.at(ii).at(0,j,k);
	   }
	   d->rate_eqtn.at(k,0).at(0,0,j)-=gsum;

	   for (int i=1; i<11; i++)
	   {
	     gsum = 0;
	     for (int ii=0; ii<10; ii++)
	     {
	       gsum+=d->g.at(ii).at(i,j,k);
	     }
	     d->rate_eqtn.at(k,0).subcube(span(i),span(i),span(j))-=gsum;
	   }

           
           for (int i=0; i<9; i++)
	   {
	     for (int ii=0; ii<11; ii++)
	     {
	       d->rate_eqtn.at(k,0).at(ii,11+i,j)+=+d->g.at(i,0).at(ii,j,k); //ii 11+i j k
	     }
	   }
	   d->rate_eqtn.at(k,0).at(1,10,j)=0;
	   

//===================
//  SOLVE
//==================  
	   vec z = zeros<vec>(21);
	   z.at(20)=1;
	     
//	   cout << "matrix for solving:" << endl;
//	   cout << d->rate_eqtn.at(k,0).slice(j).t();
	   solve(sol,d->rate_eqtn.at(k,0).slice(j).t(),z);
	   d->Nv.slice(k).col(j)= sol;  //check this to be sure the array is filling in the right dimension

//           cout << "d->Nv.slice(k).col(j):  " << d->Nv.slice(k).col(j) << endl;
           Nv_index=where(d->Nv.slice(k).col(j), [] (double datum) {return datum < 0;});
	   if (Nv_index.at(0) != -1) 
           {
             for (int jj=0; jj<Nv_index.n_elem; jj++) 
             {
               d->Nv.at(Nv_index[jj],j,k)=0;
             }
           }
//           cout << "End of solve loop" << endl;

 	 }
       }

//}
//skip_fluorcalc:
      if (coll_loop==0) {
        Nv_coll=d->Nv;
        tot_col_fluor = totalDim(d->Nv*7.55858e12,2).t();      //check all uses of accu/total in teh code!!!!!! dimension specification!
        tot_col_fluor_back=tot_col_fluor;
        d->rate_eqtn=d->rate_eqtn2;
        cerr << "Beginning pass without collisions..." << endl;
      }
      if (coll_loop==1) {
        Nv_nocoll=d->Nv;
        tot_col_fluor_nocoll = totalDim(d->Nv*7.55858e12,2).t();
      }
    }


//=========================================================================
// Angle of Incidence Correction (tweak for each star--use input file!)
//========================================================================
  cerr << "Correcting for angle of incidence..." << endl;
  vec m_disk=1.5*(5.59647e-10)*sqrt(T_rot_fl / Mstar)%pow(rdisk,0.5);
  vec m_uv = (5.59647e-10)*sqrt(T_rot_fl/Mstar)%sqrt(rdisk);

  vec phi = -atan((m_uv-m_disk)/(1+m_uv%m_disk));
  phi.at(0)=3.1415926535897/2;                         //standardize pi!

  auto xi = tot_col_fluor.n_cols-1;
  auto yi = tot_col_fluor.n_rows-1;

  for (int j=0; j<yi; j++)                    //looks like I need to transpose indices on everything involving tot_col_fluor 
  {
    tot_col_fluor.row(j)=sin(phi.at(j))*tot_col_fluor.row(j);
    tot_col_fluor_nocoll.row(j)=sin(phi.at(j))*tot_col_fluor_nocoll.row(j);
  }
  mat tot_col=tot_col_fluor(span::all,span(0,9));
  mat tot_col_coll=tot_col;

  double r = dist/1.496e13;
  rdisk=rdisk/1.496e13;            //convert to AU



//===================================================================================================
//  MOLECULAR DATA
//==================================================================================================
  
  cerr << "Beginning the processing of molecular data..." << endl;

  cube N_12CO_vj = zeros<cube>(tot_col.n_cols-1,120,tot_col.n_rows);                 //flipping these indices...
  cube N_13CO_vj = zeros<cube>(3,120,tot_col.n_rows);
  cube N_C18O_vj = zeros<cube>(1,120,tot_col.n_rows);


  //cerr << "tot_col rows,cols:  " << tot_col.n_rows << " " << tot_col.n_cols << endl;
  rowvec Jup;
  rowvec Jdn;
  rowvec wvn;
  rowvec EinA;
cerr << "First..." << endl;
  for (int i=0; i<steps; i++)
  {

    cerr << "molecular steps loop i=" << i << endl;

    //=======================
    //   12CO
    //======================

    cerr << "12CO" << endl;

    for (int j=0; j< tot_col.n_cols-1; j++)
    {
      Jup = X12CO(span(j),span(3),span::all);          //because the indices are reversed in c++, do all these first/second indices flipped.
      Jdn = X12CO(span(j),span(1),span::all);
      wvn = X12CO(span(j),span(6),span::all);
      EinA= X12CO(span(j),span(4),span::all);
cerr << "tot_col_flupr.at(i,j+1):  " << tot_col_fluor.at(i,j+1) << endl;
cerr << "tot_col_fluor.at(j+1,i):  " << tot_col_fluor.at(j+1,i)<< endl;
cin.get();
      N_12CO_vj.slice(i).row(j)=tot_col_fluor.at(i,j+1) * (2*Jup+1) % exp(-hck*Bv12[j]*Jup % (Jup+1)/T_rot_fl[i])/(T_rot_fl[i]/(hck*Bv12[j])) 
+ tot_col_coll.at(i,j+1)*(2*Jup+1) % exp(-hck*Bv12[j] * Jup % (Jup+1)/T_rot_cl[i])/(T_rot_cl[i]/(hck*Bv12[j]));

//cerr << N_12CO_vj.slice(i) << endl;
//cin.get();
    }
  cerr << "13CO" << endl;
  //====================
  //  13CO
  //====================
  for (int j=0; j<3; j++)
  {
    Jup = X13CO(span(j),span(3),span::all);
    Jdn = X13CO(span(j),span(1),span::all);
    wvn = X13CO(span(j),span(6),span::all);
    EinA= X13CO(span(j),span(4),span::all);

    N_13CO_vj.slice(i).row(j)=(tot_col_fluor_nocoll.at(i,j+1)/X12CO_13CO_fl)*(2*Jup+1) % exp(-hck*Bv13[j] * Jup % (Jup+1) / T_rot_fl[i]) / (T_rot_fl[i]/(hck*Bv13[j]));

  }
  //====================
  //  C18O
  //===================

  Jup = XC18O(span(0),span(3),span::all);
  Jdn = XC18O(span(0),span(1),span::all);
  wvn = XC18O(span(0),span(6),span::all);
  EinA= XC18O(span(0),span(4),span::all);

  N_C18O_vj.slice(i).row(0)=(tot_col_fluor_nocoll.at(i,1)/X12CO_C18O_fl)*(2*Jup+1) % exp(-hck*Bv18[0]*Jup%(Jup+1)/T_rot_fl[i]) / (T_rot_fl[i]/(hck*Bv18[0]));

  }

  cerr << "Ended steps loops.";
  double freq_size=log10(f_f/f_i)/log10(1+v/(3*c));
  cerr << "1" << endl;
  vec freq = zeros<vec>(freq_size+1);
  freq.at(0)=f_i;
  vec vfreq = freq;
  vfreq.fill(0);

  cerr << "2"<< endl;
  cerr << "freq.n_elem" << freq.n_elem << endl;
  cerr << "freq_size" << freq_size << endl;
  for (int i=1; i<freq_size; i++)
  {
    freq.at(i)=f_i*(pow(1+v/(3*c),i));
  }
  cerr << "3" << endl;
  vfreq=(freq(round(freq.n_elem/2))-freq)*c/freq(round(freq.n_elem/2));                      //check
  cerr << 4 << endl;

  ivec freq_index=where(freq, [] (double datum) {return (datum >= f_i) && (datum <= f_f);});
  cerr << 5 << endl;
  //define annuli

  vec annuli = zeros<vec>(rdisk.n_elem);
  cerr << 6 << endl;
  for (int i=0; i<steps-1; i++)
  {
    annuli.at(i)=3.1415926535897*(pow(rdisk.at(i+1),2)-pow(rdisk.at(i),2));
  }
  cerr << 7 << endl;
  annuli(steps-1)=3.1415926535897*(pow(rdisk.at(steps-1)+1,2)-pow(rdisk.at(steps-1),2));
  cerr << 8 << endl;
  mat stick_spec_12CO = zeros<mat>(freq_size+1,steps);
  mat stick_spec_13CO = stick_spec_12CO;
  mat stick_spec_C18O = stick_spec_12CO;
  mat stick_spec_tot  = stick_spec_12CO;
  cerr << 8 << endl;
  double A0;
  double A1;
  double A2;

  cerr << "Beginning molecular processing over annuli:" << endl;
  //======================
  //  X12CO
  //=====================

  for (int i=0; i<steps; i++)  //Loop over annuli
  {
double btotcomp=b_tot.at(i)/cexp;
    cerr << "annulus i=" << i << endl;
    //==============================
    //  X12CO
    //=============================
    int upbound=tot_col.row(0).n_elem-1;
    for (int j=0; j<upbound; j++)  //vibrational levels
    {
      Jup = X12CO(span(j),span(3),span::all);
      Jdn = X12CO(span(j),span(1),span::all);
      wvn = X12CO(span(j),span(6),span::all);
      EinA= X12CO(span(j),span(4),span::all);
      for (int k=0; k<X12CO.n_slices; k++) 
      {
        A1=wvn.at(k);
	if ( (A1 >= f_i) && (A1 <= f_f) )
	{

//cerr << "i,j,k:  " << i << " " << j << " " << k << endl;
	  A0=N_12CO_vj.at(j,k,i)*hc*A1*EinA.at(k);                         //should the first two indices be reversed here?
	  A2=btotcomp*A1;
	  stick_spec_12CO.col(i)+=(A0/(rpi*A2)) * exp (-pow(((A1-freq)/A2),2));
cerr << "N12CO_vj.at(j,k,i):  " << N_12CO_vj.at(j,k,i) << endl;;
cerr << "A1:  "<< A1 << endl;
cerr << "A0:  "<< A0 << endl;
cerr << "A2:  " << A2 << endl;
cerr <<"calc coeff:  " << (A0/(rpi*A2)) << endl;
cin.get();
cerr << "exponent:  " << exp(-pow(((A1-freq)/A2),2)) << endl;;
cin.get();
//cin.get();;
	}
      }
    }
    //==============================
    //  X13CO
    //==============================
    
    for (int j=0; j<3; j++)   //Loop over vibrational levels--check for-loop bounds in some earlier integral loops... might need to be adjusted!
    {

      Jup = X13CO(span(j),span(3),span::all);
      Jdn = X13CO(span(j),span(1),span::all);
      wvn = X13CO(span(j),span(6),span::all);
      EinA= X13CO(span(j),span(4),span::all);

      for (int k=0; k<X13CO.n_slices; k++)
      {
          A1=wvn.at(k);
          if ((A1 >= f_i) && (A1 <= f_f)) 
          {
	    A0=N_13CO_vj.at(j,k,i)*hc*A1*EinA.at(k);                         //should the first two indices be reversed here?
	    A2=btotcomp*A1;
	    stick_spec_13CO.col(i)+=(A0/(rpi*A2)) * exp (pow(-((A1-freq)/A2),2));
          }
      }    

    }
    //=============================
    //  XC180
    //=============================
    
     
      Jup = XC18O(span(0),span(3),span::all);
      Jdn = XC18O(span(0),span(1),span::all);
      wvn = XC18O(span(0),span(6),span::all);
      EinA= XC18O(span(0),span(4),span::all);

      for (int k=0; k<XC18O.n_slices; k++)
      {
          A1=wvn.at(k);
          if ((A1 >= f_i) && (A1 <= f_f))
          {
	    A0=N_13CO_vj.at(0,k,i)*hc*A1*EinA.at(k);                         //should the first two indices be reversed here?
	    A2=btotcomp*A1;
	    stick_spec_C18O.col(i)+=(A0/(rpi*A2)) * exp (pow(-((A1-freq)/A2),2));
          }
      }
    stick_spec_tot.row(i)=stick_spec_12CO.row(i)+stick_spec_13CO.row(i)+stick_spec_C18O.row(i);  //Note:  I have the notation backwards... since IDL is backwards, (*,i) is col(i).  Fix all of these!
  }
cerr << stick_spec_tot << endl;
cerr << "^ stick_spec_tot" << endl;
//cin.get();
  cerr << "Left loop!"  << endl;
  annuli.at(0)=1e28;
  mat iten_tot=stick_spec_tot;
  iten_tot.row(0)=iten_tot.row(0)*2.5;
  mat Lum_tot = iten_tot;

  for (int i=0; i<steps; i++)
  {
    Lum_tot.row(i)=iten_tot.row(i)*annuli.at(i);
  }

  mat Flux_tot = Lum_tot/(4*3.1415926535897*pow(data->stardist,2));


//===================================
//  SCRIPT 4 of 4
//====================================
 cerr << "Script 4!" << endl;
  double r_in=rdisk.at(0);
  double r_out=rdisk.at(rdisk.n_elem-1);

  double dv=1;

  mat flux_tot_slit=Flux_tot;

//account for slit loss
cerr << "Slit loss.." << endl;
  vec slit_loss=100.17*pow(rdisk,-1.260);
  double disk_Beyond_slit = 63;

  for (int i=0; i<rdisk.n_elem; i++)
  {
    if (rdisk.at(i) > disk_Beyond_slit)
    {
      flux_tot_slit.col(i)=flux_tot_slit.col(i)*slit_loss.at(i);
    }
  }

  auto n_rings=rdisk.n_elem;

  vec dr = rdisk;

  for (int i=0; i<n_rings; i++)
  {
    if (rdisk(i) < 0.1)    {
      dr(i)=0.01;
    } else {
    dr(i)=1.0;
    }
  }
cerr << "vmax" << endl;

  vec vmax = zeros<vec>(n_rings);
  vmax=round(sqrt(887.2*Mstar/rdisk));

//construct frequency scale!

//  vec vfreq = zeros<vec>(freq.n_elem);
//  vfreq=(freq(round(freq.n_elem/2))-freq)*c/freq(round(freq.n_elem/2));

//  freq_index=where(freq, [] (double datum) {return (freq >= f_i) && (freq <= f_f);});     //these have already been done before.. redundant?

  auto total_spec=iten_tot;
  auto vel_spec=total_spec.col(0);
cerr << "END" << endl;

//=======================

  int num_lines=22;

  vec v_line5=(2032.3528-freq)*2.9979e5/2032.3528;
  vec v_line4=(2034.4083-freq)*2.9979e5/2034.4083;
  vec v_line1=(2033.4174-freq)*2.9979e5/2033.4174;
  vec v_line2=(2033.1423-freq)*2.9979e5/2033.1423;
  vec v_line3=(2032.8258-freq)*2.9979e5/2032.8258;
  vec v_line0=(2030.1586-freq)*2.9979e5/2030.1586;

  vec fco_list = {2034.9211, 2034.7209, 2034.1352,2034.0473,2032.8701, 2032.8096, 2032.5616, 2032.2011, 2030.5160, 2030.3076, 2029.6559, 2029.4227, 2029.4297, 2029.3679, 2029.2362, 2029.1276};
  
  mat v_line_arr = zeros<mat>(20,freq.n_elem);
cerr << "Loop.." << endl;
  for (int k=0; k<16; k++)
  {
    v_line_arr.row(k)=((fco_list.at(k)-freq)*2.9979e5/fco_list.at(k)).t();
  }
cerr << "Loopdone." << endl;
//lowj lines--skipped!

  mat v_line = zeros<mat>(num_lines,v_line0.n_elem);
  v_line.row(0)=v_line0.t();
  v_line.row(1)=v_line1.t();
  v_line.row(2)=v_line2.t();
  v_line.row(3)=v_line3.t();
  v_line.row(4)=v_line4.t();
  v_line.row(5)=v_line5.t();

  for (int k=0; k<16; k++) 
  {
    v_line.row(6+k)=v_line_arr.row(k);
  }

//determine grid size here!
  mat grid = zeros<mat>(15000,26);

  field <ivec> v_line_indices(22,1);
  field<vec> iten_lines(22,1);
cerr << v_line.row(0) << endl;
  v_line_indices(0 ,0) = whererow(v_line.row(0),  [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(1 ,0) = whererow(v_line.row(1),  [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(2 ,0) = whererow(v_line.row(2),  [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(3 ,0) = whererow(v_line.row(3),  [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(4 ,0) = whererow(v_line.row(4),  [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(5 ,0) = whererow(v_line.row(5),  [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(6 ,0) = whererow(v_line.row(6),  [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(7 ,0) = whererow(v_line.row(7),  [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(8 ,0) = whererow(v_line.row(8),  [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(9 ,0) = whererow(v_line.row(9),  [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(10,0) = whererow(v_line.row(10), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(11,0) = whererow(v_line.row(11), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(12,0) = whererow(v_line.row(12), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(13,0) = whererow(v_line.row(13), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(14,0) = whererow(v_line.row(14), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(15,0) = whererow(v_line.row(15), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(16,0) = whererow(v_line.row(16), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(17,0) = whererow(v_line.row(17), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(18,0) = whererow(v_line.row(18), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(19,0) = whererow(v_line.row(19), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(20,0) = whererow(v_line.row(20), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(21,0) = whererow(v_line.row(21), [] (double datum) {return ((datum > -15) && (datum < 15));}); 


  for (int i=0; i<22; i++)
  {
    cerr << "v_line_indices(" <<i<< ",0).n_elem:  " << v_line_indices(i,0).n_elem << endl;
  }

  for (int i=0; i<22; i++)
  {
    int itencols=iten_tot.n_cols;
    int v_line_num=v_line_indices(i,0).n_elem;
    mat temp = zeros<mat>(v_line_indices(i,0).n_elem,itencols);
    for (int j=0; j<v_line_num; j++)
    {
      temp.row(j)=iten_tot.row(v_line_indices(i,0).at(j));
cerr << "temp.row(" << j << "):  "  << temp.row(j) << endl;;
    }
//cin.get();
    iten_lines(i,0)=arma::sum(temp,0).t()*5.65e-3;

  }
//cin.get();
/*
  vec iten_line0=totalDim(iten_tot.row(v_line_index0),1)*5.65e-3;
  vec iten_line1=totalDim(iten_tot.row(v_line_index1),1)*5.65e-3;
  vec iten_line2=totalDim(iten_tot.row(v_line_index2),1)*5.65e-3;
  vec iten_line3=totalDim(iten_tot.row(v_line_index3),1)*5.65e-3;
  vec iten_line4=totalDim(iten_tot.row(v_line_index4),1)*5.65e-3;
  vec iten_line5=totalDim(iten_tot.row(v_line_index5),1)*5.65e-3;
  vec iten_line6=totalDim(iten_tot.row(v_line_index6),1)*5.65e-3;
  vec iten_line7=totalDim(iten_tot.row(v_line_index7),1)*5.65e-3;
  vec iten_line8=totalDim(iten_tot.row(v_line_index8),1)*5.65e-3;
  vec iten_line9=totalDim(iten_tot.row(v_line_index9),1)*5.65e-3;
  vec iten_line10=totalDim(iten_tot.row(v_line_index10),1)*5.65e-3;
  vec iten_line11=totalDim(iten_tot.row(v_line_index11),1)*5.65e-3;
  vec iten_line12=totalDim(iten_tot.row(v_line_index12),1)*5.65e-3;
  vec iten_line13=totalDim(iten_tot.row(v_line_index13),1)*5.65e-3;
  vec iten_line14=totalDim(iten_tot.row(v_line_index14),1)*5.65e-3;
  vec iten_line15=totalDim(iten_tot.row(v_line_index15),1)*5.65e-3;
  vec iten_line16=totalDim(iten_tot.row(v_line_index16),1)*5.65e-3;
  vec iten_line17=totalDim(iten_tot.row(v_line_index17),1)*5.65e-3;
  vec iten_line18=totalDim(iten_tot.row(v_line_index18),1)*5.65e-3;
  vec iten_line19=totalDim(iten_tot.row(v_line_index19),1)*5.65e-3;
  vec iten_line20=totalDim(iten_tot.row(v_line_index20),1)*5.65e-3;
  vec iten_line21=totalDim(iten_tot.row(v_line_index21),1)*5.65e-3;
*/
  mat gs = zeros<mat>(11,2);

  for (int i=0; i<11; i++)
  {
    gs.at(i,0)=i-5; //<<-5 << -4 << -3 << -2<<-1<<0<<1<<2<<3<<4<<5<<endr;
  }
  gs.col(1)=exp(-pow(gs.col(0),2)/pow((12/1.665),2))/6.1967;

  int grid_ptr=0;
//==========================================
//  RINGS LOOP                     
//==========================================
cerr << "n_rings:  " << n_rings << endl;
cerr << "steps:  " << steps  << endl;


cerr << "iten_lines(0,0).n_elem:  " << iten_lines(0,0).n_elem <<  endl;
for (int i=0; i<22; i++)
{
  cerr << "iten_lines(" << i << "):  " << iten_lines(i,0) << endl;
}
//cin.get();
  for (int j=0; j<n_rings; j++)
  {
cerr << j << endl;
    double n_seg=4*round(vmax.at(j))/dv;
cerr << "n_seg" << n_seg << endl;

//cin.get();
    vec vseg=zeros<vec>(n_seg);
{ 
    int z=-1;

    while (1)
    {
      z++;
      vseg.at(z)=vmax.at(j)-z;
      if (vseg.at(z) > -vmax.at(j)) break;
    }

    while (1)
    {
      z++;
      vseg.at(z)=vseg.at(z-1)+1;
      if (vseg.at(z) > -vmax.at(j)) break;
    }
}

    vec phase = zeros<vec>(n_seg);
    phase(0)=0;

    for (int ii=0; ii<n_seg/2; ii++)
    {
      phase.at(ii)=acos(vseg.at(ii)/vseg.at(0));
    }
    for (int ii=n_seg/2+1; ii <n_seg; ii++)
    {
      phase.at(ii)=2*3.1415926535897-acos(vseg(ii)/vseg(0));
    }
    vec dphase=zeros<vec>(n_seg);
    
    for (int i=1; i<n_seg-1; i++)
    {
      dphase.at(i)=(phase.at(i+1)-phase.at(i))/2+(phase.at(i)-phase.at(i-1))/2;
    }
    dphase.at(0)=phase.at(1);
    dphase.at(n_seg-1)=dphase.at(1);

    vec area = zeros<vec>(n_seg);
    if (j != n_rings-1)
    {
      for (int i=0; i< n_seg; i++)
      {
	area.at(i)=dphase.at(i)*(rdisk.at(j)+dr.at(j)/2)*dr.at(j);
        vel_spec=vel_spec+interpol(total_spec.col(j)*area.at(i),freq+vseg.at(i)*sin(inc)*freq/c,freq);
        grid.at(grid_ptr,0)=rdisk.at(j);
        grid.at(grid_ptr,1)=vseg.at(i);
        grid.at(grid_ptr,2)=phase.at(i);
        grid.at(grid_ptr,3)=area.at(i)*2.25;
        for (int ii=4; ii<22; ii++) grid.at(grid_ptr,ii)=iten_lines(ii,0)(j);
        grid_ptr++;
      } 
    }
    else
    {
      for (int i=0; i<n_seg; i++)
      {
        area.at(i)=dphase.at(i)*(rdisk.at(j)+dr.at(j)/2)*dr.at(j);
        vel_spec=vel_spec+interpol(total_spec.col(j)*area.at(i),freq+vseg.at(i)*sin(inc)*freq/c,freq);
        grid.at(grid_ptr,0)=rdisk.at(j);
        grid.at(grid_ptr,1)=vseg.at(i);
        grid.at(grid_ptr,2)=phase.at(i);
        grid.at(grid_ptr,3)=area.at(i)*2.25;
        for (int ii=4; ii<22; ii++) grid.at(grid_ptr,ii)=iten_lines(ii,0)(j);
        grid_ptr++;

      } 
    }
    total_spec.col(j)=vel_spec;

  }
cerr << "OUT OF GRID LOOP" << endl;
  double v_planet=5;
  double r_planet=13;
  double phase_planet=53*3.1415926535897/180;
  double planet_intens=((2*6.62-27*pow(2.9979,2)*pow(2030,3)/(exp(6.626e-27*2.9979e10*2030/(1.38e-16*1e3))-1)));
  double planet_size=0;

  for (int i=0; i<11;  i++)
  {
cerr << i << endl;
    grid.at(grid_ptr,0)=r_planet;
    grid.at(grid_ptr,1)=v_planet+gs(i,0);
    grid.at(grid_ptr,2)=phase_planet;
    grid.at(grid_ptr,3)=planet_size;
    for (int j=4; j<22; j++)
    {
      grid.at(grid_ptr,j)=0;
    }
    grid.at(grid_ptr,9)=planet_intens*gs(10,1);
    grid_ptr++;
  }

  cerr << "grid_ptr:  " << grid_ptr << endl; 
  vec centroid = freq;
cout << "grid:  " << grid << endl;
  centroid.fill(0);
  
  double omega=55*3.1415926535897/180;
  double Lc=5.13-23*2.9979247e10*4*3.1415926535897*pow((103*3.08),2)*(.05/1.16); // continuum luminosity
  ivec index1;

  int maxloop=2*max(grid.row(1));

  for (int j=0; j<maxloop; j++)
  {
    int index1n=index1.n_elem;
    cerr << "j=" << j << endl;
    field <ivec> indexn(22,1);
    index1=whererow(grid.row(1), [&] (double datum) {return ((datum <= (max(grid.row(1)-j))) && ( datum > ( max(grid.row(1))-(j+1)) ) );});
   if (index1.at(0)!=-1)
   {

     for (int k=0; k<22; k++)
     {
       vec temp2=zeros<vec>(index1.n_elem);
       for (int i=0; i<index1n; i++)  {temp2.at(i)=grid(1,index1.at(i));}
       double mean=arma::mean(temp2);
       indexn(k,0)=whererow(v_line.row(k), [&] (double datum) {return (datum <= (mean+0.5));});
     } 

    
     vec phi2=zeros<vec>(index1.n_elem);
     for (int i=0; i<index1.n_elem; i++)
      {
	phi2.at(i)=grid(2,index1.at(i));
      }

      vec rp=zeros<vec>(index1.n_elem);
      for (int i=0; i<index1.n_elem; i++)
      {
	rp.at(i)=grid(0,index1.at(i))*sqrt(pow(cos(phi2.at(i)),2)+pow(sin(phi2.at(i)),2)*pow(cos(inc),2));
      }

      vec theta=zeros<vec>(index1.n_elem);
      for (int i=0; i<phi2.n_elem; i++)
      {
	if (phi2.at(i) < 3.1415926535897) 
	{
	  theta.at(i)=acos( cos(phi2.at(i)) / sqrt(cos(pow(phi2.at(i),2)))) + pow(sin(phi2.at(i)),2)*pow(cos(inc),2);
	}
	else  
	{
	  theta.at(i)=2*3.1415926535897 - acos( cos(phi2.at(i)) / sqrt(cos(pow(phi2.at(i),2)) + pow(sin(phi2.at(i)),2)*pow(cos(inc),2)));
	}
      }
      vec deltay=rp%cos(theta+omega);
      
      for (int i=4;i<26;i++)
      {
        int siz=indexn.at(i).n_elem;
        for (int k=0;k<siz;k++)
        {
          vec tempv=zeros<vec>(index1n);
          for (int q=0;q<index1n;q++)
          {
            tempv.at(q)=grid(i,index1.at(q))*grid(3,index1.at(q));
          }
          centroid.at(indexn(i,0).at(k))=(arma::sum(tempv*deltay))/(arma::sum(tempv+Lc));
        }
      }
    }
  }
  for (int j=0; j<n_rings; j++)
  {
    total_spec.col(j)=total_spec.col(j)*arma::sum(flux_tot_slit.col(j))/arma::sum(total_spec.col(j));
  }
  vec final_spec=arma::sum(total_spec,1);
  vec inst_prof=final_spec;
  inst_prof.fill(0);
  inst_prof=exp(-pow(vfreq,2)/pow(inst_res/1.665,2));
//shift inst_prof?
  inst_prof=inst_prof/arma::sum(inst_prof);
cout << "inst_prof:  "  << inst_prof;
cout << "final_spec: "<< final_spec << endl;
cerr << "FFT..." << endl;
  cx_vec conv_spec=ifft(fft(final_spec)%fft(inst_prof))/2;
  cx_vec cent_conv=ifft(fft(centroid)%fft(inst_prof));
cout << "conv_spec:  " << conv_spec << endl;
cout << "cent_conv:  " <<  cent_conv << endl;
  delete(d);
cin.get();
  return 0;
}








int FitData::runTrials() {
  cerr << "Running trials...";
  //start with best chi-by-eye fit and run one trial
  //then each time, reset parameters to a value from the matrix and rerun.

  //set parameters for best-chi;
  layers=300;
  disk_in=13.;
  dist=1.496e13*disk_in;
  disk_out=100.0;
  v_turb=3e5;
  T_rot0_fl=2500;             //check this.... is t_rot_fl0?
  T_rot_alpha_fl=0.25;
  T_rot0_cl=T_rot0_fl;
  T_rot_alpha_cl=T_rot_alpha_fl;
  rel_lum=20;
  
  this->runTrial();
 cerr << "=============================" << endl << "END TRIAL ONE" << endl << "=========================" << endl;

  for(int i=0; i<this->numGuesses;i++) {
    //set parameters here
  cerr << "Aux trial number " << i << " begin now" << endl;
    layers=randData[0][i];
    disk_in=randData[1][i];
    dist=1.496e13*disk_in;
    disk_out=randData[2][i];
    v_turb=randData[3][i];
    T_rot0_fl=randData[4][i];
    T_rot_alpha_fl=randData[5][i];
    T_rot0_cl=T_rot0_fl;
    T_rot_alpha_cl=T_rot_alpha_fl;
    rel_lum=randData[6][i];
    this->runTrial();
  }

  return 0;
}

FitData::FitData(int numGuesses)
{
  //class variables

  FitData::numGuesses=numGuesses;


  //read in data from files
  cerr << "Reading in files...";

  std::ifstream fin;

  fin.open("inc/EinA.txt");
  for ( auto& entry : einA) {
     fin >> entry;
  }
  fin.close(); 

  cout << "EinA:  " << einA << endl;

  fin.open("inc/lambda.txt");
  for ( auto& entry : lam_ang) {
     fin >> entry;
  }
  fin.close();


  fin.open("inc/wavenum.txt");
  for ( auto& entry : wavenum) {
     fin >> entry;
  }
  fin.close();
 
  cout << "wavenum:  " << wavenum << endl;

  fin.open("inc/HD100546_luminosity.txt");
  for (auto& entry : HD100546_luminosity) {
     fin >> entry;
  }
  fin.close();

  cout << "HD100546_luminosity:  " << endl;
  cout << HD100546_luminosity << endl;
//============================
//     MOLECULAR DATA
//===========================

/*  fin.open("CO_molecdat/X12CO");
  for (int i=0; i<119; i++) {
    for (int j=0; j<7; j++) {
      for (int k=0; k<10; k++) {
        fin >> X12CO(k,j,i);
      }
    }
  }
  fin.close();

  cout << "X12CO:  " << X12CO << endl;

  fin.open("CO_molecdat/X13CO");
  for (int i=0; i<119; i++) {
    for (int j=0; j<7; j++) {
      for (int k=0; k<3; k++) {
         fin >> X13CO(k,j,i);
      }
    }
  }
  fin.close();*/

  mat temp = zeros<mat>(7,10);
  fin.open("CO_molecdat/X12CO");

  for (int k=0; k<120; k++)
  {
    for (int i=0; i<7; i++)
    {
      for (int j=0; j<10; j++)
      {
        fin >> temp(i,j);
      }
    }
    X12CO.slice(k)=temp.t();
  }

  fin.close();

//cerr << X12CO << endl;
//cin.get();

  fin.open("CO_molecdat/X13CO");
  temp = zeros<mat>(7,3);

  for (int k=0; k<120; k++)
  {
    for (int i=0; i<7; i++)
    {
      for (int j=0; j<3; j++)
      {
        fin >> temp(i,j);
      }
    }
    X13CO.slice(k)=temp.t();
  }

  fin.close();

  fin.open("CO_molecdat/XC18O");
  temp = zeros<mat>(7,1);
  for (int k=0; k<120; k++)
  {
    for (int i=0; i<7; i++)
    {
      fin >> temp(i,0);
    }
    XC18O.slice(k)=temp.t();
  }

  fin.close();



/*//transpose these to get the write dimensinos--eventually transpose instead in host files
  HD100546_luminosity=HD100546_luminosity.t();
  wavenum=wavenum.t();
  einA=einA.t();
  lam_ang=lam_ang.t();

  cout << "HD100546_luminosity:  " << HD100546_luminosity << endl;
  cout << "wavenum:  " << wavenum << endl;
  cout << "EinA:  " << einA << endl;
  cout << "lam_ang:  " << lam_ang << endl;
 
//figure out how to deal with these... might just have to access them in a different order, or reverse the cube slicing?
  cout << "X12CO:  " << X12CO << endl;
  cout << "X13CO:  " << X13CO << endl;*/



//initialize arrays with random data
  cerr << "Generating arrays...";
  this->randData[0]=FillArray(900, 100);
  this->randData[1]=FillArray(15, 6);
  this->randData[2]=FillArray(76, 50);
  this->randData[3]=FillArray(500000, 100000);      //check to be sure all of these first arguments don't need to be offset by +/- 1
  this->randData[4]=FillArray(3500,1000);
  this->randData[5]=FillArray(50,20);
  this->randData[6]=FillArray(99, 1);
cerr << "Ended randdata" << endl;
  for (int i=0; i < numGuesses; i++) {
    this->randData[5][i]=this->randData[5][i]/100;
  };
cerr << "Printing randdata..." << endl;

 //print arrays for debugging purposes
  for (int c=0; c<7; c++)
  {
    cout << "randData[" << c << "]:  ";
    for (int i=0; i<numGuesses; i++) {
      cout << this->randData[c][i] << ' ' << endl ;
    }
  };
 
 /*=================
 *
 * Collision Data                                      EVERYTHING AFTER LAYERS ----> INTO COLLISIONS!  CollisionData.cpp will need to be revisited...
 *
 *==================*/
  cerr << "Preparing collision data." << endl;
  fAX = (3.038/2.03)*einA/(wavenum % wavenum);   //be sure this is right!
  fXA = 2*fAX;
  this->runTrials();
}

FitData::~FitData()
{
  for(int i=0; i<=7;i++)
  {
    delete this->randData[i];
  }
  delete this->randData;
}

int createInput() 
{
  

  return 1;
} 

int readInput(string inpFile)
{

  return 1;
}

int main(int argc, char* argv[]) 
{
  cout << "Begin!" << endl;
  data = new FitData(atoi(argv[1]));
  data->runTrials();
  cerr << "Trials ran!" << endl;
  delete data;
  return 1;
}
