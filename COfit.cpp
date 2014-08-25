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

//        add input file
//        fix g-factors (?)
//        fix grid



using namespace std;
using namespace idlarma;FitData* data;

constexpr double FitData::vib_einA[];
constexpr double FitData::Bv12[];
constexpr double FitData::Bv13[];

constexpr int

int FitData::numGuesses;

//returns an array of size numGuesses filled with random numbers on the interval (offset, modulus+offset)
//used to generate the matrix of guess data
double* FitData::FillArray(int modulus, int offset)
{
  double* array;
  array  = new double[numGuesses];
  for(int i=0; i<numGuesses;i++) {
    array[i]=rand() % modulus + offset;
  }
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
      for (auto k=0; k<10; k++) {
        for (auto q=0; q<11; q++)
        {
          d->rate_eqtn(i,0)(k+11,q,j)=einA.at(k,q);//k+11, q,j,i,   k ,q  //why no zero here
        }
      }
    }
  }

  for (auto j=0; j<steps; j++) 
  {

    for (int z=0; z<9; z++)
    {
      for (int q=0; q<layers; q++) 
      {
        d->rate_eqtn.at(j,0).at(z+1,z+1,q)=-vib_einA[z];  
        d->rate_eqtn.at(j,0).at(z+2,z+1,q)=vib_einA[z+1]; 
      }
    }
   
    d->rate_eqtn.at(j,0)(span(1),span(0),span::all).fill(vib_einA[0]);   
    d->rate_eqtn.at(j,0)(span(10),span(10),span::all).fill(-vib_einA[9]); 

    for (auto i=0; i<9; i++)
    {
      d->rate_eqtn.at(j,0)(span(i+11),span(i+11),span::all).fill(-sum(einA.row(i)));
    }
  }

  for (auto j=0; j<steps; j++)
  {
    for (int i=0; i<8; i++)
    {
      d->rate_eqtn.at(j,0).subcube(span(i+1),span(i+1),span::all).fill(-vib_einA[i]);
      d->rate_eqtn.at(j,0).subcube(span(i+2),span(i+1),span::all).fill(vib_einA[i+1]);
    } 

    d->rate_eqtn.at(j,0).subcube(span(1),span(0),span::all).fill(vib_einA[0]);
    d->rate_eqtn.at(j,0).subcube(span(9),span(9),span::all).fill(-vib_einA[8]);
  }

  d->rate_eqtn2=d->rate_eqtn;   //Clone colldata for the second pass without collisions
   

  T_rot_fl = T_rot0_fl *pow((1.5E13/rdisk),T_rot_alpha_fl);
  
  ivec T_rot_index=where(T_rot_fl, [] (double datum) {return datum >= 3500;});

  if (T_rot_index.at(0) != -1) {
    for (auto elem : T_rot_index )
    {
      T_rot_fl[elem]=3500;
    }
  }


  T_rot_cl=T_rot0_cl * pow((1.496E13/rdisk), T_rot_alpha_cl);
  H_den=H_den0 * pow((1.496E13/rdisk),H_den_alpha);
  T_rot_index = where(T_rot_cl, [] (double datum) {return datum >= 3500;});
  if (T_rot_index.at(0) != -1) {
    for (auto elem : T_rot_index)
    {
      T_rot_cl[elem]=3500;
    }
  }
    vec b_tot = sqrt(2*1.38e-16*6.02e23*T_rot_fl/28 + pow(v_turb,2));
   fvec b_tot2 = sqrt(2*1.38e-16*6.02e23*conv_to<fvec>::from(T_rot_fl)/28 + pow(v_turb,2));

if (rel_lum <= 1e-3) {cerr << "REL LUM TRIGGERED " << endl; cin.get();}

    vec tau = linspace<vec>(0,19.99,2000); 
    vec F_tau=zeros<vec>(2000);

    //==============CALCULATE DFDT BY INTEGRAATION USING TRAPEZOID RULE================
   
    float sum;
    float n = 1000;
    float a = 0;
    float b = 5;
    float delta=(b-a)/n;

    for (int i=0; i<2000; i++) 
    {
      sum = 0.5 * ( (1 - exp(-tau(i)*exp(-pow(a,2)))) + (1 - exp(-tau(i)*exp(-pow(b,2)))) );
      for (int j=1; j<n; j++)  
      {
        float x=a+delta*j;
        sum=sum+(1-exp(-tau(i)*exp(-pow(x,2))));
      }
      sum = sum*delta;
      F_tau(i)=sum;
    }

    
    vec dFdt=deriv2(tau,F_tau);
      fcube Nv_coll;
      fmat tot_col_fluor;
      fmat tot_col_fluor_back;                      //move, perhaps, to colldata?

      fcube Nv_nocoll;
      fmat tot_col_fluor_nocoll;


      d->Nv = zeros<fcube>(21,layers,steps);


//==========================================================
//   BIG COLLISION LOOP
//=========================================================
    for (int coll_loop=0; coll_loop<2; coll_loop++) 
    {
      mat U;
      mat s;
      mat v;
       
      vec sol;
      ivec Nv_index;
d->Nv.fill(0);
d->dwdn.fill(0);
d->tau_0.fill(0);
for (int zz=0; zz<10; zz++) {d->g.at(zz,0).fill(0);}

      double gsum;

      for (int k=0; k<steps; k++) 
      {

        mat Fuv = HD100546_luminosity / ( 4*datum::pi*pow(rdisk[k],2) );

       
       //=============COLLISIONS=============

       if (coll_loop == 0)
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
	      d->tau_0.subcube(span(i),span::all,span(j))=7.55858e12 * 0.02654 * conv_to<fvec>::from(arma::sum(d->Nv.slice(k).submat(span(0,11),span(0,j)),1)) % conv_to<fvec>::from(fXA.submat(span(i),span::all).t()) % conv_to<fvec>::from(lam_ang.submat(span(i),span::all).t())*1e-8/(sqrt(datum::pi)*b_tot2[k]);

	     }
	   }

           auto max_iter=d->tau_0.slice(j).row(0).n_elem;
           

	   for (auto ii=0; ii<max_iter; ii++)
	   {

	     if (d->tau_0.at(0,ii,j) < 20)
	     {
               double rnd=round(d->tau_0.at(0,ii,j)*10)/10.0;
	       ivec dFdt_0_index=where(tau, [rnd] (double elem) {return (abs(elem - rnd) < 0.001);});
	       auto count = dFdt_0_index.size();
	       if (dFdt_0_index.at(0) != -1)
	       {
	         d->dFdt_0.subcube(span::all,span(ii),span(j)).fill(dFdt.at(dFdt_0_index.at(0)));    //weird line!!!!!
	       }
             
	     }

	     else 
	     {
	       d->dFdt_0.slice(j).col(ii).fill((1/(2*d->tau_0.at(0,ii,j)*sqrt(log(d->tau_0.at(0,ii,j))))));    //(span::all,span(i),span(j))=1/(2*tau_0.at(0,ii,j))*sqrt(alog(tau_0.at(0,ii,j)));   CHECK THIS WITH SOURCE

	     }

	   }

	   d->dwdn.slice(j)=d->dFdt_0.slice(j)*.02654*2.0 % (conv_to<fmat>::from(lam_ang)*1e-4) % (conv_to<fmat>::from(lam_ang)*1e-8) % conv_to<fmat>::from(fXA)/(sqrt(datum::pi)*c*1e5);
	   for (int ii=0; ii<10; ii++) 
	   {
	     d->g.at(ii,0).slice(k).col(j) = (d->dwdn.slice(j).row(ii) * datum::pi  % conv_to<frowvec>::from(Fuv.row(ii)) / (hc * conv_to<frowvec>::from(wavenum.row(ii)))).t();  //check this to be sure the constants are filling in right...
           }

//=======================
//  Add in G-factors
//======================

	   float gsum=0;

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
	       d->rate_eqtn.at(k,0).at(ii,11+i,j)+=d->g.at(i,0).at(ii,j,k); //ii 11+i j k
	     }
	   }
	   d->rate_eqtn.at(k,0).at(1,10,j)=0;
	   
//===================
//  SOLVE
//==================  
	   vec z = zeros<vec>(21);
	   z.at(20)=1;
	     

	   solve(sol,d->rate_eqtn.at(k,0).slice(j).t(),z);
if (k>0){
fmat temp = zeros<fmat>(10,12);
for (int zz=0; zz<10; zz++)
{temp.row(zz)=d->g.at(zz,0).slice(k).col(j).t();}

	   d->Nv.slice(k).col(j)= conv_to<fvec>::from(sol);  //check this to be sure the array is filling in the right dimension

           Nv_index=where(d->Nv.slice(k).col(j), [] (double datum) {return datum < 0;});
	   if (Nv_index.at(0) != -1) 
           {
             for (int jj=0; jj<Nv_index.n_elem; jj++) 
             {
               d->Nv.at(Nv_index[jj],j,k)=0;
             }
           }
 	 }
       }

skip_fluorcalc:

      if (coll_loop==0) {
        Nv_coll=d->Nv;
        tot_col_fluor = totalDimf(d->Nv*7.55858e12,2).t();      //check all uses of accu/total in teh code!!!!!! dimension specification!
        tot_col_fluor_back=tot_col_fluor;
        d->rate_eqtn=d->rate_eqtn2;
      }
      if (coll_loop==1) {
        Nv_nocoll=d->Nv;
        tot_col_fluor_nocoll = totalDimf(d->Nv*7.55858e12,2).t();
      }
    }

//=========================================================================
// Angle of Incidence Correction (tweak for each star--use input file!)
//========================================================================
  vec m_disk=1.5*(5.59647e-10)*sqrt(T_rot_fl / Mstar)%pow(rdisk,0.5);
  vec m_uv = (5.59647e-10)*sqrt(T_rot_fl/Mstar)%sqrt(rdisk);

  vec phi = -atan((m_uv-m_disk)/(1+m_uv%m_disk));
  phi.at(0)=datum::pi/2;                         //standardize pi!
  
  auto xi = tot_col_fluor.n_cols-1;
  auto yi = tot_col_fluor.n_rows-1;

  for (int j=1; j<yi+1; j++)                    //looks like I need to transpose indices on everything involving tot_col_fluor 
  {
    tot_col_fluor.row(j)=sin(phi.at(j))*tot_col_fluor.row(j);
    tot_col_fluor_nocoll.row(j)=sin(phi.at(j))*tot_col_fluor_nocoll.row(j);
  }

  fmat tot_col=tot_col_fluor(span::all,span(0,9));
  fmat tot_col_coll=tot_col;

  double r = dist/1.496e13;
  rdisk=rdisk/1.496e13;            //convert to AU




//===================================================================================================
//  MOLECULAR DATA
//==================================================================================================
  

  cube N_12CO_vj = zeros<cube>(tot_col.n_cols-1,120,tot_col.n_rows);                 //flipping these indices...
  cube N_13CO_vj = zeros<cube>(3,120,tot_col.n_rows);
  cube N_C18O_vj = zeros<cube>(1,120,tot_col.n_rows);


  rowvec Jup;
  rowvec Jdn;
  rowvec wvn;
  rowvec EinA;

  for (int i=0; i<steps; i++)
  {


    //=======================
    //   12CO
    //======================


    for (int j=0; j< tot_col.n_cols-1; j++)
    {
      Jup = X12CO(span(j),span(3),span::all);          //because the indices are reversed in c++, do all these first/second indices flipped.
      Jdn = X12CO(span(j),span(1),span::all);
      wvn = X12CO(span(j),span(6),span::all);
      EinA= X12CO(span(j),span(4),span::all);

      N_12CO_vj.slice(i).row(j)=tot_col_fluor.at(i,j+1) * (2*Jup+1) % exp(-hck*Bv12[j]*Jup % (Jup+1)/T_rot_fl[i])/(T_rot_fl[i]/(hck*Bv12[j])) 
+ tot_col_coll.at(i,j+1)*(2*Jup+1) % exp(-hck*Bv12[j] * Jup % (Jup+1)/T_rot_cl[i])/(T_rot_cl[i]/(hck*Bv12[j]));

    }
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

  int freq_size=floor(log10(f_f/f_i)/log10(1+v/(3*c)));
  vec freq = zeros<vec>(freq_size);
  freq.at(0)=f_i;
  vec vfreq = freq;
  vfreq.fill(0);

  for (int i=1; i<freq_size; i++)
  {
    freq.at(i)=f_i*(pow(1+v/(3*c),i));
  }
  vfreq=(freq(round(freq.n_elem/2))-freq)*c/freq(round(freq.n_elem/2));                      //check

  ivec freq_index=where(freq, [] (double datum) {return (datum >= f_i) && (datum <= f_f);});
  //define annuli

  vec annuli = zeros<vec>(rdisk.n_elem);
  for (int i=0; i<steps-1; i++)
  {
    annuli.at(i)=datum::pi*(pow(rdisk.at(i+1),2)-pow(rdisk.at(i),2));
  }

  annuli(steps-1)=datum::pi*(pow(rdisk.at(steps-1)+1,2)-pow(rdisk.at(steps-1),2));

  annuli=annuli*2.25e26;
  mat stick_spec_12CO = zeros<mat>(freq_size,steps);
  mat stick_spec_13CO = stick_spec_12CO;
  mat stick_spec_C18O = stick_spec_12CO;
  mat stick_spec_tot  = stick_spec_12CO;

  double A0;
  double A1;
  double A2;

  //======================
  //  X12CO
  //=====================
  for (int i=0; i<steps; i++)  //Loop over annuli
  {

    double btotcomp=b_tot.at(i)/cexp;
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

	  A0=N_12CO_vj.at(j,k,i)*hc*A1*EinA.at(k);                         //should the first two indices be reversed here?
	  A2=btotcomp*A1;
	  stick_spec_12CO.col(i)+=(A0/(rpi*A2)) * exp (-pow(((A1-freq)/A2),2));

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
	    stick_spec_13CO.col(i)+=(A0/(rpi*A2)) * exp(-pow(((A1-freq)/A2),2));

          }
      }    
    }
    //=============================
    //  XC180
    //============================ 
     
      Jup = XC18O(span(0),span(3),span::all);
      Jdn = XC18O(span(0),span(1),span::all);
      wvn = XC18O(span(0),span(6),span::all);
      EinA= XC18O(span(0),span(4),span::all);

      for (int k=0; k<XC18O.n_slices; k++)
      {
          A1=wvn.at(k);
          if ((A1 >= f_i) && (A1 <= f_f))
          {
	    A0=N_1C18O_vj.at(0,k,i)*hc*A1*EinA.at(k);                         //should the first two indices be reversed here?
	    A2=btotcomp*A1;
	    stick_spec_C18O.col(i)+=(A0/(rpi*A2)) * exp(-pow(((A1-freq)/A2),2));
          }
      }



    stick_spec_tot.col(i)=stick_spec_12CO.col(i)+stick_spec_13CO.col(i)+stick_spec_C18O.col(i);  //Note:  I have the notation backwards... since IDL is backwards, (*,i) is col(i).  Fix all of these!

  }
  
  annuli.at(0)=1e28;
  mat iten_tot=stick_spec_tot;
  iten_tot.col(0)=iten_tot.col(0)*2.5;
  mat Lum_tot = iten_tot;
  for (int i=0; i<steps; i++)
  {
    Lum_tot.col(i)=iten_tot.col(i)*annuli.at(i);
  }

  mat Flux_tot = Lum_tot/(4*datum::pi*pow(data->stardist,2));


//===================================
//  SCRIPT 4 of 4
//====================================
  double r_in=rdisk.at(0);
  double r_out=rdisk.at(rdisk.n_elem-1);

  double dv=1;

  mat flux_tot_slit=Flux_tot;

//account for slit loss
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
    } else 
      {  if (rdisk(i) <1.0) {
        dr(i)=0.1;
      } else
        { dr(i) = 1.0;
        }
      }
  }
  vec vmax = zeros<vec>(n_rings);
  vmax=round(sqrt(887.2*Mstar/rdisk));

//construct frequency scale!

//  vec vfreq = zeros<vec>(freq.n_elem);
//  vfreq=(freq(round(freq.n_elem/2))-freq)*c/freq(round(freq.n_elem/2));

//  freq_index=where(freq, [] (double datum) {return (freq >= f_i) && (freq <= f_f);});     //these have already been done before.. redundant?

  mat total_spec=iten_tot;
  vec vel_spec=total_spec.col(0);

//=======================

  int num_lines=22;

  vec v_line5=(2032.3528-freq)*2.9979e5/2032.3528;
  vec v_line4=(2034.4083-freq)*2.9979e5/2034.4083;
  vec v_line1=(2033.4174-freq)*2.9979e5/2033.4174;
  vec v_line2=(2033.1423-freq)*2.9979e5/2033.1423;
  vec v_line3=(2032.8258-freq)*2.9979e5/2032.8258;
  vec v_line0=(2030.1586-freq)*2.9979e5/2030.1586;

  vec fco_list = {2034.9211, 2034.7209, 2034.1352,2034.0473,2032.8701, 2032.8096, 2032.5616, 2032.2011, 2030.5160, 2030.3076, 2029.6559, 2029.4427, 2029.4297, 2029.3679, 2029.2362, 2029.1276};
  
  mat v_line_arr = zeros<mat>(20,freq.n_elem);

  for (int k=0; k<16; k++)
  {
    v_line_arr.row(k)=((fco_list.at(k)-freq)*2.9979e5/fco_list.at(k)).t();
  }
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
  mat grid = zeros<mat>(1,26);

  field <ivec> v_line_indices(22,1);
  field<vec> iten_lines(22,1);
  v_line_indices(0 ,0) = whererow(v_line.row(0 ), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(1 ,0) = whererow(v_line.row(1 ), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(2 ,0) = whererow(v_line.row(2 ), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(3 ,0) = whererow(v_line.row(3 ), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(4 ,0) = whererow(v_line.row(4 ), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(5 ,0) = whererow(v_line.row(5 ), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(6 ,0) = whererow(v_line.row(6 ), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(7 ,0) = whererow(v_line.row(7 ), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(8 ,0) = whererow(v_line.row(8 ), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
  v_line_indices(9 ,0) = whererow(v_line.row(9 ), [] (double datum) {return ((datum > -15) && (datum < 15));}); 
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


// ITEN LINES

  for (int i=0; i<22; i++)
  {
    int itencols=iten_tot.n_cols;
    int v_line_num=v_line_indices(i,0).n_elem;
    mat temp = zeros<mat>(v_line_indices(i,0).n_elem,itencols);
    for (int j=0; j<v_line_num; j++)
    {
      temp.row(j)=iten_tot.row(v_line_indices(i,0).at(j));
    }
    iten_lines(i,0)=arma::sum(temp,0).t()*5.65e-3;
;
  }
// gs

  mat gs = zeros<mat>(11,2);

  for (int i=0; i<11; i++)
  {
    gs.at(i,0)=i-5; //<<-5 << -4 << -3 << -2<<-1<<0<<1<<2<<3<<4<<5<<endr;
  }
  gs.col(1)=exp(-pow(gs.col(0),2)/pow((12/1.665),2))/6.1967;
  int grid_ptr=1;
//==========================================
//  RINGS LOOP                     
//==========================================

  for (int j=0; j<n_rings; j++)
  {
    double n_seg=4*round(vmax.at(j))/dv;
    
    grid.resize(grid.n_rows+n_seg,26);

    vec vseg=zeros<vec>(n_seg);

    int z=-1;

    while (1)
    {
      z++;
      vseg.at(z)=vmax.at(j)-z;
      if (vseg.at(z) <= -vmax.at(j)) break;
    }
    while (1)
    {
      z++;
      vseg.at(z)=vseg.at(z-1)+1;
      if (vseg.at(z) == vmax.at(j)-1) break;
    }

    vec phase = zeros<vec>(n_seg);
    phase(0)=0;

    for (int ii=0; ii<n_seg/2; ii++)
    {
      phase.at(ii)=acos(vseg.at(ii)/vseg.at(0));
    }
   
    for (int ii=n_seg/2; ii <n_seg; ii++)
    {
      phase.at(ii)=2*datum::pi-acos(vseg(ii)/vseg(0));
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
        grid.at(grid_ptr,3)=area.at(i)*2.25e26;
        for (int ii=4; ii<26; ii++) { grid.at(grid_ptr,ii)=iten_lines(ii-4,0)(j);}
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
        grid.at(grid_ptr,3)=area.at(i)*2.25e26;
        for (int ii=4; ii<26; ii++) grid.at(grid_ptr,ii)=iten_lines(ii-4,0)(j);
        grid_ptr++;
      } 
    }

    total_spec.col(j)=vel_spec;

  }
  ivec index_grid=where(grid.col(0), [&] (double datum) {return ((datum <= rdisk.at(0)) && (datum > 0));});

  vec  grid_tmp=zeros<vec>(index_grid.n_elem);
  for (int i=0; i<index_grid.n_elem; i++)
  { 
    grid_tmp.at(i)=grid(index_grid.at(i),2);
  }
 ivec index_grid2=where(grid_tmp,[] (double datum) {return (datum > arma::datum::pi);});
 for (int i=0; i<index_grid2.n_elem; i++)
 {
   grid.at(index_grid.at(index_grid2.at(i)),3)=0;
 }
  double v_planet=5;
  double r_planet=13;
  double phase_planet=53*datum::pi/180;
  double planet_intens=((2*6.62e-27*pow(2.9979e10,2)*pow(2030,3)/(exp(6.626e-27*2.9979e10*2030/(1.38e-16*1e3))-1)));
  double planet_size=0;
  grid.resize(grid.n_rows+11,26);
  for (int i=0; i<11;  i++)
  {
    grid.at(grid_ptr,0)=r_planet;
    grid.at(grid_ptr,1)=v_planet+gs(i,0);
    grid.at(grid_ptr,2)=phase_planet;
    grid.at(grid_ptr,3)=planet_size;
    for (int j=4; j<22; j++)
    {
      grid.at(grid_ptr,j)=0;
    }
    grid.at(grid_ptr,9)=planet_intens*gs(i,1);
    grid_ptr++;
  }

  vec centroid = freq;
  centroid.fill(0);
  
  double omega=55*datum::pi/180;
  double Lc=   8.3837025e28;   //pow(5.13,-23*2.9979247e10*4*3.1415926535897*pow((103*3.08),2)*(.05/1.16)); // continuum luminosity
  ivec index1;
  
  vec gridrow=grid.col(1);
  double maxloop=2*gridrow.max()+1;

//======================================
//   CENTROID LOOP
//=====================================
  for (int j=0; j<maxloop; j++)
  {
    field <ivec> indexn(22,1);
    index1=where(grid.col(1), [&] (double datum) {return ((datum <= (max(grid.col(1)))-j) && ( datum > ( max(grid.col(1))-(j+1)) ) );});
   int index1n=index1.n_elem;
   if (index1.at(0)!=-1)
   {

//=================
//   i0-i21!
//=================

       vec temp2=zeros<vec>(index1n);

       for (int i=0; i<index1n; i++)  {temp2.at(i)=grid(index1.at(i),1);}

     double mean=arma::mean(temp2);

     for (int k=0; k<22; k++)
     {
       indexn(k,0)=whererow(v_line.row(k), [&] (double datum) {return (datum > (mean-0.5) && datum <= (mean+0.5));});
     } 

     vec phi2=zeros<vec>(index1n);
     for (int i=0; i<index1n; i++)
     {
       phi2.at(i)=grid(index1.at(i),2);
     }

      vec rp=zeros<vec>(index1n);
      for (int i=0; i<index1n; i++)
      {
	rp.at(i)=grid(index1.at(i),0)*sqrt(pow(cos(phi2.at(i)),2)+pow(sin(phi2.at(i)),2)*pow(cos(inc),2));
      }

      vec theta=zeros<vec>(index1n);
      for (int i=0; i<phi2.n_elem; i++)
      {
	if (phi2.at(i) < datum::pi) 
	{
	  theta.at(i)=acos( cos(phi2.at(i)) / sqrt(pow(cos(phi2.at(i)),2) + pow(sin(phi2.at(i)),2)*pow(cos(inc),2)));
	}
	else  
	{
	  theta.at(i)=2*datum::pi - acos( cos(phi2.at(i)) / sqrt(pow(cos(phi2.at(i)),2) + pow(sin(phi2.at(i)),2)*pow(cos(inc),2)));
	}
      }
      vec deltay=rp%cos(theta+omega);
      //first loop--look over each grid index 
        //second loop-loop over each indexn
        for (int z=0; z<22; z++) 
        {
          int siz=indexn.at(z).n_elem;
          vec tempv=zeros<vec>(index1n);
          for (int q=0; q<index1n; q++)
	  {
	      tempv.at(q)=grid(index1.at(q),z+4)*grid(index1.at(q),3);
	  }
	  for (int k=0;k<siz;k++)
	  {
	    centroid.at(indexn(z,0).at(k))=(arma::sum(tempv%deltay))/(arma::sum(tempv)+Lc);
	  }
        }
    }
  }
ivec notzero=where(centroid, [] (double datum) {return (datum != 0);});
vec final1_spec=arma::sum(total_spec,1);
  for (int j=0; j<n_rings; j++)
  {
    total_spec.col(j)=total_spec.col(j)*arma::sum(flux_tot_slit.col(j))/arma::sum(total_spec.col(j));
  }
  vec final_spec=arma::sum(total_spec,1);
  vec inst_prof=final_spec;
  inst_prof.fill(0);
  inst_prof=exp(-pow(vfreq,2)/pow(inst_res/1.665,2));
  inst_prof=shift(inst_prof,round(vfreq.n_elem)/2);
inst_prof=inst_prof/arma::sum(inst_prof);

  cx_vec conv_spec=ifft(fft(final_spec)%fft(inst_prof))/2;
  cx_vec cent_conv=ifft(fft(centroid)%fft(inst_prof));
cent_conv=cent_conv*as_scalar(arma::sum(abs(centroid)))/as_scalar(arma::sum(abs(cent_conv)));
conv_spec=conv_spec*as_scalar(arma::accu(flux_tot_slit))/as_scalar(arma::sum(conv_spec));
  delete(d);
cin.get();
  return 0;
}








int FitData::runTrials() {
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
  cerr << "Aux trial number " << i+1 << " begin now" << endl;
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
 

  fin.open("inc/HD100546_luminosity.txt");
  for (auto& entry : HD100546_luminosity) {
     fin >> entry;
  }
  fin.close();


//============================
//     MOLECULAR DATA
//===========================


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






//initialize arrays with random data
  cerr << "Generating arrays...";
  this->randData[0]=FillArray(900, 100);
  this->randData[1]=FillArray(15, 6);
  this->randData[2]=FillArray(76, 50);
  this->randData[3]=FillArray(500000, 100000);      //check to be sure all of these first arguments don't need to be offset by +/- 1
  this->randData[4]=FillArray(3500,1000);
  this->randData[5]=FillArray(50,20);
  this->randData[6]=FillArray(99, 1);
  for (int i=0; i < numGuesses; i++) {
    this->randData[5][i]=this->randData[5][i]/100;
  };

 //print arrays for debugging purposes
  for (int c=0; c<7; c++)
  {
    for (int i=0; i<numGuesses; i++) {
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
  data = new FitData(atoi(argv[1]));
  data->runTrials();
  delete data;
  return 1;
}
