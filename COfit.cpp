#include "COfit.h"
#include "CollData.cpp"
#include <iostream>
#include <cstdlib>

#include <mpi.h>

using namespace std;

/*constexpr double FitData::vib_einA[];
constexpr double FitData::Bv12[];
constexpr double FitData::Bv13[];
vec FitData::freq;*/

constexpr double CollData::fco_list[];

//int FitData::numGuesses;

//returns an array of size numGuesses filled with random numbers on the interval (offset, modulus+offset)
//used to generate the matrix of guess data
double* FitData::FillArray(int modulus, int offset)
{
  double* array;
  array  = new double[numGuesses-2];
  for(int i=0; i<numGuesses-2;i++) {
    array[i]=rand() % modulus + offset;
  }
  return array;
}

int FitData::runTrial(double layers, double disk_in, double disk_out, double v_turb, double T_rot0_fl, double T_rot_alpha_fl, double rel_lum) {

//============================================
// MAIN MODEL-FITTING FUNCTION
// -This is the part to make parallel
//============================================

//========================================
// Set up CollData!
// The CollData object contains all of the variables specific to this trial!
// CollData::CollData sets up static rate equation varibles, v_line indices, and another static varibles
//=========================================

  double dist = 1.496e13*disk_in;
  double T_rot0_cl=T_rot0_fl;
  double T_rot_alpha_cl=T_rot_alpha_fl;

  CollData* d = new CollData(layers,disk_in,disk_out,v_turb,T_rot0_fl,T_rot_alpha_fl,rel_lum);


//============================
//Divide into annulli
//=============================

/*  d->r=disk_in;  

  double r_index_a = 0;
  double r_index_b = 0;
  double r_index_c = 0;

  d->lum = HD100546_luminosity*rel_lum;
  d->ranulli.push_back(disk_in);
  if (d->r<0.1)
  {
    while ((d->ranulli.back() < 1) && (d->ranulli.back() < disk_out))
    {
      r_index_a=r_index_a+1;
      d->ranulli.push_back(disk_in+0.01*r_index_a);
    }
  }
  double maxra = d->ranulli.back();
  
  if ((maxra < 1) && (maxra >= .1))
  {
    while ((d->ranulli.back() < 1.0) || (d->ranulli.back() < disk_out))
    {
     r_index_b++;
     d->ranulli.push_back(maxra+0.1*r_index_b);
    }
  }

  double maxrb=d->ranulli.back();
  if ((maxrb <= disk_out) && (maxrb >= 1.0))
  {
    while  (d->ranulli.back() < disk_out)
    {
      r_index_c++;
      d->ranulli.push_back(maxrb+1.0*r_index_c);
    }
  }

  steps=d->ranulli.size();  //check that this size function is correct
  for (int i=0; i < steps; i++ )
  {
    d->rdisk.at(i)=d->ranulli.at(i)*1.496E13;
  }*/
//benchmark here!  Make sure these results are the same from the IDL code and the C++ code!
 
//=====================
//  Set up temperature profile
// ====================
  d->T_rot_fl = T_rot0_fl *pow((1.5E13/d->rdisk),T_rot_alpha_fl);
  
  ivec T_rot_index=where(d->T_rot_fl, [] (double datum) {return datum >= 3500;});

  if (T_rot_index.at(0) != -1) {
    for (auto elem : T_rot_index )
    {
      d->T_rot_fl[elem]=3500;
    }
  }


  d->T_rot_cl=T_rot0_cl * pow((1.496E13/d->rdisk), T_rot_alpha_cl);
  d->H_den=H_den0 * pow((1.496E13/d->rdisk),H_den_alpha);
  T_rot_index = where(d->T_rot_cl, [] (double datum) {return datum >= 3500;});
  if (T_rot_index.at(0) != -1) {
    for (auto elem : T_rot_index)
    {
      d->T_rot_cl[elem]=3500;
    }
  }
    vec b_tot = sqrt(2*1.38e-16*6.02e23*d->T_rot_fl/28 + pow(v_turb,2));
   fvec b_tot2 = sqrt(2*1.38e-16*6.02e23*conv_to<fvec>::from(d->T_rot_fl)/28 + pow(v_turb,2));

if (rel_lum <= 1e-3) {cerr << "REL LUM TRIGGERED " << endl; cin.get();}


    //==============CALCULATE DFDT BY INTEGRAATION USING TRAPEZOID RULE================
    for (int i=0; i<2000; i++) 
    {
      d->sum = 0.5 * ( (1 - exp(-d->tau(i)*exp(-pow(d->a,2)))) + (1 - exp(-d->tau(i)*exp(-pow(d->b,2)))) );
      for (int j=1; j<d->n; j++)  
      {
        float x=d->a+d->delta*j;
        d->sum=d->sum+(1-exp(-d->tau(i)*exp(-pow(x,2))));
      }
      d->sum = d->sum*d->delta;
      d->F_tau.at(i)=d->sum;
    }
      vec dFdt=deriv2(d->tau,d->F_tau);
      d->Nv = zeros<fcube>(21,layers,d->steps);


//==========================================================
//   BIG COLLISION LOOP
//=========================================================
    for (int coll_loop=0; coll_loop<2; coll_loop++) 
    {
      vec sol;
      ivec Nv_index;
      d->Nv.fill(0);
      d->dwdn.fill(0);
      d->tau_0.fill(0);
      for (int zz=0; zz<10; zz++) {d->g.at(zz,0).fill(0);}

      double gsum;

      for (int k=0; k<d->steps; k++) 
      {

        mat Fuv = d->lum / ( 4*datum::pi*pow(d->rdisk[k],2) );

       
       //=============COLLISIONS=============

       if (coll_loop == 0)
       {
	 auto k_HCO_dn = (7.57e-15*d->T_rot_cl[k]/(1-exp(-3084/d->T_rot_cl[k])))*d->H_den[k];
	 auto k_HCO_up=k_HCO_dn*exp(-3084/d->T_rot_cl[k]);

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
	       ivec dFdt_0_index=where(d->tau, [rnd] (double elem) {return (abs(elem - rnd) < 0.001);});
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
        d->Nv_coll=d->Nv;
        d->tot_col_fluor = totalDimf(d->Nv*7.55858e12,2).t();      //check all uses of accu/total in teh code!!!!!! dimension specification!
        d->tot_col_fluor_back=d->tot_col_fluor;
        d->rate_eqtn=d->rate_eqtn2;
      }
      if (coll_loop==1) {
        d->Nv_nocoll=d->Nv;
        d->tot_col_fluor_nocoll = totalDimf(d->Nv*7.55858e12,2).t();
      }
    }
//=========================================================================
// Angle of Incidence Correction (tweak for each star--use input file!)
//========================================================================
  d->m_disk=1.5*(5.59647e-10)*sqrt(d->T_rot_fl / Mstar)%pow(d->rdisk,0.5);
  d->m_uv = (5.59647e-10)*sqrt(d->T_rot_fl/Mstar)%sqrt(d->rdisk);

  d->phi = -atan((d->m_uv-d->m_disk)/(1+d->m_uv%d->m_disk));
  d->phi.at(0)=datum::pi/2;                         //standardize pi!
  
  auto xi = d->tot_col_fluor.n_cols-1;
  auto yi = d->tot_col_fluor.n_rows-1;

  for (int j=1; j<yi+1; j++)                    //looks like I need to transpose indices on everything involving tot_col_fluor 
  {
    d->tot_col_fluor.row(j)=sin(d->phi.at(j))*d->tot_col_fluor.row(j);
    d->tot_col_fluor_nocoll.row(j)=sin(d->phi.at(j))*d->tot_col_fluor_nocoll.row(j);
  }

  fmat tot_col=d->tot_col_fluor(span::all,span(0,9));
  fmat tot_col_coll=tot_col;

  double r = dist/1.496e13;
  d->rdisk=d->rdisk/1.496e13;            //convert to AU




//===================================================================================================
//  MOLECULAR DATA
//==================================================================================================
  

  d->N_12CO_vj = zeros<cube>(tot_col.n_cols-1,120,tot_col.n_rows); 
  d->N_13CO_vj = zeros<cube>(3,120,tot_col.n_rows);
  d->N_C18O_vj = zeros<cube>(1,120,tot_col.n_rows);

  for (int i=0; i<d->steps; i++)
  {

    //=======================
    //   12CO
    //======================


    for (int j=0; j< tot_col.n_cols-1; j++)
    {
      d->Jup = X12CO(span(j),span(3),span::all);          //because the indices are reversed in c++, do all these first/second indices flipped.
      d->Jdn = X12CO(span(j),span(1),span::all);
      d->wvn = X12CO(span(j),span(6),span::all);
      d->EinA= X12CO(span(j),span(4),span::all);

      d->N_12CO_vj.slice(i).row(j)=d->tot_col_fluor.at(i,j+1) * (2*d->Jup+1) % exp(-hck*Bv12[j]*d->Jup % (d->Jup+1)/d->T_rot_fl[i])/(d->T_rot_fl[i]/(hck*Bv12[j])) 
+ tot_col_coll.at(i,j+1)*(2*d->Jup+1) % exp(-hck*Bv12[j] * d->Jup % (d->Jup+1)/d->T_rot_cl[i])/(d->T_rot_cl[i]/(hck*Bv12[j]));

    }
  //====================
  //  13CO
  //====================
  for (int j=0; j<3; j++)
  {
    d->Jup = X13CO(span(j),span(3),span::all);
    d->Jdn = X13CO(span(j),span(1),span::all);
    d->wvn = X13CO(span(j),span(6),span::all);
    d->EinA= X13CO(span(j),span(4),span::all);

    d->N_13CO_vj.slice(i).row(j)=(d->tot_col_fluor_nocoll.at(i,j+1)/X12CO_13CO_fl)*(2*d->Jup+1) % exp(-hck*Bv13[j] * d->Jup % (d->Jup+1) / d->T_rot_fl[i]) / (d->T_rot_fl[i]/(hck*Bv13[j]));


 }

  //====================
  //  C18O
  //===================

  d->Jup = XC18O(span(0),span(3),span::all);
  d->Jdn = XC18O(span(0),span(1),span::all);
  d->wvn = XC18O(span(0),span(6),span::all);
  d->EinA= XC18O(span(0),span(4),span::all);

  d->N_C18O_vj.slice(i).row(0)=(d->tot_col_fluor_nocoll.at(i,1)/X12CO_C18O_fl)*(2*d->Jup+1) % exp(-hck*Bv18[0]*d->Jup%(d->Jup+1)/d->T_rot_fl[i]) / (d->T_rot_fl[i]/(hck*Bv18[0]));

}


  d->annuli = zeros<vec>(d->rdisk.n_elem);
  for (int i=0; i<d->steps-1; i++)
  {
    d->annuli.at(i)=datum::pi*(pow(d->rdisk.at(i+1),2)-pow(d->rdisk.at(i),2));
  }

  d->annuli(d->steps-1)=datum::pi*(pow(d->rdisk.at(d->steps-1)+1,2)-pow(d->rdisk.at(d->steps-1),2));

  d->annuli=d->annuli*2.25e26;

  //======================
  //  X12CO
  //=====================
  for (int i=0; i<d->steps; i++)  //Loop over annuli
  {

    double btotcomp=b_tot.at(i)/cexp;
    //==============================
    //  X12CO
    //=============================
    int upbound=tot_col.row(0).n_elem-1;
    for (int j=0; j<upbound; j++)  //vibrational levels
    {
      d->Jup = X12CO(span(j),span(3),span::all);
      d->Jdn = X12CO(span(j),span(1),span::all);
      d->wvn = X12CO(span(j),span(6),span::all);
      d->EinA= X12CO(span(j),span(4),span::all);
      for (int k=0; k<X12CO.n_slices; k++) 
      {
        d->A1=d->wvn.at(k);
	if ( (d->A1 >= f_i) && (d->A1 <= f_f) )
	{

	  d->A0=d->N_12CO_vj.at(j,k,i)*hc*d->A1*d->EinA.at(k);                         //should the first two indices be reversed here?
	  d->A2=btotcomp*d->A1;
	  d->stick_spec_12CO.col(i)+=(d->A0/(rpi*d->A2)) * exp (-pow(((d->A1-freq)/d->A2),2));
	}
      }
    }
    //==============================
    //  X13CO
    //==============================
    
    for (int j=0; j<3; j++)   //Loop over vibrational levels--check for-loop bounds in some earlier integral loops... might need to be adjusted!
    {

      d->Jup = X13CO(span(j),span(3),span::all);
      d->Jdn = X13CO(span(j),span(1),span::all);
      d->wvn = X13CO(span(j),span(6),span::all);
      d->EinA= X13CO(span(j),span(4),span::all);

      for (int k=0; k<X13CO.n_slices; k++)
      {
          d->A1=d->wvn.at(k);
          if ((d->A1 >= f_i) && (d->A1 <= f_f)) 
          {
	    d->A0=d->N_13CO_vj.at(j,k,i)*hc*d->A1*d->EinA.at(k);                         //should the first two indices be reversed here?
	    d->A2=btotcomp*d->A1;
	    d->stick_spec_13CO.col(i)+=(d->A0/(rpi*d->A2)) * exp(-pow(((d->A1-freq)/d->A2),2));

          }
      }    
    }
    //=============================
    //  XC180
    //============================ 
     
      d->Jup = XC18O(span(0),span(3),span::all);
      d->Jdn = XC18O(span(0),span(1),span::all);
      d->wvn = XC18O(span(0),span(6),span::all);
      d->EinA= XC18O(span(0),span(4),span::all);

      for (int k=0; k<XC18O.n_slices; k++)
      {
          d->A1=d->wvn.at(k);
          if ((d->A1 >= f_i) && (d->A1 <= f_f))
          {
	    d->A0=d->N_C18O_vj.at(0,k,i)*hc*d->A1*d->EinA.at(k);                         //should the first two indices be reversed here?
	    d->A2=btotcomp*d->A1;
	    d->stick_spec_C18O.col(i)+=(d->A0/(rpi*d->A2)) * exp(-pow(((d->A1-freq)/d->A2),2));
          }
      }



    d->stick_spec_tot.col(i)=d->stick_spec_12CO.col(i)+d->stick_spec_13CO.col(i)+d->stick_spec_C18O.col(i);  //Note:  I have the notation backwards... since IDL is backwards, (*,i) is col(i).  Fix all of these!
  }
  
  d->annuli.at(0)=1e28;
  mat iten_tot=d->stick_spec_tot;
  iten_tot.col(0)=iten_tot.col(0)*2.5;
  mat Lum_tot = iten_tot;
  for (int i=0; i<d->steps; i++)
  {
    Lum_tot.col(i)=iten_tot.col(i)*d->annuli.at(i);
  }

  mat Flux_tot = Lum_tot/(4*datum::pi*pow(data->stardist,2));


//===================================
//  SCRIPT 4 of 4
//====================================

  double r_in=d->rdisk.at(0);
  double r_out=d->rdisk.at(d->rdisk.n_elem-1);

  double dv=1;

  mat flux_tot_slit=Flux_tot;

//account for slit loss
  vec slit_loss=100.17*pow(d->rdisk,-1.260);
  double disk_Beyond_slit = 63;

  for (int i=0; i<d->rdisk.n_elem; i++)
  {
    if (d->rdisk.at(i) > disk_Beyond_slit)
    {
      flux_tot_slit.col(i)=flux_tot_slit.col(i)*slit_loss.at(i);
    }
  }

  auto n_rings=d->rdisk.n_elem;

  vec dr = d->rdisk;

  for (int i=0; i<n_rings; i++)
  {
    if (d->rdisk(i) < 0.1)    {
      dr(i)=0.01;
    } else 
      {  if (d->rdisk(i) <1.0) {
        dr(i)=0.1;
      } else
        { dr(i) = 1.0;
        }
      }
  }

  vec vmax=round(sqrt(887.2*Mstar/d->rdisk));

  mat total_spec=iten_tot;
  vec vel_spec=total_spec.col(0);

//==========V-LINES==========
/*
  for (int k=0; k<22; k++)
  {
    d->v_line.row(k)=((d->fco_list[k]-freq)*2.9979e5/d->fco_list[k]).t();
  }

  for (int i=0; i<22; i++) d->v_line_indices(i,0) = whererow(d->v_line.row(i), [] (double datum) {return ((datum > -15) && (datum < 15));});
*/
//=========ITEN LINES=======

  for (int i=0; i<22; i++)
  {
    int itencols=iten_tot.n_cols;
    int v_line_num=d->v_line_indices(i,0).n_elem;
    mat temp = zeros<mat>(d->v_line_indices(i,0).n_elem,itencols);
    for (int j=0; j<v_line_num; j++)
    {
      temp.row(j)=iten_tot.row(d->v_line_indices(i,0).at(j));
    }
    d->iten_lines(i,0)=arma::sum(temp,0).t()*5.65e-3;
  }

  int grid_ptr=1;
//==========================================
//  RINGS LOOP                     
//==========================================

  for (int j=0; j<n_rings; j++)
  {
    double n_seg=4*round(vmax.at(j))/dv;
    
    d->grid.resize(d->grid.n_rows+n_seg,26);

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
	area.at(i)=dphase.at(i)*(d->rdisk.at(j)+dr.at(j)/2)*dr.at(j);

        vel_spec=vel_spec+interpol(total_spec.col(j)*area.at(i),freq+vseg.at(i)*sin(inc)*freq/c,freq);
        d->grid.at(grid_ptr,0)=d->rdisk.at(j);
        d->grid.at(grid_ptr,1)=vseg.at(i);
        d->grid.at(grid_ptr,2)=phase.at(i);
        d->grid.at(grid_ptr,3)=area.at(i)*2.25e26;
        for (int ii=4; ii<26; ii++) { d->grid.at(grid_ptr,ii)=d->iten_lines(ii-4,0)(j);}
        grid_ptr++;
      } 
    }
    else
    {
      for (int i=0; i<n_seg; i++)
      {
        area.at(i)=dphase.at(i)*(d->rdisk.at(j)+dr.at(j)/2)*dr.at(j);
        vel_spec=vel_spec+interpol(total_spec.col(j)*area.at(i),freq+vseg.at(i)*sin(inc)*freq/c,freq);
        d->grid.at(grid_ptr,0)=d->rdisk.at(j);
        d->grid.at(grid_ptr,1)=vseg.at(i);
        d->grid.at(grid_ptr,2)=phase.at(i);
        d->grid.at(grid_ptr,3)=area.at(i)*2.25e26;
        for (int ii=4; ii<26; ii++) d->grid.at(grid_ptr,ii)=d->iten_lines(ii-4,0)(j);
        grid_ptr++;
      } 
    }

    total_spec.col(j)=vel_spec;

  }
  ivec index_grid=where(d->grid.col(0), [&] (double datum) {return ((datum <= d->rdisk.at(0)) && (datum > 0));});

  vec  grid_tmp=zeros<vec>(index_grid.n_elem);
  for (int i=0; i<index_grid.n_elem; i++)
  { 
    grid_tmp.at(i)=d->grid(index_grid.at(i),2);
  }
 ivec index_grid2=where(grid_tmp,[] (double datum) {return (datum > arma::datum::pi);});
 for (int i=0; i<index_grid2.n_elem; i++)
 {
   d->grid.at(index_grid.at(index_grid2.at(i)),3)=0;
 }
  d->grid.resize(d->grid.n_rows+11,26);
  for (int i=0; i<11;  i++)
  {
    d->grid.at(grid_ptr,0)=d->r_planet;
    d->grid.at(grid_ptr,1)=d->v_planet+d->gs(i,0);
    d->grid.at(grid_ptr,2)=d->phase_planet;
    d->grid.at(grid_ptr,3)=d->planet_size;
    for (int j=4; j<22; j++)
    {
      d->grid.at(grid_ptr,j)=0;
    }
    d->grid.at(grid_ptr,9)=d->planet_intens*d->gs(i,1);
    grid_ptr++;
  }

  //pow(5.13,-23*2.9979247e10*4*3.1415926535897*pow((103*3.08),2)*(.05/1.16)); // continuum luminosity
  
  vec gridrow=d->grid.col(1);
  double maxloop=2*gridrow.max()+1;

//======================================
//   CENTROID LOOP
//=====================================
  for (int j=0; j<maxloop; j++)
  {
   // d->indexn(j,0).fill(0);
    d->index1=where(d->grid.col(1), [&] (double datum) {return ((datum <= (max(d->grid.col(1)))-j) && ( datum > ( max(d->grid.col(1))-(j+1)) ) );});
   int index1n=d->index1.n_elem;
   if (d->index1.at(0)!=-1)
   {
//=================
//   i0-i21!
//=================

       vec temp2=zeros<vec>(index1n);

       for (int i=0; i<index1n; i++)  {temp2.at(i)=d->grid(d->index1.at(i),1);}

     double mean=arma::mean(temp2);

     for (int k=0; k<22; k++)
     {
       d->indexn(k,0)=whererow(d->v_line.row(k), [&] (double datum) {return (datum > (mean-0.5) && datum <= (mean+0.5));});
     } 

     vec phi2=zeros<vec>(index1n);
     for (int i=0; i<index1n; i++)
     {
       phi2.at(i)=d->grid(d->index1.at(i),2);
     }

      vec rp=zeros<vec>(index1n);
      for (int i=0; i<index1n; i++)
      {
	rp.at(i)=d->grid(d->index1.at(i),0)*sqrt(pow(cos(phi2.at(i)),2)+pow(sin(phi2.at(i)),2)*pow(cos(inc),2));
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
      vec deltay=rp%cos(theta+d->omega);
        for (int z=0; z<22; z++) 
        {
          int siz=d->indexn.at(z).n_elem;
          vec tempv=zeros<vec>(index1n);
          for (int q=0; q<index1n; q++)
	  {
	      tempv.at(q)=d->grid(d->index1.at(q),z+4)*d->grid(d->index1.at(q),3);
	  }
	  for (int k=0;k<siz;k++)
	  {
	    d->centroid.at(d->indexn(z,0).at(k))=(arma::sum(tempv%deltay))/(arma::sum(tempv)+d->Lc);
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
  inst_prof=shift(inst_prof,round(vfreq.n_elem)/2);
  inst_prof=inst_prof/arma::sum(inst_prof);

  cx_vec conv_spec=ifft(fft(final_spec)%fft(inst_prof))/2;
  cx_vec cent_conv=ifft(fft(d->centroid)%fft(inst_prof));

  cent_conv=cent_conv*as_scalar(arma::sum(abs(d->centroid)))/as_scalar(arma::sum(abs(cent_conv)));
  conv_spec=conv_spec*as_scalar(arma::accu(flux_tot_slit))/as_scalar(arma::sum(conv_spec));

  ofstream fout;
  fout.open(folderpath+"/cent_conv");
  fout << real(cent_conv);
  fout.close();
  fout.open(folderpath+"/conv_spec");
  fout << real(conv_spec);
  fout.close();
cerr << "substracting..." << endl;
  ivec indexdiff=where(r1big, [] (double datum) {return (datum != -9999);});
  int indexsiz=indexdiff.n_elem;
  d->diff = ((r1big-1)*1.5e-12 - interpol(real(conv_spec),freq-freq*5/2.9979e5,f1big));

  for (int i=0; i<indexsiz; i++) d->diff.at(indexdiff.at(i)) = 0;
  
  double chisq = arma::sum(abs(d->diff));

  cerr << chisq << endl;
  cin.get();  
  delete(d);

  return 0;
}

int FitData::runTrials() {
  //start with best chi-by-eye fit and run one trial
  //then each time, reset parameters to a value from the matrix and rerun.
  //set parameters for best-chi;


/*  int rank;
  int numtasks;
  int rc;

  if (rc != MPI_SUCCESS) 
  {
    cerr << "Error initializing MPI environment." << endl;
    MPI_Abort(MPI_COMM,WORLD, rc);
  }*/
  

  for(int i=0; i<this->numGuesses-1;i++) {
    double receivedMessage[2];
    for (int j=0; j<numtasks; j++)
    {
    if (i==0)
    { 

      double layers=300;
      double disk_in=13.;
      //double dist=1.496e13*disk_in;
      double disk_out=100.0;
      double v_turb=3e5;
      double T_rot0_fl=2500;             //check this.... is t_rot_fl0?
      double T_rot_alpha_fl=0.25;
      //double T_rot0_cl=T_rot0_fl;
      //double T_rot_alpha_cl=T_rot_alpha_fl;
      double rel_lum=20;
    }
    else 
    {
      cerr << "Aux trial number " << i << " begin now" << endl;
      layers=randData[0][i-1];
      disk_in=randData[1][i-1];
      //dist=1.496e13*disk_in;
      disk_out=randData[2][i-1];
      v_turb=randData[3][i-1];
      T_rot0_fl=randData[4][i-1];
      T_rot_alpha_fl=randData[5][i-1];
      //T_rot0_cl=T_rot0_fl;
      //T_rot_alpha_cl=T_rot_alpha_fl;
      rel_lum=randData[6][i-1];
    }
    isSent(i)=1;
    this->runTrial(layers,disk_in,disk_out,v_turb,T_rot0_fl,T_rot_alpha_fl,rel_lum);
    }
   // RECEIVE MPI HERE

//   MPI_Recv(&receivedMessage,2, MPI_DOUBLE,)
// receive MPI conv_spec cent_conv here if difference is best

  return 0;
}

FitData::FitData(int numGuesses, string folder)
{
  //class variables

  wavenum             = zeros<mat>(10,12);
  einA                = zeros<mat>(10,12);
  lam_ang             = zeros<mat>(10,12);
  HD100546_luminosity = zeros<mat>(10,12);

  FitData::numGuesses=numGuesses;
  folder=path=folder;
  //read in data from files

  FitData::readInput(folderpath+"/input");

  std::ifstream fin;

  fin.open("ratedat/EinA.txt");
  for ( auto& entry : einA) {
     fin >> entry;
  }
  fin.close(); 


  fin.open("ratedat/lambda.txt");
  for ( auto& entry : lam_ang) {
     fin >> entry;
  }
  fin.close();


  fin.open("ratedat/wavenum.txt");
  for ( auto& entry : wavenum) {
     fin >> entry;
  }
  fin.close();
 

  fin.open("ratedat/HD100546_luminosity.txt");
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

  fin.open("HD100546/r1big");
  r1big=zeros<vec>(4096);
  for (int z=0; z<4096; z++)
  {
    fin >> r1big.at(z);
  }
  fin.close();


  fin.open("HD100546/f1big");
  f1big=zeros<vec>(4096);
  for (int z=0; z<4096; z++)
  {
    fin >> f1big.at(z);
  }
  fin.close();

  freq = zeros<vec>(freq_size);
  freq.at(0)=f_i;
  vfreq = zeros<vec>(freq_size);

  for (int i=1; i<freq_size; i++)
  {
    freq.at(i)=f_i*(pow(1+v/(3*c),i));
  }
  vfreq=(freq(round(freq.n_elem/2))-freq)*c/freq(round(freq.n_elem/2));  

  isSent = new double(numGuesses);
  this->randData[0]=FillArray(900, 100);
  this->randData[1]=FillArray(15, 6);
  this->randData[2]=FillArray(76, 50);
  this->randData[3]=FillArray(500000, 100000);      //check to be sure all of these first arguments don't need to be offset by +/- 1
  this->randData[4]=FillArray(3500,1000);
  this->randData[5]=FillArray(50,20);
  this->randData[6]=FillArray(99, 1);
  for (int i=0; i < numGuesses-1; i++) {
    this->randData[5][i]=this->randData[5][i]/100;
    isSent(i)=false;
  };

 /*=================
 *
 * Collision Data                                      EVERYTHING AFTER LAYERS ----> INTO COLLISIONS!  CollisionData.cpp will need to be revisited...
 *
 *==================*/
  fAX = (3.038/2.03)*einA/(wavenum % wavenum);   //be sure this is right!
  fXA = 2*fAX;
  this->runTrials();
}

FitData::~FitData()
{
  for(int i=0; i<7;i++)
  {
    delete[] this->randData[i];
  }
  delete[] this->randData;
}

int createInput() 
{
  

  return 1;
} 

int extractValue(string sin, string varname, &double var)
{
  string split1, split2;
  if (sin[0]=='#') return 0;
  int len = sin.length();
  int sep;

  for (int i=0; i<len; i++)
  {
    sep=0;
    if (sin[i]=='=')
    {
      if (sep != 0)
      {
        cerr << "Invalid input file." << endl;
        abort();
      }
      sep=i;
    }
    split1=sin.substr(0,sep);
    split2=sin.substr(sep-1,len-sep);
    for (int i=0; i<5; i++)
    {
      if (split1==varname) var=atod(split2);
    }
}

int FitData::readInput(string inpFile)
{
  string inputStrings[6] = {"numguesses","inc","mass","dist","Lc"};
  double inputVars[6]= {&numGuesses,&inc,&Mstar,&stardist,&Lc};

  ifstream fin;
  fin.open(inpFile);
  string sin;

  while (getline(fin,sin))
  {
    //exit immediately if the line is commented out
    if (sin[0]=='#') return 0;

    int len = sin.length();
    int sep;
    string split1,split2;

    for (int i=0; i<len; i++)
    {
      sep=0;
      if (sin[i]=='=');
      {
        if (sep != 0)
        {
          cerr << "Invalid input file!" << endl;
          abort();
        }
        sep=i;
      }
    }

    split1=sin.substr(0,sep);
    split2=sin.substr(sep-1,len-sep);
    for (int i=0; i<5; i++)
    {
      if (split1==inputvars[i]) {inputVars[i]=atod(split2);  cerr << inputStrings[i] << ": " << inputVars[i] << endl;}
    }
  }
}

int main(int argc, char* argv[]) 
{
  data = new FitData(atoi(argv[1]), "HD100546");
  delete data;
  return 1;
}
