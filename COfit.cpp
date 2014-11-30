#include "COfit.h"
#include "CollData.cpp"
#include <iostream>
#include <cstdlib>
#include "time.h"

using namespace std;

#include "mpi.h"

constexpr double CollData::fco_list[];

//int FitData::numGuesses;

//returns an array of size numGuesses filled with random numbers on the interval (offset, modulus+offset)
//used to generate the matrix of guess data
double* FitData::FillArray(int min, int max)
{
  double* array;
  array  = new double[numGuesses];
  for(int i=0; i<numGuesses;i++) {
    array[i]=rand() % (max-min) + min;
  }
  return array;
}

int FitData::dbm(string message)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if (rank==0 && doDebug==1)  cout << message << endl;
}

int FitData::runTrial(double layers, double disk_in, double disk_out, double v_turb, double T_rot0_fl, double T_rot_alpha_fl, double rel_lum, int locali) {
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
         //                  (BIG INTERIOR LOOP!)
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

	   d->Nv.slice(k).col(j) = conv_to<fvec>::from(sol);  

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
        d->tot_col_fluor = totalDimf(d->Nv*7.55858e12,2).t();
        d->rate_eqtn=d->rate_eqtn2;
        d->rate_eqtn2.reset();
      }
      if (coll_loop==1) {
        d->tot_col_fluor_nocoll = totalDimf(d->Nv*7.55858e12,2).t();
      }
    }

    d->Nv.reset();
    d->rate_eqtn.reset();
    d->g.reset();
    d->tau_0.reset();
    d->dwdn.reset();
    d->dFdt_0.reset();

//=========================================================================
// Angle of Incidence Correction (tweak for each star--use input file!)
//========================================================================
  d->m_disk=1.5*(5.59647e-10)*sqrt(d->T_rot_fl / Mstar)%pow(d->rdisk,0.5);
  d->m_uv = (5.59647e-10)*sqrt(d->T_rot_fl/Mstar)%sqrt(d->rdisk);

  d->phi = -atan((d->m_uv-d->m_disk)/(1+d->m_uv%d->m_disk));
  d->phi.at(0)=datum::pi/2;      
  
  auto xi = d->tot_col_fluor.n_cols-1;
  auto yi = d->tot_col_fluor.n_rows-1;

  for (int j=1; j<yi+1; j++)    
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
	  d->stick_spec_tot.col(i)+=(d->A0/(rpi*d->A2)) * exp (-pow(((d->A1-freq)/d->A2),2));
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
	    d->stick_spec_tot.col(i)+=(d->A0/(rpi*d->A2)) * exp(-pow(((d->A1-freq)/d->A2),2));

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
	    d->stick_spec_tot.col(i)+=(d->A0/(rpi*d->A2)) * exp(-pow(((d->A1-freq)/d->A2),2));
          }
      }
  }
 
 //Calculates total luminosity, then total flux through slit
  d->annuli.at(0)=1e28;
  d->stick_spec_tot.col(0)=d->stick_spec_tot.col(0)*2.5;
  mat flux_tot_slit = d->stick_spec_tot;
  for (int i=0; i<d->steps; i++)
  {
    flux_tot_slit.col(i)=flux_tot_slit.col(i)*d->annuli.at(i);
  }

  flux_tot_slit=flux_tot_slit/(4*datum::pi*pow(data->stardist,2));
//===================================
//  SCRIPT 4 of 4
//====================================

  double r_in=d->rdisk.at(0);
  double r_out=d->rdisk.at(d->rdisk.n_elem-1);

  double dv=1;

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

  vec vel_spec=d->stick_spec_tot.col(0);

//=========ITEN LINES=======

  for (int i=0; i<22; i++)
  {
    int itencols=d->stick_spec_tot.n_cols;
    int v_line_num=d->v_line_indices(i,0).n_elem;
    mat temp = zeros<mat>(d->v_line_indices(i,0).n_elem,itencols);
    for (int j=0; j<v_line_num; j++)
    {
      temp.row(j)=d->stick_spec_tot.row(d->v_line_indices(i,0).at(j));
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

        vel_spec=vel_spec+interpol(d->stick_spec_tot.col(j)*area.at(i),freq+vseg.at(i)*sin(inc)*freq/c,freq);
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
        vel_spec=vel_spec+interpol(d->stick_spec_tot.col(j)*area.at(i),freq+vseg.at(i)*sin(inc)*freq/c,freq);
        d->grid.at(grid_ptr,0)=d->rdisk.at(j);
        d->grid.at(grid_ptr,1)=vseg.at(i);
        d->grid.at(grid_ptr,2)=phase.at(i);
        d->grid.at(grid_ptr,3)=area.at(i)*2.25e26;
        for (int ii=4; ii<26; ii++) d->grid.at(grid_ptr,ii)=d->iten_lines(ii-4,0)(j);
        grid_ptr++;
      } 
    }

    d->stick_spec_tot.col(j)=vel_spec;

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
	    d->centroid.at(d->indexn(z,0).at(k))=(arma::sum(tempv%deltay))/(arma::sum(tempv)+Lc);
	  }
        }

    }
  }
  for (int j=0; j<n_rings; j++)
  {
    d->stick_spec_tot.col(j)=d->stick_spec_tot.col(j)*arma::sum(flux_tot_slit.col(j))/arma::sum(d->stick_spec_tot.col(j));
  }
  vec final_spec=arma::sum(d->stick_spec_tot,1);

  vec inst_prof=final_spec;
  inst_prof.fill(0);
  inst_prof=exp(-pow(vfreq,2)/pow(inst_res/1.665,2));
  inst_prof=shift(inst_prof,round(vfreq.n_elem)/2);
  inst_prof=inst_prof/arma::sum(inst_prof);

  cx_vec conv_spec=ifft(fft(final_spec)%fft(inst_prof))/2;
  cx_vec cent_conv=ifft(fft(d->centroid)%fft(inst_prof));

  cent_conv=cent_conv*as_scalar(arma::sum(abs(d->centroid)))/as_scalar(arma::sum(abs(cent_conv)));
  conv_spec=conv_spec*as_scalar(arma::accu(flux_tot_slit))/as_scalar(arma::sum(conv_spec));
   
  if (round(locali) == 0)
  {
    ofstream fout;
    fout.open(folderpath+"/cent_conv_0");
    fout << real(cent_conv);
    fout.close();
    fout.open(folderpath+"/conv_spec_0");
    fout << real(conv_spec);
    fout.close();

  }
  
  ivec indexdiff=where(r1big, [] (double datum) {return (datum == -9999);});
  int indexsiz=indexdiff.n_elem;
  d->diff = (r1big - interpol(real(conv_spec),freq-freq*5/2.9979e5,f1big));
  for (int i=0; i<indexsiz; i++) d->diff.at(indexdiff.at(i)) = 0;
  double chisq = ( (arma::sum(pow(abs(d->diff.subvec(0,1023)),2)/pow(segment1,2)) + arma::sum(pow(abs(d->diff.subvec(1024,2047)),2)/pow(segment2,2)) + arma::sum(pow(abs(d->diff.subvec(2048,3071)),2)/pow(segment3,2)) + arma::sum(pow(abs(d->diff.subvec(3072,4095)),2)/pow(segment4,2) ))/4)/(4096-indexsiz);

  //===========================
  //===  MPI Communication  ===
  //===========================

  int rank,numtasks;
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  double sendMessage[2];

  sendMessage[0]=locali;
  sendMessage[1]=chisq;

  if (rank!=0)
  {
    MPI_Send(&sendMessage,2,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
  }
 
  if (rank==0)
  {
    local_i=locali;
    local_chisq=chisq;
  } 

  //==========================  

  delete(d);
  return 0;
}

int FitData::runTrials() 
{
  
  int rank;
  int numtasks;
  int rc;
  int I=0;
  int quit=0;
  int locali;

  MPI_Init(NULL,NULL);
  if (rc != MPI_SUCCESS) 
  {
    cerr << "Error initializing MPI environment." << endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  MPI_Status Stat;
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);




  //devote procesor #1 as a scheduler
  if (rank==0) 
  {


 
    ifstream fin(folderpath+"/input");
    vector<string> inp;
    string line;
    while (getline(fin,line)) inp.push_back(line);

    ofstream fout(folderpath+"/output"+to_string(static_cast<int>(fileCount)));
    fout << "======================================================================" << endl;
    fout << "======================= COSYN GENERATED OUTPUT =======================" << endl;
    fout << "======================================================================" << endl;
    fout << endl;
    fout << "COSyn " << folderpath << ": " << numGuesses << " trials" << endl; 
    fout << endl;
    fout << "Input echo: " << endl;
    fout << endl;
    for (int i=0; i<inp.size(); i++) fout << "       " <<  inp.at(i) << endl;
    fout << endl;
    fout << endl;
    fout << "======================================================================" << endl;;
    fout << "i chisq layers disk_in disk_out v_turb T_rot0_fl T_rot_alpha_fl rel_lum" << endl;;

    /*//////////////////////////////////////////////////////////////////////
    //  Two step process:
    //  
    //  (1) First, send out trial info to each processor SYNCHRONOUSLY;
    //      check to see if number of processes is greater than the number
    //      of trials, and terminate once enough communications are made.
    //
    //  (2) Second, begin asynchronous transmission.  Begin looking for
    //     incoming messages, and respond with another set of data to
    //      model.  Repeat this until the maximum number of trials
    //      have been run.
    //
    *///////////////////////////////////////////////////////////////////////
    double recvMsg[2];
    int sendMsg[2];

    finchivec=zeros<vec>(numGuesses);

    if (I >= numGuesses) quit=1;
    if (numGuesses+1<numtasks) 
    {
      numtasks=numGuesses;
      cerr << "WARNING:  number of proceses > number of guesses.  Downsizing." << endl;
      quit=1;
    }
    sendMsg[1]=quit;
    for (int i=1; i<numtasks; i++)
    {
      sendMsg[0]=I;
      MPI_Send(&sendMsg,2,MPI_INT,i,0,MPI_COMM_WORLD);
      I++;
    }

    

    double tag;
    double srcRank;
    double recvCount=0;

    /*///////////////////////////////////////////////////////////////////
    //
    //                      ASYNCHRONOUS LOOP
    //
    //  Important to keep in mind there are two different breakpoints;
    //  the point at which the asynch loop is finished sending data
    //  and the point at which it is done receiving.  When it finishes
    //  sending, it should send a quit=1 flag to all other processes,
    //  but it should continue receiving until it has received data from
    //  numGuesses proceses.
    //
    *////////////////////////////////////////////////////////////////////

    int index;
    while (recvCount<numGuesses)
    {
      MPI_Recv(&recvMsg,2,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&Stat);
      srcRank=Stat.MPI_TAG;

      recvCount++;

      finchivec.at(recvMsg[0])=recvMsg[1];

      cerr << "****    TRIAL " << recvMsg[0] << " COMPLETE    **** chisq " <<endl;

    index = static_cast<int>(recvMsg[0]);
    if (recvMsg[0]==0) 
    {
     fout << 0 << " " << finchivec.at(0) << " " << layers_0 << " " << disk_in_0 << " " << disk_out_0 << " " << v_turb_0 << " " << T_rot0_fl_0 << " " << T_rot_alpha_fl_0 << " " << rel_lum_0 << endl;
    } 
    else {
        fout << index << " "  << finchivec.at(index) << " " << randData[0][index] << " "  << randData[1][index] << " "  << randData[2][index] << " " << randData[3][index] << " " << randData[4][index] << " " << randData[5][index] << " " << randData [6][index] << endl;
    }
      if (I >=numGuesses) quit=1;
      sendMsg[0]=I;
      sendMsg[1]=quit;
      MPI_Send(&sendMsg,2,MPI_INT,srcRank,0,MPI_COMM_WORLD);
      
      I++;
    }
    fout.close();
  }
  

  // END MASTER PROCESS WRAPPER


  //  SLAVE PROCESS WRAPPER
  //  fairly simply, all other processes used as slaves

  else
  { 

    double layers;
    double disk_in;
    //double dist=1.496e13*disk_in;
    double disk_out;
    double v_turb;
    double T_rot0_fl;             //check this.... is t_rot_fl0?
    double T_rot_alpha_fl;
    //double T_rot0_cl=T_rot0_fl;
    //double T_rot_alpha_cl=T_rot_alpha_fl;
    double rel_lum;

    int recvMsg[2];
    
    while(42)
    {
      MPI_Recv(&recvMsg,2,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&Stat);
      locali=recvMsg[0];
      quit = recvMsg[1]; 
      if (quit) break;

      /////////////////////////////////////////////////////////////
      //
      //                         SELECT DATA
      //
      //  If this is the first job being sent out, then use the best
      //  chi-by-eye data given in the input file.  Otherwise, use
      //  a randomly-generated data point from the span read in from
      //  the input file.
      //
      //////////////////////////////////////////////////////////////
      if (locali==0)
      { 
	layers=layers_0;
	disk_in=disk_in_0;
	//double dist=1.496e13*disk_in;
	disk_out=disk_out_0;
	v_turb=v_turb_0;
	T_rot0_fl=T_rot0_fl_0;             //check this.... is t_rot_fl0?
	T_rot_alpha_fl=T_rot_alpha_fl_0;
	//double T_rot0_cl=T_rot0_fl;
	//double T_rot_alpha_cl=T_rot_alpha_fl;
	rel_lum=rel_lum_0;
      }
      else 
      {
	layers=randData[0][locali];
	disk_in=randData[1][locali];
	//dist=1.496e13*disk_in;
	disk_out=randData[2][locali];
	v_turb=randData[3][locali];
	T_rot0_fl=randData[4][locali];
	T_rot_alpha_fl=randData[5][locali];
	//T_rot0_cl=T_rot0_fl;
	//T_rot_alpha_cl=T_rot_alpha_fl;
	rel_lum=randData[6][locali];
      }

      //Run a trial with this data; MPI_Send will be called within this trial to return data to the master proces

      cout << "runTrial " << layers << " " << disk_in << " " << disk_out << " " << v_turb << " " << T_rot0_fl << " " << T_rot_alpha_fl << " " << rel_lum << " " << locali << endl;
      this->runTrial(layers,disk_in,disk_out,v_turb,T_rot0_fl,T_rot_alpha_fl,rel_lum,locali);
        //this->runTrial(layers_0,disk_in_0,disk_out_0,v_turb_0,T_rot0_fl_0,T_rot_alpha_fl_0,rel_lum_0,locali);
    }
// receive MPI conv_spec cent_conv here if difference is best
  }
  if (rank==0)
  {
    uword mindex;
    double min;
    min = finchivec.min(mindex);
    ofstream fout(folderpath+"/output"+to_string(static_cast<int>(fileCount)), fstream::app|fstream::out);
    fout << endl;
    fout << "Minimum:  i=" << mindex << ", chi^2=" << min << endl;
    fout.close();
  }

  MPI_Finalize();
  return 0;


}

int FitData::runTrialsSlave() {
  //start with best chi-by-eye fit and run one trial
  //then each time, reset parameters to a value from the matrix and rerun.
  //set parameters for best-chi;
  int rank;
  int numtasks;
  int rc;

  MPI_Init(NULL,NULL);
  if (rc != MPI_SUCCESS) 
  {
    cerr << "Error initializing MPI environment." << endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  MPI_Status Stat;
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  double layers;
  double disk_in;
  //double dist=1.496e13*disk_in;
  double disk_out;
  double v_turb;
  double T_rot0_fl;             //check this.... is t_rot_fl0?
  double T_rot_alpha_fl;
  //double T_rot0_cl=T_rot0_fl;
  //double T_rot_alpha_cl=T_rot_alpha_fl;
  double rel_lum;

  int locali;

  if (rank==0) finchivec=zeros<vec>(numGuesses);
    for (int i=0; i<this->numGuesses;) {
    double receivedMessage[2];

    int sent=0;
    int exitWhile=0;
    int quit=0;
    int* useProcess = new int[numtasks];
    int useThis=1;
    for (int k=0; k<numtasks; k++)  useProcess[k]=1;
    if (rank==0)
    { 

      int* useProcess = new int[numtasks];

      while (sent<numtasks)
      {
        int doLoop=1;

        if (i>=numGuesses) 
        {
          cerr << "Terminated!  Maximum number of trials reached!  Process " << sent << endl;
          for (int k=sent; k<numtasks; k++) 
          {
            useProcess[k]=0;
            MPI_Send(&useProcess[k],1,MPI_INT,k,useProcess[k], MPI_COMM_WORLD);
          }
          break;
        }
        else
        {
          useProcess[sent]=1;
        }
       
        //for each core, find the first i in randData that has not been assigned a processor and send that i to that process
        for (int j=0;doLoop==1;j++)
        {
          if (!isSent[j])
          {

            isSent[j]=1;
            i++;
            if (sent==0)
            {
              locali=j;
              doLoop=0;
            }
            else
            {
              MPI_Send(&j,1,MPI_INT,sent,useProcess[sent], MPI_COMM_WORLD);
              doLoop=0;
            }

            sent++;
          }

          if (j>numGuesses) {cerr << "ERROR in guess communication! i = " << i << endl;}
        }
      if (exitWhile) break;
      }
    }

    if (rank==0)  useThis=1;     

    if (rank!=0)
    {
      MPI_Recv(&locali,1,MPI_INT,0,MPI_ANY_TAG, MPI_COMM_WORLD, &Stat);
      useThis=Stat.MPI_TAG;
    }


    if (locali==0)
    { 

      layers=layers_0;
      disk_in=disk_in_0;
      //double dist=1.496e13*disk_in;
      disk_out=disk_out_0;
      v_turb=v_turb_0;
      T_rot0_fl=T_rot0_fl_0;             //check this.... is t_rot_fl0?
      T_rot_alpha_fl=T_rot_alpha_fl_0;
      //double T_rot0_cl=T_rot0_fl;
      //double T_rot_alpha_cl=T_rot_alpha_fl;
      rel_lum=rel_lum_0;
    }
    else 
    {
      layers=randData[0][locali-1];
      disk_in=randData[1][locali-1];
      //dist=1.496e13*disk_in;
      disk_out=randData[2][locali-1];
      v_turb=randData[3][locali-1];
      T_rot0_fl=randData[4][locali-1];
      T_rot_alpha_fl=randData[5][locali-1];
      //T_rot0_cl=T_rot0_fl;
      //T_rot_alpha_cl=T_rot_alpha_fl;
      rel_lum=randData[6][locali-1];
    }
    if (useThis) 
    {
      this->runTrial(layers,disk_in,disk_out,v_turb,T_rot0_fl,T_rot_alpha_fl,rel_lum,locali);
    }
    // RECEIVE MPI HERE
    if (rank==0) 
    {
      finchivec.at(local_i)=local_chisq;
      cerr << "****    TRIAL " << locali << " COMPLETE    ****" << endl;
      for (int j=1; j<sent; j++)
      {
        MPI_Recv(&receivedMessage,2, MPI_DOUBLE,j,2,MPI_COMM_WORLD, &Stat);
        finchivec.at(receivedMessage[0])=receivedMessage[1];


        cerr << "****    TRIAL " << receivedMessage[0] << " COMPLETE    ****" << endl;
      }
    }
 
      if (rank==0)
      { 
        if (i>=numGuesses)
        {
          quit=1; 
          for (int j=1; j<numtasks; j++) MPI_Send(&quit,1,MPI_INT,j,3,MPI_COMM_WORLD);
        } else 
        {
          quit=0;
          for (int j=1; j<numtasks; j++) MPI_Send(&quit,1,MPI_INT,j,3,MPI_COMM_WORLD);
        }
      }
      else
      {
        MPI_Recv(&quit,1,MPI_INT,0,3,MPI_COMM_WORLD,&Stat);
      }

    if (quit) break;
// receive MPI conv_spec cent_conv here if difference is best
  }

  if (rank==0)
  {

    ifstream fin(folderpath+"/input");
    vector<string> inp;
    string line;
    while (getline(fin,line)) inp.push_back(line);

    ofstream fout(folderpath+"/output"+to_string(fileCount));
    fout << "======================================================================" << endl;
    fout << "======================= COSYN GENERATED OUTPUT =======================" << endl;
    fout << "======================================================================" << endl;
    fout << endl;
    fout << "COSyn " << folderpath << ": " << numGuesses << " trials" << endl; 
    fout << endl;
    fout << "Input echo: " << endl;
    for (int i=0; i<inp.size(); i++) fout << inp.at(i) << endl;
    fout << endl;
    fout << endl;
    fout << "======================================================================" << endl;;
    fout << "i chisq layers disk_in disk_out v_turb T_rot0_fl T_rot_alpha_fl rel_lum" << endl;;
    for (int i=0; i<numGuesses; i++)
    {
      fout << i << " "  << finchivec.at(i) << " " << randData[0][i] << " "  << randData[1][i] << " "  << randData[2][i] << " " << randData[3][i] << " " << randData[4][i] << " " << randData[5][i] << " " << randData [6][i] << endl;
    }
    fout.close();

    uword mindex;
    double min;
    min = finchivec.min(mindex);
  }

  MPI_Finalize();
  return 0;

}

FitData::FitData(string folder)
{
  //class variables

  wavenum             = zeros<mat>(10,12);
  einA                = zeros<mat>(10,12);
  lam_ang             = zeros<mat>(10,12);
  HD100546_luminosity = zeros<mat>(10,12);

  //FitData::numGuesses=numGuesses;
  folderpath=folder;
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

  fin.open(folderpath+"/r1big");
  r1big=zeros<vec>(4096);
  for (int z=0; z<4096; z++)
  {
    fin >> r1big.at(z);
  }
  fin.close();


  fin.open(folderpath+"/f1big");
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

  isSent = new bool[numGuesses];

  this->randData[0]=FillArray(layers_min, layers_max);
  this->randData[1]=FillArray(disk_in_min, disk_in_max);
  this->randData[2]=FillArray(disk_out_min,disk_out_max);
  this->randData[3]=FillArray(v_turb_min, v_turb_max);     
  this->randData[4]=FillArray(T_rot0_fl_min,T_rot0_fl_max);
  this->randData[5]=FillArray(T_rot_alpha_fl_min,T_rot_alpha_fl_max);
  this->randData[6]=FillArray(rel_lum_min, rel_lum_max);
  for (int i=0; i < numGuesses; i++) {
    this->randData[5][i]=this->randData[5][i]/1000;
    isSent[i]=0;
  }
  isSent[numGuesses-1]=0;

 /*=================
 *
 * Collision Data                                      EVERYTHING AFTER LAYERS ----> INTO COLLISIONS!  CollisionData.cpp will need to be revisited...
 *
 *==================*/
  fXA = 2*(3.038/2.03)*einA/(wavenum % wavenum);   //be sure this is right!
  this->runTrials();
}

FitData::~FitData()
{
  for(int i=0; i<7;i++)
  {
    delete[] this->randData[i];
  }
  delete[] this->randData;
  delete[] this->isSent;
}

int createInput() 
{
  return 1;
} 

int FitData::extractValue(string sin, string varname, double& var)
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
    for (int i=0; i<6; i++)
    {
        if (split1==varname) var=stod(split2);
    }
    if (split1=="numGuesses") {this->numGuesses=atoi(split2.c_str());}
  }
}

int FitData::readInput(string inpFile)
{
  ifstream fin;
  fin.open(inpFile);
  string sin;

  while (getline(fin,sin))
  {
    //exit immediately if the line is commented out or is an empty string
    if (sin[0]=='#' || sin.empty()) continue;

    int len = sin.length();
    int sep=-1;
    string split1,split2;
    for (int i=0; i<len; i++)
    {
      if (sin[i] =='=')
      {
        if (sep != -1)
        {
          cerr << "Invalid input file!" << endl;
          abort();
        }
        sep=i;
      }
    }

    if (sep==0)
    {
      cerr << "Invalid input file!  Please preface comments with '#'." << endl;
      abort();
    }
    split1=sin.substr(0,sep);
    split2=sin.substr(sep+1,len-sep);
    for (int i=0; i<inputs; i++)
    {
      if (split1==inputStrings[i]) 
      {
        *inputVars[i]=stod(split2);
        if (split1=="T_rot_alpha_fl_min" || split1== "T_rot_alpha_fl_max") *inputVars[i]=*inputVars[i]*1000;
      }
    }
    if (split1=="numGuesses") this->numGuesses=atoi(split2.c_str());
  }

}

int main(int argc, char* argv[]) 
{
  data = new FitData("HD100546");
  delete data;
  return 0;
}
