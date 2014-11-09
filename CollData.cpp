#include "COfit.h"

#ifndef COLLDATA_H
#define COLLDATA_H


#ifndef ARMA
#define ARMA_DONT_USE_WRAPPER
#include "include/armadillo"
using namespace arma;
#endif

class CollData {

public:
 int steps;
 float A0;
 float A1;
 float A2;
 double r;
 float sum;

 double dist;

 static constexpr double n=1000;
 static constexpr double a=0;
 static constexpr double b=5;

 static constexpr double delta=0.005; // (b-a)/n
 static constexpr double omega = 55*datum::pi/180;
 static constexpr double Lc=8.3837025e28;  //pow(5.13,-23*2.9979247e10*4*3.1415926535897*pow((103*3.08),2)*(.05/1.16)); // continuum luminosity
 
 //planet parameters
 static constexpr double v_planet=5;
 static constexpr double r_planet=13;
 static constexpr double phase_planet=53*datum::pi/180;
 static constexpr double planet_intens=((2*6.62e-27*pow(2.9979e10,2)*pow(2030,3)/(exp(6.626e-27*2.9979e10*2030/(1.38e-16*1e3))-1)));
 static constexpr double planet_size=0;

 static constexpr int num_lines=22;
 static constexpr double fco_list[] = {2032.3528,2034.4083,2033.4174,2033.1423,2032.8258,2030.1586,2034.9211, 2034.7209, 2034.1352,2034.0473,2032.8701, 2032.8096, 2032.5616, 2032.2011, 2030.5160, 2030.3076, 2029.6559, 2029.4427, 2029.4297, 2029.3679, 2029.2362, 2029.1276};

 std::vector<double> ranulli;

 vec rdisk;
 vec phi;
 vec m_disk;
 vec m_uv;
 vec T_rot_fl;
 vec T_rot_cl;
 vec H_den;
 vec annuli;
 vec tau; 
 vec F_tau;
 vec centroid;
 vec diff;

 rowvec Jup;
 rowvec Jdn;
 rowvec wvn;
 rowvec EinA;

 ivec index1;

 fmat tot_col_fluor;
 fmat tot_col_fluor_nocoll;
 mat stick_spec_tot;
 mat grid;
 mat gs;
 mat v_line;
 mat lum;

 fcube dFdt_0;
 fcube tau_0;
 fcube dwdn;
 fcube Nv;
 cube N_12CO_vj;
 cube N_13CO_vj;
 cube N_C18O_vj;
 
 field <cube> rate_eqtn;
 field <cube> rate_eqtn2;
 field <fcube> g;
 field <ivec> v_line_indices;
 field<ivec>indexn;
 field <vec> iten_lines;

  CollData(double layers, double disk_in, double disk_out, double v_turb, double T_rot0_fl, double T_rot_alpha_fl, double rel_lum)  {
    double dist=1.496e13*disk_in;

    dFdt_0 = zeros<fcube>(10,12,layers);
    tau_0  = zeros<fcube>(10,12,layers);
    dwdn   = zeros<fcube>(10,12,layers);

    tau = linspace<vec>(0,19.99,2000); 
    F_tau = zeros<vec>(2000);
   
    g=field<fcube>(10,1);       
    v_line=zeros<mat>(num_lines,data->freq_size);
    v_line_indices=field<ivec>(22,1);
    indexn = field<ivec>(22,1);
    iten_lines=field<vec>(22,1);
    diff = zeros<vec>(4096);   
    gs = zeros<mat>(11,2);
    grid=zeros<mat>(1,26);   

  //=================================
  //
  //    DIVIDE INTO ANNULI
  //    CALCULATE STEPS
  //
  //================================                               

    r=disk_in;  

    double r_index_a = 0;
    double r_index_b = 0;
    double r_index_c = 0;

    lum = data->HD100546_luminosity*rel_lum;
    ranulli.push_back(disk_in);
    if (r<0.1)
    {
      while ((ranulli.back() < 1) && (ranulli.back() < disk_out))
      {
	r_index_a=r_index_a+1;
	ranulli.push_back(disk_in+0.01*r_index_a);
      }
    }
    double maxra = ranulli.back();

    if ((maxra < 1) && (maxra >= .1))
    {
      while ((ranulli.back() < 1.0) || (ranulli.back() < disk_out))
      {
       r_index_b++;
       ranulli.push_back(maxra+0.1*r_index_b);
      }
    }

    double maxrb=ranulli.back();
    if ((maxrb <= disk_out) && (maxrb >= 1.0))
    {
      while  (ranulli.back() < disk_out)
      {
	r_index_c++;
	ranulli.push_back(maxrb+1.0*r_index_c);
      }
    }
    steps=ranulli.size();  //check that this size function is correct
    rdisk=zeros<vec>(steps);

    for (int i=0; i < steps; i++ )
    {
      rdisk.at(i)=ranulli.at(i)*1.496E13;
    }

    Nv     = zeros<fcube>(21,layers,steps);

    rate_eqtn = field<cube>(steps,1);       
    rate_eqtn2 = field<cube>(steps,1);

    for (int i=0; i<10; i++)                                                                                                             
    {                                                                                                                                    
      g.at(i,0)=zeros<fcube>(12,layers,steps);                                                                                            
    }         
    for (int i=0; i<steps; i++)                                                                           
    {                                                                                                                                    
      rate_eqtn.at(i,0)=zeros<cube>(21,21,layers);                                                                                    
      rate_eqtn.at(i,0)(span::all,span(20),span::all).fill(1);                                                                           
    }                                                                                                                                    
    for (auto i=0; i<steps;i++) {
      for (auto j=0; j<layers;j++) {
	for (auto k=0; k<10; k++) {
	  for (auto q=0; q<11; q++)
	  {
	    rate_eqtn(i,0)(k+11,q,j)=data->einA.at(k,q);//k+11, q,j,i,   k ,q  //why no zero here
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
	  rate_eqtn.at(j,0).at(z+1,z+1,q)=-data->vib_einA[z];  
	  rate_eqtn.at(j,0).at(z+2,z+1,q)=data->vib_einA[z+1]; 
	}
      }
     
      rate_eqtn.at(j,0)(span(1),span(0),span::all).fill(data->vib_einA[0]);   
      rate_eqtn.at(j,0)(span(10),span(10),span::all).fill(-data->vib_einA[9]); 

      for (auto i=0; i<9; i++)
      {
	rate_eqtn.at(j,0)(span(i+11),span(i+11),span::all).fill(-arma::sum(data->einA.row(i)));
      }
    }

    for (auto j=0; j<steps; j++)
    {
      for (int i=0; i<8; i++)
      {
	rate_eqtn.at(j,0).subcube(span(i+1),span(i+1),span::all).fill(-data->vib_einA[i]);
	rate_eqtn.at(j,0).subcube(span(i+2),span(i+1),span::all).fill(data->vib_einA[i+1]);
      } 

      rate_eqtn.at(j,0).subcube(span(1),span(0),span::all).fill(data->vib_einA[0]);
      rate_eqtn.at(j,0).subcube(span(9),span(9),span::all).fill(-data->vib_einA[8]);
    }

    rate_eqtn2=rate_eqtn;

    for (int i=0; i<11; i++)
    {
      gs.at(i,0)=i-5; //<<-5 << -4 << -3 << -2<<-1<<0<<1<<2<<3<<4<<5<<endr;
    }
    gs.col(1)=exp(-pow(gs.col(0),2)/pow((12/1.665),2))/6.1967;
    stick_spec_tot = zeros<mat>(data->freq_size,steps);
    centroid = zeros<vec>(data->freq_size);
    for (int k=0; k<22; k++)
    {
      v_line.row(k)=((fco_list[k]-data->freq)*2.9979e5/fco_list[k]).t();
    }
    for (int i=0; i<22; i++) v_line_indices(i,0) = whererow(v_line.row(i), [] (double datum) {return ((datum > -15) && (datum < 15));});
  }
};

#endif
