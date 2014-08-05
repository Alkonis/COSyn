#define _USE_MATH_DEFINES

#include <vector>
#include <iostream>
#include <math.h>
#include <fstream>
#include <armadillo>
#include "lib/idlarma.h"
#include "CollData.cpp"

//namespace stdmath = math;

using namespace std;
using namespace arma;

class FitData{

 private:

  static int numGuesses;

  double* FillArray(int modulus, int offset);


  //parameters constant for all runs

  static constexpr double f_i = 1995.;
  static constexpr double f_f = 2179.;
  static constexpr double mass = 2.4;
  static constexpr double H_den0 =2.5e10;
  static constexpr double H_den_alpha=.15;
  static constexpr double inc = 40 * 3.14159265358979323 / 180;
  static constexpr double X12CO_13CO_fl = 65./30;
  static constexpr double X12CO_C18O_fl = 550./16.25;
  static constexpr double X12CO_13CO_cl = 65.;
  static constexpr double X12CO_C18O_cl = 560.;

  static constexpr double Mstar=2.4;
  static constexpr double stardist=3.1781680e20;
  static constexpr double inst_res=6.0;

  static constexpr double c=2.997924562e5;
  static constexpr double cexp=29979245620;
  static constexpr double hc=  6.626196e-27*2.997924562e10;
  static constexpr double hck=(6.626196e-27*2.997924562e10)/(1.380622e-16);
  static constexpr double cer=8.85282e-13;
  static constexpr double rpi=1.7724538509055159;
  
  static constexpr double vib_einA[10]={34.60,67.68,98.40,126.99,153.59,178.31,201.35,223.10,244.15,265.21};

  static constexpr double v = 2.5;
  static constexpr double kb = 1.380622E-16;
  static constexpr double B = 1.9225;       

//Rotational constants
  static constexpr double B10   = 1.93124;
  static constexpr double B21   = 1.91376;
  static constexpr double B32   = 1.89628;
  static constexpr double B43   = 1.87880;
  static constexpr double B54   = 1.86132;
  static constexpr double B65   = 1.84384;
  static constexpr double B76   = 1.82636;
  static constexpr double B87   = 1.80888;
  static constexpr double B98   = 1.79140;
  static constexpr double B109  = 1.77392;
  static constexpr double B1110 = 1.76644;
  static constexpr double B1310 = 1.847;
  static constexpr double B1321 = 1.830;
  static constexpr double B1332 = 1.813;
  static constexpr double B1810 = 1.840;

  static constexpr double Bv12[11] = {B10,B21,B32,B43,B54,B65,B76,B87,B98,B109,B1110};
  static constexpr double Bv13[3]  = {B1310,B1321,B1332};
  static constexpr double Bv18[1]  = {B1810};

//script 2 parameters
  static constexpr double disk_beyond_slit=63;

  //These are the three arrays read in from the data files
  //float* einA[12];          //change these to doubles?
  //float* lam_ang[12];
  //float* wavenum[12];
  vec T_rot_fl;
  vec T_rot_cl;
  vec H_den;
 
  mat wavenum             = zeros<mat>(10,12);
  mat einA                = zeros<mat>(10,12);
  mat lam_ang             = zeros<mat>(10,12);
  mat HD100546_luminosity = zeros<mat>(10,12);

//Molecular data
  cube X12CO = zeros<cube>(10,7,120);
  cube X13CO = zeros<cube>(3,7,120);
  cube XC18O = zeros<cube>(1,7,120);                         //get this done!

 //Collision data

  mat fAX;
  mat fXA; 

  cube dFdt_0;
  cube tau_0; 
  cube dwdn;
  //cube Nv;
  
  field<cube> rate_eqtn;  //copy this over to colldata on the start of each coll_loop! do this for others that need reiteration as well
  field<cube> g;


//maybe make these fixed

//variable parameters; these are changed with each iteration of the loop
  double layers;
  double disk_in;
  double dist;
  double disk_out;
  double v_turb;
  double T_rot0_fl;
  double T_rot_alpha_fl;
  double T_rot0_cl;
  double T_rot_alpha_cl;
  double rel_lum;

//radial anulli
  double r_index_a=0;
  double r_index_b=0;
  double r_index_c=0;
  double rdisk_index;
  double steps;
  double r;
  double rb[2];
  double maxra;
  double maxrb;

 public:

  double** randData = new double*[7];
  
  std::vector<double> ranulli;
  vec rdisk;    //rdisk may be unnecessary; see if ranulli can replace entirely
  
  FitData(int numGuesses);
  ~FitData();

  double* FillArrays(int modulus, int offset);
  int runCollisions(bool doCols);
  int runTrial();
  int runTrials();

};
