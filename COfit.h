#define _USE_MATH_DEFINES

#ifndef COFIT_H
#define COFIT_H

#include <vector>
#include <iostream>
#include <math.h>
#include <fstream>

#ifndef ARMA
#define ARMA_DONT_USE_WRAPPER
#define ARMA
#define ARMA_NO_DEBUG
//#include "include/OpenBLAS/cblas.h"


#include "include/armadillo"
using namespace arma;
#endif



#ifndef IDLARMA
#define IDLARMA

#include "lib/idlarma.h"
using namespace idlarma;

#endif

//namespace stdmath = math;

using namespace std;

class FitData{

 private:

  //random data generation function
  double* FillArray(int modulus, int offset);
  double* BlankArray(double val);


  //Used by rank 0 MPI
  int local_i;
  int minChi;
  double local_chisq;
  double doMPI;

  //final vector containing chi for each i in numguesses 
  vec finchivec;

  //=================================
  //====  input file parameters  ====
  //=================================
  
  //star parameters
  static double Mstar;
  static double stardist;
  static double inst_res;
  static double f_i;
  static double f_f;
  static unsigned int numGuesses;
  static double Lc;  //pow(5.13,-23*2.9979247e10*4*3.1415926535897*pow((103*3.08),2)*(.05/1.16)); // continuum luminosity

  static double order1;
  static double order2;
  static double order3;
  static double order4;

  static double doDebug;

  //chi-by-eye parameters
  static double disk_in_0;
  static double disk_out_0;
  static double v_turb_0;
  static double T_rot0_fl_0;
  static double T_rot_alpha_fl_0;
  static double rel_lum_0;
  static double layers_0;
  static double inclination_0;
  static double H_den0_0;
  static double H_den_alpha_0;
  static double X12CO_13CO_fl_0;
  static double X12CO_C18O_fl_0;
  static double X12CO_13CO_cl_0;
  static double X12CO_C18O_cl_0;        

  //range parameters
  static double disk_in_min;
  static double disk_in_max;
  static double disk_out_min;
  static double disk_out_max;
  static double v_turb_min;
  static double v_turb_max;
  static double T_rot0_fl_min;
  static double T_rot0_fl_max;
  static double T_rot_alpha_fl_min;
  static double T_rot_alpha_fl_max;
  static double rel_lum_min;
  static double rel_lum_max;
  static double layers_min;
  static double layers_max;
   
  static double inclination_min;
  static double inclination_max;
  static double H_den0_min;
  static double H_den0_max;
  static double H_den_alpha_min;
  static double H_den_alpha_max;

  static double X12CO_13CO_fl_min;
  static double X12CO_13CO_fl_max;
  static double X12CO_C18O_fl_min;
  static double X12CO_C18O_fl_max;
  static double X12CO_13CO_cl_min;
  static double X12CO_13CO_cl_max;
  static double X12CO_C18O_cl_min;
  static double X12CO_C18O_cl_max;

  static double fileCount;

  static double orders;
  static double orderpix;
  
  static double model_layers;
  static double model_disk_in;
  static double model_disk_out;
  static double model_v_turb;
  static double model_T_rot0_fl;
  static double model_T_rot_alpha_fl;
  static double model_rel_lum;
  static double model_inclination;
  static double model_H_den0;
  static double model_H_den_alpha;
  static double model_X12CO_13CO_fl;
  static double model_X12CO_C18O_fl;
  static double model_X12CO_13CO_cl;
  static double model_X12CO_C18O_cl;

  static double doModel[14];

  //setup variables for reading input

  static constexpr unsigned int inputs = 71;

  static constexpr double* inputVars[inputs] = 
  {
    &inclination_0,
    &Mstar,
    &stardist,
    &Lc,
    &inst_res,
    &doDebug,
    &disk_in_0,
    &disk_out_0,
    &v_turb_0,
    &T_rot0_fl_0,
    &T_rot_alpha_fl_0,
    &rel_lum_0,
    &layers_0,
    &inclination_0,
    &H_den0_0,
    &H_den_alpha_0,
    &X12CO_13CO_fl_0,
    &X12CO_C18O_fl_0,
    &X12CO_13CO_cl_0,
    &X12CO_C18O_cl_0,

    &disk_in_min,   
    &disk_in_max,
    &disk_out_min,
    &disk_out_max,
    &v_turb_min,
    &v_turb_max,
    &T_rot0_fl_min,
    &T_rot0_fl_max,
    &T_rot_alpha_fl_min,
    &T_rot_alpha_fl_max,
    &rel_lum_min,
    &rel_lum_max,
    &layers_min,
    &layers_max,

    &inclination_min,
    &inclination_max,
    &H_den0_min,
    &H_den0_max,
    &H_den_alpha_min,
    &H_den_alpha_max,
			 
    &X12CO_13CO_fl_min,
    &X12CO_13CO_fl_max,
    &X12CO_C18O_fl_min,
    &X12CO_C18O_fl_max,
    &X12CO_13CO_cl_min,
    &X12CO_13CO_cl_max,
    &X12CO_C18O_cl_min,
    &X12CO_C18O_cl_max,

    &order1,
    &order2,
    &order3,
    &order4,
    &fileCount,
    &orders,
    &orderpix,
    &f_f,
    &f_i,
    &model_layers,
    &model_disk_in,
    &model_disk_out,
    &model_v_turb,
    &model_T_rot0_fl,
    &model_T_rot_alpha_fl,
    &model_rel_lum,
    &model_inclination,
    &model_H_den0,
    &model_H_den_alpha,
    &model_X12CO_13CO_fl,
    &model_X12CO_C18O_fl,
    &model_X12CO_13CO_cl,
    &model_X12CO_C18O_cl
  };

  string inputStrings[inputs] =
  {
    "inclination",
    "Mstar",
    "stardist",
    "Lc",
    "inst_res",
    "doDebug",
    "disk_in_0",
    "disk_out_0",
    "v_turb_0",
    "T_rot0_fl_0",
    "T_rot_alpha_fl_0",
    "rel_lum_0",
    "layers_0",
    
    "inclination_0",
    "H_den0_0",
    "H_den_alpha_0",

    "X12CO_13CO_fl_0",
    "X12CO_C18O_fl_0",
    "X12CO_13CO_cl_0",
    "X12CO_C18O_cl_0",

    "disk_in_min",
    "disk_in_max",
    "disk_out_min",
    "disk_out_max",
    "v_turb_min",
    "v_turb_max",
    "T_rot0_fl_min",
    "T_rot0_fl_max",
    "T_rot_alpha_fl_min",
    "T_rot_alpha_fl_max",
    "rel_lum_min",
    "rel_lum_max",
    "layers_min",
    "layers_max",

    "inclination_min",
    "inclination_max",
    "H_den0_min",
    "H_den0_max",
    "H_den_alpha_min",
    "H_den_alpha_max",
 
    "X12CO_13CO_fl_min",
    "X12CO_13CO_fl_max",
    "X12CO_C18O_fl_min",
    "X12CO_C18O_fl_max",
    "X12CO_13CO_cl_min",
    "X12CO_13CO_cl_max",
    "X12CO_C18O_cl_min",
    "X12CO_C18O_cl_max",

    "order1",
    "order2",
    "order3",
    "order4",
    "fileCount",
    "orders",
    "orderpix",
    "f_f",
    "f_i",

    "model_layers",
    "model_disk_in",
    "model_disk_out",
    "model_v_turb",
    "model_T_rot0_fl",
    "model_T_rot_alpha_fl",
    "model_rel_lum",
    "model_inclination",
    "model_H_den0",
    "model_H_den_alpha",
    "model_X12CO_13CO_fl",
    "model_X12CO_C18O_fl",
    "model_X12CO_13CO_cl",
    "model_X12CO_C18O_cl"

  };  



  //constant constants
  static constexpr double c=2.997924562e5;
  static constexpr double cexp=29979245620;
  static constexpr double hc=  6.626196e-27*2.997924562e10;
  //static constexpr double hck=(6.626196e-27*2.997924562e10)/(1.380622e-16);
  static constexpr double hck=1.43800;
  static constexpr double cer=8.85282e-13;
  static constexpr double rpi=1.7724538509055159;
  
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

  vec vfreq;

private:
//Molecular data
  cube X12CO = zeros<cube>(10,7,120);
  cube X13CO = zeros<cube>(3,7,120);
  cube XC18O = zeros<cube>(1,7,120);                         

  mat fXA; 

  



//variable parameters; these are changed with each iteration of the loop
/*  double layers;
  double dist;
  double stepsn&;
  double v_turb;
  double T_rot0_fl;
  double T_rot_alpha_fl;
  double T_rot0_cl;
  double T_rot_alpha_cl;
  double rel_lum;
  double disk_out;
  double disk_in;*/

public: 
  static mat wavenum            ;
  static mat einA               ;
  static mat lam_ang            ;
  static mat HD100546_luminosity;

  static vec freq;
  static vec r1big;
  static vec f1big;

  static double freq_size;  

  static constexpr double vib_einA[10]={34.60,67.68,98.40,126.99,153.59,178.31,201.35,223.10,244.15,265.21};

  //double** randData = new double*[7];
  bool* isSent;

  string folderpath;

  FitData(string folder);
  ~FitData();

  double* FillArrays(int modulus, int offset);
  int runCollisions(bool doCols);
 

int runTrial(double layers, double disk_in, double disk_out, double v_turb, double T_rot0_fl, double T_rot_alpha_fl, double rel_lum, double inclination, double H_den0, double H_den_alpha, double X12CO_13CO_fl, double X12CO_C18O_fl, double X12CO_13CO_cl, double X12CO_C18O_cl,  int locali);

  int runTrials();
  int runTrialsSlave();

  int readInput(string inpFile);
  int extractValue(string sin, string varname, double& var);

  int dbm(string message);

  string dirname;
};

//==============================================
//==external declarations of static parameters==
//=============================================

unsigned int FitData::numGuesses;

//tabular data
mat FitData::wavenum;
mat FitData::einA;
mat FitData::lam_ang;
mat FitData::HD100546_luminosity;

//input star parameters
double FitData::Mstar;
double FitData::stardist;
double FitData::inst_res;
double FitData::Lc; 
double FitData::doDebug;

//chi-by-eye parameters
double FitData::disk_in_0;
double FitData::disk_out_0;
double FitData::v_turb_0;
double FitData::T_rot0_fl_0;
double FitData::T_rot_alpha_fl_0;
double FitData::rel_lum_0;
double FitData::layers_0;

double FitData::inclination_0;
double FitData::H_den0_0;
double FitData::H_den_alpha_0;
double FitData::X12CO_13CO_fl_0;
double FitData::X12CO_C18O_fl_0;
double FitData::X12CO_13CO_cl_0;
double FitData::X12CO_C18O_cl_0;      

//range parameters
double FitData::disk_in_min;
double FitData::disk_in_max;
double FitData::disk_out_min;
double FitData::disk_out_max;
double FitData::v_turb_min;
double FitData::v_turb_max;
double FitData::T_rot0_fl_min;
double FitData::T_rot0_fl_max;
double FitData::T_rot_alpha_fl_min;
double FitData::T_rot_alpha_fl_max;
double FitData::rel_lum_min;
double FitData::rel_lum_max;
double FitData::layers_min;
double FitData::layers_max;   

double FitData::inclination_min;
double FitData::inclination_max;
double FitData::H_den0_min;
double FitData::H_den0_max;
double FitData::H_den_alpha_min;
double FitData::H_den_alpha_max;
                        
double FitData::X12CO_13CO_fl_min;
double FitData::X12CO_13CO_fl_max;
double FitData::X12CO_C18O_fl_min;
double FitData::X12CO_C18O_fl_max;
double FitData::X12CO_13CO_cl_min;
double FitData::X12CO_13CO_cl_max;
double FitData::X12CO_C18O_cl_min;
double FitData::X12CO_C18O_cl_max;


double FitData::order1;
double FitData::order2;
double FitData::order3;
double FitData::order4;

double FitData::fileCount;

static double fileCount;

double FitData::f_f;
double FitData::f_i;
double FitData::freq_size;
double FitData::orders;
double FitData::orderpix;

double FitData::model_layers         = 1;
double FitData::model_disk_in        = 1;
double FitData::model_disk_out       = 1;
double FitData::model_v_turb         = 1;
double FitData::model_T_rot0_fl      = 1;
double FitData::model_T_rot_alpha_fl = 1;
double FitData::model_rel_lum        = 1;

double FitData::model_inclination=0;
double FitData::model_H_den0=0;
double FitData::model_H_den_alpha=0;
double FitData::model_X12CO_13CO_fl=0;
double FitData::model_X12CO_C18O_fl=0;
double FitData::model_X12CO_13CO_cl=0;
double FitData::model_X12CO_C18O_cl=0;

double FitData::doModel[14];


//vibrational constants
constexpr double FitData::vib_einA[];
constexpr double FitData::Bv12[];
constexpr double FitData::Bv13[];

//input file parameters
constexpr unsigned int FitData::inputs;
constexpr double* FitData::inputVars[];

//frequency scale
vec FitData::freq;

//input data to be fit
vec FitData::r1big;
vec FitData::f1big;

//global FitData object
FitData* data;
#endif

