#include "COfit.h"
#include <iostream>
#include <cstdlib>

//optimizations so far:  read in files ONCE in FitData() rather than in each iteration
//                       removal of redundant declarations
//

//What to do:



//To do:  Fix variable names
//        Benchmark!
//        Keep going!!!!!!!!!
//        Woohoo!!!

using namespace std;
using namespace idlarma;

//FitData* Data;

constexpr double FitData::vib_einA[];


int FitData::numGuesses;

//returns an array of size numGuesses filled with random numbers on the interval (offset, modulus+offset)
double* FitData::FillArray(int modulus, int offset)
{
  double* array = new double[numGuesses];
  for(int i=0; i<numGuesses;i++) {
    array[i]=rand() % modulus + offset;
  }
  return array;
}



int FitData::runCollisions(bool doCols){             //declare F_tau!   Get matrix from cube view?!?!?!?!!?!?!@?!?!?!?!?!? 
  if (rel_lum > 0.001)   
  {
    vec b_tot = sqrt(2*1.38e-16*6.02e23*T_rot_fl/28 + pow(v_turb,2));
    vec tau = linspace<vec>(0,20,2000);   //check valueis here   idl = findgen(2e3)/1e2) === findgen(2000/100)  == [0, 20], N=2000
    vec F_tau=zeros<vec>(2000);
 
    //==============INTEGRATE FUNCTION USING TRAPEZOID RULE================
    
    double sum = 0;
    auto n = 1000;
    auto a = 0;
    auto b = 5;

    auto delta=(b-a)/n;
    
    for (int i=0; i<2000; i++) 
    {
      sum=0.5 * ( (1 - exp(-tau(i)*exp(-a^2))) + (1 - exp(-tau(i)*exp(-b^2))) );
      for (int j=0; j<n; j++)  
      {
        auto x=a+delta*j;
        sum=sum+(1-exp(-tau(i)*exp(-x^2)));
      }
      sum = sum*delta;
      F_tau(i)=sum;
    }
    
    vec dFdt=deriv(tau,F_tau);
    
    //=============COLLISIONS=============
    for (int k=0; k<steps; k++)
    {
      for (int z=0; z<8; z++)
      {
        //rate_eqtn.at(i+1,0).subcube(span(i+1),span::all,span(k)).fill(-vib_einA[i]);
       // rate_eqtn.at(z+1,0).subcube(span(z+1),span::all,span(k))=-vib_einA[z];
        rate_eqtn.at(z+1,0).slice(k).col(z+1).fill(-vib_einA[z]);
        
        rate_eqtn.at(z+2,0).slice(k).col(z+1).fill(vib_einA[z+1]);
        //rate_eqtn.at(z+2,0).subcube(span(z+1),span::all,span(k))=vib_einA[z+1];
      }
      rate_eqtn.at(1,0).subcube(span(0),span::all,span(k))=vib_einA[0];
      rate_eqtn.at(9,0).subcube(span(9),span::all,span(k))=-vib_einA[8];

     if (doCols == 0)
     {
       auto k_HCO_dn = (7.57e-15*T_rot_cl[k]/(1-exp(-3084/T_rot_cl[k])))*H_den[k];
       auto k_HCO_up=k_HCO_dn*exp(-3084/T_rot_cl[k]);

       rate_eqtn.at(0,0).subcube(span(0),span::all,span(k))=rate_eqtn.at(0,0).subcube(span(0),span::all,span(k))-k_HCO_up;
       rate_eqtn.at(1,0).subcube(span(0),span::all,span(k))=rate_eqtn.at(1,0).subcube(span(0),span::all,span(k));
       
       for (int i=1; i<9; i++) 
       {
         rate_eqtn.at(i-1,0).subcube(span(i),span::all,span(k)) = rate_eqtn.at(i-1,0).subcube(span(i),span::all,span(k))+k_HCO_up;
         rate_eqtn.at(i  ,0).subcube(span(i),span::all,span(k)) = rate_eqtn.at(i  ,0).subcube(span(i),span::all,span(k))-k_HCO_dn-k_HCO_up;
         rate_eqtn.at(i+1,0).subcube(span(i),span::all,span(k)) = rate_eqtn.at(i+1,0).subcube(span(i),span::all,span(k))+k_HCO_dn;
       }
     
    
       //===================ULTRAVIOLET PUMPING========================
       
       auto Fuv=HD100546_luminosity/(4*3.1415926535897)*rdisk[k];
       
       for (int j=0; j<layers; j++)
       {
         if (j>0)
         {
           for (int i=0; i<10; i++)
           {
            // tau_0.subcube(span(i),span::all,span(j))=sum(Nv.subcube(span(0,11),span(0,j),span(k))) * 7.55858e12 * 0.02654*fXA.submat(span(i),span::all)*lam_ang.submat(span(i),span::all)*1e-8/(sqrt(3.1415926535897)*b_tot[k]);   //Check this!
            tau_0.subcube(span(i),span::all,span(j))=accu(Nv.slice(k).submat(span(0,11),span(0,j))) *  7.55858e12 * 0.02654*fXA.submat(span(i),span::all)*lam_ang.submat(span(i),span::all)*1e-8/(sqrt(3.1415926535897)*b_tot[k]);
           }
         }
         
         for (int ii=0; ii<tau_0.slice(j).col(0).n_elem; ii++)
         {
           if (tau_0.at(0,ii,j) < 200)
           {
             ivec dFdt_0_index=where(tau, [&] (double elem) {return elem == round(tau_0.at(0,ii,j)*10)/10;});
             auto count = dFdt_0_index.size();
             
             if (count != 0)
             {
                 ///dFdt_0.subcube(span::all,span(ii),span(j))=dFdt[dFdt_index];    weird line!!!!!
             }

           }
           else 
           {
             dFdt_0.slice(j).row(ii)=(1/(2*tau_0.at(0,ii,j))*sqrt(log(tau_0.at(0,ii,j))));    //(span::all,span(i),span(j))=1/(2*tau_0.at(0,ii,j))*sqrt(alog(tau_0.at(0,ii,j)));
           }
         }
         //dWdN.subcube(span::all,span::all,span(j))=dfdt_0.subcube(span::all,span::all,span(j))*.02654*2.0 % (lam_ang*1e-4) % (lam_ang*1e-8) % fXA/(sqrt(3.1415926535897)*c*1e5);
         for (int ii=0; ii<10; ii++) 
         {
          // g.at(ii,0).subcube(span::all,span(j),span(k))=dwdn.subcube(span::all,span::all,span(j))*3.1415926535897 % Fuv / (hc * wavenum);
          g.at(ii,0).slice(k).row(j) = dwdn.slice(j) * 3.1415926535897 % Fuv / (hc * wavenum);
         }         
         //add in g-factors:
         
         double gsum=0;
         
         for (int ii=0; ii<10; ii++)
         { 
           gsum+=g.at(ii).at(0,j,k);
         }
         rate_eqtn.at(0,0).at(0,j,k)=rate_eqtn.at(0,0).at(0,j,k)-gsum;
         
         for (int i=0; i<10; i++)
         {
           gsum = 0;
           for (int ii=0; ii<10; ii++)
           {
             gsum+=g.at(ii).at(i,j,k);
           }
           rate_eqtn.at(i,0).subcube(span(i),span(j),span(k))=rate_eqtn.at(i,0).subcube(span(i),span(j),span(k)) - gsum;
         }

         for (int i=0; i<8; i++)
         {
           for (int ii=0; i<11; i++)
           {
             rate_eqtn.at(ii,0).at(11+i,j,k)=rate_eqtn.at(ii,0).at(ii+i,j,k)+g.at(i,0).at(ii,j,k);
           }
         }
         rate_eqtn.at(1,0).at(10,j,k)=0;
         
         mat rateToSolve= zeros<mat>(21,21);
         for (int i=0; i<21; i++)
         {
           //rateToSolve.submat(span(i),span::all)=rate_eqtn.at(i).subcube(span::all,span(j),span(k));
           rateToSolve.col(i)=rate_eqtn.at(i).slice(k).row(j).t();
         }
        // vec rateSVD= svd(rateToSolve);
        // vec rateSolution=solve(rateToSolve);
       }
     }
   }
  }
  return 0;
}

int FitData::runTrial() {
//fit model, compare, etc...
//ranulli is the port of the IDL variable r--fix current r variable?


//Define variables
/*HD100546_luminosity = HD100546_luminosity*L;  //Luminosity scaling
mat fAX = (3.038/2.03)*EinA/wavenum^2
mat fXA = 2*fAX*/

//============================
//Divide into annulli
//=============================
  cerr << "Diving into annulli...";
  ranulli.push_back(r); 
  if (r<0.1)
  {
    while ((ranulli.back() >= 1) && (ranulli.back() >= disk_out))
    {
      r_index_a=r_index_a+1;
      ranulli.push_back(disk_in+0.01*r_index_a);
    }
  }

  maxra = ranulli.back();
 
  if ((maxra < 1) && (maxra >= .1))
  {
    while ((ranulli.back() >= 1.0) || (ranulli.back() >= disk_out))
    {
     r_index_b++;
     ranulli.push_back(maxra+0.1*r_index_b);
    }
  }
  maxrb=ranulli.back();
 
  if ((maxrb <= disk_out) && (maxrb >= 1.0))
  {
    while  (ranulli.back() <= disk_out)
    {
      r_index_c++;
      ranulli.push_back(maxrb+1.0*r_index_c);
    }
  }
 
  steps=ranulli.size();  //check that this size function is correct
  
  vec rdisk= zeros<vec>(steps);
 
  for (int i=0; i < steps; i++ )  {
    rdisk.at(i)=ranulli.at(i)*1.496E13;
  }

//benchmark here!  Make sure these results are the same from the IDL code and the C++ code!
  
  T_rot_fl = T_rot0_fl *pow((1.5E13/rdisk),T_rot_alpha_fl);
  ivec T_rot_index=where(T_rot_fl, [] (double datum) {return datum >= 3500;});

  int T_rot_cnt = T_rot_index.size();

  if (T_rot_cnt != 0) {
    for (auto elem : T_rot_index )
    {
      T_rot_fl[elem]=3500;
    }
  }

  T_rot_cl=T_rot0_cl * pow((1.496E13/rdisk), T_rot_alpha_cl);
  H_den=H_den0 * pow((1.496E13/rdisk),H_den_alpha);
  T_rot_index = where(T_rot_cl, [] (double datum) {return datum > 3500;});
  T_rot_cnt=T_rot_index.size();

  if (T_rot_cnt != 0) {
    for (auto elem : T_rot_index)
    {
      T_rot_cl[elem]=3500;
    }
  }

  cerr << "Starting collisions...";
  runCollisions(0);
  runCollisions(1);
  
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

  for(int i=0; i<this->numGuesses;i++) {
    //set parameters here
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

  fin.open("EinA.txt");
  for ( auto& entry : einA) {
     fin >> entry;
  }
  fin.close(); 

  fin.open("lambda.txt");
  for ( auto& entry : lam_ang) {
     fin >> entry;
  }
  fin.close();

  fin.open("wavenum.txt");
  for ( auto& entry : wavenum) {
     fin >> entry;
  }
  fin.close();

  fin.open("HD100546_luminosity.txt");
  for (auto& entry : HD100546_luminosity) {
     fin >> entry;
  }
  fin.close();

  fin.open("CO_molecdat/X12CO");
  for (int i=0; i<119; i++) {
    for (int j=0; j<7; j++) {
      for (int k=0; k<10; k++) {
        fin >> X12CO(k,j,i);
      }
    }
  }
  fin.close();

  fin.open("CO_molecdat/X13CO");
  for (int i=0; i<119; i++) {
    for (int j=0; j<7; j++) {
      for (int k=0; k<3; k++) {
         fin >> X13CO(k,j,i);
      }
    }
  }
  fin.close();



//initialize arrays with random data
  cerr << "Generating arrays...";
  this->randData[0]=FillArray(900, 100);
  this->randData[1]=FillArray(20, 5);
  this->randData[2]=FillArray(75, 50);
  this->randData[3]=FillArray(500000, 50000);
  this->randData[4]=FillArray(3500,1000);
  this->randData[5]=FillArray(101,0);
  this->randData[6]=FillArray(99, 1);
  for (int i=0; i <= numGuesses; i++) {
    this->randData[5][i]=this->randData[5][i]/100;
  };
  cerr << "Done.";
cin.get();
 //print arrays for debugging purposes
  for (int c=0; c<7; c++)
  {
    cout << "randData[" << c << "]:  ";
    for (int i=0; i<numGuesses; i++) {
      cout << this->randData[c][i] << ' ' ;
    }
    cin.get();
  };
 
 /*=================
 *
 * Collision Data                                      EVERYTHING AFTER LAYERS ----> INTO COLLISIONS!  CollisionData.cpp will need to be revisited...
 *
 *==================*/
  cerr << "Preparing collision data." << endl;
  HD100546_luminosity = HD100546_luminosity*rel_lum;
  
  fAX = (3.038/2.03)*einA/(wavenum % wavenum);   //be sure this is right!
  fXA = 2*fAX;

  cerr << "Creating cubes." << endl;
  cerr << "Layers,steps:" << endl;
  cerr << layers;
  cerr << steps;
  dFdt_0 = zeros<cube>(10,12,layers);
  tau_0  = zeros<cube>(10,12,layers);
  dwdn   = zeros<cube>(10,12,layers);
  Nv     = zeros<cube>(21,layers,steps);
  
  cerr << "Creating rate equation." << endl;
  rate_eqtn = field<cube>(21,1);
  g=field<cube>(10,1);
  cerr << "Test1." << endl;

  for (int i=0; i<10; i++)
  {
    g.at(i,0)=zeros<cube>(12,layers,steps);
  }
  cerr <<"Made g." << endl;
  for (int i=0; i<21; i++) 
  {
    cout << i << endl;
    rate_eqtn.at(i,0)=zeros<cube>(21,layers,steps);
    cout << rate_eqtn.at(i,0);
    rate_eqtn.at(i,0)(span(20),span::all,span::all).fill(1);
  }
  cerr << "Test." << endl;
  for (auto i=0; i<steps;i++) {
    for (auto j=0; i<layers;j++) {
      for (auto k=0; k<9;k++) {
        for (auto q=0; q<10; q++)
        {
          rate_eqtn.at(k+11,0)(q,j,k)=einA(k,q);
        }
      }
    }
  }
  
/*  for (int k=0; k<steps; k++)
  {
    for (int i=0; i<8 i++)
    {
      rate_eqtn.at(i+1,0).subcube(span(i+1),span::all,span(k))=-vib_einA(i);
      rate_eqtn.at(i+2,0).subcube(span(i+1,span::all,span(k))=vib_einA(i+1);

    }

  }*/
  
  for (auto i=0; i<9; i++)
  {
    //rate_eqtn.at(i+1,0)(span(i+1),span::all,span::all) = vib_einA[i];
    //rate_eqtn.at(i+2,0)(span(i+1),span::all,span::all) = vib_einA[i+1];   //vib_EinA is not defined in this scope--check this and move
  }

  rate_eqtn.at(1,0)(span(0),span::all,span::all)=vib_einA[0];   //vib_EinA will need to be some kind of Arma class for this to work... I will need to test this
  rate_eqtn.at(10,0)(span(10),span::all,span::all)=vib_einA[9];  //look up vib_EinA
  
  for (auto i=0; i<9; i++) 
  {
    rate_eqtn.at(i+11,0)(span(i+1),span::all,span::all) = sum(einA.row(i));
  } 

  //check these threshold-sets ^^^^^^^
  cerr << "Done.";
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

int main(int argc, char* argv[]) 
{
  cout <<  "Starting...";
  FitData* data = new FitData(atoi(argv[1]));
  data->runTrials();
  delete data;
  return 1;
}
