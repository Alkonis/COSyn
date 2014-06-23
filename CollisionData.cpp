#include <fstream>

class CollisionData{


private:
  int layers;
  int steps;

public:
 /* double*** tau_0;    //these must be dynamically created according to the number of layers/steps we have... see inside RunTrial()!
  double*** dfdt_0;
  double*** dwdn;
  double*** Nv;
  double**** g;
  double**** rate_eqtn;
*/
  //=====MOVED TO COFIT.cpp
  /*cube dfdt_0;
  cube tau_0;
  cube dwdn;
  cube Nv;
  field<cube> rate_eqtn;
  field<cube> g;*/
//Sets up arrays for collisional data generation
  CollisionData(int layers, int steps) {

 //   this->layers=layers;                          //may be able to optimize layers and steps... declare elsewhere
  //  this->steps=steps;


//arrays of form [10,12,layers]  -- check to make sure the 10 and 12 aren't backwards, here
    /*dfdt_0 = new double**[layers];
    dwdn = new double**[layers];
    tau_0 = new double**[layers];

    for(int i=0; i < layers; i++) {
      dfdt_0[i] = new double*[12];
      dwdn[i]=new double*[12];
      tau_0[i] = new double*[12];
      for( int j=0; j < 12; j++) {
        dfdt_0[i][j] = new double[10];
        dwdn[i][j] = new double[10];
        tau_0[i][j] = new double[10];
      }
    }*/


//==============MOVED TO COFIT.H====== 
  /*dfdt_0 = zeros<cube>(10,12,layers);
  tau_0  = zeros<cube>(10,12,layers);
  dwdn   = zeros<cube>(10,12,layers);
  Nv     = zeros<cube>(21,layers,steps);
  
  rate_eqtn = field<cube>(10,1);
  g=field<cube>(10,1);

  for (int i=0; i<10; i++) {
    rate_eqtn.at(i,0).subcube(span(20),span::all,span:all).fill(1);
  }

  for (auto i=0; i<steps;i++) {
    for (auto j=0; i<layers;j++) {
      for (auto k=0; k<9;k++) {
        for (auto q=0; q<10; q++) {
        {
          rate_eqtn.at(k+11,0)(q,j,k)=EinA(k,q)
        }
      }
    }
  }

  for (auto i=0; i<9; i++) {
                //two declarations right in front of each other in IDL.  Investigate maybe
    rate_eqtn.at(i+2,0)(span(i+1),span::all,span::all)=vib_EinA(i+1)   //vib_EinA is not defined in this scope--check this and move!
  }
  rate_eqtn.at(1,0)(span(0),span::all,span::all)=vib_EinA(0)   //vib_EinA will need to be some kind of Arma class for this to work... I will need to test this
  rate_eqtn.at(10,0)(span(10),span::all,span::all)=vib_EinA(9)  //look up vib_EinA

  for (auto i=0; i<9; i++) {
    rate_eqtn.at(i+11,0)(span(i+1),span::all,span::all) = sum(EinA.row(i)) 
  } 
   */

  T_rot_fl=T_rot0_fl*(1.5e13/rdisk)^T_rot_alpha_fl;
  T_rot_index=where(T_rot_cl, [] (double datum) {return datum > 3.5E3;};
  T_rot_cnt=T_rot_index.n_elem;
  
  if (T_rot_cnt > 0) 
  {
    for (auto index : T_rot_cnt) 
    {
      T_rot_fl(index)=3.5E3;
    }
  }
 
  T_rot_cl=T_r0t0_cl * (1.496e13/rdisk)^T_rot_alpha_cl;
  H_den_H_den0 * (1.496E13/rdisk)^H_den_alpha;

  T_rot_index=0;



 
/*  for (int i=0; i<10; i++) {
    g(i)=zeros<cube>(12,layers,steps);
    rate_eqtn=zeros<cube>(12,layers,steps)
  }

//array of form [21, layers, steps]
    Nv = new double**[steps];
    for(int i=0; i<steps; i++) {
      Nv[i] = new double*[layers];
      for(int j=0; j < layers; j++) {
        Nv[i][j]=new double[21];
      }
    }

//array of form [*,*,layers,steps]
    g= new double***[steps];
    rate_eqtn = new double***[steps];
    for (int i=0; i<steps; i++) {
      g[i] = new double**[layers];
      rate_eqtn[i] = new double**[layers];
      for(int j=0; j<layers; j++) {
        g[i][j] = new double*[12];
        rate_eqtn[i][j] = new double*[21];
        for(int k=0; k<12; k++) {
          g[i][j][k] = new double[10];
        }
        for (int k=0; k<21; k++) {
          rate_eqtn[i][j][k] = new double[21];
          for (int l=0; l<21; l++) {
            if (k==20) {rate_eqtn[i][j][k][l]= 1;}
            else {rate_eqtn[i][j][k][l]=0;}
          }
        }
      }
    }



    for (int i=0; i<steps-1; i++) {
        for (int j=0; i<layers-1; j++) {
            for (int k=11; k<21;k++) {
                for (int l=0; l<11; l++) {
                       //rate_eqtn[k][l]

                }
                }

            }
            }
        
*/


  }

  ~CollisionData() {

  /*  for(int i=0; i < layers; i++) {
      for( int j=0; j < 10; j++) {
        delete dfdt_0[i][j];
        delete dwdn[i][j];
        delete tau_0[i][j];
      }
      delete dfdt_0[i];
      delete dwdn[i];
      delete tau_0[i];
    }
    delete dfdt_0;
    delete dwdn;
    delete tau_0;

    for(int i=0; i<steps; i++) {
      for(int j=0; j < layers; j++) {
        delete Nv[i][j];
      }
      delete Nv[i];
    }
    delete Nv;

    for(int i=0; i<steps; i++) {
      for(int j=0; j<layers; j++) {
        for(int k=0; k<12; k++) {
          delete g[i][j][k];
        }
        for(int k=0; k<21; k++) {
          delete rate_eqtn[i][j][k];
        }
        delete g[i][j];
        delete rate_eqtn[i][j];
      }
      delete g[i];
      delete rate_eqtn[i];
    }
    delete g;
    delete rate_eqtn;

	}
*/
};
