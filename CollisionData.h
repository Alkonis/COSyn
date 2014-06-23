#include <fstream>

class CollisionData{


private:
  int layers;
  int steps;

public:
  double*** tau_0;    //these must be dynamically created according to the number of layers/steps we have... see inside RunTrial()!
  double*** dfdt_0;
  double*** dwdn;
  double*** Nv;
  double**** g;
  double**** rate_eqtn;

//Sets up arrays for collisional data generation
  CollisionData(int layers, int steps) {

    this->layers=layers;                          //may be able to optimize layers and steps... declare elsewhere
    this->steps=steps;


//arrays of form [10,12,layers]  -- check to make sure the 10 and 12 aren't backwards, here
    dfdt_0 = new double**[layers];
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


/*for k=0,steps-1 do begin
	for j=0,layers-1 do begin
		rate_eqtn(11:20,0:10,j,k)=EinA(0:9,0:10)
	endfor
endfor

		FOR i=0,8 DO BEGIN
			rate_eqtn(i+1,i+1,*,*)=-vib_EinA(i)
			rate_eqtn(i+2,i+1,*,*)=vib_EinA(i+1)
		ENDFOR
		rate_eqtn(1,0,*,*)=vib_EinA(0)
		rate_eqtn(10,10,*,*)=-vib_EinA(9)

		FOR i=0,8 do begin
			rate_eqtn(i+11,i+11,*,*)=-TOTAL(EinA(i,*),2)
		ENDFOR*/

    for (int i=0; i<steps-1; i++) {
        for (int j=0; i<layers-1; j++) {
            for (int k=11; k<21;k++) {
                for (int l=0; l<11; l++) {
                        rate_eqtn[k][l]

                }
                }

            }
            }
        }
    }



  }

  ~CollisionData() {

    for(int i=0; i < layers; i++) {
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

};
