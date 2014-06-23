#include <iostream>
#include <stdio.h>
#include <armadillo>
#include <fstream>

using namespace arma;
using namespace std;

int main() {
  ifstream fin;

  fin.open("CO_molecdat/X12CO");
  cube X12CO = zeros<cube>(10,7,120);
  for (int i=0; i<119; i++) {
    for (int j=0; j<7; j++) {
      for (int k=0; k<10; k++) {
        fin >> X12CO(k,j,i);
      }
    }
  }
  fin.close();
  
  cube X13CO = zeros<cube>(3,7,120);
  fin.open("CO_molecdat/X13CO");
  for (int i=0; i<119; i++) {
    for (int j=0; j<7; j++) {
      for (int k=0; k<3; k++) {
         fin >> X13CO(k,j,i);
      }
    }
  }
  fin.close();

  cout << X13CO;
  return 0;
}
