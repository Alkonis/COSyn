#include<stdio.h>
#include<iostream>
#include<armadillo>

using namespace std;
using namespace arma;

int main()  {
  field<cube> F(2,1);
  for (int i=0; i<2; i++) {
    F(i,0) = randi<cube>(5,5,2);
  }
  field<cube> F2 = F;
  F2(1,0)=F2(1,0)*5;
  cout << "F:" << endl;
  cout << F;
  cout << "F2:" << endl;
  cout << F2;
  return 0;
}
