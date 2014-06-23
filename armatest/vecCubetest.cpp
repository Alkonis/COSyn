#include<iostream>
#include<stdio.h>
#include<armadillo>

using namespace std;
using namespace arma;



int main() {

  field<cube> data(1000,1);

  for (int i=0; i<1000; i++)
  {
    data.at(i,0)=randi<cube>(100,100,100);
  cout << i << endl;
  }
//  cout << data;
  return 0;
}



