#include<stdio.h>
#include<iostream>
#include<armadillo>
#include "idlwhere.h"

using namespace std;
using namespace arma;


int main()  {
  vec x = {1,2,3,4,5};
  vec y = {2,4,6,8,10};

  vec z = idlwhere::deriv(x,y);
  cout << z;
  return 0;
}
