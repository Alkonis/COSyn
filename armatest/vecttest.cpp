#include<iostream>
#include<stdio.h>
#include<armadillo>

using namespace std;
using namespace arma;

int main() {
  mat gg = zeros<mat>(3,5);
  rowvec vectest = ones<rowvec>(5);

  gg(1,span::all)=vectest;

  cout << gg;
  return 0;
}
