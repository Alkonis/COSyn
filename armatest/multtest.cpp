#include <iostream>
#include<armadillo>
#include<stdio.h>

using namespace arma;

int main()  {
  mat A = mat(5,5);
  A.fill(123);
 
  mat B = mat(5,5);
  B.fill(10);

  mat C = A % B;
  cout << C;

  mat D = A * B;
  cout << D;

  return 0;
}
