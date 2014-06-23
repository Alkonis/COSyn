#include<iostream>
#include<stdio.h>
#include<armadillo>

using namespace arma;
using namespace std;


int main()  {
  mat A = randi<mat>(10,12);
  mat B = randi<mat>(10,12);
  cout << A;
  cout << B;
  mat C = A % B;
  
  cout << C;
  
  mat D = randi<mat>(100000,100000);  
  cout<<"Done.";
  auto count=0;
  auto count2=0;
  for (auto elem : D) 
  {
    count++;
    if (count>999) 
    {
      count2++;
      cout << count2*1000 << endl;;
      count=0;
    }

  }
  
  
  return 0;
}
