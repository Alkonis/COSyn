#include <iostream>
#include <stdio.h>
#include <armadillo>

using namespace arma;
using namespace std;


int main() {
  field<cube> data(10);
  auto j=0;
  for (auto i=0; i<10; i++) {
    data(i)=zeros<cube>(2,3,4);
  /*  for (auto& entry : data(i))
    {
     entry=j;
     j++;
    }*/
    cube sub = data.at(i);
//    cout << sub;
    cube subsub = sub(span(1),span::all,span::all);
    subsub.fill(123);                                                                                 
    cout << subsub;
    sub(span(1),span::all,span::all).fill(123);
    cout<<sub;
  }
// cout << "======================";
// cout << data;

  return 0;
}



