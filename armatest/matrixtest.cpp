#include <iostream>
#include <stdio.h>
#include <armadillo>

using namespace arma;
using namespace std;


int main() {

  mat testmat = zeros<mat>(10,12);
  ifstream fin;
  fin.open("lambda.txt");

  for ( auto& entry : testmat) {
       fin >> entry;
  }
  cout << testmat;
  fin.close();

  ofstream fout;
  fout.open("lamba2.txt");
  
auto i=0;
  for (auto& entry : testmat) {
      fout << entry;  
      i++;
      if(i==10){
      fout << endl;
      i=0;}
      else{ fout << " ";;} 
  }

  fout.close();

  return 0;

}




