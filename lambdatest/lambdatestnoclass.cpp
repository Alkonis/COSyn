#include<iostream>
#include<stdio.h>
#include<armadillo>
#include<vector>

using namespace arma;
using namespace std;



//*  WHERECLASS
//   C++ port of IDL WHERE() class for armadillo row/column vectors using lambda functions
//
//   USAGE:  idlwhere::where(input_data, [] (double datum) {return (input_condition)});
//      

namespace idlwhere {

  template<typename cond>
  
  ivec where(vec data, cond lambda) {

    int count=0;
    int outputPointer=0;
    ivec output = zeros<ivec>(data.size());

    for (auto elem : data) 
    {
      if (lambda(elem)) {
        output.at(outputPointer)=count;
        outputPointer++; 
      }
      count++;
    }

    output.resize(outputPointer);

    return output;
  }

  template<typename cond>

  ivec where(rowvec data, cond lambda) {
    
    int count=0;
    int outputPointer=0;
    ivec output = zeros<ivec>(data.size());
    
    for (auto elem : data)
    {
      if (lambda(elem)) {
        output.at(outputPointer)=count;
        outputPointer++;
      }
      count++;
    }
    
    output.resize(outputPointer);
    
    return output;
  }

  template<typename cond>
  ivec where(std::vector<double> data, cond lambda) {
    
    int count=0;
    int outputPointer=0;
    ivec output = zeros<ivec>(data.size());

    for (auto elem: data) 
    {
      if (lambda(elem)) {
        output.at(outputPointer)=count;
        outputPointer++;
      }
      count++;
    }

    output.resize(outputPointer);

    return output;
  }

}

/*EXAMPLE PROGRAM

using namespace idlwhere;

int main() {
  vec data = randi<vec>(100,distr_param(100,200));
  auto indices=where(data, [] (double datum) {return datum > 175;} );
  cout << "Indices (should be 0-100)";
  cout << indices;
  cout << "Data value: (should be >175)";
  for( auto i : indices) {
    cout << data[i] << endl;
  }
  cout << "Data:";
  for (int i=0; i<100; i++ ) {cout << data.at(i) << " ";}
  return 0;
}
*/
