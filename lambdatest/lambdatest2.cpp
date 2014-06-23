#include<iostream>
#include<stdio.h>
#include<armadillo>

using namespace arma;
using namespace std;



//*  WHERECLASS
//   C++ port of IDL WHERE() class for armadillo row/column vectors using lambda functions
//
//   USAGE:  where.wvec(input_vector, [] (double datum) {return (input_condition)});
//      OR:  where.wrow(input_row, ... )
//  
//
//

class whereclass {
  public:
 
  
  vec wvec(vec data, bool lambda) {

    int count=0;
    vec output = zeros<vec>(data.size());

    for (auto elem : data) 
    {
      if (lambda(elem)) {
        output.at(count)=elem;
        count++; 
      }
    }

    output.resize(count);

    return output;
  }

  vec wrow(rowvec data, bool lambda) {
    
    int count=0;
    vec output = zeros<vec>(data.size());
    
    for (auto elem : data)
    {
      if (lambda(elem)) {
        output.at(count)=elem;
        count++;
      }
    }
    
    output.resize(count);
    
    return output;
  }

};

whereclass where;


int main() {
  vec data = randi<vec>(100,distr_param(100,200));
  auto indices=where.wvec(data, [] (double datum) {return datum > 175;} );
  cout << "Indices";
  cout << indices;
  cout << "Data value:";
  for( auto i : indices) {
    cout << data[i] << endl;;
  }
  return 0;
}
