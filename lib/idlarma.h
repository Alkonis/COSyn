#include<iostream>
#include<stdio.h>
//#include<armadillo>
#include<vector>




//*  WHERECLASS
//   C++ port of IDL WHERE() class for armadillo row/column vectors using lambda functions
//
//   USAGE:  idlwhere::where(input_data, [] (double datum) {return (input_condition)});
//      

namespace idlarma {

  using namespace std;
  using namespace arma;

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

  vec deriv(vec xin, vec yin) {
    int inputSize=xin.size();
    vec output = zeros<vec>(inputSize);
    
    double x0, x1, x2;
    double y0, y1, y2;
    
    for (int i=1; i<inputSize-1; i++)
    {
      x0=xin.at(i-1);  x1=xin.at(i);  x2=xin.at(i+1);
      y0=yin.at(i-1);  y1=yin.at(i);  y2=yin.at(i+1);
      output.at(i)= y0*(x1-x2)/( (x0-x1)*(x0-x2) ) +
                    y1*(1/(x1-x2)-1/(x0-x1)) -
                    y2*(x0-x1)/((x0-x2)*(x1-x2)); 
    }

    output.at(0) = yin.at(0)*(xin.at(0)-xin.at(1)+xin.at(0)-xin.at(2))/((xin.at(0)-xin.at(1))*(xin.at(0)-xin.at(2))) -
                   yin.at(1)*(xin.at(0)-xin.at(2))/((xin.at(0)-xin.at(1))*(xin.at(1)-xin.at(2))) +
                   yin.at(2)*(xin.at(0)-xin.at(1))/((xin.at(0)-xin.at(2))*(xin.at(1)-xin.at(2)));

    x0=xin.at(inputSize-3);   x1=xin.at(inputSize-2); x2=xin.at(inputSize-1);
    y0=yin.at(inputSize-3);   y1=yin.at(inputSize-2); y2=yin.at(inputSize-1);
    
    output.at(inputSize-1) = - y0*(x1-x2)/((x0-x1)*(x0-x2))
                             + y1*(x0-x2)/((x0-x1)*(x1-x2))
                             - y2*(x0-x2+x1-x2)/((x0-x2)*(x1-x2));
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
