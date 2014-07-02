#include<iostream>
#include<stdio.h>
#include<armadillo>
#include<vector>




//*  IDLARMA
//  
//  C++ implentation of useful IDL functions wrapped through the Armadillo linear algebra libary.


//   idlarma::where
//   C++ port of IDL WHERE() function for armadillo row/column vectors using lambda functions
//   USAGE:  idlwhere::where(input_data, [] (double datum) {return (input_condition)});

//   idlarma::deriv
//   C++ port of the IDL deriv function
//   USAGE:  deriv(y,x), computes dy/dx


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


    // Resize to the array of the size needed to support the indices.
    // If no matches were found, then it will return an ivec of length 1 with value -1.
    // Since -1 is an invalid index, this can be used to detect a match failure.    
    if (outputPointer!=0) 
    {
      output.resize(outputPointer);
    }
    else
    {
      output.resize(1);
      output.at(0)=-1;
    }

    return output;
  }

  template<typename cond>

  ivec whererow(rowvec data, cond lambda) {
    
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
    
   
    // Resize to the array of the size needed to support the indices.
    // If no matches were found, then it will return an ivec of length 1 with value -1.
    // Since -1 is an invalid index, this can be used to detect a match failure.    
    if (outputPointer!=0) 
    {
      output.resize(outputPointer);
    }
    else
    {
      output.resize(1);
      output.at(0)=-1;
    }
 
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



    // Resize to the array of the size needed to support the indices.
    // If no matches were found, then it will return an ivec of length 1 with value -1.
    // Since -1 is an invalid index, this can be used to detect a match failure.    
    if (outputPointer!=0) 
    {
      output.resize(outputPointer);
    }
    else
    {
      output.resize(1);
      output.at(0)=-1;
    }

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


  mat totalDim(cube input, int dim)
  {

    mat output;
    
    int i_max;
    int j_max;

    if (dim==1) {
      i_max=input.n_cols;
      j_max=input.n_slices;
    } else if (dim==2) {
      i_max=input.n_rows;
      j_max=input.n_slices;
    } else if (dim==3) {
      i_max=input.n_rows;
      j_max=input.n_cols;
    } else {
      output << 0 << 0 << endr
             << 0 << 0 << endr;
      return output;
    }

    output=zeros<mat>(i_max,j_max);
   
    for (int i=0; i<i_max; i++)
    {
      for (int j=0; j<j_max; j++)
      {
        if (dim==1) {
          output(i,j)=accu(input(span::all,span(i),span(j)));
        } else if (dim==2) {
	  output(i,j)=accu(input(span(i),span::all,span(j)));
        } else if (dim==3) {
          output(i,j)=accu(input(span(i),span(j),span::all));
        }
      }
    }
    return output;
  }
}

/*EXAMPLE PROGRAM--WHERE

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

/*EXAMPLE PROGRAM--TOTALDIM */


/*int main() {
  cube data = zeros<cube>(2,2,2);
 
  data.slice(0) << 1 << 3 << endr
                << 5 << 7 << endr;
  data.slice(1) << 2 << 4 << endr
                << 6 << 8 << endr;

  cout << data << endl;

  mat tot = totalDim(data, 2);
  
  cout << tot << endl;

  return 0;

}*/
