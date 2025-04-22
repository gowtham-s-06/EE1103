/*

            EE1103 Assignment - 9 : Newton and Lagrange Interpolation
            
Developer : EE24B017
Date : 4th October 2024
Purpose : To find the value of a function at a particular point using Lagrange and Newton interpolation.
Generate the Newton and lagrange polynomial n(x) and l(x) and hence test the code by using xo = 0.8.
Input(s) : polynomial_order xstart xend xo
Outputs(s): f(xo) n(xo) l(xo)
*/

//Including the necessary header files
#include <stdio.h>
#include <stdlib.h>


//Defining a function to generate the lagrange polynomial l(x) and hence calculate the value of the function at xx value.
//Here n refers to number of points (ie order+1)
float Lagrange(float* x, float* y, int n, float xx){
  float sum = 0;
  float product;

  for (int i=0;i<n;i++){
    product = y[i];    //Initialising the product term

    for (int j=0;j<n;j++){
      if (i!=j){
        product = product*(xx-x[j])/(x[i]-x[j]);
      }
    }

    sum = sum + product;  //Adding the product
  }

  return sum;
}

//Function to calculate f(xo) values using Newton Interpolation
//Here n refers to number of points (ie order+1)
float Newton(float* x, float* y, int n, float xx){
  float div_diff[n][n];     //Initializing a 2d array to store divided difference values
  
  //Initizlizing first column
  for (int i=0;i<n;i++){
    div_diff[i][0] = y[i];
  }

  //Computing divided difference
  for (int j=1;j<n;j++){
    for (int i=0;i<n-j;i++){
      div_diff[i][j] = (div_diff[i+1][j-1]-div_diff[i][j-1]) / (x[i+j]-x[i]);
    }
  }

  //Computing actual value
  float val = div_diff[0][0];
  float term = 1;

  for (int j=1;j<n;j++){
    term = term * (xx-x[j-1]);
    val = val + div_diff[0][j] * term;
  }

  return val;
}

float function(float x){
  return 1/(1 + 25*x*x);
}

int main(int argc, char *argv[]){

  if (argc != 5){
    printf("Usage: polynomial_order xstart xend x0\n");
    return 0;
  }

  //Taking in input values
  int order = atoi(argv[1]);
  float xstart = atof(argv[2]);
  float xend = atof(argv[3]);
  float xo = atof(argv[4]);
  
  //Checking the values
  if (xstart >= xend) {
    printf("xstart should be less than xend\n");
    return 0;
  }
  if (order < 1) {
    printf("Order should be a positive integer\n");
    return 0;
  }

  //Generating the datapoints
  float xvalues[order+1];
  float yvalues[order+1];
  float step = (xend - xstart)/order;
  
  for (int i=0;i<=order;i++){
    xvalues[i] = xstart + step * i;
    yvalues[i] = function(xvalues[i]);
  }
  
  //Calculating the Interpolated and actual values
  float lng = Lagrange(xvalues,yvalues,order+1,xo);
  float nwt = Newton(xvalues,yvalues,order+1,xo);
  float val = function(xo);
  
  //Printing the output values
  printf("%.4f %.4f %.4f\n",val,nwt,lng);
  return 0;

}
