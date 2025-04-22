/*

            EE1103 Assignment - 2 : Series Truncation
            
Developer : EE24B017,EE24B069
Date : 23rd August 2024
Purpose : To calculate the value of sinx using a truncated series upto a a given number of terms and hence 
calculate the error between the truncated value of sinx and the actual value of sinx obtained using the
math.h library
Input(s) : N (the number of terms upto which sinx is to be calculated) , x (the value for which the sin() function must be evaluated)
Outputs(s): n (the value obtained by evaluating the truncated series of sinx), e (the error between the actual value of sinx obtained
using the sin function in the math.h library and the calculated value of sinx calculated upto N terms as a truncated series)
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//using unsigned long long int because int is not sufficient after 12!
unsigned long long int factorial(unsigned long long int n){
  unsigned long long int fact=1;
  if (n<1) {
    return 0;
  }
  else if (n>1) {
    fact=n*factorial(n-1);
  }
  else if (n==1) {
    return fact;
  }
}
int main(int argc, char **argv) {
  if (argc<=2) {
    printf("Please input atleast 2 arguments\n");   //To confirm atleast 2 arguments are given
    return 0; 
  }
  int N=atoi(argv[1]);
  float x=atof(argv[2]);
  double sinx=0;       //To store all the decimals of sinx
  int l=1;             //To change the sign of the additional term to be added in sinx
  for (int i=1; i<=N; i++){
    int m=i*2-1;       //To calculate the exponent and factorial
    sinx+=l*pow(x,m)/factorial(m);
    l*=(-1);
  }
  float e=fabs(sin(x)-sinx);       //Calculating the absolute error between truncated value and actual value
  printf("%.3f,%.3f\n",sinx,e);    //Truncating to 3 decimals
  return 0;
}  
    
  
