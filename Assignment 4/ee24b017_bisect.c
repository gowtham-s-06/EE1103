/*

            EE1103 Assignment - 4a : Bracketing Methods - Bisection Method
            
Developer : EE24B017
Date : 29th August 2024
Purpose : To find the root of the given function f(x) = -25 + 82x - 90 x^2 + 44 x^3 - 8 x^4 + 0.7 x^5 using the bisection method.
Input(s) : xl (lower bound) , xu (upper bound) , epsilon (the error percentage)
Outputs(s): ans ( the value of the root as calculated by the bisection method)
*/

//Including the necessary header files

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


float f(float x) {             //Declaring a function f to calculate the value of f(x) = -25 + 82x - 90 x^2 + 44 x^3 - 8 x^4 + 0.7 x^5
  return -25 + 82 * x - 90 * pow(x, 2) + 44 * pow(x, 3) - 8 * pow(x, 4) + 0.7 * pow(x, 5);
}


float bisect(float xl,float xu, float epsilon) {            //Declaring a function bisect to perform the bisection method

  float midold,mid,check,error=100.0;          //Declaring the required variables and setting error to 100%
  int i = 1;
  
  if (f(xl)*f(xu) > 0) {                                  //Checking if the root lies in the current bound between xl and xu
    printf("The root does not lie in the current bound between xl and xu. Please try again using a different bound\n");
    exit(1);
  }
  
  mid = (xl+xu)/2.0;      //Initializing mid value
  check = epsilon/100.0;   //Initializing a check value against which f(mid) will be compared to check for equality with 0
  
  
  while (error > epsilon) {
  
    if (i > 1) {                    //Creating a if condition to update the error value for each iteration
      midold = mid;
      mid = (xl+xu)/2.0;
      error = fabs(midold-mid)*100/midold;
    }
    
    if (fabs(f(mid)) < check)           //Checking if the f(mid) value is close to zero within prescribed limits
      return mid;
      
    else if (f(xl)*f(mid) < 0)          //Checking if root lies between xl and mid
      xu=mid;
      
    else                                //else means root lies between mid and xu
      xl=mid;
      
    i++;                                //Incrementing the iteration count variable
  }
  float ans = mid;
  printf("ee24b017, %.2f, %d\n",ans,i); //Prints the root value that is within the prescribed error limits
  return 0.0;                          
}

int main(int argc, char **argv) {
  if (argc != 4) {                    //Checking if correct number of arguments are entered
    printf("Please enter lower bound (xl), upper bound (xu) and error percentage (epsilon)");
    return 0;
  }
  
  float xl = atof(argv[1]);           //Getting the inputs from argv
  float xu = atof(argv[2]);
  float epsilon = atof(argv[3]);
  
  bisect(xl,xu,epsilon);           //Calling the bisect function
                                  //Printing the outuput inside the function itself
  
  return 0;
}
