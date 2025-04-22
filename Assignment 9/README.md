Convert the pseudo code to c-code for both Newton and Lagrangian interpolation  
Prob 18.26: f(x) = 1/(1+25 x^2)  
Generate the 4th order Lagrange polynomial using 5 equidistant points for x in the interval [-1,1]  
Find the value of f(x0) using a 4th order Newton interpolating polynomial. Assume x0=0.8 to test your code  
Generalize the above example to obtain the interpolating functions n(x) and l(x), using the Newton and Lagrangian methods and estimate the value from each at x=x0.  
  
Filename: rollnum_interpolation.c  
Command line inputs: polynomial_order xstart xend x0  
Output: print the values f(x0) n(x0) and l(x0)  
