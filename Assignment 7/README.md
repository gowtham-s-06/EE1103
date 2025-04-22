Generate normal distributions for noise eps_T, eps_a, n  
Use a randn() or similar function or write a function to implement the Box - Muller transform   
Generate a time series of Lorentzians  s(t) separated by T with width proportional to a  
 (a^2/((t-m*T)^2+a^2), where m = 1 to M  
Define a structure for each peak, with its location, amplitude and width as attributes.  
Add "normal" noise to the location and the width of each Lorentzian  
a -> a + eps_a  
mT -> mT + eps_T  
Add "normal" noise n(t) to the amplitude   
s(t) + n(t)  
Identify the locations of the peaks and estimate  
The average time <T> and the average width <a> and their standard deviations  
Compare it to the values of T and a that you used  
After your code has been tested, read from a datafile and process it  
Program: groupname_timeseries.c  
Inputs:  Check the number of inputs and take as input  either of these  
M, T, a  as argv inputs to simulate a data file  
any input data file  from this folder  
Outputs: avg(T) avg(a) stdev(T) stdev(a)  
