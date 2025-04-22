Read data from any datafile in this folder (save the timestamps in a txt or csv file).  The first column of the datafile has the timestamps at which a single photon detector registers a click. We are expecting clicks in a window of 3 ns, but the centre of the window is unknown.  The expt is repeated every 32 ns  
  
Use module-32 to determine the timestamps within a 32ns window  
Identify the centre of a 3ns window where the total counts is a maximum  
Identify 3 time bins of 1ns each within the 3ns window and estimate the counts C1, D1 and C2  
Estimate the BER (across the full data file) as BER1=D1/(C1+D1+C2)  
Estimate the visibility as V1=1-D1/(C1+C2)  
Now, assume that we have a 100ps guard bands between two consequtive time bins, where the timestamps are likely to be erroneous. Throw away the timestamps within the guard bands, and estimate BER2 and Visibility2  
  
Filename: group_ber.c  
Input: datafile  
Output: group, BER1, Visibility1, BER2, Visibility2  
When uploading your submission, use the comments box to comment on the utility of the guard band, and whether 100ps was a good enough choice or if we should increase/decrease it to maximize V and minimize BER.  

