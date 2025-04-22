Part 1:   

Use the Hamming distance as a measure to find the most probable location of the codeword (length M) within a block of binary values (length N).  
  
Use pointers to the arrays of integers  
Use malloc() and free() to allocate and free memory of size M and N integers  
Use srand() and rand() while generating the codeword and block of binary integers (bits).  
The starting location in the block is 1 (not zero).  
Filename: rollnum_hamming.c  
Inputs: M N  
Output: rollnum,  location, hamming_dist. (comma separated)  
  
Part2:  

Fix a code word of size M, and a message of size 2M. This becomes a transmitted packet. At the receiver, you will have a block of N, with the packet at some random location (between 0 and N-3M). Simulate this by inserting the code+msg at a random location within N bits.  Now repeat this process P times. Calculate the efficiency of transmitting a message from sender to receiver (this will be 1 without any noise)  
  
Filename: rollnum_packet.c  
Inputs: M N P  
Outputs: rollnum, M, N, P, efficiency  
