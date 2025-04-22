/*

            EE1103 Assignment - 5 : Random Bit Generation
            
Developer : EE24B017
Date : 3rd September 2024
Purpose : To find the most probable location of the codeword (length M) within a block of binary values (length N) using the Hamming distance.
Input(s) : M N
Outputs(s): rollnum location hamming_dist
*/

//Including the necessary header files
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

int main( int argc, char **argv){
  if (argc!=3){             //Checking if the correct number of arguments are entered
    printf("Please input correct M and N values\n");
    return 0;
  }
  
  int M = atoi(argv[1]);     //Taking input for the length
  int N = atoi(argv[2]);
  
  if (M>N) {                //Checking if M > N
    printf("M should be less than or equal to N\n");
    return 0;
  }
  
  int *block = (int *)malloc(N*sizeof(int));    //Allocating memory for the block and codeword
  int *code = (int *)malloc(M*sizeof(int));
  
  
  srand(time(NULL));               //Setting seed for random number generator
  
  
  int i,j,k;
  
  for (i=0;i<N;i++)      //Generating block
    *(block+i)=rand()%2;
    
  for (j=0;j<M;j++)      //Generating codeword
    *(code+j)=rand()%2;
    
  int store = M;         //Initialising minimum hamming distance as pos with value as M
  int pos;
  
  for (k=0;k <= N-M;k++){  
  
    int ham = 0;        //Initialise ham as the measure of hamming distance of current selection
    int n = 0;
    
    while (n < M){
      if (*(block+k+n) != *(code+n))
        ham += 1;        //Calculating hamming distance 
      n++;
    }
    
    if (ham < store){
      pos = k;            //If the minimum hamming distance is greater than ham, we reassign minimum hamming distance
      store = ham;
    }
    
    if (ham==0){          //If the hamming distance is zero, we can terminate the program
      printf("ee24b017,%d,%d\n",pos+1,ham);
      free(block);
      free(code);
      return 0;
    }
  }
  
  printf("ee24b017,%d,%d\n",pos+1,store); //Printing the output
  free(code);
  free(block);
}
