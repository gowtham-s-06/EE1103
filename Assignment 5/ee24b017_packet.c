/*
  EE1103 Assignment - 4 : Packet Insertion
  Developer : EE24B017
  Date : 4th September 2024
  Purpose : To simulate insertion of a codeword in a block of binary values and compute the efficiency of transmission.
  Input(s) : M N P
  Output(s): rollnum, M, N, P, efficiency
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv) {
    if (argc != 4) {  // Checking if the correct number of arguments are entered
        printf("Please input correct M,N and P values\n");
        return 0;
    }

    int M = atoi(argv[1]);     // Length of the codeword
    int N = atoi(argv[2]);     // Length of the block
    int P = atoi(argv[3]);     // Number of codewords

    if (M > N) {  // Checking if M > N
        printf("M should be less than or equal to N\n");
        return 0;
    }

    int *block = (int *)malloc(N * sizeof(int));    // Allocating memory for the block and codeword
    int *code = (int *)malloc(M * sizeof(int));

    srand(1);                                       // Setting seed for random number generator

                                                    // Generate random block
    for (int i = 0; i < N; i++) {
        block[i] = rand() % 2;
    }

                                            // Generate random codeword
    for (int j = 0; j < M; j++)
        code[j] = rand() % 2;
  
   
    int length = N / P;  // Size of each sub-block
    
    // Insert codeword into random positions within each sub-block
    for (int i = 0; i < P; i++) {
        int start = i * length;                      // Start position of the current sub-block
        int max_position = length - 3 * M;           // Ensure that there is at least 2M space from the end of the sub-block
        int position = rand() % (max_position + 1);  // Random position within the allowed range
        int insert_pos = start + position;           // Calculate the actual position within the block

        // Insert the codeword at the calculated position
        for (int j = 0; j < M; j++) {
            block[insert_pos + j] = code[j];
        }
    }

    int counter = 0;                                 // Count the occurrences of the codeword in the block
    for (int k = 0; k <= N - M; k++) {
        int ham = 0;  // Measure of Hamming distance
        int n = 0;

        while (n < M) {     // Calculate Hamming distance
            if (block[k + n] != code[n])
                ham += 1;
            n++;
        }

        if (ham == 0) {  // If Hamming distance is zero, we found a codeword
            counter++;
        }
    }

    float efficiency = counter / P ;
    printf("ee24b017,%d,%d,%d,%.2f\n", M, N, P, efficiency); // Printing the output


    free(code);
    free(block); 

    return 0;
}

