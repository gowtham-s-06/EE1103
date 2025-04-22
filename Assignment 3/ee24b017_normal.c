/*

            EE1103 Assignment - 3 : Normal Distributions
            
Developer : EE24B017
Date : 23rd August 2024
Purpose : To generate a set of random numbers from a normal distribution and confirm
the number of values that lie withing 1,2 and 3 standard deviations of the mean.
Input(s) : N (the number of terms) , x (the mean of the distribution) , y (the standard deviation of the distribution)
Outputs(s): n1 (the number of values within 1 standard deviation of mean), n2 (the number of values within 2 standard deviations of mean), n3 (the number of values within 3 standard deviations of mean)
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Please input 3 values N, mean, standard deviation\n"); // Check if the correct number of arguments is provided
        return 0;
    }
    int N = atoi(argv[1]);
    double mean = atof(argv[2]);
    double stdev = atof(argv[3]); //Reading the inputs
    if (N >= 10001) {
        printf("N should be less than 10001\n"); // Making sure N is less than 10001
        return 0;
    }
    double values[N];
    int within_1_stdev = 0, within_2_stdev = 0, within_3_stdev = 0;
    
    srand((unsigned int)time(NULL));  // Seed the random number generator

    for (int i = 0; i < N; i += 2) {
        double u1 = rand() / (RAND_MAX + 1.0);  // Generate N random numbers using Box-Muller Transformation
        double u2 = rand() / (RAND_MAX + 1.0);
        double z0 = sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);  // Box-Muller transformation
        double z1 = sqrt(-2.0 * log(u1)) * sin(2 * M_PI * u2);
        values[i] = z0 * stdev + mean;  // Scale and shift to get the desired mean and stdev
        if (i + 1 < N) {
            values[i + 1] = z1 * stdev + mean;
        }
    }
    // Calculate the number of values within 1, 2, and 3 standard deviations of the mean
    for (int i = 0; i < N; i++) {
        double d = fabs(values[i] - mean);
        if (d <= stdev) within_1_stdev++;
        if (d <= 2 * stdev) within_2_stdev++;
        if (d <= 3 * stdev) within_3_stdev++;
    }
    printf("%d,%d,%d\n", within_1_stdev , within_2_stdev , within_3_stdev);
    return 0;
}

