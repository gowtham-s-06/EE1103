/*
            EE1103 Assignment - 8 : Estimating BER          
Developer : Group G (EE24B017,EE24B058,EE24B068) 
Date : 25th September 2024
Purpose : 1) To read a csv file and then use modulo 32 to get the timestamps and frequencies.
          2) To find the 3 ns windows where total counts is maximum and find c1,d1,c2.
          3) To find Bit Error Rate BER1 = D1/(C1+D1+C2) and Visibility V1 = 1-D1/(C1+C2)
          4) To implement 100 ps guardbands and recalculate BER2,V2.
Input(s) : <filename>.csv
Outputs(s) : group, BER1, Visibility1, BER2, Visibility2
Used ChatGPT for help
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAX_LINE_LENGTH 10000
//Defining maximum length of line to read from txt file

//Function to calculate BER and Visibility without using guard bands
float* calculate_BER_Vis(int ws, int we, float * arr){
  FILE *f = fopen("output.txt","r");

  int bin1 = (ws + 1000);      //Defining bins
  int bin2 = (we - 1000);

  int c1=0;                    //Initializing counters of clicks
  int c2=0;
  int d1=0;

  int value=0;
  int freq=0;
  char line[100];

  while (fgets(line, sizeof(line), f)) {
        // Tokenize the line to get two columns
        char *token1 = strtok(line, " ");
        char *token2 = strtok(NULL, " ");

        if (token1 != NULL && token2 != NULL ) {
            value = atoi(token1);
            freq = atoi(token2);
        }

        if (ws <= value && value < bin1)   //Counting frequencies in each bin
          c1+=freq;
        else if (bin1 <= value && value < bin2)
          d1+=freq;
        else if (bin2 <= value && value <= we)
          c2+=freq;
    }

  //Calculating BER
  float BER = (float) d1/(d1+c1+c2);
  float Vis = 1 - (float) d1/(c1+c2);
  arr[0] = BER;
  arr[1] = Vis;
  fclose(f);
  return arr;

}


//Function to calculate BER and Visibility with guard bands in place
float* calculate_BER_Vis_with_guards(int ws, int we, float * arr){
  FILE *f = fopen("output.txt","r");

  int bin1 = (ws + 950);      //Defining bins with guard bands
  int bin2 = (we - 950);

  int c1=0;                   //Initializing counters of clicks
  int c2=0;
  int d1=0;

  int value=0;
  int freq=0;
  char line[100];

  while (fgets(line, sizeof(line), f)) {
        // Tokenize the line to get two columns
        char *token1 = strtok(line, " ");
        char *token2 = strtok(NULL, " ");

        if (token1 != NULL && token2 != NULL ) {
            value = atoi(token1);
            freq = atoi(token2);
        }

        if (ws <= value && value < bin1)      //Counting frequencies in each bin
          c1+=freq;
        else if ((bin1 + 100) <= value && value < (bin2-100))
          d1+=freq;
        else if (bin2 <= value && value <= we)
          c2+=freq;
    }

  float BER = (float) d1/(d1+c1+c2);          //Calculating BER and Visibility
  float Vis = 1 - (float) d1/(c1+c2);


  arr[0] = BER;
  arr[1] = Vis;
  fclose(f);
  return arr;

}


int compare(const void *a, const void *b) {
    return (*(int *)a - *(int *)b);
}


int main(int argc, char *argv[]) 
{
    if (argc != 2)
    {
        printf("Check the number of inputs");
    }
    else
    {
    FILE *input_file = fopen(argv[1], "r");
    FILE *output_file = fopen("output.txt", "w+");

    int l=0;
    char line[MAX_LINE_LENGTH];

    while (fgets(line, sizeof(line), input_file)) 
        {
        // Extract the first column
            char *token = strtok(line, ",");
            if (token != NULL) 
              {
                l++;
              }
         }
         
    int arr[l];
    int arr1[32000];
    int arr2[32000];
    int a=0;
    int A;
    int b=0;
    int w=3000;
    int S[32000];


    rewind(input_file);

    while (fgets(line, sizeof(line), input_file)) 
    {
      // Extract the first column
      char *token = strtok(line, ",");
      if (token != NULL) 
        {
      // Convert token to integer and apply modulus
          // Convert to int and apply modulus
          arr[a]=(atoi(token) % 32000);
        }
        a++;
    }

    // Sort the array
    qsort(arr, l, sizeof(int), compare);

    for(int i=0;i<a;i++)
      {             //Creating a file with frequencies of each timestamp
        int c=0;
        for(int j=0;j<a;j++)
          {
          if(arr[i]==arr[j])
              {
                c++;
              } 
          }
        if (i == 0 || arr[i] != arr[i - 1]) 
        {
          fprintf(output_file,"%d %d\n",arr[i],c);
          if(i!=(a-1))
          {
            for( int h=(arr[i]+1);h<arr[i+1];h++)
            {
                fprintf(output_file,"%d %d\n",h,0);
            }
          }
          if(i==(a-1))
          {
            for(int h=arr[i]+1;h<32000;h++)
            {
                fprintf(output_file,"%d %d\n",h,0);
            }
          }
            
        }  

      }
      
    rewind(output_file);
    while (fgets(line, sizeof(line), output_file)) {
        // Tokenize the line to get two columns
        char *token1 = strtok(line, " ");
        char *token2 = strtok(NULL, " ");
        if (token1 != NULL && token2 != NULL ) {
            arr1[b] = atoi(token1);
            arr2[b] = atoi(token2);
            //printf("%d %d\n",arr1[b],arr2[b]);
            b++;
        }
    }

    while(w<=32000)
        {
          int s=0;
          for(int u=(w-3000);u<w;u++)
            {
              s=s+arr2[u];
            }
          S[w-3000]=s;
          w++;
        }

    int max_sum=S[0];
    for(int i=0;i<29001;i++)
    {
      if(S[i]>max_sum)
        {
          max_sum=S[i];
          A=i;
        }
    }
    //printf("The window is:%d to %d\n",A,A+3000); 
    fclose(input_file);
    fclose(output_file);
    float *output = (float *)malloc(2 * sizeof(float));
    
    output = calculate_BER_Vis(A, A+3000,output);

    float *output1 = (float *)malloc(2 * sizeof(float));
    output1 = calculate_BER_Vis_with_guards(A, A+3000,output1);

    printf("Group G, %.3f, %.3f, %.3f, %.3f\n",output[0],output[1],output1[0],output1[1]);
    free(output);
    free(output1);
    remove("output.txt");
    return 0;
  } 
}

