/*

            EE1103 Assignment 12 - Ising Model
            
Developer : EE24B017,EE24B069
Date : 27th October 2024
Purpose :   To simulate the ising model of magnetism and implement the Metropolis-Hastings algorithm to either reach the ferromagnetic or antiferromagnetic case
Input(s) : N J T_norm
Outputs(s): rollnum N J T maxiter totdelE total_m
*/

//Including necessary libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//Defining required macros
#define T_NORM_START 0.5
#define T_INCREMENT 0.1
#define T_NORM_END 4.0
#define MAX_ITERATIONS 100
#define TRANSIENT_ITERATIONS 10
#define no_of_times_to_do_calculations_for_same_t_norm_different_lattices 50
#define acceptable_ratio_of_stddev_by_mean 0.01
#define acceptable_stddev 5
#define probablity_of_flipping_for_delta_E_equal_to_zero 0.25

//Defining cell as a structure
typedef struct Cell
{
    int value;
    struct Cell* up_ptr;
    struct Cell* down_ptr;
    struct Cell* right_ptr;
    struct Cell* left_ptr;
}cell;

//Defining a lattice
typedef struct Lattice
{
    cell* start_cell_ptr;
    int no_of_rows;
    int no_of_columns;
}lattice;


//Initialising all functions
void rows_connector(int no_of_rows, int no_of_columns, cell* previous_row_start_cell_ptr, cell* row_start_cell_ptr, cell* start_cell_ptr, int row_number);
cell* cells_linked_matrix_gen(int no_of_rows, int no_of_columns);
void cells_free(cell *start_cell_ptr, int no_of_rows, int no_of_columns);
void cells_display(cell* start_cell_ptr, int no_of_rows, int no_of_columns);
void cells_ptr_display(cell* start_cell_ptr, int no_of_rows, int no_of_columns);
void lattice_random_intializer(lattice lattice1);
float cell_energy_calculator(cell* cell_ptr, float j_value);
float lattice_energy_calculator(lattice lattice1, float j_value);
void lattice_copy(lattice* lattice_ptr1, lattice* lattice_ptr2);
float lattice_update(lattice* lattice_ptr, lattice* temp_lattice_ptr, float j_value, double t_norm);
void mean_stddev_calculator(int* array_of_five_values, double* mean_ptr, double* stddev_ptr);
double curie_temp_calculator(FILE *fptr,int no_of_points, int no_of_rows_in_lattice);
int magnetisation(lattice lattice1);

void main(int argc, char **argv)
{   //Checking for correct number of arguments
    if (argc!=4){
        printf("Usage ./ee24b017_ee24b069_ising N J T");
        return;
    }
    //Reading the variables
    int N = atoi(argv[1]);
    float j_val=atof(argv[2]);
    float t_nor=atof(argv[3]);

    srand(time(NULL));

    //Initialising the lattice
    lattice latt;
    latt.no_of_rows = N;
    latt.no_of_columns = N;
    //Initialising the 2nd lattice
    lattice temp_latt;
    //Generating the lattices
    latt.start_cell_ptr = cells_linked_matrix_gen(latt.no_of_rows, latt.no_of_columns);
    temp_latt.start_cell_ptr = cells_linked_matrix_gen(latt.no_of_rows, latt.no_of_columns);
    temp_latt.no_of_rows = latt.no_of_rows;
    temp_latt.no_of_columns = latt.no_of_columns;

    //Creating variables to store the energy values for a convergernce criterion
    int array_to_store_previous_5_values_of_total_energy[5];
    double mean_of_previous_5_values_of_total_energy = 0.0;
    double stddev_of_previous_5_values_of_total_energy = 0.0;
    int total_energy;
    double avg_total_energy = 0.0;
    double array_to_store_previous_5_values_of_avg_total_energy[5];
    double mean_of_previous_5_values_of_avg_total_energy = 0.0;
    double stddev_of_previous_5_values_of_avg_total_energy = 0.0;
    int maxiter=0;
    int mag=0;
    float avg_mag = 0.0;
    //We run the loop for a fixed number of times to average out the values and avoid outliers
    for(int j=0; j<no_of_times_to_do_calculations_for_same_t_norm_different_lattices; j++)
    {
        lattice_random_intializer(latt);
        //ignoring the  transient iterations
        for(int i=0; i<TRANSIENT_ITERATIONS; i++)
        {
            //cells_display(latt.start_cell_ptr, latt.no_of_rows, latt.no_of_columns);
            total_energy = lattice_update(&latt, &temp_latt, j_val, t_nor);
        }

        //Running the actual iterations
        for(int i=0;i<300;i++)
        {
            total_energy = lattice_update(&latt, &temp_latt, j_val, t_nor);
            //Defining a convergence criterion
            if(i>4)
            {
                mean_stddev_calculator(array_to_store_previous_5_values_of_total_energy, &mean_of_previous_5_values_of_total_energy, &stddev_of_previous_5_values_of_total_energy);
                //printf("%lf %lf\n", mean_of_previous_5_values_of_total_energy, stddev_of_previous_5_values_of_total_energy);
                if((fabs(stddev_of_previous_5_values_of_total_energy/mean_of_previous_5_values_of_total_energy) < 0.0001) || stddev_of_previous_5_values_of_total_energy<1)
                {
                    //Exiting the loop once the st dev of last 5 entries becomes small enough 
                    maxiter=i;
                    //totdelE = total_energy;
                    //cells_display(latt.start_cell_ptr,latt.no_of_rows,latt.no_of_columns);
                    break;
                }
                mag=magnetisation(latt);
                for(int k=4; k>0 ;k--)
                {
                    //Storing the previous 5 energy values
                    array_to_store_previous_5_values_of_total_energy[k] = array_to_store_previous_5_values_of_total_energy[k-1];
                    //printf("%d\n", array_to_store_previous_5_values_of_total_energy[k]);
                }
                array_to_store_previous_5_values_of_total_energy[0] = total_energy;
            }
            else
            {
                //If array doesn't have all 5 elements
                array_to_store_previous_5_values_of_total_energy[4-i] = total_energy;
            }

        }
        avg_total_energy += mean_of_previous_5_values_of_total_energy;
        avg_mag += (float)mag;
    }
    if (maxiter==0){
        //In case the loop never breaks, it means it has run for MAX_ITERATIONS
        maxiter=MAX_ITERATIONS;
    }
    //printf("%f\n",totdelE);
    avg_total_energy = avg_total_energy/(double)no_of_times_to_do_calculations_for_same_t_norm_different_lattices;
    avg_mag = avg_mag/(float)no_of_times_to_do_calculations_for_same_t_norm_different_lattices;
    printf("ee24b017_ee24b069 %d %.3f %.3f %.3d %.4f %.3f\n",N,j_val,t_nor,maxiter,0.0001,avg_mag);
    //Freeing the lattices
    cells_free(latt.start_cell_ptr, latt.no_of_rows, latt.no_of_columns);
    cells_free(temp_latt.start_cell_ptr, temp_latt.no_of_rows, temp_latt.no_of_columns);


    //The following part of the code is for generating the plot of M vs T norm
    //Initialising the lattice
    lattice lattice1;
    lattice1.no_of_rows = 16;
    lattice1.no_of_columns = 16;
    int j_value = 1;
    //Initialising the 2nd lattice
    lattice temp_lattice;
    //Generating the lattices
    lattice1.start_cell_ptr = cells_linked_matrix_gen(lattice1.no_of_rows, lattice1.no_of_columns);
    temp_lattice.start_cell_ptr = cells_linked_matrix_gen(lattice1.no_of_rows, lattice1.no_of_columns);
    temp_lattice.no_of_rows = lattice1.no_of_rows;
    temp_lattice.no_of_columns = lattice1.no_of_columns;
    
    //Creating variables to store the energy values for a convergernce criterion
    float avg_magnetisation = 0.0;
    array_to_store_previous_5_values_of_total_energy[5];
    mean_of_previous_5_values_of_total_energy = 0.0;
    stddev_of_previous_5_values_of_total_energy = 0.0;
    total_energy;
    avg_total_energy = 0.0;
    array_to_store_previous_5_values_of_avg_total_energy[5];
    mean_of_previous_5_values_of_avg_total_energy = 0.0;
    stddev_of_previous_5_values_of_avg_total_energy = 0.0;
    //Writing data points to this file
    FILE* fileptr = fopen("ee17_ee69_display.txt", "w");
    //cells_ptr_display(lattice1.start_cell_ptr, lattice1.no_of_rows, lattice1.no_of_columns);
    for(double t_norm=T_NORM_START; t_norm<=T_NORM_END; t_norm+=T_INCREMENT)
    {    //We run the loop for a fixed number of times to average out the values and avoid outliers
        for(int j=0; j<no_of_times_to_do_calculations_for_same_t_norm_different_lattices; j++)
        {
            lattice_random_intializer(lattice1);
            //ignoring the  transient iterations
            for(int i=0; i<TRANSIENT_ITERATIONS; i++)
            {
                //cells_display(lattice1.start_cell_ptr, lattice1.no_of_rows, lattice1.no_of_columns);
                total_energy = lattice_update(&lattice1, &temp_lattice, j_value, t_norm);
            }
            //Running the actual iterations
            for(int i=0;i<MAX_ITERATIONS;i++)
            {
                total_energy = lattice_update(&lattice1, &temp_lattice, j_value, t_norm);
                //Defining a convergence criterion
                if(i>4)
                {
                    mean_stddev_calculator(array_to_store_previous_5_values_of_total_energy, &mean_of_previous_5_values_of_total_energy, &stddev_of_previous_5_values_of_total_energy);
                    //printf("%lf %lf\n", mean_of_previous_5_values_of_total_energy, stddev_of_previous_5_values_of_total_energy);
                    if((fabs(stddev_of_previous_5_values_of_total_energy/mean_of_previous_5_values_of_total_energy) < acceptable_ratio_of_stddev_by_mean) || stddev_of_previous_5_values_of_total_energy<acceptable_stddev)
                    {
                        //Exiting the loop once the st dev of last 5 entries becomes small enough 
                        break;
                    }
                    for(int k=4; k>0 ;k--)
                    {
                        //Storing the previous 5 energy values
                        array_to_store_previous_5_values_of_total_energy[k] = array_to_store_previous_5_values_of_total_energy[k-1];
                        //printf("%d\n", array_to_store_previous_5_values_of_total_energy[k]);
                    }
                    array_to_store_previous_5_values_of_total_energy[0] = total_energy;
                }
                else
                {   //If array doesn't have all 5 elements
                    array_to_store_previous_5_values_of_total_energy[4-i] = total_energy;
                }

            }
            avg_total_energy += mean_of_previous_5_values_of_total_energy;
        }
        avg_total_energy = avg_total_energy/(double)no_of_times_to_do_calculations_for_same_t_norm_different_lattices;
        
        //cells_display(lattice1.start_cell_ptr, lattice1.no_of_rows, lattice1.no_of_columns);
        fprintf(fileptr,"%f %lf\n", t_norm, -0.25*avg_total_energy/(4*16*16)+0.5);
    }
    fclose(fileptr);

    //Reopening the data file to find the curie temperature
    FILE* fptr = fopen("ee17_ee69_display.txt", "r");
    double curie_temp = curie_temp_calculator(fptr, (int)ceil((T_NORM_END-T_NORM_START)/T_INCREMENT), lattice1.no_of_rows);
    fclose(fptr);

    //Calling gnuplot inside c to plot the points
    FILE* gplot = popen("gnuplot -persist", "w");
    fprintf(gplot, "set terminal jpeg\n");
    fprintf(gplot, "set output 'ee24b017_ee24b069_ising.jpg'\n");
    fprintf(gplot, "set ylabel 'relative magnetisation(calculated using energies)'\n");
    fprintf(gplot, "set xlabel 'normalised temperature'\n");
    fprintf(gplot, "plot 'ee17_ee69_display.txt' with points title sprintf('Curie Temperature : %.2lf')\n" , curie_temp);
    pclose(gplot);
    remove("ee17_ee69_display.txt");
    cells_free(lattice1.start_cell_ptr, lattice1.no_of_rows, lattice1.no_of_columns);
    cells_free(temp_lattice.start_cell_ptr, temp_lattice.no_of_rows, temp_lattice.no_of_columns);
    return;
}

//Function to actually link all the cells and ensure periodic boundary conditions are implemented
void rows_connector(int no_of_rows, int no_of_columns, cell* previous_row_start_cell_ptr, cell* row_start_cell_ptr, cell* start_cell_ptr, int row_number)
{
    //Initialising various ptrs to connect adjacent rows
    cell* start_row_traverse_ptr = start_cell_ptr;
    cell* previous_row_traverse_ptr = previous_row_start_cell_ptr;
    cell* current_row_traverse_ptr = row_start_cell_ptr;
    for(int i=0; i < no_of_columns; i++)
    {
        //Linking up and down ptrs
        previous_row_traverse_ptr->down_ptr = current_row_traverse_ptr;
        current_row_traverse_ptr->up_ptr = previous_row_traverse_ptr;
        if (row_number == no_of_rows)
        {
            //Linking last row starting ptr with first row starting ptr
            current_row_traverse_ptr->down_ptr = start_row_traverse_ptr;
            start_row_traverse_ptr->up_ptr = current_row_traverse_ptr;
            start_row_traverse_ptr = (start_row_traverse_ptr->right_ptr);
            //This ensures that the next column elements are similarly linked up
        }
        //Shifting to next column
        previous_row_traverse_ptr = (previous_row_traverse_ptr->right_ptr);
        current_row_traverse_ptr = (current_row_traverse_ptr->right_ptr);
    }
    return;
}

//Function to generate a linked matrix and dynamically allocate space accordingly
cell* cells_linked_matrix_gen(int no_of_rows, int no_of_columns)
{
    cell* current_cell_ptr;
    cell* previous_cell_ptr;
    cell* start_cell_ptr;
    cell* previous_row_start_cell_ptr;

    //Going through the rows
    for(int i=0; i<no_of_rows; i++)
    {
        cell* row_start_cell_ptr = (cell*)malloc(sizeof(cell));
        current_cell_ptr = row_start_cell_ptr;
        //Going through each element of a given row
        for(int j=0; j<no_of_columns-1; j++)
        {
            cell* next_cell_ptr = (cell*)malloc(sizeof(cell));
            //Linking the cells horizontally
            current_cell_ptr->right_ptr = next_cell_ptr;
            //For the last column, we need to circularly link it to the first cell
            if(j == (no_of_columns-2))
            {
                row_start_cell_ptr->left_ptr = next_cell_ptr;
                next_cell_ptr->right_ptr = row_start_cell_ptr;
                next_cell_ptr->left_ptr = current_cell_ptr;
            }
            //Link the left ptrs
            if (j != 0)
            {
                current_cell_ptr->left_ptr = previous_cell_ptr;
            }
            //Proceeding to the next cell
            previous_cell_ptr = current_cell_ptr;
            current_cell_ptr = next_cell_ptr;
        }
        //Initialising the first element of each row
        if (i==0)
        {
            start_cell_ptr = row_start_cell_ptr;
        }
        else
        {
            //This is for horizontally linking the matrix
            rows_connector(no_of_rows, no_of_columns, previous_row_start_cell_ptr, row_start_cell_ptr, start_cell_ptr, i+1);
        }
        previous_row_start_cell_ptr = row_start_cell_ptr;
    }
    return start_cell_ptr;
}

//Function to free up the memory allocated to a matrix
void cells_free(cell *start_cell_ptr, int no_of_rows, int no_of_columns)
{
    cell* next_cell_ptr;
    cell* current_cell_ptr = start_cell_ptr->right_ptr;
    cell* row_start_cell_ptr = start_cell_ptr;
    cell* next_row_start_cell_ptr;
    for(int i=0; i<no_of_rows; i++)
    {
        for(int j=0; j<no_of_columns-1; j++)
        {
            next_cell_ptr = current_cell_ptr->right_ptr;
            //Removing the space allocated to each cell individually
            free(current_cell_ptr);
            current_cell_ptr = next_cell_ptr;
        }
        //Freeing the start cell as well
        free(row_start_cell_ptr);
        //Iterating forward, except for the last row
        if(i != no_of_rows-1)
        {
            next_row_start_cell_ptr = (row_start_cell_ptr->down_ptr);
            current_cell_ptr = (next_row_start_cell_ptr->right_ptr);
            row_start_cell_ptr = next_row_start_cell_ptr;
        }
    }
    return;
}
//Function to print the 2D matrix
void cells_display(cell* start_cell_ptr, int no_of_rows, int no_of_columns)
{
    cell* row_start_cell_ptr = start_cell_ptr;
    cell* current_cell_ptr = start_cell_ptr;
    for(int i=0; i<no_of_rows; i++)
    {
        for(int j=0; j<no_of_columns; j++)
        {
            //We use # for value =1 and ' ' for value = -1
            if ((current_cell_ptr->value) == 1)
            {
                printf("# ");
            }
            if ((current_cell_ptr->value) == -1)
            {
                printf("  ");
            }
            //Going forward to the next ptr horizontally
            current_cell_ptr = current_cell_ptr->right_ptr;
        }
        //printing newline as we progress to the next row
        printf("\n");
        current_cell_ptr = row_start_cell_ptr->down_ptr;
        row_start_cell_ptr = current_cell_ptr;
    }
    return;
}

//Function to check if the linking is done correctly
void cells_ptr_display(cell* start_cell_ptr, int no_of_rows, int no_of_columns)
{
    cell* row_start_cell_ptr = start_cell_ptr;
    cell* current_cell_ptr = start_cell_ptr;
    for(int i=0; i<no_of_rows; i++)
    {
        for(int j=0; j<no_of_columns; j++)
        {
            //Printing all neighbours
            printf("( ");
            printf("C:%3.0lu ", ((unsigned long int)current_cell_ptr)%1000);
            printf("R:%3.0lu ", ((unsigned long int)current_cell_ptr->right_ptr)%1000);
            printf("L:%3.0lu ", ((unsigned long int)current_cell_ptr->left_ptr)%1000);
            printf("U:%3.0lu ", ((unsigned long int)current_cell_ptr->up_ptr)%1000);
            printf("D:%3.0lu ", ((unsigned long int)current_cell_ptr->down_ptr)%1000);
            printf(") ");
            current_cell_ptr = current_cell_ptr->right_ptr;
        }
        printf("\n");
        current_cell_ptr = row_start_cell_ptr->down_ptr;
        row_start_cell_ptr = current_cell_ptr;
    }
    return;
}

//Function to assign random values to a pre initialized matrix 
void lattice_random_intializer(lattice lattice1)
{
    cell* row_start_cell_ptr = lattice1.start_cell_ptr;
    cell* current_cell_ptr = lattice1.start_cell_ptr;
    //Going through the rows
    for(int i=0; i<lattice1.no_of_rows; i++)
    {
        //Going through the columns
        for(int j=0; j<lattice1.no_of_columns; j++)
        {
            //Random number generator
            current_cell_ptr->value = (rand()%2)*2-1;
            current_cell_ptr = current_cell_ptr->right_ptr;
        }
        //Iterating to the next row
        current_cell_ptr = row_start_cell_ptr->down_ptr;
        row_start_cell_ptr = current_cell_ptr;
    }
    return;
}

//Function to calculate energy of a individual cell
float cell_energy_calculator(cell* cell_ptr,float j_value)
{
    float energy = ((cell_ptr->up_ptr)->value) + ((cell_ptr->down_ptr)->value) + ((cell_ptr->right_ptr)->value) + ((cell_ptr->left_ptr)->value);
    energy = -1*j_value*energy*(cell_ptr->value);
    return energy;
}

//Function to copy the new lattice(lattice_ptr1) over to the old lattice(lattice_ptr2)
void lattice_copy(lattice* lattice_ptr1, lattice* lattice_ptr2)
{
    //Initialising 2 sets of ptrs for the 2 lattices
    cell* row_start_cell_ptr2 = lattice_ptr2->start_cell_ptr;
    cell* row_start_cell_ptr1 = lattice_ptr1->start_cell_ptr;
    cell* current_cell_ptr2 = lattice_ptr2->start_cell_ptr;
    cell* current_cell_ptr1 = lattice_ptr1->start_cell_ptr;
    //Going through rows
    for(int i=0; i<lattice_ptr2->no_of_rows; i++)
    {
        //Going through columns
        for(int j=0; j<lattice_ptr2->no_of_columns; j++)
        {
            //Copying the values
            current_cell_ptr1->value = current_cell_ptr2->value;
            //Moving to next element to the right
            current_cell_ptr1 = current_cell_ptr1->right_ptr;
            current_cell_ptr2 = current_cell_ptr2->right_ptr;
        }
        //Going to the next row
        current_cell_ptr1 = row_start_cell_ptr1->down_ptr;
        current_cell_ptr2 = row_start_cell_ptr2->down_ptr;
        row_start_cell_ptr1 = current_cell_ptr1;
        row_start_cell_ptr2 = current_cell_ptr2;
    }
}

//Function to calculate the energy of the lattice as a whole
float lattice_energy_calculator(lattice lattice1, float j_value)
{
    cell* row_start_cell_ptr = lattice1.start_cell_ptr;
    cell* current_cell_ptr = lattice1.start_cell_ptr;
    float total_energy = 0;
    //Going through rows
    for(int i=0; i<lattice1.no_of_rows; i++)
    {
        //Going through each element
        for(int j=0; j<lattice1.no_of_columns; j++)
        {
            total_energy += cell_energy_calculator(current_cell_ptr, j_value);
            //Going to the next cell on the right
            current_cell_ptr = current_cell_ptr->right_ptr;
        }
        //Going to the next row
        current_cell_ptr = row_start_cell_ptr->down_ptr;
        row_start_cell_ptr = current_cell_ptr;
    }
    return total_energy;
}

//Function to update the lattice at each iteration
float lattice_update(lattice* lattice_ptr, lattice* temp_lattice_ptr, float j_value, double t_norm)
{
    //Initialising the energy variables
    float total_energy = 0;
    float energy_of_current_cell;
    //Copying the old lattice into temp_lattice_ptr
    lattice_copy(temp_lattice_ptr, lattice_ptr);
    lattice temp_lattice = (*temp_lattice_ptr);
    cell* row_start_cell_ptr_in_temp_lattice = temp_lattice.start_cell_ptr;
    cell* current_cell_ptr_in_temp_lattice = temp_lattice.start_cell_ptr;
    cell* current_cell_ptr_in_original_lattice = lattice_ptr->start_cell_ptr;
    cell* row_start_cell_ptr_in_original_lattice = lattice_ptr->start_cell_ptr;
    //Going through rows
    for(int i=0; i<temp_lattice.no_of_rows; i++)
    {
        //Going through each element
        for(int j=0; j<temp_lattice.no_of_columns; j++)
        {
            //Calculating the energy of the current cell
            energy_of_current_cell = (2)*cell_energy_calculator(current_cell_ptr_in_original_lattice, j_value);
            //If energy difference is greater than 1, do flip
            if (energy_of_current_cell>1)
            {
                current_cell_ptr_in_temp_lattice->value *= -1;
            }
            //If it is zero, flip with probability of 0.2
            else if(energy_of_current_cell == 0)
            {
                if ((((double)rand()/((double)RAND_MAX)) < probablity_of_flipping_for_delta_E_equal_to_zero))
                {
                    current_cell_ptr_in_temp_lattice->value *= -1;
                }
            }
            //If it is negative, flip with probability of exp(-energy/Tnorm)
            else if(((double)rand()/(double)RAND_MAX) < exp((((double)energy_of_current_cell)/t_norm)))
            {
               current_cell_ptr_in_temp_lattice->value *= -1;
            }
            //Updating the total energy
            total_energy += energy_of_current_cell;
            current_cell_ptr_in_temp_lattice = current_cell_ptr_in_temp_lattice->right_ptr;
            current_cell_ptr_in_original_lattice = current_cell_ptr_in_original_lattice->right_ptr;
        }
        current_cell_ptr_in_temp_lattice = row_start_cell_ptr_in_temp_lattice->down_ptr;
        current_cell_ptr_in_original_lattice = row_start_cell_ptr_in_original_lattice->down_ptr;
        row_start_cell_ptr_in_temp_lattice = current_cell_ptr_in_temp_lattice;
        row_start_cell_ptr_in_original_lattice = current_cell_ptr_in_original_lattice;
    }
    //total_energy=lattice_energy_calculator(*lattice_ptr,j_value);
    lattice_copy(lattice_ptr, temp_lattice_ptr);
    return total_energy;
}

//Function to calculate the mean and standard deviation of a given array
void mean_stddev_calculator(int* array_of_five_values, double* mean_ptr, double* stddev_ptr)
{
    double sum_of_values = 0.0;
    double sum_of_squares_of_values = 0.0;
    for(int i=0; i<5; i++)
    {
        sum_of_values += (double)*(array_of_five_values+i);
        sum_of_squares_of_values += (double)(*(array_of_five_values+i))*(*(array_of_five_values+i));
    }
    *(mean_ptr) = sum_of_values/5.0;
    *(stddev_ptr) = sqrt(((sum_of_squares_of_values)/5.0) - (*mean_ptr)*(*mean_ptr));
}

//Function to calculate the curie temperature
double curie_temp_calculator(FILE* fptr, int no_of_points, int no_of_rows_in_lattice)
{
    double tnorm_values[no_of_points];
    double relative_magnetisation[no_of_points];
    int index_of_rel_magn_min = 0;
    for(int i=0; i<no_of_points; i++)
    {
        fscanf(fptr, "%lf %lf\n", tnorm_values+i, relative_magnetisation+i);
        //printf("%lf %lf\n", tnorm_values[i], relative_magnetisation[i]);
        if(*(relative_magnetisation+i)<*(relative_magnetisation+index_of_rel_magn_min))
        {
            index_of_rel_magn_min = i;
        }
    }
    //Finding the minimum value of the energy and using a certain threshold to get to the critical value
    double check = *(relative_magnetisation+index_of_rel_magn_min) + 0.05;
    for(int i=0; i<no_of_points; i++)
    {
        if(check > *(relative_magnetisation+i))
        {
            return (*(tnorm_values+i));
        }
    }
}


int magnetisation(lattice lattice1)
{
    cell* row_start_cell_ptr = lattice1.start_cell_ptr;
    cell* current_cell_ptr = lattice1.start_cell_ptr;
    float total_m = 0;
    //Going through rows
    for(int i=0; i<lattice1.no_of_rows; i++)
    {
        //Going through each element
        for(int j=0; j<lattice1.no_of_columns; j++)
        {
            total_m += current_cell_ptr->value;
            //Going to the next cell on the right
            current_cell_ptr = current_cell_ptr->right_ptr;
        }
        //Going to the next row
        current_cell_ptr = row_start_cell_ptr->down_ptr;
        row_start_cell_ptr = current_cell_ptr;
    }
    return abs(total_m);
}

//End of program