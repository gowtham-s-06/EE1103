/*

            EE1103 Assignment - 11 : Damped Gyromagnetic Switching
            
Developer : EE24B017
Date : 17th October 2024
Purpose : To solve the Landau Lifshitz Gilbert Equation in spherical and cartesian coordinate systems using Heun's Method and hence plot the trajectory and rms error between the two methods.
Input(s) : theta_start theta_stop alpha delta_t
Outputs(s): alpha R^2 switch_time_ns

WARNING : Use time step of less than 5e-12 or else the program will overshoot the target.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846
#define GAMMA 1.76e11   // Gyromagnetic ratio in rad/(s·T)

// Define the parameters of the system
double H = 1.0;         // Applied field in Tesla
double Hk = 0.0;        // Anisotropy field in Tesla
double alpha;           // Damping constant (will be set from argv)

// Function to compute d(theta)/dt
double dtheta_dt(double theta) {
    return (GAMMA * alpha * (H * sin(theta) - Hk * sin(theta) * cos(theta))) / (1 + alpha * alpha);
}

// Function to compute d(phi)/dt
double dphi_dt(double theta, double dtheta) {
    if (sin(theta) == 0) {
        return 0.0;  // Avoid division by zero when theta is 0 or pi
    }
    return -dtheta / (alpha * sin(theta));
}

// Convert from spherical to Cartesian coordinates
void spherical_to_cartesian(double theta, double phi, double *Mx, double *My, double *Mz) {
    *Mx = sin(theta) * cos(phi);
    *My = sin(theta) * sin(phi);
    *Mz = cos(theta);
}

// Heun's method to estimate theta(t) and phi(t)
void heun_step(double *theta, double *phi, double delta_t) {
    // Predict theta using Heun's method
    double k1_theta = dtheta_dt(*theta);                        // Slope at the beginning of the interval for theta
    double theta_predict = *theta + delta_t * k1_theta;         // Predictor for theta
    double k2_theta = dtheta_dt(theta_predict);                 // Slope at the end of the interval for theta
    *theta = *theta + delta_t * 0.5 * (k1_theta + k2_theta);    // Heun's corrected theta

    // Predict phi using Heun's method based on theta
    double k1_phi = dphi_dt(*theta, k1_theta);                  // Slope at the beginning of the interval for phi
    double phi_predict = *phi + delta_t * k1_phi;               // Predictor for phi
    double k2_phi = dphi_dt(theta_predict, k2_theta);           // Slope at the end of the interval for phi
    *phi = *phi + delta_t * 0.5 * (k1_phi + k2_phi);            // Heun's corrected phi
}

// Function to compute the cross product of two vectors
void cross_product(double v1[3], double v2[3], double result[3]) {
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

// Function to normalize a vector
void normalize(double v[3]) {
    double magnitude = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] /= magnitude;
    v[1] /= magnitude;
    v[2] /= magnitude;
}

// Function to compute dM/dt using the LLG equation
void dM_dt(double M[3], double dM[3]) {
    double H_eff[3] = {0.0, 0.0, -H}; // Effective field along -z
    double M_cross_H[3], M_cross_M_cross_H[3];
    
    // Compute M x H
    cross_product(M, H_eff, M_cross_H);
    
    // Compute M x (M x H)
    cross_product(M, M_cross_H, M_cross_M_cross_H);
    
    // Compute dM/dt
    for (int i = 0; i < 3; i++) {
        dM[i] = -GAMMA * M_cross_H[i] - alpha * GAMMA * M_cross_M_cross_H[i];
    }
}

// Heun's method to estimate Mx(t), My(t), and Mz(t)
void heun_step_Cartesian(double M[3], double delta_t) {
    double k1[3], k2[3], M_predict[3];

    // Slope at the beginning of the interval
    dM_dt(M, k1);

    // Predictor step: M_predict = M + delta_t * k1
    for (int i = 0; i < 3; i++) {
        M_predict[i] = M[i] + delta_t * k1[i];
    }

    // Slope at the end of the interval
    dM_dt(M_predict, k2);

    // Corrector step: M_new = M + delta_t * 0.5 * (k1 + k2)
    for (int i = 0; i < 3; i++) {
        M[i] = M[i] + delta_t * 0.5 * (k1[i] + k2[i]);
    }

    // Normalize M to stay on the unit sphere
    normalize(M);
}
//Function to plot the trajectory obtained using spherical coordinates
void plot_trajectory() {
    FILE *gnuplot = popen("gnuplot -persistent", "w");
    if (gnuplot == NULL) {
        printf("Error: Could not open Gnuplot.\n");
        return;
    }
    char *img="ee24b017_trajectory.jpg";
    // Gnuplot commands to plot the trajectory on the unit sphere
    fprintf(gnuplot, "set terminal jpeg size 1000,1000\n");
    fprintf(gnuplot, "set output '%s'\n",img);
    fprintf(gnuplot, "set xrange [-1:1]\n");
    fprintf(gnuplot, "set yrange [-1:1]\n");
    fprintf(gnuplot, "set zrange [-1:1]\n");
    fprintf(gnuplot, "set xlabel 'Mx(X Axis)'\n");
    fprintf(gnuplot, "set ylabel 'My(Y Axis)'\n");
    fprintf(gnuplot, "set zlabel 'Mz(Z Axis)'\n");
    fprintf(gnuplot, "set title 'Magnetization Trajectory on Unit Sphere'\n");
    fprintf(gnuplot, "splot 'trajectory.dat' using 1:2:3 with lines title 'Trajectory'\n");

    fflush(gnuplot);
    pclose(gnuplot);
}

//Function to calculate euclidean distance between two points. This function is used in calculation og RMS error.
double distance_sq(double x1, double y1, double z1, double x2, double y2, double z2){
    double d1 = pow((x1-x2),2);
    double d2 = pow((y1-y2),2);
    double d3 = pow((z1-z2),2);
    return d1+d2+d3;
}


//Following functions are meant for the optional part 2 to plot switching time vs alpha
// Function to compute d(theta)/dt for different values of alpha
double dtheta_dt_2(double theta, double alpha) {
    return (GAMMA * alpha * (H * sin(theta) - Hk * sin(theta) * cos(theta))) / (1 + alpha * alpha);
}

// Heun's method to estimate theta(t) alone
void heun_step_2(double *theta, double delta_t, double alpha) {
    // Predict theta using Heun's method
    double k1_theta = dtheta_dt_2(*theta, alpha);             // Slope at the beginning of the interval for theta
    double theta_predict = *theta + delta_t * k1_theta;     // Predictor for theta
    double k2_theta = dtheta_dt_2(theta_predict, alpha);      // Slope at the end of the interval for theta
    *theta = *theta + delta_t * 0.5 * (k1_theta + k2_theta); // Heun's corrected theta
}
// Function to calculate switching time for a given alpha
double calculate_switching_time(double theta_start, double theta_stop, double delta_t, double alpha) {
    double theta = theta_start;
    double t = 0.0;  // Start time

    // Time integration loop until theta exceeds theta_stop
    while (theta <= theta_stop) {
        heun_step_2(&theta, delta_t, alpha);
        t += delta_t;
    }

    // Return total time in nanoseconds
    return t * 1e9;
}

// Function to generate data for plotting switching time vs. alpha
void plot_switching_time_vs_alpha(double theta_start_deg, double theta_stop_deg, double delta_t) {
    FILE *data_file = fopen("switching_time_vs_alpha.dat", "w");
    if (data_file == NULL) {
        printf("Error: Could not open file for writing.\n");
        return;
    }

    double theta_start = theta_start_deg * PI / 180.0;  // Convert degrees to radians
    double theta_stop = theta_stop_deg * PI / 180.0;    // Convert degrees to radians

    // Vary alpha between 0.01 and 0.2
    for (double alpha = 0.01; alpha <= 0.2; alpha += 0.01) {
        // Calculate switching time for each alpha
        double switching_time_ns = calculate_switching_time(theta_start, theta_stop, delta_t, alpha);

        // Write alpha and switching time to file
        fprintf(data_file, "%.6f %.6f\n", alpha, switching_time_ns);
    }

    fclose(data_file);

    // Plot the data using Gnuplot
    FILE *gnuplot = popen("gnuplot -persistent", "w");
    if (gnuplot == NULL) {
        printf("Error: Could not open Gnuplot.\n");
        return;
    }

    char *img="ee24b017_switching_time_vs_alpha.jpg";
    fprintf(gnuplot, "set terminal jpeg size 1000,1000\n");
    fprintf(gnuplot, "set output '%s'\n",img);
    fprintf(gnuplot, "set xlabel 'Alpha (Damping Constant)'\n");
    fprintf(gnuplot, "set ylabel 'Switching Time (ns)'\n");
    fprintf(gnuplot, "set title 'Switching Time vs. Alpha'\n");
    fprintf(gnuplot, "plot 'switching_time_vs_alpha.dat' using 1:2 with linespoints title 'Switching Time'\n");

    fflush(gnuplot);
    pclose(gnuplot);
}


int main(int argc, char *argv[]) {
    // Ensure proper usage
    if (argc != 5) {
        printf("Usage: %s theta_start theta_stop alpha delta_t\n", argv[0]);
        return 1;
    }

    // Convert input alpha from command-line argument
    alpha = atof(argv[3]);

    // Initial conditions
    double theta_start = atof(argv[1]) * PI / 180; // Initial angle theta 
    double theta_stop = atof(argv[2]) * PI /180 ;    // Stopping condition when theta exceeds this value
    double delta_t = atof(argv[4]);          // Time step in seconds
    double phi_start = 0.0;          // Initial angle phi (assumed to be 0)
    if (delta_t>5e-12){
        printf("Enter a time step of less than 5e-12");
    }
    if (theta_start>theta_stop){
        printf("theta_start should be less than theta_stop");
    }
    double theta = theta_start;
    double phi = phi_start;

    // Time parameters
    double t = 0.0; // Initial time
    // Variables to store the trajectory in Cartesian coordinates
    double Mx, My, Mz;

    // Open file to save the trajectory data
    FILE *data_file = fopen("trajectory.dat", "w");
    if (data_file == NULL) {
        printf("Error: Could not open trajectory.dat for writing.\n");
        return 1;
    }
    // Print the initial conditions
    //printf("Using alpha = %.6f\n", alpha);
    //printf("t (s)\t\ttheta (rad)\tphi (rad)\n");
    //printf("%.3e\t%.6f\t%.6f\n", t, theta, phi);

    // Time integration loop
    while (theta <= theta_stop) {

        // Convert the current theta and phi to Cartesian coordinates
        spherical_to_cartesian(theta, phi, &Mx, &My, &Mz);

        // Save the current coordinates to the data file
        fprintf(data_file, "%.6f %.6f %.6f\n", Mx, My, Mz);

        // Perform a Heun step to update theta and phi
        heun_step(&theta, &phi, delta_t);

        // Increment time
        t += delta_t;

        // Output the current time, theta, and phi values
        //printf("%.3e\t%.6f\t%.6f\n", t, theta, phi);
    }

    //printf("Simulation stopped. Theta reached %.6f radians.\n", theta);
    fclose(data_file);

    // Call Gnuplot to plot the trajectory
    plot_trajectory();

    //Below part of the program repeats the exercise using Cartesian coordinates
    // Initialize magnetization vector M on the unit sphere (polar angle theta_start)
    double M[3] = {sin(theta_start), 0.0, cos(theta_start)};

    // Variables to store the theta angle in polar and Cartesian methods
    double theta_cartesian=theta_start;
    theta = theta_start;
    phi = phi_start;
    int i=0;
    // Print initial conditions
    //printf("Using alpha = %.6f\n", alpha);
    //printf("t (s)\t\tMx\t\tMy\t\tMz\t\tTheta_Cartesian (rad)\tR² Difference\n");
    double total_time;
    double sum_sq;
    
    double x,y,z,dist;
    // Time integration loop
    while (theta_cartesian <= theta_stop) {
        // Perform a Heun step in Cartesian coordinates
        heun_step_Cartesian(M, delta_t);

        // Compute theta in Cartesian coordinates: theta = acos(Mz)
        theta_cartesian = acos(M[2]);
        heun_step(&theta, &phi, delta_t);

        spherical_to_cartesian(theta, phi, &x, &y, &z);
        dist=distance_sq(M[0],M[1],M[2],x,y,z);
        sum_sq+=dist;
        // Print the time, Cartesian coordinates, and R² difference
        //printf("%.3e\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", total_time, M[0], M[1], M[2], theta_cartesian, R_squared);

        // Increment total time
        total_time += delta_t;
        i+=1;
    }

    // Convert total time to nanoseconds for output
    total_time *= 1e9;
    t *= 1e9;
;   double avg = (total_time+t)/2.0;
    if (i!=0){
        double rms = pow(sum_sq/i,0.5);
        printf("%f %f %f\n",alpha,rms,avg);
    }
    else{
        printf("Theta start is less than theta stop");
    }
    // Call the function to plot switching time vs. alpha
    plot_switching_time_vs_alpha(theta_start, theta_stop, delta_t);
    remove("trajectory.dat");
    remove("switching_time_vs_alpha.dat");
    return 0;

}
