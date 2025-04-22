/*
  EE1103 Assignment - 6 : Generating Fractals
  Developer : EE24B017
  Date : 16th September 2024
  Purpose : To generate a Mandelbrot set zoom sequence, with black points corresponding to numbers that are outside the set, using structures to 
            represent Complex Numbers
  Input(s) : xmin xmax ymin ymax maxiter xres <out.jpg>
  Output(s): A jpg file with the name as <out.jpg>
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Complex number struct
typedef struct {
    double real;
    double imag;
} Complex;

// Parameters for Mandelbrot set
typedef struct {
    double xmin, xmax, ymin, ymax;
    int maxiter, height, width, xres, yres;
} Parameters;

// Function to calculate the Mandelbrot iteration 
Complex MandelbrotFormula(Complex z, Complex c) {
    Complex new_z;
    new_z.real = z.real * z.real - z.imag * z.imag + c.real;
    new_z.imag = 2 * z.real * z.imag + c.imag;
    return new_z;
}

// Function to calculate the modulus squared of a complex number
double Sq_Of_Mod_Of_Complex(Complex c) {
    return c.real * c.real + c.imag * c.imag;
}

// Function to generate the Mandelbrot set and save points to a file
void Generate_File_With_Points(Parameters p, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error: Could not open file for writing.\n");
        return;
    }

    // Step sizes
    double xstep = (p.xmax - p.xmin) / p.xres;
    double ystep = (p.ymax - p.ymin) / p.yres;

    for (int y = 0; y < p.yres; y++) {
        for (int x = 0; x < p.xres; x++) {
            Complex c;
            c.real = p.xmin + x * xstep;
            c.imag = p.ymin + y * ystep;

            Complex z = {0.0, 0.0};  // Initialize z to 0
            int iter = 0;

            // Mandelbrot iteration
            while (Sq_Of_Mod_Of_Complex(z) < 4.0 && iter < p.maxiter) {
                z = MandelbrotFormula(z, c);
                iter++;
            }

            // Output points: 1 for inside, 0 for outside
            if (iter == p.maxiter) {
                fprintf(fp, "1 ");
            } else {
                fprintf(fp, "0 ");
            }
        }
        fprintf(fp, "\n");  // New line after each row for matrix format
    }
    fclose(fp);
}

// Function to plot the Mandelbrot set using gnuplot and export as a JPG
void plotMandelbrotSet(const char *img, const char *data, Parameters p) {
    FILE *gp = popen("gnuplot -persistent", "w");
    if (gp == NULL) {
        printf("Error: Could not open gnuplot.\n");
        return;
    }

    // Export plot as a JPEG
    fprintf(gp, "set terminal jpeg size 800,600\n");
    fprintf(gp, "set output '%s'\n",img);
    fprintf(gp, "set title 'Mandelbrot Set'\n");
    fprintf(gp, "set palette defined (0 'black', 1 'white')\n");
    fprintf(gp, "unset colorbox\n");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set size ratio -1\n");

    // Use pixel-based plotting, but customize axis labels to reflect real-world coordinates
    // Label the x-axis with real coordinates
    fprintf(gp, "set xtics ('%f' 0, '%f' %d, '%f' %d)\n", 
            p.xmin, p.xmax / 2, p.xres / 2, p.xmax, p.xres);

    // Label the y-axis with real coordinates
    fprintf(gp, "set ytics ('%f' 0, '%f' %d, '%f' %d)\n", 
            p.ymin, p.ymax / 2, p.yres / 2, p.ymax, p.yres);

    // Plot the data without custom axis ranges, allowing pixel mapping
    fprintf(gp, "plot '%s' matrix with image\n", data);
    pclose(gp);
}


int main(int argc, char **argv) {

    //Taking input and checking for correct number of arguments
    if (argc != 8 ){
      printf("Usage :./ee24b017_mandelbrot <xmin> <xmax> <ymin> <ymax> <maxiter> <xres> <out.jpg>\n");
      return 1;
    }
    Parameters p;
    p.xmin = atof(argv[1]);
    p.xmax = atof(argv[2]);
    p.ymin = atof(argv[3]);
    p.ymax = atof(argv[4]);
    p.width = 800;
    p.height = 600;
    p.xres = atoi(argv[6]);
    p.yres = (int)(p.xres * (p.ymax - p.ymin) / (p.xmax - p.xmin));
    p.maxiter = atoi(argv[5]);

    const char *data_filename = "ee24b017_mandelbrot_data.txt";

    // Generate the Mandelbrot set and save to file
    Generate_File_With_Points(p, data_filename);
    const char *img_name = argv[7];
    // Plot the Mandelbrot set using gnuplot and export it as a JPG
    plotMandelbrotSet(img_name,data_filename,p);
    // Removing the data file with the points.
    remove("ee24b017_mandelbrot_data.txt");
    return 0;
}

