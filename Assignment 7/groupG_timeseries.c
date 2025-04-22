/*
            EE1103 Assignment - 7 : Time Series Analysis          
Developer : Group G (EE24B017,EE24B058,EE24B068) 
Date : 19th September 2024
Purpose : 1) To generate a lorentzian time series, add noise to it, then use a peak detection algorithm to get the average width and time gap between each peak.
          2) After successfully testing the code, it should be used to analyze a given CSV file.
Input(s) : M(Number of peaks), T(Time period between peaks), a(FWHM) OR <filename>.csv
Outputs(s): avg(T) avg(a) stdev(T) stdev(a)
Used ChatGPT for help
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
typedef struct {
    double location;   // Peak location (in time units)
    double amplitude;  // Peak amplitude (value at the peak)
    double width;      // Peak width (FWHM)
} Peak;
double generate_normal_noise() {
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}
double calculate_total_time(int M, double T) {
    return M * T * 1.01;
}
// Generate Lorentzian peaks with correct FWHM and time intervals
void generate_lorentzian_series(double *time_series, double total_time, int N, int M, double T, double a) {
    double dt = total_time / N; // Time increment per index
    for (int i = 0; i < N; i++) {
        time_series[i] = 0.0; // Initialize time series
    }
    for (int m = 1; m <= M; m++) {
        double peak_time = m * T; // Correct peak time location
        for (int t = 0; t < N; t++) {
            double time = t * dt; // Convert index to time
            double diff = time - peak_time;
            time_series[t] += (a * a) / (diff * diff + a * a); // Lorentzian equation
        }
    }
}
void add_noise_to_lorentzian(double *time_series, double total_time, int N, int M, double T, double a, double noise_stddev) {
    double dt = total_time / N; // Time increment per index
    for (int m = 1; m <= M; m++) {
        double noisy_peak_time = m * T + generate_normal_noise() * noise_stddev * T;
        double noisy_width = (a + generate_normal_noise() * noise_stddev * a);
        for (int t = 0; t < N; t++) {
            double time = t * dt;
            double diff = time - noisy_peak_time;
            time_series[t] += (noisy_width * noisy_width) / (diff * diff + noisy_width * noisy_width);
        }
    }
}
// Smoothing with a simple moving average filter (window size = 10)
void smooth_time_series(double *time_series, int N) {
    double *smoothed_series = (double *)malloc(N * sizeof(double));
    for (int i = 5; i < N - 5; i++) {
        smoothed_series[i] = 0.0;
        for (int j = -5; j <= 5; j++) {
            smoothed_series[i] += time_series[i + j];
        }
        smoothed_series[i] /= 11.0;
    }
    for (int i = 0; i < 5; i++) smoothed_series[i] = time_series[i];
    for (int i = N - 5; i < N; i++) smoothed_series[i] = time_series[i];
    // Copy smoothed values back
    for (int i = 0; i < N; i++) {
        time_series[i] = smoothed_series[i];
    }
    free(smoothed_series);
}
double index_to_time(int index, double total_time, int N) {
    return total_time * index / N;
}
int detect_peaks(double *time_series, int N, Peak *peaks, int M, double threshold, double total_time, int min_distance) { 
  int num_peaks = 0;
  int i = 5; // Start at index 5 to avoid edge issues

  while (i < N - 5)
  {
    int candidate_peak_idx = -1;
    double candidate_peak_value = -INFINITY;
    for (int j = i; j < i + min_distance && j < N - 5; j++)
    {
      if (time_series[j] > time_series[j - 1] &&
      time_series[j] > time_series[j + 1] &&
      time_series[j] > threshold)
      {
        if (time_series[j] > candidate_peak_value)
        {
          candidate_peak_value = time_series[j];
          candidate_peak_idx = j;
        }
      }
    }
    if (candidate_peak_idx != -1)
    {
      peaks[num_peaks].location = index_to_time(candidate_peak_idx,total_time,N);
      peaks[num_peaks].amplitude = candidate_peak_value;
      double half_max = candidate_peak_value / 2.0;// Estimate width at half maximum
      int left = candidate_peak_idx, right = candidate_peak_idx;
      while (left > 0 && time_series[left] > half_max)
      left--;
      while (right < N && time_series[right] > half_max)
      right++;
      peaks[num_peaks].width = (index_to_time(right, total_time, N) - index_to_time(left, total_time, N));
      num_peaks++;
      i = candidate_peak_idx + min_distance;
    }
    else
    {
      i++;
    }
  }
  return num_peaks;
}
// Smoothing file content with a simple moving average filter (window size = 21)
void f_smooth_time_series(double *time_series, int N)
  {
  double *smoothed_series = (double *)malloc(N * sizeof(double));
  int window_size = 10; // Half-window size for the moving average (window = 21)
  // Smooth the main part of the series (where full window size can be applied)
  for (int i = window_size; i < N - window_size; i++)
  {
    smoothed_series[i] = 0.0;
    for (int j = -window_size; j <= window_size; j++)
    {
      smoothed_series[i] += time_series[i + j];
    }
    smoothed_series[i] /= (2 * window_size + 1);
  }
  for (int i = 0; i < window_size; i++)
  {
    smoothed_series[i] = 0.0;
    int actual_window_size = i + 1; // Smaller window near the edges
    for (int j = -i; j <= i; j++)
    {
      smoothed_series[i] += time_series[i + j];
    }
    smoothed_series[i] /= (2 * actual_window_size - 1);
  }
  for (int i = N - window_size; i < N; i++)
    {
    smoothed_series[i] = 0.0;
    int actual_window_size = N - i; // Smaller window near the edges
    for (int j = -(N - i - 1); j <= (N - i - 1); j++)
    {
      smoothed_series[i] += time_series[i + j];
    }
    smoothed_series[i] /= (2 * actual_window_size - 1);
    }
  for (int i = 0; i < N; i++)
  {
  time_series[i] = smoothed_series[i];
  }
  free(smoothed_series);
  }
int file_detect_peaks(double *time_series, int N, Peak *peaks, double threshold, int min_distance)
{
    int num_peaks = 0;
    int i = 5; // Start at index 5 to avoid edge issues
    while (i < N - 5)
    {
        int candidate_peak_idx = -1;
        double candidate_peak_value = -INFINITY;
        for (int j = i; j < i + min_distance && j < N - 5; j++)
        {
            if (time_series[j] > time_series[j - 1] &&
                time_series[j] > time_series[j + 1] &&
                time_series[j] > threshold)
            {
                if (time_series[j] > candidate_peak_value)
                {
                    candidate_peak_value = time_series[j];
                    candidate_peak_idx = j;
                }
            }
        }
        if (candidate_peak_idx != -1)
        {
            peaks[num_peaks].location = candidate_peak_idx;
            peaks[num_peaks].amplitude = candidate_peak_value;
            double half_max = candidate_peak_value / 2.0;
            int left = candidate_peak_idx, right = candidate_peak_idx;
            while (left > 0 && time_series[left] > half_max)
                left--;
            while (right < N && time_series[right] > half_max)
                right++;
            peaks[num_peaks].width = right - left;
            num_peaks++;
            i = candidate_peak_idx + min_distance;
        }
        else
        {
            i++;
        }
    }
    return num_peaks;
}
void calculate_statistics(double *data, int size, double *mean, double *stddev) {
    double sum = 0.0, sum_sq = 0.0;

    for (int i = 0; i < size; i++) {
        sum += data[i];
        sum_sq += data[i] * data[i];
    }
    *mean = sum / size;
    *stddev = sqrt((sum_sq / size) - (*mean) * (*mean));
}
// Uncomment this Function to plot generated data using gnuplot
/*
void plot_data(double *time_series, int N, Peak *peaks, int num_peaks, double total_time) {
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe == NULL) {
        fprintf(stderr, "Could not open gnuplot.\n");
        return;
    }
    fprintf(gnuplotPipe, "set title 'Lorentzian Time Series with Peaks'\n");
    fprintf(gnuplotPipe, "set xlabel 'Time'\n");
    fprintf(gnuplotPipe, "set ylabel 'Amplitude'\n");
    fprintf(gnuplotPipe, "plot '-' with lines title 'Time Series', '-' with points title 'Detected Peaks'\n");
    // Plot time series
    for (int i = 0; i < N; i++) {
        fprintf(gnuplotPipe, "%lf %lf\n", index_to_time(i, total_time, N), time_series[i]);
    }
    fprintf(gnuplotPipe, "e\n");
    // Plot detected peaks
    for (int i = 0; i < num_peaks; i++) {
        fprintf(gnuplotPipe, "%lf %lf\n", peaks[i].location, peaks[i].amplitude);
    }
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
}
*/
//Below is the function to plot the given input data
/*
void Plot_data(double *time_series, int N, Peak *peaks, int num_peaks)
{
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe == NULL)
    {
        fprintf(stderr, "Could not open gnuplot.\n");
        return;
    }

    // Set up plot settings
    fprintf(gnuplotPipe, "set title 'Lorentzian Time Series with Peaks'\n");
    fprintf(gnuplotPipe, "set xlabel 'Time'\n");
    fprintf(gnuplotPipe, "set ylabel 'Amplitude'\n");
    fprintf(gnuplotPipe, "plot '-' with dots title 'Time Series', '-' with points title 'Detected Peaks'\n");

    // Plot time series
    for (int i = 0; i < N; i++)
    {
        fprintf(gnuplotPipe, "%d %lf\n", i, time_series[i]);
    }
    fprintf(gnuplotPipe, "e\n");

    // Plot detected peaks
    for (int i = 0; i < num_peaks; i++)
    {
        fprintf(gnuplotPipe, "%lf %lf\n", peaks[i].location, peaks[i].amplitude);
    }
    fprintf(gnuplotPipe, "e\n");

    // Close gnuplot pipe
    pclose(gnuplotPipe);
}
*/
int read_csv_file(const char *filename, double *time_series, int max_rows, int b)
{
    FILE *file = fopen(filename, "r"); // Open the file in read mode
    if (file == NULL)
    {
        printf("Could not open file.\n");
        return -1;
    }
    char line[100]; // Buffer to store each line
    int first_column;
    double second_column;
    fgets(line, sizeof(line), file);// Read and ignore the first line (header)
    int j = 0;
    while (fgets(line, sizeof(line), file) != NULL && j < b * max_rows)
    {
        j++;
    }
    int i = 0;
    while (fgets(line, sizeof(line), file) != NULL && i < max_rows)
    {
        if (sscanf(line, "%d,%lf", &first_column, &second_column) == 2) // Read the first and second columns, assuming the format is: int, double, and then a third column
        {
            time_series[i] = second_column;
            i++; // Increment i only if a value was successfully read
        }
    }
    
    fclose(file);
    return i;
}
int main(int argc, char *argv[]) {
    srand(time(NULL));
    int M = 1000;
    double T = 100.0;
    double a = 5.0;
    double noise_stddev_T_a = 0.07;
    double noise_stddev_amp = 0.05;
    if (argc == 4) {    //For generating data
        M = atoi(argv[1]);
        T = atof(argv[2]);
        a = atof(argv[3]);
        double total_time = calculate_total_time(M, T);
        int N = (int)(total_time * 10); // Increased resolution
        double *time_series = (double *)malloc(N * sizeof(double));
        Peak *peaks = (Peak *)malloc( 2 * M * sizeof(Peak));
        generate_lorentzian_series(time_series, total_time, N, M, T, a);
        add_noise_to_lorentzian(time_series, total_time, N, M, T, a, noise_stddev_T_a);
        smooth_time_series(time_series, N);
        double max_amplitude = 0;
        for (int i = 0; i < N; i++) {
            if (time_series[i] > max_amplitude) max_amplitude = time_series[i];
        }
        double threshold = 0.4 * max_amplitude;
        int num_peaks = detect_peaks(time_series, N, peaks, M, threshold, total_time, 19 * T / 20);
        if (num_peaks > 1) {
            double *peak_intervals = (double *)malloc((num_peaks - 1) * sizeof(double));
            double *peak_widths = (double *)malloc(num_peaks * sizeof(double));

            for (int i = 0; i < num_peaks - 1; i++) {
                peak_intervals[i] = peaks[i + 1].location - peaks[i].location;
            }
            for (int i = 0; i < num_peaks; i++) {
                peak_widths[i] = peaks[i].width;
            }
            double avg_T, avg_a, stddev_T, stddev_a;
            calculate_statistics(peak_intervals, num_peaks - 1, &avg_T, &stddev_T);
            calculate_statistics(peak_widths, num_peaks, &avg_a, &stddev_a);
            printf("%.2f %.2f %.2f %.2f\n", avg_T, avg_a/2.0,  stddev_T, stddev_a/2.0);
            free(peak_intervals);
            free(peak_widths);
        } else {
            printf("Not enough peaks detected.\n");
        }
        /* Use this to plot data 
        plot_data(time_series, N, peaks, num_peaks, total_time);
        */
        free(time_series);
        free(peaks);

        return 0;
    }
    else if (argc == 2)
{
    float countTimeAvg;
    float countTimeStd;
    float countaAvg;
    float countaStd;
    float *timeavgArray = (float *)(malloc(10 * sizeof(float)));
    float *timestdArray = (float *)(malloc(10 * sizeof(float)));
    float *aavgArray = (float *)(malloc(10 * sizeof(float)));
    float *astdArray = (float *)(malloc(10 * sizeof(float)));
    int Max_Rows = 100000;
    if (argc < 2)
    {
        printf("Please provide the input CSV file.\n");
        return 1; // Exit if no input file is given
    }
    char *file = argv[1];
    srand(time(NULL));
    for (int b = 0; b < 10; b++)
    {
        double *time_series = (double *)malloc(Max_Rows * sizeof(double));
        int N = read_csv_file(file, time_series, Max_Rows, b);
        if (N == -1)
            return 1;
        Peak *peaks = (Peak *)malloc(N * sizeof(Peak));
        f_smooth_time_series(time_series, N);
        double max_amplitude = 0;
        for (int i = 0; i < N; i++)
        {
            if (time_series[i] > max_amplitude)
                max_amplitude = time_series[i];
        }
        double threshold = 0.4 * max_amplitude;
        int num_peaks = file_detect_peaks(time_series, N, peaks, threshold, 500);
        if (num_peaks > 1)
        {
            double *peak_intervals = (double *)malloc((num_peaks - 1) * sizeof(double));
            double *peak_widths = (double *)malloc(num_peaks * sizeof(double));
            for (int i = 0; i < num_peaks - 1; i++)
            {
                peak_intervals[i] = peaks[i + 1].location - peaks[i].location;
            }
            for (int i = 0; i < num_peaks; i++)
            {
                peak_widths[i] = peaks[i].width;
            }
            double avg_T, avg_a, stddev_T, stddev_a;
            calculate_statistics(peak_intervals, num_peaks - 1, &avg_T, &stddev_T);
            calculate_statistics(peak_widths, num_peaks, &avg_a, &stddev_a);
            timeavgArray[b] = avg_T;
            timestdArray[b] = stddev_T;
            aavgArray[b] = avg_a / 2.0;
            astdArray[b] = stddev_a / 2.0;
            free(peak_intervals);
            free(peak_widths);
        }
        else
        {
            printf("Not enough peaks detected.\n");
        }
        //Plot the data
        // Plot_data(time_series, N, peaks, num_peaks);
        free(time_series);
        time_series = NULL;
        free(peaks);
        peaks = NULL;
    }
    countTimeAvg = 0;
    for (int i = 0; i < 10; i++)
    {
        countTimeAvg += timeavgArray[i];
    }
    countTimeStd = 0;
    for (int i = 0; i < 10; i++)
    {
        countTimeStd += timestdArray[i];
    }
    countaAvg = 0;
    for (int i = 0; i < 10; i++)
    {
        countaAvg += aavgArray[i];
    }
    countaStd = 0;
    for (int i = 0; i < 10; i++)
    {
        countaStd += astdArray[i];
    }
    printf("%0.2f %0.2f %0.2f %0.2f\n", countTimeAvg / 10.0, countaAvg / 10.0, countTimeStd / 10.0, countaStd / 10.0);
    free(timeavgArray);
    free(timestdArray);
    free(aavgArray);
    free(astdArray);
    return 0;
}
else
{
    printf("Usage: M(Number of peaks) T(Time period between peaks) a(FWHM) \nOR\n<filename>.csv\n");
}
}
