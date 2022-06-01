// This file declares and defines various functions that are used
// throughout the text.  Updated 1-10-2020.

// These are some standard C function libraries.
#include <stdlib.h>  // "standard library"
#include <math.h>    // various math functions, such as exp()
#include <stdio.h>   // various "input/output" functions
#include <time.h>    // functions for timing computations

// Function definitions are found below.

// - Mersenne Twister RNG.
double MTUniform (unsigned int);

// - Histogram functions.
void Histogram (double, double, double, int, int);
void DiscreteHistogram (int, int, int, int);
void NormalHistogram (double, int, int);
void ExponentialHistogram (double, int, int);
void UniformHistogram (double, int, int);

// - Functions related to standard normals.
double Psi (double);
double PsiInv (double);

// - Miscellanceous functions.
void    Pause (void);
void    Exit (void);
double  Time (void);
int     Equal (double, double, double);
int     GetInteger (char *);
double  GetDouble (char *);
double *List (int);

// - Linear algebra functions.
double **Array (int, int);
void     Show (double **);
void     Write (double **, FILE *);
double **Transpose (double **);
double **Multiply (double **, double **);
double **Invert (double **);
double **Copy (double **);
double **Identity (int);
double **Cholesky (double **);
double **ScalarMultiple (double, double **);
double **Add (double **, double **);
int      Rows (double **);
int      Columns (double **);
void     Free (double **);
double **Evector (double **, double **);
double **MeanZero (double **);
double **Covariance (double **);
double **Correlation (double **);
double **QRalgorithm (double **);
double **GramSchmidt (double **);

// - Black-Scholes call option functions.
double BlackScholes (double, double, double, double, double);
double ImpliedVol (double, double, double, double, double);


#define ZERO 1e-8

////////////////////////////////////////////////////////////////////////////////
// FUNCTION DEFINITIONS FOLLOW
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// MERSENNE TWISTER
// By M. Matsumoto and T. Nishimura (1998).
// "Mersenne Twister: a 623-dimensionally equidistributed uniform pseudo-random
//   number generator".
// ACM Transactions of Modeling and Computer Simulation 8(1):3-30.
// Any coding errors introduced are my own (C.D. Howard).

// An unsigned integer is represented by 32 bits in base 2.  The largest
// unsigned integer is: 2^32 - 1 = 4,294,967,295 (base 10);
//          = ffffffff (base 16);
//          = 1111 1111 1111 1111 1111 1111 1111 1111 (base 2).
// The digits in hexadecimal (base 16) are 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, a, b,
//                                         c, d, e, f.
////////////////////////////////////////////////////////////////////////////////
double MTUniform (unsigned int seed) {

   static unsigned int X[1248], m[2], initialized = 0, k;
   unsigned int N, Y;

   // Seed the RNG when a new seed is passed or it has not yet been initialized.
   if (seed || !initialized) {
      // If no seed is specified, default is 1.
      X[0] = (seed ? seed : 1);
      // Seed X[1],X[2],...,X[623] with your favorite LCG.
      for (k = 1; k < 624; k++) {
         X[k] = 22695477 * X[k-1] + 1;
      }
      m[0] = 0; m[1] = 0x9908b0df;
      // The counter "k" is now 624.
      initialized = 1;
   }

   // Concatenate the first bit of X[k-624] with the last 31 bits of X[k-623],
   //    "|" is "bit-wise or".
   Y = (X[k-624] & 0x80000000) | (X[k-623] & 0x7fffffff);

   // Now apply an invertible linear transformation Y and bit-wise add X[k-227].
   X[k] = ((Y >> 1) ^ m[Y & 1] ^ X[k-227]);

   // Re-load X[0],X[1],...,X[623] as you go.
   X[k-624] = X[k];

   // Apply the tempering function (also an invertible linear transformation).
   N = X[k];
   N ^= (N >> 11);
   N ^= (N << 7) & 0x9d2c5680;
   N ^= (N << 15) & 0xefc60000;
   N ^= (N >> 18);

   // Increment the counter; shift vectors when appropriate.
   k ++;
   if (k == 1248) k = 624;

   // Now 0 <= N <= 4,294,967,295; scale it to be on the interval (0,1).
   return ( (N + 0.5) / 4294967296.0 );

}

////////////////////////////////////////////////////////////////////////////////
// This function waits for a user-input "enter", then CONTINUES the program.
////////////////////////////////////////////////////////////////////////////////
void Pause () {

   char input[100];

   printf ("\n");
   printf ("Hit Enter to continue program... ");
   fgets (input, 9, stdin);

   return;

}

////////////////////////////////////////////////////////////////////////////////
// This function waits for a user-input "enter", then EXITS the program.
// It prevents the window from closing up before the output can be viewed.
////////////////////////////////////////////////////////////////////////////////
void Exit () {

   char input[100];

   printf ("\n");
   printf ("Hit Enter to exit program... ");
   fgets (input, 9, stdin);

   exit (0);

}

////////////////////////////////////////////////////////////////////////////////
// This function computes elapsed computation time in seconds.
////////////////////////////////////////////////////////////////////////////////
double Time () {

   static clock_t time;
   static int initialized = 0;

   // With the first call to this function, "time" is initialized.
   if (!initialized) {
      time = clock ();
      initialized = 1;
   }

   // With subsequent calls, elapsed time since the first call is returned.
   return ((double) (clock() - time)) / CLOCKS_PER_SEC;

}

////////////////////////////////////////////////////////////////////////////////
// This function gets an integer value typed in by the user.
////////////////////////////////////////////////////////////////////////////////
int GetInteger (char *question) {

   char input[100];
   int n;

   // Print the question.
   printf (question);

   // Get the answer.
   fgets (input, 99, stdin);
   sscanf (input, "%d", &n);

   return (n);

}

////////////////////////////////////////////////////////////////////////////////
// This function gets a double precision value typed in by the user.
////////////////////////////////////////////////////////////////////////////////
double GetDouble (char *question) {

   char input[100];
   double x;

   // Print the question.
   printf (question);

   // Get the answer.
   fgets (input, 99, stdin);
   sscanf (input, "%lf", &x);

   return (x);

}

////////////////////////////////////////////////////////////////////////////////
// This function allocates array space for a list of doubles.
////////////////////////////////////////////////////////////////////////////////
double *List (int n) {

   double *a;

   a = (double *) calloc (n+1, sizeof (double));

   return a;

}

////////////////////////////////////////////////////////////////////////////////
// This function creates a histogram of randomly generated numbers. It is
// typically used with continuous (not discrete)  data.
////////////////////////////////////////////////////////////////////////////////
void   Histogram (double x,        // Item of data to be added to the histogram.
                  double a,        // Lower bound of histogram.
                  double b,        // Upper bound of histogram.
                  int n,           // Number of "buckets".
                  int done)        // 0 if not done, 1 if done.

{

   int k;
   double x0, x1, biggest;
   FILE *fp;

   static int initialized=0;
   static double dx, *freq;

   // Initialize certain variables on the first call to this function.
   if (!initialized) {

      // Allocate array space for frequencies.
      freq = (double *) calloc (n, sizeof (double));

      // Compute size of each "bucket".
      dx = (b - a) / n;

      // Variables have now been initialized.
      initialized = 1;

   }

   // If adding a new piece of data.
   if (!done) {

      // See if data is below the lower bound.
      if (x <= a) {
         freq[0] ++;
      }

      // See if data is above the upper bound.
      else if (x >= b) {
         freq[n-1] ++;
      }

      // Otherwise find the appropriate bucket.
      else {
         k = (x - a) / dx;
         freq[k] ++;
      }

   }

   // If finished, create the TeX files that generate the plot. This is the
   // final call to the histogram function.
    else {

      // Report the raw data to a TXT file for viewing without TeX.
      fp = fopen ("HistogramTxtData.txt", "w");
      fprintf (fp, "      From         To  Frequency\n");
      fprintf (fp, "========== ========== ==========\n");
      for (k = 0; k < n; k++) {
         x0 = a + k * dx;
         x1 = a + (k+1)* dx;
         fprintf (fp, "%10.5f %10.5f %10.0f\n", x0, x1, freq[k]);
      }
      fprintf (fp, "\n");
      fprintf (fp, "Data outside histogram range is put in first or last bucket.\n");
      fclose (fp);

      // Generate TeX output files.
      // Scale so that the biggest bucket has 1 (for graphing convenience).

      // First find the largest bucket value.
      biggest = 0;
      for (k = 0; k < n; k++) {
         if (freq[k] > biggest) {
            biggest = freq[k];
         }
      }

      // Now re-scale all the bucket values.
      for (k = 0; k < n; k++) {
         freq[k] /= biggest;
      }

      // Report data to the output file.
      fp = fopen ("HistogramData.txt", "w");
      fprintf (fp, "%10.5f  %10.5f\n", a, 0.0);
      for (k = 0; k < n; k++) {
         fprintf (fp, "%10.5f  %10.5f\n", a + (k+1) * dx, freq[k]);
      }
      fclose (fp);

      // Free up the freq[] array.
      free (freq);

      fp = fopen ("Histogram.tex", "w");
      fprintf (fp, "\\input pictex\\magnification=\\magstep1\\nopagenumbers\n");
      fprintf (fp, "\\centerline {\n");
      fprintf (fp, "\\beginpicture\n");
      fprintf (fp, "\\setcoordinatesystem units <%8.3f truein, 2.0 truein>\n", 4.0 / (b-a));
      fprintf (fp, "\\setplotarea x from %8.3f to %8.3f, y from  0 to 1.0\n", a, b);
      fprintf (fp, "\\axis bottom\n");
      fprintf (fp, "ticks numbered from %8.2f to %8.2f by %8.2f\n", a, b, (b-a)/10);
      fprintf (fp, "/\n");
      fprintf (fp, "\\plot %8.3f 0  %8.3f 0  %8.3f 1  %8.3f 1 /\n",
                a-.045*(b-a), a-.03*(b-a), a-.03*(b-a), a-.045*(b-a) );
      fprintf (fp, "\\put {0} [cr] at %8.3f  0\n", a-.055*(b-a));
      fprintf (fp, "\\put {%d} [cr] at %8.3f  1\n", (int) biggest, a-.055*(b-a));
      fprintf (fp, "\\sethistograms\n");
      fprintf (fp, "\\plot \"HistogramData.txt\"\n");
      fprintf (fp, "\\put {\\sl Histogram for a Continuous Random Variable} at %8.3f 1.2\n", (a+b)/2.0);
      fprintf (fp, "\\endpicture}\\vfill\\end\n");
      fclose (fp);


   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// This function creates a histogram of randomly generated numbers. It is
// typically used with discrete data, like binomially distributed data.
////////////////////////////////////////////////////////////////////////////////
void   DiscreteHistogram (int x,          // Data to be added to the histogram.
                          int a,          // Lower bound of histogram.
                          int b,          // Upper bound of histogram.
                          int done)       // 0 if not done, 1 if done.

{

   static double *freq;
   static int initialized = 0;
   int k;
   double biggest;
   FILE *fp;

   // Allocate space on the first call to this function.
   if (! initialized) {

      freq = (double *) calloc (b-a+2, sizeof (double));

      initialized = 1;

   }

   // If not yet finished, add the current item of data to the histogram.
   if (!done) {

      // Truncate if data is below the lower bound.
      if (x < a) {
         x = a;
      }

      // Truncate if data is above the upper bound.
      if (x > b) {
         x = b;
      }

      // Increment the appropriate bucket.
      freq[x-a] ++;


   }

   // If finished, create the TeX files that generate the plot. This is the
   // final call to the histogram function.
   else {

      // Report the raw data to a TXT file for viewing without TeX.
      fp = fopen ("HistogramTxtData.txt", "w");
      fprintf (fp, "         X  Frequency\n");
      fprintf (fp, "========== ==========\n");
      for (k = a; k <= b; k++) {
         fprintf (fp, "%10d %10.0f\n", k, freq[k]);
      }
      fprintf (fp, "\n");
      fprintf (fp, "Data outside histogram range is put in first or last bucket.\n");
      fclose (fp);

      fp = fopen ("Histogram.tex", "w");
      fprintf (fp, "\\input pictex\\magnification=\\magstep1\\nopagenumbers\n");
      fprintf (fp, "\\centerline {\n");
      fprintf (fp, "\\beginpicture\n");
      fprintf (fp, "\\setcoordinatesystem units <%8.3f truein, 2.0 truein>\n", 4.0 / (b-a));
      fprintf (fp, "\\setplotarea x from %d to %d, y from  0 to 1.0\n", a, b);
      fprintf (fp, "\\axis bottom\n");
      fprintf (fp, "ticks numbered from %d to %d by 1\n", a, b);
      fprintf (fp, "/\n");
      fprintf (fp, "\\linethickness=%ftruein\n", 1.0/(b-a));
      fprintf (fp, "\\input HistogramData.txt\n");
      fprintf (fp, "\\put {\\sl Histogram for a Discrete Random Variable} at %8.3f 1.2\n", (a+b)/2.0);
      fprintf (fp, "\\endpicture}\\vfill\\end\n");
      fclose (fp);

      // Scale so that the biggest bucket has 1 (for graphing convenience).

      // First find the largest bucket.
      biggest = 0;
      for (k = 0; k <= b-a; k++) {
         if (freq[k] > biggest) {
            biggest = freq[k];
         }
      }

      // Now scale all the buckets.
      for (k = 0; k <= b-a; k++) {
         freq[k] /= biggest;
      }

      // Report data to the output file.
      fp = fopen ("HistogramData.txt", "w");
      for (k = 0; k <= b-a; k++) {
         fprintf (fp, "\\putrule from %d 0 to %d %f \n", a+k, a+k, freq[k]);
      }
      fclose (fp);

      // Free up the freq[] array.
      free (freq);

   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// This histogram function is specially designed for data that should be
// Normal(0,1) in distribution.  The function computes a properly scaled
// standard normal density function for comparison to the data.
////////////////////////////////////////////////////////////////////////////////
void NormalHistogram (double x,          // Data to be added to the histogram.
                      int n,             // Number of "buckets".
                      int done)          // 0 if not done, 1 if done.

{

   int k;
   double z, x0, x1, area;
   FILE *fp;

   static int N, initialized=0;
   static double a, b, dx, *freq;

   // Initialize certain variables on the first call to this function.
   if (!initialized) {

      // Allocate array space for frequencies.
      freq = (double *) calloc (n, sizeof (double));

      // Hardwire upper and lower bounds for the standard normal.
      a = -5; b = 5;

      // Compute size of each "bucket".
      dx = (b - a) / n;

      // Data counter.
      N = 0;

      // Variables have now been initialized.
      initialized = 1;

   }

   // If adding a new piece of data.
   if (!done) {

      // See if data is below the lower bound.
      if (x <= a) {
         freq[0] ++;
      }

      // See if data is above the upper bound.
      else if (x >= b) {
         freq[n-1] ++;
      }

      // Otherwise find the appropriate bucket.
      else {
         k = (x - a) / dx;
         freq[k] ++;
      }

      // Increment the data counter.
      N ++;

   }

   // ...otherwise, when finished, create the TeX output files for viewing.
   else {

      // Report the raw data to a TXT file for viewing without TeX.
      fp = fopen ("HistogramTxtData.txt", "w");
      fprintf (fp, "                                   Expected\n");
      fprintf (fp, "      From         To  Frequency  Frequency\n");
      fprintf (fp, "========== ========== ========== ==========\n");
      for (k = 0; k < n; k++) {
         x0 = a + k * dx;
         x1 = a + (k+1)* dx;
         fprintf (fp, "%10.5f %10.5f %10.0f %10.0f\n",
                       x0, x1, freq[k], N * (Psi(x1) - Psi(x0)));
      }
      fprintf (fp, "\n");
      fprintf (fp, "Data outside histogram range is put in first or last bucket.\n");
      fclose (fp);

      // Create the main TeX file.
      fp = fopen ("Histogram.tex", "w");
      fprintf (fp, "\\input pictex\\magnification=\\magstep1\\nopagenumbers\n");
      fprintf (fp, "\\centerline {\n");
      fprintf (fp, "\\beginpicture\n");
      fprintf (fp, "\\setcoordinatesystem units <0.45 truein, 1.5 truein>\n", 5.0 / (b-a));
      fprintf (fp, "\\setplotarea x from -5 to 5, y from  0 to 1.0\n");
      fprintf (fp, "\\axis bottom\n");
      fprintf (fp, "ticks numbered from -5 to 5 by 1\n");
      fprintf (fp, "/\n");
      fprintf (fp, "\\plot \"Normal.txt\"\n");
      fprintf (fp, "\\plot -5.2 0 -5.1 0 -5.1 1  -5.2 1 /\n");
      fprintf (fp, "\\put {0} [cr] at -5.3 0\n");
      fprintf (fp, "\\put {$\\frac{1}{\\sqrt{2\\pi}}$} [cr] at -5.3 1\n");
      fprintf (fp, "\\put {\\sl Standard Normal Histogram} at 0 1.2\n");
      fprintf (fp, "\\sethistograms\n");
      fprintf (fp, "\\plot \"HistogramData.txt\"\n");
      fprintf (fp, "\\endpicture}\\vfill\\end\n");
      fclose (fp);

      // Scale so that the area is sqrt(2pi) (for comparing to the Normal(0,1) density).
      area = 0;
      for (k = 0; k < n; k++) {
         area += freq[k] * dx;
      }
      area /= sqrt(2*3.14159);
      for (k = 0; k < n; k++) {
         freq[k] /= area;
      }

      // Report the histogram to an output file.
      fp = fopen ("HistogramData.txt", "w");
      fprintf (fp, "%10.5f  %10.5f\n", a, 0.0);
      for (k = 0; k < n; k++) {
         fprintf (fp, "%10.5f  %10.5f\n", a + (k+1)*dx, freq[k]);
      }
      fclose (fp);

      // Report the appropriately scaled normal density function for comparison.
      fp = fopen ("Normal.txt", "w");
      for (z = -5; z <= 5; z += .01) {
         fprintf (fp, "%8.3f %8.4f\n", z, exp(-z*z/2));
      }
      fclose (fp);

      // Restore to initial conditions.
      initialized = 0;
      free (freq);

   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// This function creates a histogram of the variable x and compares it to
// the density function for a mean 1 exponential.
////////////////////////////////////////////////////////////////////////////////
void ExponentialHistogram (double x,          // data to be added to the histogram
                           int n,             // number of "buckets"
                           int done)          // 1 when done, 0 otherwise

{

   int k;
   double z, x0, x1, area;
   FILE *fp;

   static int N=0, initialized=0;
   static double a, b, dx, *freq;

   // Initialize certain variables on the first call to this function.
   if (!initialized) {

      // Allocate array space for frequencies.
      freq = (double *) calloc (n, sizeof (double));

      // Hardwire upper and lower bounds for the mean 1 exponential distribution.
      // The cutoff of 5 is arbitrary but it captures most of the action.
      a = 0; b = 5;

      // Compute size of each "bucket".
      dx = (b - a) / n;

      // Variables have now been initialized.
      initialized = 1;

   }

   // If adding a new piece of data...
   if (!done) {

      // See if data is below the lower bound.
      if (x <= a) {
         freq[0] ++;
      }

      // See if data is above the upper bound.
      else if (x >= b) {
         freq[n-1] ++;
      }

      // Otherwise find the appropriate bucket.
      else {
         k = (x - a) / dx;
         freq[k] ++;
      }

      // Increment the data counter.
      N ++;

      // ...otherwise, when finished, create the TeX output files for viewing.
   } else {

      // Report the raw data to a TXT file for viewing without TeX.
      fp = fopen ("HistogramTxtData.txt", "w");
      fprintf (fp, "                                   Expected\n");
      fprintf (fp, "      From         To  Frequency  Frequency\n");
      fprintf (fp, "========== ========== ========== ==========\n");
      for (k = 0; k < n; k++) {
         x0 = a + k * dx;
         x1 = a + (k+1)* dx;
         fprintf (fp, "%10.5f %10.5f %10.0f %10.0f\n",
                       x0, x1, freq[k], N * (exp(-x0) - exp(-x1)));
      }
      fprintf (fp, "\n");
      fprintf (fp, "Data outside histogram range is put in last bucket.\n");
      fclose (fp);

      // Create the main TeX file.
      fp = fopen ("Histogram.tex", "w");
      fprintf (fp, "\\input pictex\\magnification=\\magstep1\\nopagenumbers\n");
      fprintf (fp, "\\centerline {\n");
      fprintf (fp, "\\beginpicture\n");
      fprintf (fp, "\\setcoordinatesystem units <1 truein, 1.5 truein>\n");
      fprintf (fp, "\\setplotarea x from 0 to 5, y from  0 to 1.0\n");
      fprintf (fp, "\\axis bottom\n");
      fprintf (fp, "ticks numbered from 0 to 5 by 1\n");
      fprintf (fp, "/\n");
      fprintf (fp, "\\put {$\\leftarrow$} at 5.1 %8.3f\n", exp(-b) * n / (b-a));
      fprintf (fp, "\\plot \"Exponential.txt\"\n");
      fprintf (fp, "\\plot -.3 0 -.2 0 -.2 1  -.3 1 /\n");
      fprintf (fp, "\\put {0} [cr] at -.35 0\n");
      fprintf (fp, "\\put {1} [cr] at -.35 1\n");
      fprintf (fp, "\\put {\\sl Mean-1 Exponential Histogram} at 2.5 1.2\n");
      fprintf (fp, "\\sethistograms\n");
      fprintf (fp, "\\plot \"HistogramData.txt\"\n");
      fprintf (fp, "\\endpicture}\\vfill\\end\n");
      fclose (fp);

      // Scale frequencies so that the area is 1 (for comparing to the Exp(1) density).
      area = 0;
      for (k = 0; k < n; k++) {
         area += freq[k] * (b-a) / n;
      }
      for (k = 0; k < n; k++) {
         freq[k] /= area;
      }

      // Report histogram data.
      fp = fopen ("HistogramData.txt", "w");
      fprintf (fp, "%10.5f  %10.5f\n", a, 0.0);
      for (k = 0; k < n; k++) {
         fprintf (fp, "%10.5f  %10.5f\n",
                       a + (k + 1) * ((b - a) / n), freq[k]);
      }
      fclose (fp);

      // Report the exp(1) density function.
      fp = fopen ("exponential.txt", "w");
      for (z = 0; z <= b; z += .01) {
         fprintf (fp, "%8.3f %8.4f\n", z, exp(-z));
      }
      fclose (fp);

   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// This program creates a histogram of the data x, and plots the Uniform[0,1]
// density function for comparison.
////////////////////////////////////////////////////////////////////////////////
void UniformHistogram (double x,         // data to be added to the histogram
                      int n,             // number of "buckets"
                      int done)          // 1 when done, 0 otherwise
{

   int k;
   double x0, x1, area;
   FILE *fp;

   static int N=0, initialized=0;
   static double a, b, dx, *freq;

   // Initialize certain variables on the first call to this function.
   if (!initialized) {

      // Allocate array space for frequencies.
      freq = (double *) calloc (n, sizeof (double));

      // Hardwire upper and lower bounds for the Uniform[0,1].
      a = 0; b = 1;

      // Compute size of each "bucket".
      dx = (b - a) / n;

      // Variables have now been initialized.
      initialized = 1;

   }

   // If adding a new piece of data.
   if (!done) {

      // See if data is below the lower bound.
      if (x <= a) {
         freq[0] ++;
      }

      // See if data is above the upper bound.
      else if (x >= b) {
         freq[n-1] ++;
      }

      // Otherwise find the appropriate bucket.
      else {
         k = (x - a) / dx;
         freq[k] ++;
      }

      // Increment the data counter.
      N ++;

   }

   // When finished, create the output files for viewing.
   else {

      // Report the raw data to a TXT file for viewing without TeX.
      fp = fopen ("HistogramTxtData.txt", "w");
      fprintf (fp, "                                   Expected\n");
      fprintf (fp, "      From         To  Frequency  Frequency\n");
      fprintf (fp, "========== ========== ========== ==========\n");
      for (k = 0; k < n; k++) {
         x0 = a + k * dx;
         x1 = a + (k+1) * dx;
         fprintf (fp, "%10.5f %10.5f %10.0f %10.0f\n", x0, x1, freq[k], N*dx);
      }
      fclose (fp);

      // Create the main TeX file.
      fp = fopen ("Histogram.tex", "w");
      fprintf (fp, "\\input pictex\\magnification=\\magstep1\\nopagenumbers\n");
      fprintf (fp, "\\centerline {\n");
      fprintf (fp, "\\beginpicture\n");
      fprintf (fp, "\\setcoordinatesystem units <4.5 truein, 1.5 truein>\n", 5.0 / (b-a));
      fprintf (fp, "\\setplotarea x from 0 to 1, y from  0 to 1.0\n");
      fprintf (fp, "\\axis bottom\n");
      fprintf (fp, "ticks numbered from 0 to 1 by .1\n");
      fprintf (fp, "/\n");
      fprintf (fp, "\\plot 0 1  1 1 /\n");    // the Uniform[0,1] density function
      fprintf (fp, "\\plot -.06 0 -.05 0 -.05 1  -.06 1 /\n");
      fprintf (fp, "\\put {0} [cr] at -.07 0\n");
      fprintf (fp, "\\put {1} [cr] at -.07 1\n");
      fprintf (fp, "\\put {\\sl Uniform (0,1) Histogram} at 0.5 1.2\n");
      fprintf (fp, "\\sethistograms\n");
      fprintf (fp, "\\plot \"HistogramData.txt\"\n");
      fprintf (fp, "\\endpicture}\\vfill\\end\n");
      fclose (fp);

      // Scale so that the area is 1 (for comparing to the uniform[0,1] density).
      area = 0;
      for (k = 0; k < n; k++) {
         area += freq[k] * dx;
      }
      for (k = 0; k < n; k++) {
         freq[k] /= area;
      }

      // Report the histogram data.
      fp = fopen ("HistogramData.txt", "w");
      fprintf (fp, "%10.5f  %10.5f\n", a, 0.0);
      for (k = 0; k < n; k++) {
         fprintf (fp, "%10.5f  %10.5f\n", a + (k+1)*dx, freq[k]);
      }
      fclose (fp);


   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Test to see if a = b up to an allowable error of epsilon.
////////////////////////////////////////////////////////////////////////////////
int Equal (double a, double b, double epsilon) {

   int equal;

   equal = (fabs(a-b) <= epsilon);

   return (equal);

}

////////////////////////////////////////////////////////////////////////////////
// This function computes Psi(x) via its Taylor series expansion.
// for -5.5 <= x <= 5.5.  For values of x outside this
// range the computation is handled differently.
////////////////////////////////////////////////////////////////////////////////
double Psi (double x) {

   static double  OneOverSqrtTwoPi = 0.398942280401433, r[201];
   static int initialized = 0;
   double u, term;
   int n, n0;

   if (!initialized) {
      for (n = 1; n <= 200; n++) {
         r[n] = -(2.0*n-1.0)/(2.0*n+1.0)/n/2.0;
      }
      initialized = 1;
   }

   // Extreme cases (correct to within 10^-9).
   if (x < -5.5) return -(1.0 - 1.0/x/x) * (OneOverSqrtTwoPi/x) * exp(-x*x/2.0);
   if (x >  5.5) return 1.0 - (1.0 - 1.0/x/x) * (OneOverSqrtTwoPi/x) * exp(-x*x/2.0);

   // General case to at least within 10^-9.
   n0 = (int) (0.5 * x * x) + 1;
   n = 0;
   term = OneOverSqrtTwoPi * x;
   u = 0.5 + term;
   while (1) {
      n ++;
      term *= (r[n] * x * x);
      u += term;
      if (n > n0 && fabs(term) < 1e-9) break;
   }

   return (u);

}

////////////////////////////////////////////////////////////////////////////////
// This algorithm is due to Peter John Acklam. Any coding errors are my own.
// It is a very good approximation of Psi^(-1) -- according to Acklam the
//    relative error is always <= 1.15e-9.
////////////////////////////////////////////////////////////////////////////////
double PsiInv (double u) {

   static double
    A1 =  -3.969683028665376e+01,
    A2 =   2.209460984245205e+02,
    A3 =  -2.759285104469687e+02,
    A4 =   1.383577518672690e+02,
    A5 =  -3.066479806614716e+01,
    A6 =   2.506628277459239e+00,
    B1 =  -5.447609879822406e+01,
    B2 =   1.615858368580409e+02,
    B3 =  -1.556989798598866e+02,
    B4 =   6.680131188771972e+01,
    B5 =  -1.328068155288572e+01,
    C1 =  -7.784894002430293e-03,
    C2 =  -3.223964580411365e-01,
    C3 =  -2.400758277161838e+00,
    C4 =  -2.549732539343734e+00,
    C5 =   4.374664141464968e+00,
    C6 =   2.938163982698783e+00,
    D1 =   7.784695709041462e-03,
    D2 =   3.224671290700398e-01,
    D3 =   2.445134137142996e+00,
    D4 =   3.754408661907416e+00,
    P0 =  0.02425,
    P1 =  0.97575;

   double N, q, r;

   // Left tail.
   if (u < P0) {
      q = sqrt(-2*log(u));
      N = (((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
   }

   // Middle, general case.
   else if (u <= P1) {
      q = u - 0.5;
      r = q*q;
      N = (((((A1*r+A2)*r+A3)*r+A4)*r+A5)*r+A6)*q /(((((B1*r+B2)*r+B3)*r+B4)*r+B5)*r+1);
   }

   // Right tail.
   else {
      q = sqrt(-2*log(1.0-u));
      N = -(((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
   }

   return (N);

}

////////////////////////////////////////////////////////////////////////////////
// Allocate array space for an m x n array.
////////////////////////////////////////////////////////////////////////////////
double **Array (int m, int n) {

   int i;
   double **A;

   A = (double **) calloc (m+1, sizeof (double *));
   for (i = 0; i <= m; i++) {
      A[i] = (double *) calloc (n+1, sizeof (double));
   }

   // Record dimensions in the 0^th row, which is unused in matrix operations.
   A[0][0] = m;
   A[0][1] = n;

   return (A);

}

////////////////////////////////////////////////////////////////////////////////
// Multiply the matrix A by the scalar c.
////////////////////////////////////////////////////////////////////////////////
double **ScalarMultiple (double c, double **A) {

   int n, m, i, j;
   double **cA;

   m = Rows (A);
   n = Columns (A);

   cA = Array (m, n);

   for (i = 1; i <= m; i++) {
      for (j = 1; j <= n; j++) {
         cA[i][j] = c * A[i][j];
      }
   }

   return cA;

}

////////////////////////////////////////////////////////////////////////////////
// Compute the transpose of the matrix A.
////////////////////////////////////////////////////////////////////////////////
double **Transpose (double **A) {

   int i, j, n, m;
   double **At;

   m = Rows (A);
   n = Columns (A);

   At = Array (n,m);

   for (i = 1; i <= n; i++) {
      for (j = 1; j <= m; j++) {
         At[i][j] = A[j][i];
      }
   }

   return At;

}

////////////////////////////////////////////////////////////////////////////////
// Use Gaussian elimination to invert the square matrix A0.
////////////////////////////////////////////////////////////////////////////////
double **Invert (double **A0) {

   int i, j, k, n, rmax;
   double **A, **Ainv, c;

   // Make a copy so that the original matrix is not changed during G.E.
   A = Copy (A0);

   n = Rows (A);

   if (n != Columns (A)) {
      printf ("Trying to invert a non-square matrix.\n");
      Pause ();
      exit(1);
   }

   // Start with the n x n identity matrix.  This matrix will eventually hold
   //    the inverse of A.
   Ainv = Identity (n);

   // Work on column j of A.
   for (j = 1; j <= n; j++) {

      // Find the largest non-zero entry in the column on or below the diagonal.
      c = 0;
      rmax = 0;
      for (i = j; i <= n; i++) {
         if (fabs(A[i][j]) > c) {
            c = fabs (A[i][j]);
            rmax = i;
         }
      }

      // If they are all 0 the matrix is singular.
      if (c < ZERO) {
         printf ("Trying to invert a singular matrix.\n");
         Exit();
      }

      // Swap rows j and rmax in both A and Ainv.
      i = j;
      for (k = 1; k <= n; k++) {
         c = A[i][k];
         A[i][k] = A[rmax][k];
         A[rmax][k] = c;
         c = Ainv[i][k];
         Ainv[i][k] = Ainv[rmax][k];
         Ainv[rmax][k] = c;
      }

      // Scale so the pivot is 1.
      c = A[i][i];
      for (k = 1; k <= n; k++) {
         A[i][k] /= c;
         Ainv[i][k] /= c;
      }

      // Make rest of column j equal to 0. Apply same row operations to Ainv.
      for (i = 1; i <= n; i++) if (i != j) {
         c = A[i][j];
         for (k = 1; k <= n; k++) {
            A[i][k] += -c * A[j][k];
            Ainv[i][k] += -c * Ainv[j][k];
         }
      }

   }

   // We're done with A.
   Free (A);

   return Ainv;

}

////////////////////////////////////////////////////////////////////////////////
// Compute the eigenvalues for the PD matrix V. Put the result in
//    E[1][1], E[1][2], ... , E[1][n].  Do this using the QR decomposition.
//    Order them from largest to smallest.
////////////////////////////////////////////////////////////////////////////////
double **QRalgorithm (double **V) {

   int i, j, k, n, done, i_max;
   double **Q, **A, **Qt, **AQ, **QtAQ, **E, lambda_max;

   // Record the size of V for use below.
   n = Rows (V);

   // Make a working copy of V so that V's entries remain unchanged.
   A = Copy (V);


   done = 0;
   while (!done) {

      // Next iteration of the QR algorithm.
      Q    = GramSchmidt (A);
      Qt   = Transpose (Q);
      AQ   = Multiply (A,Q);
      QtAQ = Multiply (Qt,AQ);

      // Replace A with (Q transpose) * A * Q.
      Free (A);
      A = Copy (QtAQ);

      // Free up this memory for the next iteration.
      Free (Q); Free(Qt); Free(AQ); Free (QtAQ);

      // See if all off-diagonal entries are small.
      done = 1;
      for (i = 1; i < n && done; i++) {
         for (j = 1; j < i && done; j++) {
            if (fabs(A[i][j]) > ZERO) {
               done = 0;
            }
         }
      }

   }

   // Put the diagonal entries (the eigenvalues) in the row matrix E.
   // Arrange them in decreasing order of magnitude.
   E = Array (1,n);
   for (k = 1; k <= n; k++) {

      lambda_max = 0.0;
      for (i = 1; i <= n; i++) {
         if (A[i][i] > lambda_max) {
            lambda_max = A[i][i];
            i_max = i;
         }
      }

      // lambda_max is the k^th largest eigenvalue.
      E[1][k] = lambda_max;

      // Zero out this diagonal entry so it will not be found again.
      A[i_max][i_max] = 0.0;

   }

   Free (A);

   return E;

}

////////////////////////////////////////////////////////////////////////////////
// This program applies Gram-Schmidt to the column vectors of an n x n matrix
//    A and puts the result into the n x n matrix Q.
////////////////////////////////////////////////////////////////////////////////
double **GramSchmidt (double **A) {

   int i, j, k, n;
   double num, den, scale;
   double **Q;

   // Allocate memory for Q.
   n = Rows (A);
   Q = Array (n,n);

   for (k = 1; k <= n; k++) {

      // Load column k of A into column k of Q.
      for (i = 1; i <= n; i++) {
         Q[i][k] = A[i][k];
      }

      // Project column k of A onto the span of columns 1,...,k-1 of Q, and
      //    subtract that from column k of A (which is also column k of Q).
      for (j = 1; j < k; j++) {

         num = den = 0.0;
         for (i = 1; i <= n; i++) {
            num += Q[i][j] * A[i][k];
            den += Q[i][j] * Q[i][j];
         }

         scale = num / den;

         for (i = 1; i <= n; i++) {
            Q[i][k] -= scale * Q[i][j];
         }

      }

   }


   // Normalize each column.
   for (j = 1; j <= n; j++) {

      scale = 0.0;
      for (i = 1; i <= n; i++) {
         scale += Q[i][j] * Q[i][j];
      }

      if (scale < ZERO) {
         printf ("Matrix A is not of full rank in GramSchmidt().\n");
         Exit ();
      }

      // Scale is the length of column j.
      scale = sqrt(scale);

      // Normalize column j.
      for (i = 1; i <= n; i++) {
         Q[i][j] /= scale;
      }

   }

   return Q;

}

////////////////////////////////////////////////////////////////////////////////
// Compute the eigenvectors of V corresponding to the eigenvalues that are held
//   in E[1][1], E[1][2], ... , E[1][n]. This is done by solving the system
//   of equations V*x = lambda*x, where lambda is the eigenvalue in question.
//   The eigenvectors are put into the columns of the matrix Q.
////////////////////////////////////////////////////////////////////////////////
double **Evector (double **V, double **E) {

   int i, j, i0, j0, n, r, rmax, imax;
   double **M, **Q, c, Qmax;

   // Make sure that V is square and symmetric.
   n = Rows (V);
   if (n != Columns(V)) {
      printf ("In Eigenvalue function V is not square.\n");
      Pause ();
      exit (1);
   }

   for (i = 1; i <= n; i++) {
      for (j = 1; j < i; j++) {
         if (V[i][j] != V[j][i]) {
            printf ("In Eigenvalue function V is not symmetric.\n");
            Pause ();
            exit (1);
         }
      }
   }

   // Allocate working array space.
   M = Array (n,n);
   Q = Array (n,n);


   // Work on eigenvalue number r, which is E[1][r].
   for (r = 1; r <= n; r++) {

      // Subtract the eigenvalue from the diagonal of the matrix V0, put the
      //    result in V.
      for (i = 1; i <= n; i++) {
         for (j = 1; j <= n; j++) {
            M[i][j] = (i == j ? E[1][r] - V[i][j] : -V[i][j]);
         }
      }


      // Use Gaussian elimination to solve the system M*x=0.
      // Work on column j0.
      for (j0 = 1; j0 <= n; j0++) {

         i0 = j0;

         if (i0 == n) {
            break;
         }

         // Find the largest non-zero entry in the column on or below the diagonal.
         c = 0;
         rmax = 0;
         for (i = i0; i <= n; i++) {
            if (fabs(M[i][j0]) > c) {
               c = fabs (M[i][j0]);
               rmax = i;
            }
         }

         // Get out of j0 loop if no such non-zero entry.
         if (c < ZERO) {
            break;
         }

         // Swap rows i0 and rmax if different.
         if (rmax != i0) {
            for (j = 1; j <= n; j++) {
               c = M[i0][j];
               M[i0][j] = M[rmax][j];
               M[rmax][j] = c;
            }
         }


         // Scale so the pivot is 1.
         c = M[i0][j0];
         for (j = 1; j <= n; j++) {
            M[i0][j] /= c;
         }

         // Make rest of column j0 equal to 0.
         for (i = 1; i <= n; i++) if (i != i0) {
            c = M[i][j0];
            for (j = 1; j <= n; j++) {
               M[i][j] += -c * M[i0][j];
            }
         }

      }

      // Put the eigenvector into the r^th column of the matrix Q, and
      //    compute its length. See which component is largest (in abs val)
      //    to give that component a positive sign (this is arbitrary).
      c = 0;
      Qmax = 0.0;
      imax = 1;
      for (i = 1; i <= n; i++) {
         if (i < i0) {
            Q[i][r] = -M[i][j0];
         }
         else if (i == i0) {
            Q[i][r] = 1.0;
         }
         else {
            Q[i][r] = 0.0;
         }
         c += Q[i][r] * Q[i][r];
         if (fabs (Q[i][r]) > Qmax) {
            Qmax = fabs (Q[i][r]);
            imax = i;
         }
      }


      // Normalize the eigenvector making the largest component positive.
      c = (Q[imax][r] > 0 ? 1 : -1) * sqrt(c);
      for (i = 1; i <= n; i++) {
         Q[i][r] /= c;
      }

   } // Next r.

   // We're done with M.
   Free (M);

   return Q;

}

////////////////////////////////////////////////////////////////////////////////
// Compute the matrix product AB.
////////////////////////////////////////////////////////////////////////////////
double **Multiply (double **A, double **B) {

   int i, j, k, nA, mA, nB, mB;
   double **AB;

   mA = Rows (A);
   nA = Columns (A);

   mB = Rows (B);
   nB = Columns (B);

   if (nA != mB) {
      printf ("Dimensions don't match in matrix multiplication.\n");
      Pause ();
      exit(1);
   }

   // Allocate space for the product.
   AB = Array (mA, nB);

   // Compute (AB)_ij.
   for (i = 1; i <= mA; i++) {
      for (j = 1; j <= nB; j++) {
         for (k = 1; k <= nA; k++) {
            AB[i][j] += A[i][k] * B[k][j];
         }
      }
   }

   return AB;

}

////////////////////////////////////////////////////////////////////////////////
// Compute the matrix sum A+B.
////////////////////////////////////////////////////////////////////////////////
double **Add (double **A, double **B) {

   int i, j, n, m;
   double **S;

   m = Rows (A);
   n = Columns (A);

   // Check that the dimensions match.
   if ( (m != Rows (B)) || (n != Columns (B)) ) {
      printf ("Dimensions don't match in matrix addition.\n");
      Exit();
   }

   // Allocate space for the sum.
   S = Array (m, n);

   // Compute (A+B)_ij
   for (i = 1; i <= m; i++) {
      for (j = 1; j <= n; j++) {
         S[i][j] = A[i][j] + B[i][j];
      }
   }

   return S;

}

////////////////////////////////////////////////////////////////////////////////
// Compute the Cholesky decomposition V = LL^T for the symmetric PD matrix V.
////////////////////////////////////////////////////////////////////////////////
double **Cholesky (double **V) {

   int n, i, j, k;
   double **L, sum;

   // Check that the matrix V is square and symmetric.
   n = Rows (V);
    if (n != Columns(V)) {
      printf ("In Cholesky function V is not square.\n");
      Exit();
   }

   for (i = 1; i <= n; i++) {
      for (j = 1; j < i; j++) {
         if (V[i][j] != V[j][i]) {
            printf ("In Cholesky function V is not symmetric.\n");
            Exit();
         }
      }
   }

   // Allocate space for L.
   L = Array (n, n);

   // Work on row k of L.
   for (k = 1; k <= n; k++) {

      // Work on L[k][j]
      for (j = 1; j < k; j++) {
         sum = 0;
         for (i = 1; i < j; i++) {
            sum += L[k][i] * L[j][i];
         }
         L[k][j] = (V[k][j] - sum) / L[j][j];
      }

      // Now compute L[k][k].
      sum = 0.0;
      for (i = 1; i < k; i++) {
         sum += L[k][i] * L[k][i];
      }

      // At this point should have sum < V[k][k].
      if (sum >= V[k][k] + ZERO) {
         printf ("In Cholesky function V is not positive semi-definite.\n");
         Exit ();
      }

      if (Equal (V[k][k], sum, ZERO)) {
         printf ("Warning: In Cholesky function V is positive semi-definite.\n");
         L[k][k] = 0.0;
      }

      else {
         L[k][k] = sqrt (V[k][k] - sum);
      }

   }

   return L;

}

////////////////////////////////////////////////////////////////////////////////
// Print the matrix A to the screen.
////////////////////////////////////////////////////////////////////////////////
void Show (double **A) {

   int i, j, m, n;

   m = Rows (A);
   n = Columns (A);

   for (i = 1; i <= m; i++) {
      for (j = 1; j <= n; j++) {
         printf ("%9.4f", A[i][j]);
      }
      printf ("\n");
   }

   printf ("\n");
   Pause ();
   printf ("\n");

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Print the matrix A to a file.
////////////////////////////////////////////////////////////////////////////////
void Write (double **A, FILE *fp) {

   int i, j, m, n;

   m = Rows (A);
   n = Columns (A);

   for (i = 1; i <= m; i++) {
      for (j = 1; j <= n; j++) {
         fprintf (fp, "%9.4f", A[i][j]);
      }
      fprintf (fp, "\n");
   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Copy the matrix A into another matrix (C).
////////////////////////////////////////////////////////////////////////////////
double **Copy (double **A) {

   double **C;
   int i, j, m, n;

   // Allocate space for C.
   m = Rows (A);
   n = Columns (A);
   C = Array (m, n);

   // Copy A_ij into C_ij.
   for (i = 1; i <= m; i++) {
      for (j = 1; j <= n; j++) {
         C[i][j] = A[i][j];
      }
   }

   return C;

}

////////////////////////////////////////////////////////////////////////////////
// Generate the n x n identity matrix.
////////////////////////////////////////////////////////////////////////////////
double **Identity (int n) {

   double **I;
   int i, j;

   I = Array (n, n);

   for (i = 1; i <= n; i++) {
      for (j = 1; j <= n; j++) {
         I[i][j] = (i == j);
      }
   }

   return I;

}

////////////////////////////////////////////////////////////////////////////////
// Report the number of rows in the matrix A.  This is held in A[0][0].  Note
//   that row 0 of A is not used in any matrix computations.
////////////////////////////////////////////////////////////////////////////////
int Rows (double **A) {

   return (int) A[0][0];

}

////////////////////////////////////////////////////////////////////////////////
// Report the number of columns in the matrix A. Held in A[0][1].
////////////////////////////////////////////////////////////////////////////////
int Columns (double **A) {

   return (int) A[0][1];

}

////////////////////////////////////////////////////////////////////////////////
// Release allocated space for the matrix A.
////////////////////////////////////////////////////////////////////////////////
void Free (double **A) {

   int i, m;

   m = Rows (A);
   for (i = 0; i <= m; i++) {
      free (A[i]);
   }
   free (A);

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Make the columns of an array have mean 0.
////////////////////////////////////////////////////////////////////////////////
double **MeanZero (double **A) {

   double **C, **mu;
   int i, j, m, n;

   // Allocate space for C.
   m = Rows (A);
   n = Columns (A);
   C = Array (m, n);
   mu = Array (1, n);

   // Compute the means of each column.
   for (i = 1; i <= m; i++) {
      for (j = 1; j <= n; j++) {
         mu[1][j] += A[i][j] / m;
      }
   }

   for (i = 1; i <= m; i++) {
      for (j = 1; j <= n; j++) {
         C[i][j] = A[i][j] - mu[1][j];
      }
   }

   Free (mu);

   return C;

}

////////////////////////////////////////////////////////////////////////////////
// Calculate the sample covariance matrix for the data in X0.
////////////////////////////////////////////////////////////////////////////////
double **Covariance (double **X0) {

   int n;
   double **X, **XT, **C, **Cov;

   n = Rows (X0);

   X = MeanZero (X0);

   XT = Transpose(X);

   C = Multiply (XT, X);

   Cov = ScalarMultiple (1.0/(n-1.0), C);

   Free (X);
   Free (XT);
   Free (C);

   return Cov;

}

////////////////////////////////////////////////////////////////////////////////
// Calculate the sample correlation for the data in X0.
////////////////////////////////////////////////////////////////////////////////
double **Correlation (double **X0) {

   int m, n, i, j;
   double s_ii;
   double **X, **XT, **C, **Cor;

   m = Rows (X0);
   n = Columns (X0);

   X = MeanZero (X0);

   XT = Transpose(X);

   C = Multiply (XT, X);

   Cor = ScalarMultiple (1.0/(m-1.0), C);

   for (i = 1; i <= n; i++) {
      if (Cor[i][i] < ZERO) {
         printf ("Zero variance in Correlation function.\n");
         Exit ();
      }
      s_ii = sqrt(Cor[i][i]);
      for (j = 1; j <= n; j++) {
         Cor[i][j] /= s_ii;
         Cor[j][i] /= s_ii;
      }
   }

   Free (X);
   Free (XT);
   Free (C);

   return Cor;

}

////////////////////////////////////////////////////////////////////////////////
// This function computes the implied B-S volatility of a call option given its price.
// tau = time to expiration.
// s   = current stock price.
// k   = option strike.
// r   = risk-free rate as a decimal (e.g. 0.05 for 5%).
// c   = current value of the option.
////////////////////////////////////////////////////////////////////////////////
double ImpliedVol (double tau, double s, double k, double r, double c) {

   int i;
   double step = 0.10, sigma = 0.0, value, value0;

   // First make sure the price "c" is at least the zero-volatility value of
   //   the option.
   value0 = BlackScholes (tau, s, k, 0.0, r);
   if (value0 > c) {
      printf ("Zero-volatility value exceeds option price in ImpliedVol.\n");
      return -1;
   }

   // If the price is the zero-volatility value return sigma = 0;
   if (c < value0 + ZERO) {
      return 0.0;
   }

   // Now handle the general case where sigma > 0.
   value = value0;
   for (i = 1; i <= 5; i++) {

      // Increase the volatility until the B-S value of the call exceeds its price.
      while (value <= c) {

         sigma += step;
         value = BlackScholes (tau, s, k, sigma, r);

      }

      // We now have value > c; went too far so go back one step.
      sigma -= step;
      value = BlackScholes (tau, s, k, sigma, r);

      // Decrease the step size.
      step /= 10;

   }

   return sigma;

}

////////////////////////////////////////////////////////////////////////////////
// The Black-Scholes formula for a call option.  Here tau is T-t.
////////////////////////////////////////////////////////////////////////////////
double BlackScholes (double tau, double s, double k, double sigma, double r) {

   double value, d_plus, d_minus, v, d0;


   // The B-S formula cannot have tau or sigma = 0.
   v = sigma * sqrt(tau);
   if (v < ZERO) {

      // This code handles the cases: (i) tau = 0; (ii) tau > 0 but sigma = 0.
      s *= exp(r*tau);
      value = (s > k ? s-k : 0);
      value *= exp(-r*tau);

   }

   // Otherwise use Black-Scholes.
   else {

      // Compute d+ and d-.
      d0 = (log (s/k) + r*tau) / v;
      d_plus  = d0 + 0.5 * v;
      d_minus = d0 - 0.5 * v;

      // Black-Scholes formula.
      value = s * Psi(d_plus) - k * exp(-r*tau) * Psi(d_minus);

   }

   return (value);

}





