////////////////////////////////////////////////////////////////////////////////
// This program generates 5 sample stock price paths using a GARCH(1,1)
//    stochastic volatility model. Results are printed to the screen. The
//    current volatility is 35% and the long-run volatility is 30% annually.
////////////////////////////////////////////////////////////////////////////////

#include "Functions.h"

int main () {

   int i, n, N;
   double s_start, s, T, r, sigma, sigma2,  sigma2_start, mu, alpha, beta,
          gamma, dt, Z, R, U, sigma2_LT, call, tildes, PV, tildePV, Xbar,
          X2bar, error, X;

   // Specify the alpha, beta, and gamma parameters.
   n = GetInteger ("Which problem (2, 3 ,or 4) are you working on?... ");
   if (n == 2) {
      alpha = 0.0; beta = 0.0; gamma = 1.0;
   }
   else if (n == 3) {
      alpha = 0.01; beta = 0.0; gamma = 0.99;
   }
   else { // Problems 1 and 4.
      alpha = 0.01; beta = 0.10; gamma = 0.89;
   }

   // Time to expiration.
   T = 0.5;

   // Risk-free interest rate.
   r = 0.05;

   // Number of days to expiration (252 days per year), and one day in years.
   N = 252 * T;
   dt = 1.0/252.0;

   // Convert "r" to a daily risk-free rate.
   r *= dt;

   // Current stock price.
   s_start = 100.0;

   // Current stock price variance.
   sigma2_start = 0.35 * 0.35;

   // Current daily stock price variance.
   sigma2_start *= dt;

   // Annual long-term variance.
   sigma2_LT = 0.30 * 0.30;

   // Daily long-term variance.
   sigma2_LT *= dt;

   // Seed the RNG.
   MTUniform (1);

   //set estimator variables to equal 0 
    Xbar = X2bar = 0; 

   //headers 
   printf("     N       Call Option\n");
   printf(" =========   ============\n");
   // Generates at most 10 million paths.
   for (n = 1; n <= 10000000; n++) {

      // Initialize the stock price and the volatility to their starting (time 0)
      //    values.
      s = tildes = s_start;
      sigma2 = sigma2_start;

      // Now march forward in time day by day.
      for (i = 1; i <= N; i++) {

         // Compute the drift parameter that corresponds to the volatility for
         //    this period.
         mu = r - 0.5*sigma2;

         // Compute the stock price at day "i"...

         // First get a standard normal Z.
         U = MTUniform (0);
         Z = PsiInv (U);

         // Apply current volatility.
         sigma = sqrt(sigma2);
         R = sigma * Z;

         // Update the stock price & antithetic stock price
         s *= exp (mu + R);
         tildes *= exp(mu - R); 

         if (s < 0.0)
            s = 0.0;
         if (tildes <0.0)
            tildes =0.0; 

         // Update the stochastic volatility according to the GARCH(1,1) model.
         sigma2 = alpha * sigma2_LT  +  beta * R*R +  gamma * sigma2;
   
      }

      // Calculate the discounted terminal stock price, and antithetic.
         PV  = exp(-.05*T) * s;
         tildePV = exp(-.05*T) * tildes;

         // Average the original present value and the antithetic present value.
         X = (PV + tildePV) / 2.0;

         

         // Update the sample moments;
         Xbar = (Xbar * (n-1) + X) / n;
         X2bar = (X2bar * (n-1) + X*X) / n;

         //prints call option value every 100,000 simulations 
         //and checks if error tolerance is met 
          if (n % 100000 ==0 ){
            printf ("%9i    %7.3f \n", n,  Xbar);
            error = 1.96 * sqrt((X2bar-Xbar*Xbar)/n);
            if (error <= .01)
               break; 

         }

   } // Next path.


      // Report the results.
      printf ("Expected present value is %6.2f +/- %6.2f ", Xbar, error);
      printf ("with 95%% confidence.\n");
   

}


