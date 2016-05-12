#include "smacof_1.h"
#include <stdexcept>
#include <climits>
#include "stress.h"
#include "cmds.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

using namespace std;
using namespace boost::numeric::ublas;

void smacof_1(const symmetric_matrix<double, lower, column_major> &R,
              int p,
              matrix<double, column_major> &X,
              double stress_tolerance,
              unsigned int max_iterations)
{
   if (R.size1() != R.size2()) throw invalid_argument("Relation matrix R must be square.");
   const unsigned int n = R.size1(); // Number of data points

   // Initially, use Torgerson Classical MDS to get initial values for X
   cmds(R, p, X);

   // Set stress to be very large (don't bother calculating it after cmds).
   double stress = numeric_limits<double>::infinity();

   // Iteration to refine X
   symmetric_matrix<double, lower, column_major> d(n);
   for (int iteration = 0; iteration < max_iterations; ++iteration) 
   {
      // 3. Calculate distances o(X). 
      for (unsigned int i = 0; i < X.size1(); i++)
      {
         for (unsigned int j = 0; j < X.size1(); j++)
         {
            d(i, j) = 0;
            for (unsigned int k = 0; k < X.size2(); k++)
            {
               // Distance function can be non metric
               d(i, j) += pow(X(i, k) - X(j, k), 2.0);
            }
            d(i, j) = sqrt(d(i, j));
         }
      }

      // Compute the Guttman transform G of X. 
      // (This is the solution to the majorizing function)
      symmetric_matrix<double, lower, column_major> G(n); 
      for (unsigned int i = 0; i < X.size1(); i++)
      {
         double sumg = 0;
         for (unsigned int j = 0; j < X.size1(); j++)
         {
            if (d(i, j) == 0) G(i, j) = 0;
            else if (i != j)
            {
               G(i, j) = -R(i, j) / d(i, j);
            }

            sumg -= G(i, j);
         }

         G(i, i) = sumg;
      }

      X = prec_prod(G, X) / double(X.size1());

      // Stop if the change in stress is too small
      double stress_prev = stress;
      stress = calculate_stress(X, R, normalised);
      if (abs(stress_prev - stress) / stress_prev < stress_tolerance) break;
   }
}
