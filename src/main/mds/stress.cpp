#include "stress.h"
#include <stdexcept>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

using namespace std;
using namespace boost::numeric::ublas;

double calculate_stress(const matrix<double, column_major> &X,
                        const symmetric_matrix<double, lower, column_major> &R,
                        stress_type type)
{
   // Construct elementary weights
   symmetric_matrix<double, lower, column_major> W(R.size1());
   for (symmetric_matrix<double, lower, column_major>::iterator1 i = W.begin1(), i_end = W.end1(); i != i_end; ++i) fill(i.begin(), i.end(), 1.);
   return calculate_stress(X, R, W, type);
}

double calculate_stress(const matrix<double, column_major> &X,
                        const symmetric_matrix<double, lower, column_major> &R,
                        const symmetric_matrix<double, lower, column_major> &W,
                        stress_type type)
{
   if (X.size1() != R.size1()) throw invalid_argument("X must have the same number of points as R."); 
   if (W.size1() != W.size1()) throw invalid_argument("W must have the same dimensions as R."); 
   if (W.size2() != W.size2()) throw invalid_argument("W must have the same dimensions as R."); 

   double stress = 0;
   double normalised_denominator = 0;
   double kruskal_denominator = 0;

   for (symmetric_matrix<double, lower, column_major>::const_iterator1 i_itr = R.begin1(), i_end = R.end1(); i_itr != i_end; ++i_itr)
   {
      for (symmetric_matrix<double, lower, column_major>::const_iterator2 j_itr = i_itr.begin(), j_end = i_itr.end(); j_itr != j_end; ++j_itr)
      {
         long i = i_itr.index1(), j = j_itr.index2();
         if (i > j)
         {
            // Actual distance
            double d2 = 0;
            for (matrix<double, column_major>::const_iterator2 k_itr = X.begin2(), k_end = X.end2(); k_itr != k_end; ++k_itr)
            {
               long k = k_itr.index2();
               d2 += pow(X(i, k) - X(j, k), 2.);
            }

            stress += pow(sqrt(d2) - R(i, j), 2.0) * W(i, j);
            normalised_denominator += pow(R(i, j), 2.0) * W(i, j);
            kruskal_denominator += d2 * W(i, j);
         }
      }
   }

   if (type == normalised) stress /= normalised_denominator;
   if (type == kruskal) stress /= kruskal_denominator;
   if (type == kruskal) stress = sqrt(stress);

   return stress;
}





