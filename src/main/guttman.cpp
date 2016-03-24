#include "guttman.h"
#include "euclidean_distance.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

using namespace boost::numeric::ublas;

void guttman_transform(const matrix<double, column_major> &X,
                       const symmetric_matrix<double, lower, column_major> &R,
                       const symmetric_matrix<double, lower, column_major> &W,
                       matrix<double, column_major> &G)
{
   const unsigned int n = R.size1(); // Number of data points

   // Calculate actual distances in feature matrix X
   symmetric_matrix<double, lower, column_major> d(n);
   euclidean_distance(X, d);

   // Compute the weighted Guttman transform G of X. 
   // (This is the solution to the majorizing function)
   for (unsigned int i = 0; i < n; ++i) G(i, i) = 0.;
   for (matrix<double, column_major>::const_iterator1 i_itr = G.begin1(), i_end = G.end1(); i_itr != i_end; ++i_itr)
   {
      for (matrix<double, column_major>::const_iterator2 j_itr = i_itr.begin(), j_end = i_itr.end(); j_itr != j_end; ++j_itr)
      {
         long i = i_itr.index1(), j = j_itr.index2();
         if (i > j) // skip the diagonal, only consider lower triangle
         {
            double g = d(i, j) != 0 ? R(i, j) * W(i, j) / d(i, j) : 0;
            g -= W(i, j);
            G(i, j) = -g;
            G(i, i) += g;
            G(j, j) += g;
         }
      }
   }
}
