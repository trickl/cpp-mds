#include "wmds.h"
#include <stdexcept>
#include <climits>
#include <mds/smds.h>
#include <mds/stress.h>
#include <mds/guttman.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/bindings/atlas/cblas.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_symmetric.hpp> 

using namespace std;
using namespace boost::numeric::bindings;
using namespace boost::numeric::ublas;

void wmds(const symmetric_matrix<double, lower, column_major> &R,
          const symmetric_matrix<double, lower, column_major> &W,
          matrix<double,  column_major> &X,
          unsigned int max_p,
          double variance_tolerance,
          double fit_tolerance,
          unsigned int max_iterations)
{
   if (R.size1() != W.size1()) throw invalid_argument("Weights W must have same dimensions are R.");
   unsigned int n = R.size1(); // Number of data points

   // Initially, use metric MDS to get initial values for X
   smds(R, X, max_p, variance_tolerance);

   // Calculate V 
   boost::numeric::ublas::vector<double> Vdiag(n);

   // Set diagonal to zero
   for (unsigned int i = 0; i < n; ++i) Vdiag(i) = 0.;
   for (symmetric_matrix<double, lower, column_major>::const_iterator1 i_itr = W.begin1(), i_end = W.end1(); i_itr != i_end; ++i_itr)
   {
      for (symmetric_matrix<double, lower, column_major>::const_iterator2 j_itr = i_itr.begin(), j_end = i_itr.end(); j_itr != j_end; ++j_itr)
      {
         unsigned int i = i_itr.index1(), j = j_itr.index2();
         if (i > j) // skip the diagonal, only consider lower triangle
         {
            double w = W(i, j);
            Vdiag(i) += w;
            Vdiag(j) += w;
         }
      }
   }

   // Iteration to refine X
   for (unsigned int iteration = 0; iteration < max_iterations; ++iteration) 
   {
      // DMA iterative improvement
      matrix<double, column_major> G(n, n);
      guttman_transform(X, R, W, G);
      matrix<double, column_major> dX(n, X.size2());
      atlas::symm(symmetric_adaptor<matrix<double, column_major>, lower>(G), X, dX);

      double max_dX = 0;
      double max_X = 0;
      for (unsigned int i = 0; i < n; ++i)
      {
         matrix_row<matrix<double,  column_major> > dX_row(dX, i);
         dX_row /= (Vdiag(i) * 2);

         matrix_row<matrix<double,  column_major> > X_row(X, i);
         X_row += dX_row;

         max_dX = max(max_dX, *max_element(dX_row.begin(), dX_row.end()));
         max_X  = max(max_X,  *max_element(X_row.begin(), X_row.end()));
      }

      // Stop if max(dX) is too small
      // cout << "Tolerance: " << fit_tolerance << ", max_dX / max_X: " << max_dX / max_X << endl;
      if (max_dX / max_X < fit_tolerance) break;
   }
}
