#include "cmds.h"
#include <stdexcept>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>

using namespace std;
using namespace boost::numeric::ublas;
using namespace boost::numeric::bindings;

void cmds(const symmetric_matrix<double, lower, column_major> &R,
          unsigned int p,
          matrix<double, column_major> &X)
{
   if (R.size1() != R.size2()) throw invalid_argument("Relation matrix R must be square.");

   // First need to construct double centred distance matrix B of scalar products
   const int n = R.size1();

   // Need d^2 matrix where d is a Euclidean distance
   symmetric_matrix<double, lower, column_major> D = R;
   for (long i = 0; i < n; ++i)
   {
      for (long j = 0; j < n; ++j)
      {
        D(i, j) = pow(R(i, j), 2.0);
      } 
   }
 
   // Double centring, so origin is centroid of all points
   symmetric_matrix<double, lower, column_major> J(n, n); // Double centering matrix
   identity_matrix<double> I(n);  
   for (long i = 0; i < n; ++i)
   {
      for (long j = 0; j < n; ++j)
      {
         J(i, j) = I(i, j) - (1.0 / double(n));
      }
   }
   matrix<double, column_major> B = -0.5 * prec_prod(J, prec_prod<matrix<double, column_major> >(D, J));

   // Use svd to find B = ELE'
   matrix<double, column_major> E(B.size1(), B.size1()), ET(B.size2(), B.size2());
   boost::numeric::ublas::vector<double> L(B.size1());
   lapack::gesvd('S', 'S', B, L, E, ET);

   // Can now use B to solve for the co-ordinates X = EL^0.5 (XX' = B)
   X.resize(R.size1(), p);
   for (long i = 0; i < p; i++)
   {
      matrix_column<matrix<double, column_major> >(X, i) = matrix_column<matrix<double, column_major> >(E, i) * sqrt(L(i));
   }
}
