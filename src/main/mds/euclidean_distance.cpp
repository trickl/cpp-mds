#include "euclidean_distance.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

using namespace std;
using namespace boost::numeric::ublas;

void euclidean_distance(const matrix<double, column_major> &X,
                        symmetric_matrix<double, lower, column_major> &R)
{
   for (matrix<double, column_major>::const_iterator1 i_itr = X.begin1(), i_end = X.end1(); i_itr != i_end; ++i_itr)
   {
      for (matrix<double, column_major>::const_iterator1 j_itr = i_itr + 1, j_end = X.end1(); j_itr != j_end; ++j_itr)
      {
         long i = i_itr.index1(), j = j_itr.index1();
         double d2 = 0;
         for (matrix<double, column_major>::const_iterator2 k_itr = j_itr.begin(), k_end = j_itr.end(); k_itr != k_end; ++k_itr)
         {
            long k = k_itr.index2();
            d2 += pow(X(i, k) - X(j, k), 2.);
         }
      
         R(i, j) = sqrt(d2);
      }
   }
}

