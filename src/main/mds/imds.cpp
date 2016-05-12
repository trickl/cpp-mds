#include "imds.h"
#include <mds/wmds.h>
//#include <mds/dgels.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/lapack/gesv.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>

using namespace std;
using namespace boost;
using namespace boost::numeric::bindings;
using namespace boost::numeric::ublas;
using boost::numeric::ublas::range;

imds::imds(unsigned int p, 
          unsigned int overlap,
          double fit_tolerance, 
          unsigned int max_iterations)
   : _p(p),
     _overlap(overlap),
     _fit_tolerance(fit_tolerance),
     _max_iterations(max_iterations)
{
}

imds::imds(const imds& rhs)
   : _p(rhs._p),
     _overlap(rhs._overlap),
     _fit_tolerance(rhs._fit_tolerance),
     _max_iterations(rhs._max_iterations),
     _features(rhs._features)
{
}

// Run MDS on the submatrix for this iteration
void imds::operator()(const symmetric_matrix<double, lower, column_major>& dissimilarities,
                      const symmetric_matrix<double, lower, column_major>& weights,
                      matrix<double, column_major>& features,
                      matrix<double, column_major>& affine_mapping)
{  
   // Variance tolerance is zero, as dimensionality is fixed
   wmds(dissimilarities, weights, features, _p, 0, _fit_tolerance, _max_iterations);

   // affine_mapping must be square (no dimensionality change)
   if (_features.size1() == 0)
   {
      affine_mapping = identity_matrix<double>(_p + 1, _p);
      _features = zero_matrix<double>(_overlap, _p + 1);
   }
   else
   {
      // Find features for new common points
      affine_mapping = matrix_range<matrix<double, column_major> >(features, range(0, _overlap), range(0, _p));

      // Solve the generalised least squares problem
      // old mds * A = new mds
      matrix<double, column_major> old_features = _features;
      //dgels(_features, affine_mapping);      
   }

   // Update the stored features and keys
   matrix_range<matrix<double, column_major> >(_features, range(0, _overlap), range(0, _p)) = matrix_range<matrix<double, column_major> >(features, range(features.size1() - _overlap, features.size1()), range(0, _p));
   matrix_column<matrix<double, column_major> >(_features, _p) = scalar_vector<double>(_overlap, 1.);
}


