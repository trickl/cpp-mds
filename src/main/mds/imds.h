#ifndef _PMDS_H
#define _PMDS_H

#include <boost/numeric/ublas/fwd.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// Incremental MDS, designed for large data sets
// Like local mds, but done in an incremental fashion
// as opposed to batch

class imds
{
public:
   imds(unsigned int p,
        unsigned int overlap,
        double fit_tolerance,
        unsigned int max_iterations);

   imds(const imds& rhs);

   void operator()(const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major>& dissimilarities,
                   const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major>& weights,
                   boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major>& features,
                   boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major>& affine_mapping);
   
private:
   unsigned int _p;
   unsigned int _overlap;
   double _fit_tolerance;
   unsigned int _max_iterations;
   boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> _features;
};

#endif
