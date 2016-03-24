#ifndef _STRESS_H
#define _STRESS_H

#include <boost/numeric/ublas/fwd.hpp>

// Calculate a stress measure between projected points in a reduced dimensional space and 
// a dissimilarity matrix
enum stress_type
{
   raw,
   normalised,
   kruskal // (stress-1) 
};

double calculate_stress(const boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> &X,
                        const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &R,
                        stress_type type);

double calculate_stress(const boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> &X,
                        const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &R,
                        const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &W,
                        stress_type type);

#endif
