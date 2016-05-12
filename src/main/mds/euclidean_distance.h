#ifndef _EUCLIDEAN_DISTANCE_H
#define _EUCLIDEAN_DISTANCE_H

#include <boost/numeric/ublas/fwd.hpp>

// Convert a feature matrix into a Euclidean distance dissimilarity matrix
void euclidean_distance(const boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> &X,
                        boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &R);

#endif
