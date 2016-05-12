#ifndef _SMACOF_H
#define _SMACOF_H

#include <boost/numeric/ublas/fwd.hpp>

// de Leeuw SMACOF-1 (1977) combines metric/nonmetric mds/wmds
// Scaling by MAjorizing a COmplicated Function
// R is a n x n relational matrix of dissimilarities (assumed to be Euclidean distances)
// W is a n x n weight matrix for each of the dissimlarities
// p is the dimensionality of the target space 
// X are the projected points
// stress_tolerance is the smallest ratio change in stress that continues the iterations
// max_iterations is the maximum number of iterations
void smacof_1(const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &R,
              int p,
              boost::numeric::ublas::matrix<double,  boost::numeric::ublas::column_major> &X,
              double stress_tolerance = 0.0001,
              unsigned int max_iterations = 100);

#endif
