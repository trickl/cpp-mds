#ifndef _WMDS_H
#define _WMDS_H

#include <boost/numeric/ublas/fwd.hpp>

// Weighted MDS with Diagonal Majorizing Algorithm
// Coded from "Multidimensional Scaling Algorithms for Large Data Sets"
// Paper written by Michael W. Trosset & Patrick J.F. Groenen
//
// R is a n x n relational matrix of dissimilarities (assumed to be Euclidean distances)
// W is a n x n weight matrix for each of the dissimlarities
// X are the projected points
// max_p is the max dimensionality of the target space
void wmds(const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &R,
          const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &W,
          boost::numeric::ublas::matrix<double,  boost::numeric::ublas::column_major> &X,
          unsigned int max_p,
          double variance_tolerance, 
          double fit_tolerance = 1e-4,
          unsigned int max_iterations = 100);

#endif
