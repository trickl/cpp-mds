#ifndef _SMDS_H
#define _SMDS_H

#include <boost/numeric/ublas/fwd.hpp>

// Sparse (classical) MDS with a large matrix eigensolver for performance
// Credit for the IETL library used for the eigensolver goes to:
// Prakash Dayal <prakash@comp-phys.org>
// Matthias Troyer <troyer@comp-phys.org>
//
// R is a n x n relational matrix of dissimilarities (assumed to be Euclidean distances)
// max_p is the max dimensionality of the target space
// p_tolerance is the amount of variance to be explained by the target space
// X are the projected points
void smds(const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &R,
          boost::numeric::ublas::matrix<double,  boost::numeric::ublas::column_major> &X,
          unsigned int max_p,
          double variance_tolerance);

#endif
