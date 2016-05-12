#ifndef _ISOMAP_H
#define _ISOMAP_H

#include <boost/numeric/ublas/fwd.hpp>

// J. B. Tenenbaum, V. de Silva, and J. C. Langford ISOMAP Algorithm
// C++ Version built using Josh Tenenbaum's MatLab version as reference
//
// Uses graph shortest-path distances to discover manifolds defined by closest points
// R is a n x n relational matrix of dissimilarities (assumed to be Euclidean distances)
// p is the dimensionality of the target space 
// k is the size of the nearest neighbourhood considered
// X are the projected points
void isomap(const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &R,
            boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &S,
            unsigned int k);

#endif
