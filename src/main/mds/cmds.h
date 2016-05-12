#ifndef _CMDS_H
#define _CMDS_H

#include <boost/numeric/ublas/fwd.hpp>

// Torgerson Classical MDS (1952)
// R is a n x n relational matrix of dissimilarities (assumed to be Euclidean distances)
// p is the dimensionality of the target space 
// X are the projected points
void cmds(const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &R,
          unsigned int p,
          boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> &X);

#endif
