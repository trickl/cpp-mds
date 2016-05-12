#ifndef _SQMDS_H
#define _SQMDS_H

#include <boost/numeric/ublas/fwd.hpp>
#include <gsl/gsl_rng.h>

// MDS algorithm from Chalmers, Morrison and Ross' 2002 paper:
// "A Hybrid Layout Algorithm for Sub-Quadratic Multidimensional Scaling"
// A force-directed model that has a fixed size neighbour of points that are measured against
// R is a n x n relational matrix of dissimilarities (assumed to be Euclidean distances)
// p is the dimensionality of the target space 
// X are the projected points
void sqmds(const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &R,
          boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> &X,
          unsigned int p,
          gsl_rng* rng);

#endif
