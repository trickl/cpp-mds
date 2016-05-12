#ifndef _SQPMDS_H
#define _SQPMDS_H

#include <boost/numeric/ublas/fwd.hpp>
#include <gsl/gsl_rng.h>

// MDS algorithm from Chalmers, Morrison and Ross' 2003 paper:
// "Improving Hybrid MDS with Pivot-Based Searching"
// A force-directed model that has a fixed size neighbour of points that are measured against
// R is a n x n relational matrix of dissimilarities (assumed to be Euclidean distances)
// p is the dimensionality of the target space 
// X are the projected points
void sqpmds(const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &R,
          boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> &X,
          unsigned int p,
          gsl_rng* rng);

#endif
