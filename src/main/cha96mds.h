#ifndef _CHA96MDS_H
#define _CHA96MDS_H

#include <boost/numeric/ublas/fwd.hpp>
#include <gsl/gsl_rng.h>

// MDS algorithm from Matthew Chalmers 1996 paper:
// "A Linear Iteration Time Layour for Visualising High-Dimensional Data"
// A force-directed model that has a fixed size neighbour of points that are measured against
// R is a n x n relational matrix of dissimilarities (assumed to be Euclidean distances)
// p is the dimensionality of the target space 
// X are the projected points
void cha96mds(const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &R,
          boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> &X,
          unsigned int p,
          gsl_rng* rng,
          unsigned int v_max = 5,
          unsigned int s_max = 10,
          bool use_planar_force = true,
          unsigned int max_iterations = 0);
          

#endif
