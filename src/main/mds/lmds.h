#ifndef _LMDS_H
#define _LMDS_H

#include <boost/numeric/ublas/fwd.hpp>

// Local MDS, designed for large data sets
// Coded from "A Fast Approximation to MDS"
// Paper by Tynia Yang1, Jinze Liu1, Leonard McMillan1, and Wei Wang1
// University of Chapel Hill at North Carolina, Chapel Hill NC 27599, USA
//
// R is a n x n relational matrix of dissimilarities (assumed to be Euclidean distances)
// W is a n x n weight matrix for each of the dissimlarities
// p is the dimensionality of the target space
// X are the projected points
// s is the number of submatrices to divide R 
void lmds(const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &R,
          const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &W,
          boost::numeric::ublas::matrix<double> &X, 
          unsigned int p,
          double variance_tolerance);

#endif
