#ifndef _RKISOMAP_H
#define _RKISOMAP_H

#include <boost/numeric/ublas/fwd.hpp>

// Robust Kernel ISOMAP Algorithm
// Implementation based on:
// Heeyoul Choi & Seungjin Choi 2007
// Department of Computer Science, Pohang University of Science and Technology,
// San 31 Hyoja-dong, Nam-gu, Pohang 790-784, Korea  
// This adds a method of removing critical outliers for robustness
// and a generalisation property by treating the Isomap like a Mercer kernel machine
//
// Uses graph shortest-path distances to discover manifolds defined by closest points
// R is a n x n relational matrix of dissimilarities (assumed to be Euclidean distances)
// p is the dimensionality of the target space 
// k is the size of the nearest neighbourhood considered
// X are the projected points
// m is the p-norm distance measure for the graph, 1 = Manhattan, 2 = Euclidean
// flow_interolance prevents too much flow going through a single node, the lower the intolerance, the more likely for short-circuits

void rkisomap(const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &R,
              boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &S,
              unsigned int k,
              double m = 2.,
              double flow_intolerance = 0.01);

#endif
