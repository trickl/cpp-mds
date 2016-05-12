#ifndef _FDMDS_H
#define _FDMDS_H

#include <boost/numeric/ublas/matrix.hpp>
#include <gsl/gsl_rng.h>

// Force-directed mds
// Like classical mds, scales as O(N^3). Models points as a damped spring system, stable state is lowest kinetic energy
// R is a n x n relational matrix of dissimilarities (assumed to be Euclidean distances)
// p is the dimensionality of the target space 
// X are the projected points
class fdmds
{
public:
   fdmds(gsl_rng* rng,
         unsigned int n,
         unsigned int p = 2,
         bool use_planar_force = true,
         double dampening = 0.5, // Affects number of iterations before system comes to rest
         double force_cap = 0.5, // Should be relative to size of the final solution
         double spring_constant = 0.5 // How many iterations required to correct a distance
         );
   void process_relations(const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &R, unsigned int iteration = 1);
   const boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major>& points() {return X;};

private:          
   boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> X;
   boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> V;
   gsl_rng* _rng;
   const unsigned int _n;
   const unsigned int _p;
   const bool _use_planar_force;
   const double _dampening; // Affects number of iterations before system comes to rest
   const double _force_cap; // Should be relative to size of the final solution
   const double _spring_constant; // How many iterations required to correct a distance
};

#endif
