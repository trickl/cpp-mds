#include "fdmds.h"
#include <stdexcept>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

using namespace std;
using namespace boost::numeric::ublas;

fdmds::fdmds(gsl_rng* rng, 
         unsigned int n,
         unsigned int p,
         bool use_planar_force,
         double dampening,
         double force_cap,
         double spring_constant)
  : _rng(rng), _n(n), _p(use_planar_force ? p + 1 : p),
    _use_planar_force(use_planar_force),
    _dampening(dampening),
    _force_cap(force_cap),
    _spring_constant(spring_constant)
{
   // Calculate initial positions and velocities for particles
   X = zero_matrix<double>(_n, _p);
   for (unsigned int i = 0; i < _n; ++i)
   {
      for (unsigned int j = 0; j < _p; ++j)
      {
         X(i, j) = gsl_rng_uniform(rng);
      }
   }

   V = zero_matrix<double>(_n, _p);
}

void fdmds::process_relations(const symmetric_matrix<double, lower, column_major> &R, unsigned int iteration)
{
   if (R.size1() != R.size2()) throw invalid_argument("Relation matrix R must be square.");

   // First need to construct double centred distance matrix B of scalar products
   
   double kinetic_energy = 0.;
   // Calculate the force on each particle due to all it's neighbours 
   for (unsigned int i = 0; i < _n; ++i)
   {
      boost::numeric::ublas::vector<double> net_force = zero_vector<double>(_p);
      for (unsigned int j = 0; j < _n; ++j)
      {
         if (i != j) 
         {
            // Force is proportional to difference between higher and lower dimensional distances
            boost::numeric::ublas::vector<double> dx = matrix_row<matrix<double, column_major> >(X, i) - 
                                                    matrix_row<matrix<double, column_major> >(X, j);
            double distance = sqrt(inner_prod(dx, dx));
            double force_size = _spring_constant * (R(i, j) - distance);

            net_force += dx * force_size / (distance * double(_n - 1));
         } 
      }

      // Additional force towards the plane, increases as iterations progress
      if (_use_planar_force)
      {  
         net_force(_p - 1) += -X(i, _p) * _spring_constant * double(iteration) / double(_n);
      }

      // Cap the force magnitude
      double net_force_size = sqrt(inner_prod(net_force, net_force));
      net_force = net_force * (min(net_force_size, _force_cap) / net_force_size);

      matrix_row<matrix<double, column_major> >(V, i) = matrix_row<matrix<double, column_major> >(V, i) * (1 - _dampening)
                                                           + net_force;
      matrix_row<matrix<double, column_major> >(X, i) += matrix_row<matrix<double, column_major> >(V, i);

      kinetic_energy += inner_prod(matrix_row<matrix<double, column_major> >(V, i),
                                   matrix_row<matrix<double, column_major> >(V, i));
   }

   cout << "Kinetic Energy = " << kinetic_energy << endl;
}
