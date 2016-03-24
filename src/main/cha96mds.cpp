#include "cha96mds.h"
#include <stdexcept>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

using namespace std;
using namespace boost::numeric::ublas;

void cha96mds(const symmetric_matrix<double, lower, column_major> &R,
          matrix<double, column_major> &X,
          unsigned int p,
          gsl_rng* rng,
          unsigned int v_max,
          unsigned int s_max,
          bool use_planar_force,
          unsigned int max_iterations)
{
   if (R.size1() != R.size2()) throw invalid_argument("Relation matrix R must be square.");

   // First need to construct double centred distance matrix B of scalar products
   const unsigned int n = R.size1();
   const unsigned int p_interim = use_planar_force ? p + 1 : p;

   // Calculate initial positions and velocities for particles
   if (X.size1() == 0)
   {
      X.resize(n, p_interim);
      for (unsigned int i = 0; i < n; ++i)
      {
         for (unsigned int k = 0; k < p_interim; ++k)
         {
            X(i, k) = gsl_rng_uniform(rng);
         }
      }
   }

   // Neighbour comparison matrices
   matrix<unsigned int> V = zero_matrix<unsigned int>(n, v_max);
   matrix<unsigned int> S = zero_matrix<unsigned int>(n, s_max);
   boost::numeric::ublas::vector<unsigned int> VS_i = zero_vector<unsigned int>(v_max + s_max);
   boost::numeric::ublas::vector<unsigned int> v_size = zero_vector<unsigned int>(n);

   matrix<double, column_major> Vl = zero_matrix<double>(n, p_interim);
   const double dampening = 0.3; // Affects number of iterations before system comes to rest
   const double force_cap = 0.5; // Should be relative to size of the final solution
   const double spring_constant = 0.2; // How many iterations required to correct a distance

   if (max_iterations == 0) max_iterations = n;
   for (unsigned int iteration = 0; iteration < max_iterations; ++iteration)
   {

      double kinetic_energy = 0.;
      for (unsigned int i = 0; i < n; ++i)
      {
         // Find the neighbour sets of i by stochastic sampling
         for (unsigned int s_size = 0; s_size < s_max;)
         {
            // Choose a potential neighbour for S randomly 
            unsigned int j =  gsl_rng_uniform_int(rng, n);

            if (i != j)
            {
               // Check j is not in V already
               unsigned int v_j = 0; 
               unsigned int v_max_j = 0;
               double v_max_distance = v_size(i) ? 0 : numeric_limits<double>::infinity();
               while (v_j < v_size(i) && V(i, v_j) != j)
               {
                  // While checking j is not in V, find the max distance in V.
                  if (R(i, V(i, v_j)) >= v_max_distance)
                  {
                     v_max_distance = R(i, V(i, v_j));
                     v_max_j = v_j; 
                  }
                  ++v_j;
               } 

               if (v_j == v_size(i))
               { 
                  if (v_size(i) < v_max)
                  {
                     // Insert into V instead
                     V(i, v_size(i)) = j;
                     v_size(i) = v_size(i) + 1;
                  }
                  else if (R(i, j) < v_max_distance)
                  {
                     // As V is full, push out the last element into S
                     S(i, s_size) = V(i, v_max_j);
                     s_size++;

                     V(i, v_max_j) = j;
                  }
                  else
                  {
                     // Add neighbour to S
                     S(i, s_size) = j; 
                     s_size++;
                  } 
               }
            }
         }

         // Calculate the force on each particle due to all it's VS_i 
         boost::numeric::ublas::vector<double> net_force = zero_vector<double>(p_interim);

         // Comparisons will be made between the points in V and S
         for (unsigned int v_j = 0; v_j < v_max; ++v_j) VS_i(v_j) = V(i, v_j);
         for (unsigned int s_j = 0; s_j < s_max; ++s_j) VS_i(v_max + s_j) = S(i, s_j);

         for (unsigned int k = 0; k < (v_max + s_max); ++k)
         {
            unsigned int j = VS_i(k);   

            // Force is proportional to difference between higher and lower dimensional distances
            boost::numeric::ublas::vector<double> dx = matrix_row<matrix<double, column_major> >(X, i) - 
                                                       matrix_row<matrix<double, column_major> >(X, j);
            double distance = sqrt(inner_prod(dx, dx));
            double force_size = spring_constant * (R(i, j) - distance);

            net_force += dx * force_size / (distance * double(v_max + s_max));
         }

         // Additional force towards the plane, increases as iterations progress
         if (use_planar_force)
         {
            net_force(p) += -X(i, p) * spring_constant * double(iteration) / double(max_iterations);
         } 

         // Cap the force magnitude
         double net_force_size = sqrt(inner_prod(net_force, net_force));
         net_force = net_force * (min(net_force_size, force_cap) / net_force_size);

         matrix_row<matrix<double, column_major> >(Vl, i) = matrix_row<matrix<double, column_major> >(Vl, i) * (1 - dampening)
                                                           + net_force;
         matrix_row<matrix<double, column_major> >(X, i) += matrix_row<matrix<double, column_major> >(Vl, i);

         kinetic_energy += inner_prod(matrix_row<matrix<double, column_major> >(Vl, i),
                                      matrix_row<matrix<double, column_major> >(Vl, i));
      }

      cout << "Iteration(" << iteration << "), Kinetic Energy = " << kinetic_energy << endl;
   }
 
   // Remove the extra degree of freedom
   X.resize(n, p);
}
