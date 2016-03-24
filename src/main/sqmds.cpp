#include "sqmds.h"
#include <mds/cha96mds.h>
#include <mds/stress.h>
#include <stdexcept>
#include <set>
#include <algorithm>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

using namespace std;
using namespace boost::numeric::ublas;

void sqmds(const symmetric_matrix<double, lower, column_major> &R,
          matrix<double, column_major> &X,
          unsigned int p,
          gsl_rng* rng)
          
{
   if (R.size1() != R.size2()) throw invalid_argument("Relation matrix R must be square.");

   // First need to construct double centred distance matrix B of scalar products
   const unsigned int n = R.size1();
   const unsigned int sample_size = 10 + sqrt(n);
   const unsigned int interp_size = 3 + sqrt(sample_size);
   const unsigned int v_max = 5;
   const unsigned int s_max = 10;
   const unsigned int max_search_iterations = 7;
   const double spring_constant = 1.0; // How many iterations required to correct a distance

   std::vector<unsigned int> samples; 
   set<unsigned int> sample_set; 
   set<unsigned int> interp_subset;

   // Select an O(N^1/2) sample of the original subset
   while (sample_set.size() < sample_size)
   {
      // No repeats and not the nearest neighbour x
      unsigned int j = gsl_rng_uniform_int(rng, n);
      if (sample_set.find(j) == sample_set.end()) sample_set.insert(j);
   } 
   copy(sample_set.begin(), sample_set.end(), back_inserter(samples));

   symmetric_matrix<double, lower, column_major> S(sample_size);
   for (unsigned int i = 0; i < sample_size; ++i) 
   {
      for (unsigned int j = 0; j < sample_size; ++j) 
      { 
         S(i, j) = R(samples[i], samples[j]);
      }
   } 

   // Start by running the Chalmer's 96 algorithm on the small subset
   matrix<double, column_major> XS;
   cha96mds(S, XS, p, rng, v_max, s_max);

   // Map this smaller result set into the full data set 
   X.resize(n, p); 
   for (unsigned int i = 0; i < sample_size; ++i) 
   {
      matrix_row<matrix<double, column_major> >(X, samples[i]) = matrix_row<matrix<double, column_major> >(XS, i);
   }
    
   // Now for each remaing point in the set, interpolate the position
   for (unsigned int i = 0; i < n; ++i)
   {
      if (!binary_search(samples.begin(), samples.end(), i))
      { 
         // Find the nearest neighbour to i in S (brute force)
         unsigned int x = 0; 
         double nearest_distance = numeric_limits<double>::infinity();
         for (unsigned int j = 0; j < sample_size; ++j) 
         {
            double distance = R(i, samples[j]);
            if (distance < nearest_distance) 
            {
               nearest_distance = distance;
               x = samples[j];
            }
         }

         // Select a random sample of the original subset S 
         interp_subset.clear();
         while (interp_subset.size() < interp_size)
         {
            // No repeats and not the nearest neighbour x
            unsigned int j = gsl_rng_uniform_int(rng, sample_size);
            if (j != x && interp_subset.find(samples[j]) == interp_subset.end()) interp_subset.insert(samples[j]);
         } 

         // Use a binary search to refine the place of i
         double angle = 0;
         double search_angle = M_PI / 2.;
         double radius = nearest_distance; 
         for (unsigned int search_iteration = 0; search_iteration < max_search_iterations; ++search_iteration)
         {
            // Find the stress of both configurations
            double stress = 0;
            for (unsigned int configuration = 0; configuration < 2; ++configuration)
            { 
               for (set<unsigned int>::const_iterator itr = interp_subset.begin(), end = interp_subset.end(); itr != end; ++itr) 
               {
                  unsigned int j = *itr;
                  boost::numeric::ublas::vector<double> direction_vector(2);

                  // TODO: Fix this polar search algorithm for the general case of p dimensions
                  if (p != 2) throw invalid_argument("Polar search algorithm can only handle 2 dimensions");

                  double query_angle = configuration ? angle - search_angle : angle + search_angle;
                  direction_vector[0] = radius * cos(query_angle);
                  direction_vector[1] = radius * sin(query_angle);

                  matrix_row<matrix<double, column_major> >(X, i) = matrix_row<matrix<double, column_major> >(X, x) +
                                                                 direction_vector; 
                  boost::numeric::ublas::vector<double> dx = matrix_row<matrix<double, column_major> >(X, i) -
                                                          matrix_row<matrix<double, column_major> >(X, j);
                  double distance = sqrt(inner_prod(dx, dx));
                  stress += pow(R(i, j) - distance, 2.);
               } 

               stress = -stress;
            }

            // Change the angle to the lower stress configuration 
            angle = stress < 0 ? angle + search_angle : angle - search_angle;
            search_angle /= 2;
         }
         
         // Use a force directed approach to refine i over a constant number of iterations
         boost::numeric::ublas::vector<double> net_force = zero_vector<double>(p);
         for (set<unsigned int>::const_iterator itr = interp_subset.begin(), end = interp_subset.end(); itr != end; ++itr) 
         {
            unsigned int j = *itr;
   
            // Force is proportional to difference between higher and lower dimensional distances
            boost::numeric::ublas::vector<double> dx = matrix_row<matrix<double, column_major> >(X, i) - 
                                                    matrix_row<matrix<double, column_major> >(X, j);
            double distance = sqrt(inner_prod(dx, dx));
            double force_size = spring_constant * (R(i, j) - distance);

            net_force += dx * force_size / (distance * double(interp_size));
         }

         // Refine the position of i due to the force aggregrate
         matrix_row<matrix<double, column_major> >(X, i) += net_force;
         
      } 
   }

   // Sanity check - output the stress
   cout << "Pretune - Classical Stress: " << calculate_stress(X, R, normalised) << endl; 

   // Finally run Chalmer's iterations on the output to reduce the final stress
   cha96mds(R, X, p, rng, v_max, s_max, false, sqrt(n));
}
