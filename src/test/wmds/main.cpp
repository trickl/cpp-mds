#include <iostream/manipulators/csv.h>
#include <fstream>
#include <mds/wmds.h>
#include <mds/stress.h>
#include <mds/euclidean_distance.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <ctime>
#include <iomanip>
#include <gsl/gsl_randist.h>

using namespace std;
using namespace boost::numeric::ublas;

int main(int argc, char** argv)
{
   // Load feature data from file 
   matrix<double, column_major> X;
   ifstream fin("/home/tgee/projects/c++/algo/data/inclined_plane.dat");
   fin >> csv >> X;

  // Setup the gsl random number generator
   gsl_rng_env_setup();
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
   long seed = time (NULL) *getpid();
   gsl_rng_set(rng, seed);

   try
   {
      // Convert into similarity data
      symmetric_matrix<double, lower, column_major> R(X.size1());
      euclidean_distance(X, R);
/*
      // Add some noise
      for (symmetric_matrix<double, lower, column_major>::iterator1 i_itr = R.begin1(), i_end = R.end1(); i_itr != i_end; ++i_itr)
      {
         // Weight low index points higher (less noise)
         for (symmetric_matrix<double, lower, column_major>::iterator2 j_itr = i_itr.begin(), j_end = i_itr.end(); j_itr != j_end; ++j_itr)
         {
            // Add high noise for high index points
            long i = i_itr.index1(), j = j_itr.index2();
            double noise = (log(i + 1) / log(R.size1())) * 0.4 * (0.5 - gsl_rng_uniform(rng));
            R(i, j) = R(i, j) * (1 + noise);
         }
      }
*/

      // Weights for stress measure
      symmetric_matrix<double, lower, column_major> W(R.size1());
      for (symmetric_matrix<double, lower, column_major>::iterator1 i_itr = W.begin1(), i_end = W.end1(); i_itr != i_end; ++i_itr)
      {
         // Weight low index points higher (less noise)
         for (symmetric_matrix<double, lower, column_major>::iterator2 j_itr = i_itr.begin(), j_end = i_itr.end(); j_itr != j_end; ++j_itr)
         {
            long i = i_itr.index1(), j = j_itr.index2();
 //           W(i, j) = (log(R.size1()) - log(i + 1)) / log(R.size1());
              W(i, j) = 1;
         }
      }

      matrix<double, column_major> Xscaled;
      clock_t start, end;
      start = clock();
      wmds(R, W, Xscaled, 2, 0.001, 0.00001, 100);
      end = clock();

      cout << setprecision(5) << "Time: " << end - start << " ticks" << endl;
      cout << "Classical Stress: " << calculate_stress(Xscaled, R, normalised) << endl;
      cout << "Weighted Stress: " << calculate_stress(Xscaled, R, W, normalised) << endl;

      ofstream fout("/home/tgee/projects/c++/algo/data/inclined_plane_wmds.dat");
      fout << csv << Xscaled;

   }
   catch (exception &e)
   {
      cerr << e.what() << endl;
   }
}
