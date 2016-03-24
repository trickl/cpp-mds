#include <iostream/manipulators/csv.h>
#include <fstream>
#include <mds/cmds.h>
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

      matrix<double, column_major> Xscaled;
      clock_t start, end;
      start = clock();
      cmds(R, 2, Xscaled);
      end = clock();

      cout << setprecision(5) << "Time: " << end - start << " ticks" << endl;
      cout << "Classical Stress: " << calculate_stress(Xscaled, R, normalised) << endl;
      cout << "Weighted Stress: " << calculate_stress(Xscaled, R, normalised) << endl;

      ofstream fout("/home/tgee/projects/c++/algo/data/inclined_plane_cmds.dat");
      fout << csv << Xscaled;

   }
   catch (exception &e)
   {
      cerr << e.what() << endl;
   }
}
