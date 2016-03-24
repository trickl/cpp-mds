#include <iostream/manipulators/csv.h>
#include <fstream>
#include <mds/smds.h>
#include <mds/stress.h>
#include <mds/euclidean_distance.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <ctime>
#include <iomanip>

using namespace boost::numeric::ublas;
using namespace std;

int main(int argc, char** argv)
{
   // Load feature data from file 
   ifstream fin("/home/tgee/projects/c++/algo/data/inclined_plane.dat");
   matrix<double, column_major> X;
   fin >> csv >>  X;

   try
   {
      // Convert into similarity data
      cout << "Creating distance matrix..." << endl;
      symmetric_matrix<double, lower, column_major> R(X.size1());
      euclidean_distance(X, R);

      matrix<double, column_major> Xscaled;
      clock_t start, end;
      start = clock();
      smds(R, Xscaled, 2, 1e-03);
      end = clock();

      cout << setprecision(5) << "Time: " << end - start << " ticks" << endl;
      cout << "Stress: " << calculate_stress(Xscaled, R, normalised) << endl;

      ofstream fout("/home/tgee/projects/c++/algo/data/inclined_plane_smds.dat"); 
      fout << csv << Xscaled;

   }
   catch (exception &e)
   {
      cerr << e.what() << endl;
   }
}
