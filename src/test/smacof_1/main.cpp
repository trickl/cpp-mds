#include <mds/rkisomap.h>
#include <mds/smacof_1.h>
#include <mds/lmds.h>
#include <mds/stress.h>
#include <mds/euclidean_distance.h>
#include <iostream/manipulators/csv.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <fstream>

using namespace boost::numeric::ublas;
using namespace std;

int main(int argc, char** argv)
{
   // Load feature data from file
   matrix<double, column_major> X;
   ifstream fin("/home/tgee/projects/c++/algo/data/swiss_roll.dat");
   fin >> csv >> X;

   try
   {
      // Convert into similarity data
      symmetric_matrix<double, lower, column_major> R;
      euclidean_distance(X, R);

      matrix<double, column_major> Xscaled;
      symmetric_matrix<double, lower, column_major> S;
      rkisomap(R, S, 5);
      smacof_1(S, 2, Xscaled);

      ofstream fout("/home/tgee/projects/c++/algo/data/swiss_roll_rkisomap.dat");
      fout << csv << Xscaled;
   }
   catch (exception &e)
   {
      cerr << e.what() << endl;
   }
}
