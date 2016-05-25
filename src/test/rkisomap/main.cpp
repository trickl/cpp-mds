#include <mds/rkisomap.h>
#include <mds/cmds.h>
#include <mds/stress.h>
#include <mds/euclidean_distance.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <iostream/manipulators/csv.h>
#include <fstream>

using namespace boost::numeric::ublas;
using namespace std;

int main(int argc, char** argv)
{
   // Load feature data from file 
   matrix<double, column_major> X;
   ifstream fin("swiss_roll_with_noise.dat");
   fin >> csv >> X;
   
   try
   {
      // Convert into similarity data
      symmetric_matrix<double, lower, column_major> R(X.size1());
      euclidean_distance(X, R);

      matrix<double, column_major> Xscaled;
      symmetric_matrix<double, lower, column_major> S(X.size1());
      rkisomap(R, S, 5, 2, 0.04);

      ofstream fout_tmp("swiss_roll_rkisomap_S.dat");
      fout_tmp << csv << S;

      cmds(S, 2, Xscaled);

      ofstream fout("swiss_roll_rkisomap.dat");
      fout << csv << Xscaled;
   }
   catch (exception &e)
   {
      cerr << e.what() << endl;
   }
}
