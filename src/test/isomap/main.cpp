#include <mds/isomap.h>
#include <mds/cmds.h>
#include <mds/stress.h>
#include <mds/euclidean_distance.h>
#include <iostream/manipulators/csv.h>
#include <fstream>

using namespace boost::numeric::ublas;
using namespace std;

int main(int argc, char** argv)
{
   // Load feature data from file 
   matrix<double, column_major> X;
   ifstream fin("/home/Sleeper/projects/cluster/data/swiss_roll.dat");
   fin >> csv >> X;
   
   try
   {
      // Convert into similarity data
      symmetric_matrix<double, lower, column_major> R;
      euclidean_distance(X, R);

      matrix<double, column_major> Xscaled;
      symmetric_matrix<double, lower, column_major> S;
      isomap(R, S, 5);

      ofstream fout_tmp("/home/Sleeper/projects/cluster/data/swiss_roll_isomap_S.dat");
      fout_tmp << csv << S;

      cmds(S, 2, Xscaled);

      ofstream fout("/home/Sleeper/projects/cluster/data/swiss_roll_isomap.dat");
      fout << csv << Xscaled;
   }
   catch (exception &e)
   {
      cerr << e.what() << endl;
   }
}
