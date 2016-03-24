#include <iostream/manipulators/csv.h>
#include <fstream>
#include <mds/imds.h>
#include <mds/stress.h>
#include <mds/euclidean_distance.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/timer.hpp>
#include <ctime>
#include <iomanip>
#include <functional>

using namespace std;
using namespace boost::numeric::ublas;

int main(int argc, char** argv)
{
   // Load feature data from file 
   matrix<double, column_major> X;
   ifstream fin("/home/tgee/projects/c++/algo/data/inclined_plane.dat"); 
   fin >> csv >> X;

   try
   {
      // Convert into similarity data
      cout << "Creating distance matrix..." << endl;
      symmetric_matrix<double, lower, column_major> R(X.size1());
      euclidean_distance(X, R);

      // Weights
      cout << "Creating weights matrix..." << endl;
      symmetric_matrix<double, lower, column_major> W(R.size1());
      for (symmetric_matrix<double, lower, column_major>::iterator1 i = W.begin1(), i_end = W.end1(); i != i_end; ++i) fill(i.begin(), i.end(), 1.);

      cout << "Now running imds..." << endl;
      unsigned int overlap = 10;
      unsigned int p = 2;
      imds solver(p, overlap, 1e-3, 100);

      matrix<double, column_major> features(R.size1(), p); 
      for (unsigned int i = 0, max_itr = 10; i < max_itr; ++i)
      {
         matrix<double, column_major> align, lfeatures;
         unsigned int local_width  = (int) (ceil(double(R.size1()) / double(max_itr)));
         unsigned int local_offset = (int) floor(double(i) * double(R.size1()) / double(max_itr));
         symmetric_matrix<double, lower, column_major> dissimilarities = matrix_range<symmetric_matrix<double, lower, column_major> >(R, range(i == 0 ? local_offset : local_offset - overlap, local_offset + local_width),
                                                            range(i == 0 ? local_offset : local_offset - overlap, local_offset + local_width));
         symmetric_matrix<double, lower, column_major> weights = matrix_range<symmetric_matrix<double, lower, column_major> >(W, range(i == 0 ? local_offset : local_offset - overlap, local_offset + local_width),
                                                    range(i == 0 ? local_offset : local_offset - overlap, local_offset + local_width));

         cout << "Iteration: " << i << endl;
         cout << "Calculating from " << (i == 0 ? local_offset : local_offset - overlap) << " to " << local_offset + local_width << endl;

         solver(dissimilarities, weights, lfeatures, align);

         // Transform the old points
         if (local_offset > 0)
         {
            matrix<double, column_major> old_features(local_offset, p + 1);
            matrix_range<matrix<double, column_major> >(old_features, range(0, local_offset), range(0, p)) =  matrix_range<matrix<double, column_major> >(features, range(0, local_offset), range(0, p));
            matrix_column<matrix<double, column_major> >(old_features, p) = scalar_vector<double>(local_offset, 1.);

            matrix_range<matrix<double, column_major> >(features, range(0, local_offset), range(0, p)) = prod(old_features, align);
         }

         matrix_range<matrix<double, column_major> >(features, range(i == 0 ? local_offset : local_offset - overlap, local_offset + local_width), range(0, p)) = lfeatures;
      }
      
      cout << "Stress: " << calculate_stress(features, R, W, normalised) << endl;

      ofstream fout("/home/tgee/projects/c++/algo/data/inclined_plane_imds.dat");
      fout << csv << features;

   }
   catch (exception &e)
   {
      cerr << e.what() << endl;
   }
}
