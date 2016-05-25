#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/random.hpp>
#include <boost/limits.hpp>
#include <iostream/manipulators/csv.h>
#include <ietl/interface/ublas.h>
#include <ietl/vectorspace.h>
#include <ietl/lanczos.h>
#include <ietl/iteration.h>

using namespace boost::numeric::ublas;
using namespace std;

int main(int argc, char** argv)
{
   // Load feature data from file 
   matrix<double, column_major> X(3, 3);
/*
   X(0,0) = 3;
   X(0,1) = 2;
   X(1,0) = 2;
   X(1,1) = 3;
*/
/*
   X(0,0) = 2;
   X(0,1) = 1;
   X(1,0) = 1;
   X(1,1) = 2;
*/
   X(0,0) = 2;
   X(1,1) = 3;
   X(2,2) = 9;
   X(1,2) = 4;
   X(2,1) = 4;

   // Use Iterative Eigensolver Template Library to solve for p largest eigenvectors
   typedef ietl::wrapper_vectorspace<boost::numeric::ublas::vector<double> > Vecspace;
   typedef boost::mt19937 Gen; 
   Gen mygen(123456789);
   Vecspace vec(X.size1());
   ietl::lanczos<matrix<double, column_major>, Vecspace> lanczos(X, vec);

   // Creation of an iteration object:    
   unsigned int max_iter = 20;

   std::vector<double> eigenvalues;
   std::vector<boost::numeric::ublas::vector<double> > eigenvectors;
   ietl::Info<double> info; 
   //ietl::fixed_lanczos_iteration<double> iter(max_iter);
   ietl::lanczos_iteration_nhighest<double> iter(max_iter, 2, 10e-4, 10e-4);

   try
   {
      // Run the lanczos solver
      lanczos.calculate_eigenvalues(iter,mygen);
      eigenvalues = lanczos.eigenvalues();

      for (std::vector<double>::const_reverse_iterator itr = eigenvalues.rbegin(), end = eigenvalues.rend(); itr != end; itr++)
      {
         double eigenvalue = *itr;
         cout << "Eigenvalue: " << eigenvalue << endl;
      }

      matrix<double, column_major> U = zero_matrix<double>(X.size1(), eigenvalues.size());
      lanczos.eigenvectors(eigenvalues.rbegin(), eigenvalues.rend(), back_inserter(eigenvectors), info, mygen); 
      for (unsigned int j = 0; j < eigenvalues.size(); ++j)
      {
         matrix_column<matrix<double, column_major> > column(U, j);
         column += eigenvectors[j];
      }
      cout << "Eigenvectors: " << csv << U << endl;

   }
   catch (exception &e)
   {
      cerr << e.what() << endl;
   }
}
