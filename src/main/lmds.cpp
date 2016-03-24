#include "lmds.h"
#include "wmds.h"
#include <stdexcept>
#include <cmath>
#include <vector>
#include <set>
#include <algorithm>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/bindings/lapack/gesvd.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

using namespace std;
using namespace boost::numeric::ublas;
using namespace boost::numeric::bindings;

// Heuristic to estimate time cost of performing normal mds on a matrix of certain size
double mds_time_cost(unsigned long size)
{
   // TODO: Get a decent heuristic by performing some tests
   // For the moment, lets naively assume a quadratic cost factor and a fixed constant
   return 100 + (size / 10.) * (size / 10.);
};

struct lmds_time_cost_params
{
   unsigned int p;   // Estimated dimensionality 
   unsigned long R_n; // Size of dissimiliarity matrix
};

// Heuristic to estimate time cost of performing local mds given subdivisions
double lmds_time_cost(double subdivisions, void *params)
{
    // Main factors are mds of alignment matrix and submatrices and a fixed constant
    lmds_time_cost_params* ltc_params = reinterpret_cast<lmds_time_cost_params*>(params);
    unsigned int s = (unsigned int) subdivisions;
    return 100 + mds_time_cost(2 * ltc_params->p * s) + s * mds_time_cost(ltc_params->R_n / s);
};

// First derivative of lmds time cost (to be minimized)
double df_lmds_time_cost(double subdivisions, void *params)
{
   return lmds_time_cost(subdivisions + 1, params) - lmds_time_cost(subdivisions, params);
}

void lmds(const symmetric_matrix<double, lower, column_major> &R,
          const symmetric_matrix<double, lower, column_major> &W,
          matrix<double, column_major> &X,
          unsigned int p,
          double variance_tolerance)
{
   if (R.size1() != W.size1()) throw invalid_argument("Weights W must have same dimensions are R.");
   long n = R.size1(); // Number of data soints

   // Calculate the subdivision size based on the mds cost heuristic
   const gsl_root_fsolver_type *solver_type = gsl_root_fsolver_brent;
   gsl_root_fsolver *solver = gsl_root_fsolver_alloc(solver_type);
   gsl_function solver_function;
   unsigned int subdivisions = 0;

   // Bounds are between 1 and 100,000
   lmds_time_cost_params solver_params = {p, R.size1()};
   solver_function.function  = &df_lmds_time_cost;
   solver_function.params = &solver_params;
   gsl_root_fsolver_set(solver, &solver_function, 1., 100000.); 
   double root0 = -1, root = -1;

   for (int iter = 0, max_iter = 100, solver_status = GSL_CONTINUE;
        solver_status == GSL_CONTINUE && iter < max_iter; ++iter)
   {
           solver_status = gsl_root_fsolver_iterate(solver);
           root0 = root;
           root = gsl_root_fsolver_root(solver);
           solver_status = gsl_root_test_delta(root, root0, 1., 0);

           if (solver_status == GSL_SUCCESS)
           {
               if (root > 2. && root < 10000)
               {
                  subdivisions = (unsigned int) floor(root);
                  cout << "Division cost solver settled on " << subdivisions << " subdivisions." << endl; 
               }
               else
               {
		  throw runtime_error("lmds: Division cost solver could not find sensible subdivision root().");
               } 
           } 
           else if (iter == max_iter)
           {
	       throw runtime_error("lmds: Division cost solver did not converge.");
           }
   }
   gsl_root_fsolver_free(solver);

   // Divide the dissimilarity matrix into submatrices
   set<long> align_points;
   long D_width  = (long) (ceil(double(n) / double(subdivisions)));

   for (unsigned int k = 0; k < subdivisions; ++k)
   {
      unsigned long D_offset = (long) floor( double(k) * double(n) / double(subdivisions));

      matrix_range<symmetric_matrix<double, lower, column_major> > D(const_cast<symmetric_matrix<double, lower, column_major> &>(R), range(D_offset, D_offset + D_width), range(D_offset, D_offset+ D_width));

      // Sample  p * 2 points from each submatrix
      unsigned int sample_size = p * 2;
      std::vector<long> D_points;
      for (long i = 0; i < D_width; ++i) D_points.push_back(D.start1() + i);
      random_shuffle(D_points.begin(), D_points.end());
      align_points.insert(D_points.begin(), D_points.begin() + sample_size);
   }

   // Construct alignment matrix from the sampled points
   symmetric_matrix<double, lower, column_major> M(align_points.size());
   symmetric_matrix<double, lower, column_major> WM(align_points.size());
   {
      unsigned long i = 0;
      for (set<long>::const_iterator i_itr = align_points.begin(), i_end = align_points.end(); i_itr != i_end; ++i_itr)
      {
         unsigned long j = i;
         for (set<long>::const_iterator j_itr = i_itr, j_end = i_end; j_itr != j_end; ++j_itr)
         {
            if (i_itr == j_itr)
            {
               M(i, j) = 0;
               WM(i, j) = 0;
            }
            else
            {
               M(i, j)  = R(*i_itr, *j_itr);
               WM(i, j) = W(*i_itr, *j_itr);
            }
            j++;
         }
         i++;
      }
   }

   // Run MDS on the alignment matrix to produce an alignment mapping A
   matrix<double, column_major> M_mds;
   wmds(M, WM, M_mds, p, variance_tolerance); 

   // Allocate space for the final feature matrix
   X.resize(R.size1(), p);

   // Run MDS on each sub matrix
   for (unsigned int k = 0; k < subdivisions; ++k)
   {
      cout << "Processing submatrix " << k + 1 << " of " << subdivisions << endl;
      unsigned long D_offset = (long) floor( double(k) * double(n) / double(subdivisions));
      matrix_range<symmetric_matrix<double, lower, column_major> > D(const_cast<symmetric_matrix<double, lower, column_major> &>(R), range(D_offset, D_offset + D_width), range(D_offset, D_offset+ D_width));
      matrix_range<symmetric_matrix<double, lower, column_major> > WD(const_cast<symmetric_matrix<double, lower, column_major> &>(W), range(D_offset, D_offset + D_width), range(D_offset, D_offset+ D_width));

      matrix<double, column_major> D_mds;
      wmds(D, WD, D_mds, p, variance_tolerance);
      set<long> D_align_points(align_points.upper_bound(D_offset - 1), align_points.lower_bound(D_offset + D_width));

      unsigned long M_offset = 0;
      for (set<long>::const_iterator itr = align_points.begin(), end = align_points.upper_bound(D_offset - 1); itr != end; itr++) M_offset++;

      matrix<double, column_major> M_mds_i = matrix_range<matrix<double, column_major> >(M_mds, range(M_offset, M_offset + D_align_points.size()), range(0, M_mds.size2()));
      matrix<double, column_major> D_mds_i(D_align_points.size(), p + 1);
      {
         unsigned long i = 0;
         for (set<long>::const_iterator itr = D_align_points.begin(), end = D_align_points.end(); itr != end; ++itr)
         {
            boost::numeric::ublas::vector<double> D_mds_row = matrix_row<matrix<double, column_major> >(D_mds, *itr - D.start1());
            D_mds_row.resize(p + 1);
            D_mds_row(p) = 1; // Allow for translation between D and M
            matrix_row<matrix<double, column_major> >(D_mds_i, i) = D_mds_row;
            i++;
         }
      } 

      // Use SVD to find pseudoinverse of D for calculation of affine mapping A
      matrix<double, column_major> U(D_mds_i.size1(), D_mds_i.size2());
      matrix<double, column_major> Vt(D_mds_i.size2(), D_mds_i.size2());
      boost::numeric::ublas::vector<double> S(D_mds_i.size2());
      matrix<double, column_major> D_mds_i2 = D_mds_i;
      lapack::gesvd(D_mds_i, S, U, Vt); 

      // A = (V * S^-1 * U^T) * M_mds where (U * S * V^T) is the SVD solution of D_mds
      // First calc S^-1 * U^T
      matrix<double, column_major> SUt(D_mds_i.size2(), D_mds_i.size1());
      for (unsigned long i = 0; i < D_mds_i.size2(); i++)
      {
         matrix_column<matrix<double, column_major> > U_col(U, i);
         matrix_row<matrix<double, column_major> > SUt_row(SUt, i);
         SUt_row = U_col / S(i);
      }

      // Get V from transpose Vt
      matrix<double, column_major> V(D_mds_i.size2(), D_mds_i.size2());
      for (matrix<double, column_major>::const_iterator1 i_itr = Vt.begin1(), i_end = Vt.end1(); i_itr != i_end; ++i_itr)
      {
         for (matrix<double, column_major>::const_iterator2 j_itr = i_itr.begin(), j_end = i_itr.end(); j_itr != j_end; ++j_itr)
         {
            unsigned long i = i_itr.index1(), j = j_itr.index2();
            V(j, i) =  Vt(i, j);
         }
      }

      matrix<double, column_major> Dinv(D_mds_i.size2(), D_mds_i.size1());
      opb_prod(V, SUt, Dinv);

      matrix<double, column_major> A(D_mds_i.size2(), M_mds_i.size2());
      opb_prod(Dinv, M_mds_i, A);

      // Last row is translation vector,  store it separately
      boost::numeric::ublas::vector<double> translation = matrix_row<matrix<double, column_major> >(A, p); 
      A = matrix_range<matrix<double, column_major> >(A, range(0, p), range(0, p));

      // Use alignment mapping to transform all points D(i) into aligned M space.
      for (unsigned long i = 0; i < D_mds.size1(); i++)
      {
         matrix_row<matrix<double, column_major> > D_mds_row(D_mds, i);
         matrix_row<matrix<double, column_major> > X_row(X, D.start1() + i);

         // TODO: Handle edge overlaps in submatrices by averaging
         axpy_prod(D_mds_row, A, X_row);
	 X_row += translation;
      } 
   }
}


