#include <iostream/manipulators/csv.h>
#include <fstream>
#include <mds/pmds.h>
#include <mds/stress.h>
#include <mds/euclidean_distance.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <gsl/gsl_randist.h>
#include <draw/drawing_pad.h>
#include <GL/gl.h>
#include <Fl/gl.h>
#include <Fl/Fl.H>

using namespace std;
using namespace boost::numeric::ublas;

class pmds_interface : public key_handler, public drawable
{
public:
   pmds_interface(pmds& solver,
                  symmetric_matrix<double, lower, column_major>& R,
                  drawing_pad &pad,
                  gsl_rng* rng)
       : _solver(solver), _R(R), _pad(pad), _rng(rng), _total_relations(0), _total_adjustments(0) {};

   virtual void draw() const
   {
      const matrix<double, column_major> &X = _solver.item_locations_as_matrix();
      for (unsigned int i = 0; i < X.size1(); ++i)
      {
         // Color gradient
         glColor3f(double(i) / double(X.size1()), 0.0, double(X.size1() - i) / double(X.size1()));
         glPointSize(6);
         glBegin(GL_POINTS);
         glVertex2d(X(i, 0), X(i, 1));
         glEnd();

         stringstream point_text;
         point_text << i;

         glColor3f(0, 0, 0);
         gl_font(1, 8);
         gl_draw(point_text.str().c_str(), float(X(i, 0)), float(X(i, 1)));
      }
   }

   virtual void handle_key(unsigned char key, int x, int y)
   {
      ofstream fout("/home/tgee/projects/c++/algo/data/inclined_plane_pmds.dat");

      switch(key)
      {
      case 'r':
         process_all_relations();
         _pad.draw();
         break;
      case 'e':
         adjust_all_items();
         _pad.draw();
         break;
      case 'c':
         process_random_relations(100);
         _pad.draw();
         break;
      case 'x':
         adjust_items(100);
         _pad.draw();
         break;
      case 's': 
         // Calculate stress 
         cout << "Points : " << _solver.item_locations_as_matrix().size1() << endl;
         cout << "Classical Stress: " << calculate_stress(_solver.item_locations_as_matrix(), _R, normalised) << endl;
         break;
      case 'f':
         // Write to file
         fout << csv << _solver.item_locations_as_matrix();
         break;
      default:
         break;
      }

      cout << "Adjustments: " << _total_adjustments << ", Relations:" << _total_relations << endl;
   } 

   void process_all_relations()
   {
      for (unsigned int i = 0; i < _R.size1(); ++i)
      {
        for (unsigned int j = 0; j < _R.size1(); ++j)
        {
            _solver.process_relation(i, j, 0., _R(i, j));
        } 
      } 
      _total_relations += _R.size1() * _R.size1(); 
   }

   void adjust_all_items()
   {
      for (unsigned int i = 0; i < _R.size1(); ++i)
      {
        _solver. adjust_item(i);
      }
      _total_adjustments += _R.size1(); 
   }

   void process_random_relations(unsigned int size)
   {
      for (unsigned int iteration = 0; iteration < size; ++iteration)
      {
         unsigned int i = gsl_rng_uniform_int(_rng, _R.size1());
         unsigned int j = gsl_rng_uniform_int(_rng, _R.size1());
         _solver.process_relation(i, j, 0., _R(i, j));
      }
      _total_relations += size; 
   }

   void adjust_items(unsigned int size)
   {
      for (unsigned int iteration = 0; iteration < size; ++iteration)
      {
         unsigned int i = gsl_rng_uniform_int(_rng, _R.size1());
         _solver.adjust_item(i);
      }
      _total_adjustments += size; 
   } 

private:
   pmds& _solver;
   symmetric_matrix<double, lower, column_major>& _R;
   drawing_pad& _pad;
   gsl_rng*     _rng;
   unsigned int _total_relations;
   unsigned int _total_adjustments;
};

int main(int argc, char** argv)
{
   // Load feature data from file 
   matrix<double, column_major> X;
   ifstream fin("/home/tgee/projects/c++/algo/data/inclined_plane.dat");
   fin >> csv >> X;

  // Setup the gsl random number generator
   gsl_rng_env_setup();
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
   long seed = 213;
   //long seed = time (NULL) *getpid();
   gsl_rng_set(rng, seed);

   try
   {
      // Convert into similarity data
      symmetric_matrix<double, lower, column_major> R(X.size1());
      euclidean_distance(X, R);

      pmds solver(20, 20, 0.2);
      drawing_pad pad(720, 720, 20, 20, "Progressive MDS Solver");
      pad.set_scale(-10, -10, 10, 10);
      pmds_interface pi(solver, R, pad, rng);
      pad.handle(&pi);
      pad.draw(&pi); 
      
      Fl::run();
   }
   catch (exception &e)
   {
      cerr << e.what() << endl;
   }

   gsl_rng_free(rng);
}
