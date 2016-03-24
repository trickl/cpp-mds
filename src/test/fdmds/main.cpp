#include <iostream/manipulators/csv.h>
#include <fstream>
#include <mds/fdmds.h>
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

class fdmds_interface : public key_handler, public drawable
{
public:
   fdmds_interface(fdmds& solver,
                  symmetric_matrix<double, lower, column_major>& R,
                  drawing_pad &pad)
       : _solver(solver), _R(R), _pad(pad), _iteration(0) {};

   virtual void draw()
   {
      const matrix<double, column_major> &X = _solver.points();
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
      ofstream fout("/home/tgee/projects/c++/algo/data/inclined_plane_fdmds.dat");

      switch(key)
      {
      case 'c': 
         break;
      case 'x':
         process_all_points();
         _pad.draw();
         break;
      case 's': 
         // Calculate stress 
         cout << "Classical Stress: " << calculate_stress(_solver.points(), _R, normalised) << endl;
         break;
      case 'f':
         // Write to file
         fout << csv << _solver.points();
         break;
      default:
         break;
      }
   } 

   void process_all_points()
   {
      _solver.process_relations(_R, _iteration++);
   } 

private:
   fdmds& _solver;
   symmetric_matrix<double, lower, column_major>& _R;
   drawing_pad& _pad;
   unsigned int _iteration; 
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
   long seed = time (NULL) *getpid();
   gsl_rng_set(rng, seed);

   try
   {
      // Convert into similarity data
      symmetric_matrix<double, lower, column_major> R(X.size1());
      euclidean_distance(X, R);

      fdmds solver(rng, R.size1(), 2, false, 1., 10., 1.);
      drawing_pad pad(720, 720, 20, 20, "Force-directed MDS Solver");
      pad.set_scale(0, 0, 2, 2);
      fdmds_interface fi(solver, R, pad);
      pad.handle(&fi);
      pad.draw(&fi); 
      
      Fl::run();
   }
   catch (exception &e)
   {
      cerr << e.what() << endl;
   }

   gsl_rng_free(rng);
}
