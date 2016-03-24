#include <iostream/manipulators/csv.h>
#include <fstream>
#include <mds/pdbmds.h>
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

class pdbmds_interface : public key_handler, public drawable
{
public:
   pdbmds_interface(pdbmds& solver,
                  symmetric_matrix<double, lower, column_major>& R,
                  drawing_pad &pad,
                  gsl_rng* rng)
       : _solver(solver), _R(R), _pad(pad), _rng(rng), _total_relations(0), _total_adjustments(0), _user_id(0) {};

   virtual void draw() const
   {
      const matrix<double, column_major> &X = _solver.item_locations_as_matrix();
      for (unsigned int i = 0; i < X.size1(); ++i)
      {
         // Color gradient
         glColor3f(double(i) / double(X.size1()), 0.0, double(X.size1() - i) / double(X.size1()));
         glPointSize(6);
         glBegin(GL_POINTS);

         // Debug
         if ((abs(X(i, 0)) > 100 && abs(X(i, 0)) < 10e20) ||
             (abs(X(i, 1)) > 100 && abs(X(i, 1)) < 10e20)) 
         {
            cout << "Item " << i << " has been displaced to " << X(i, 0) << ", " << X(i, 1) << endl; 
         } 

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
      ofstream fout("/home/tgee/projects/c++/algo/data/inclined_plane_pdbmds.dat");

      switch(key)
      {
      case 'c':
         process_random_ratings(50);
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
      default:
         break;
      }

      cout << "Adjustments: " << _total_adjustments << ", Relations:" << _total_relations << endl;
   } 

   void process_random_ratings(unsigned int size)
   {
      std::vector<pdbmds::rating> ratings; 
      for (unsigned int iteration = 0; iteration < size; ++iteration)
      {
         unsigned int i = gsl_rng_uniform_int(_rng, _R.size1());
         unsigned int j = gsl_rng_uniform_int(_rng, _R.size1());
         ratings.push_back(pdbmds::rating(_user_id, i, 0));
         ratings.push_back(pdbmds::rating(_user_id++, j, _R(i, j)));
      }
      _solver.queue_ratings(ratings);
      _solver.process_rating_queue(ratings.size(), true);
      _total_relations += size; 
   }

   void adjust_items(unsigned int size)
   {
      std::vector<unsigned int> ids;
      for (unsigned int iteration = 0; iteration < size; ++iteration)
      {
         unsigned int i = gsl_rng_uniform_int(_rng, _R.size1());
         ids.push_back(i);
      }
      _solver.adjust_items(ids);
      _total_adjustments += size; 
   } 

private:
   pdbmds& _solver;
   symmetric_matrix<double, lower, column_major>& _R;
   drawing_pad& _pad;
   gsl_rng*     _rng;
   unsigned int _total_relations;
   unsigned int _total_adjustments;
   unsigned int _user_id;
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

      pdbmds solver(20, 20, 2, 0.2);
      drawing_pad pad(720, 720, 20, 20, "Progressive MDS Solver");
      pad.set_scale(-10, -10, 10, 10);
      pdbmds_interface pi(solver, R, pad, rng);
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
