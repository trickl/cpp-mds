#include <fstream>
#include <mds/cpmds.h>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <draw/drawing_pad.h>
#include <GL/gl.h>
#include <Fl/gl.h>
#include <Fl/Fl.H>

using namespace std;

class cpmds_interface : public key_listener, public drawable
{
public:
   cpmds_interface(cpmds& solver,
                  drawing_pad &pad)
       : _solver(solver), _pad(pad) {};

   void draw() 
   {
      const item_table_type& items = _solver.items();
      for (unsigned int i = 0; i < items.size(); ++i)
      {
         item itm = items[i];
         // Draw a circle
         glColor4f(double(i) / double(items.size()), 0.0, double(items.size() - i) / double(items.size()), 0.5);
         glBegin(GL_POLYGON);

         unsigned int total_angle_steps = 36;
         for (unsigned int angle_step = 0; angle_step < total_angle_steps; ++angle_step)
         {
            glVertex2d(itm.location.x + itm.radius * cos(2 * M_PI * angle_step / total_angle_steps),
                       itm.location.y + itm.radius * sin(2 * M_PI * angle_step / total_angle_steps));
         }

         glEnd();

         // Draw a line to the first neighbour
         item neighbour;
         if (_solver.has_neighbour(itm.id, neighbour))
         {
            glColor4f(0, 0, 0, 0.5);
            glBegin(GL_LINES);
            glVertex2d(neighbour.location.x, neighbour.location.y);
            glVertex2d(itm.location.x, itm.location.y);
            glEnd();

            // Draw a line from this neighbour using the angle sum
            glColor3f(0.8, 0, 0);
            glBegin(GL_LINES);
            point angle_sum_line = rotate(neighbour.location, itm.angle_sum, itm.location);
            glVertex2d(itm.location.x, itm.location.y);
            glVertex2d(angle_sum_line.x, angle_sum_line.y);
            glEnd();
         }

         // Draw the text
         stringstream point_text;
         point_text << itm.id << " : " << itm.angle_sum * 180 /  M_PI;

         glColor3f(0, 0, 0);
         gl_font(1, 8);
         gl_draw(point_text.str().c_str(), float(itm.location.x), float(itm.location.y));
      }
   }

   virtual void on_key_down(unsigned char key, int x, int y)
   {
      ofstream fout("/home/tgee/projects/c++/algo/data/inclined_plane_cpmds.dat");

      switch(key)
      {
      case 'r':
         run_test();
         _pad.request_redraw();
         break;
      default:
         break;
      }
   } 

   void run_test()
   {
      // Low curvature test 
      _solver.process_relation(1, 2, 3, 3);
      _solver.process_relation(2, 3, 3, 3);
      _solver.process_relation(4, 3, 4, 4);
      _solver.process_relation(5, 1, 4, 4);
      _solver.process_relation(5, 6, 2, 2);
      _solver.process_relation(4, 6, 2, 2);
      _solver.process_relation(3, 7, 2, 2);
      _solver.process_relation(1, 8, 2, 2);
      _solver.process_relation(1, 9, 2, 2);

      // High curvature test
/*
      _solver.process_relation(1, 2, 3, 3);
      _solver.process_relation(1, 3, 3, 3);
      _solver.process_relation(1, 4, 3, 3);
      _solver.process_relation(1, 5, 3, 3);
      _solver.process_relation(1, 6, 3, 3);
      _solver.process_relation(1, 7, 3, 3);
*/
      _solver.test_update_curvature();

/*
      _solver.test_remove_item(1);
      _solver.test_whitehead_move(4, 5);
*/
   }

private:
   cpmds& _solver;
   drawing_pad& _pad;
};

int main(int argc, char** argv)
{
   try
   {
      // Convert into similarity data
      cpmds solver;
      drawing_pad pad(720, 720, 20, 20, "Progressive MDS Solver");
      pad.set_viewport(-10, -10, 10, 10);
      cpmds_interface pi(solver, pad);
      pad.add_key_listener(&pi);
      pad.draw(&pi); 
      
      Fl::lock();
      Fl::run();
   }
   catch (exception &e)
   {
      cerr << e.what() << endl;
   }
}
