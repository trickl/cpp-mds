#include <mds/lmds.h>

using namespace std;

class cout_fetch_callback : public lmds::range_callback
{
public:
   cout_fetch_callback() {};

   virtual void operator()(unsigned int start, unsigned int end)
   {
      cout << "On fetch data: " << start << " -> " <<  end << endl;
   };
};

class cout_free_callback : public lmds::range_callback
{
public:
   cout_free_callback() {};

   virtual void operator()(unsigned int start, unsigned int end)
   {
      cout << "On free data: " << start << " -> "  << end << endl;
   };
};

int main(int argc, char** argv)
{
   cout_fetch_callback on_fetch;
   cout_free_callback on_free;

   lmds::range_callback_buffer buffer(on_fetch, on_free, 1);
   buffer.on_fetch()(0, 500);
   buffer.on_free()(0, 500);

   buffer.on_fetch()(500, 1000);
   buffer.on_free()(500, 1000);
   buffer.on_fetch()(500, 1000);
   buffer.on_free()(500, 1000);

   buffer.on_fetch()(1000, 1500);
   buffer.on_free()(1000, 1500);
   buffer.on_fetch()(1000, 1500);
   buffer.on_free()(1000, 1500);
   buffer.on_fetch()(1000, 1500);
   buffer.on_free()(1000, 1500);

   buffer.on_fetch()(1500, 2000);
   buffer.on_free()(1500, 2000);
   buffer.on_fetch()(1500, 2000);
   buffer.on_free()(1500, 2000);

   buffer.on_fetch()(2000, 2500);
   buffer.on_free()(2000, 2500);
}
