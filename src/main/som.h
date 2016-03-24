#ifndef _SOM_H
#define _SOM_H

#include <list>
#include <map>

struct adjacency_face
{
   unsigned int first_itm_id;
   unsigned int second_itm_id;
   unsigned int third_itm_id;
};

// Self organising map
class som
{
public:
   som();

   void initialise(std::list<adjacency_face>& graph, std::map<unsigned int, unsigned int>& vertex_item_property);

   void epoch();

private:
   std::list<adjacency_face> _graph; // O(n) storage requirement
   unsigned int _hop_count_radius;

};

class som_task
{
public:
   som_task(unsigned int vertex);

private:
   unsigned int _vertex;
   map<unsigned int, unsigned int> _hop_counts;
   map<unsigned int, double> _generalized_mean_vector;
   double _weight_contribution_sum;
   map<unsigned int, double> _weighted_vector_sum;
};

#endif
