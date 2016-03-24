#include "som.h"

using namespace std;

som::som()
{
}


void som::initialise(list<adjacency_face>& planar_triangular_graph, map<unsigned int, unsigned int>& vertex_item_property)
   : _neighborhood_hop_count(1);
{
   // Create som_tasks for each vertex
   for (map<unsigned int, unsigned int>::const_iterator itr = vertex_item_property.begin(), end = vertex_item_property.end(); itr != end; ++itr)
   {
      som_task st(itr->first, planar_triangular_graph);
   }
}

void som::run_epoch()
{
   if (_neighborhood_hop_count > 1) --_neighborhood_hop_count;

   for (map<unsigned int, unsigned int>::const_iterator itr = vertex_item_property.begin(), end = vertex_item_property.end(); itr != end; ++itr)
   {
      st.representation_stage();
   }

   // Exchange the weighted distances of nearby items and accumulate
   for (map<unsigned int, unsigned int>::const_iterator = vertex_item_property.begin(), end = vertex_item_property.end(); itr != end; ++itr)
   {
      st.affection_stage();
   }
}

som_task::som_task(unsigned int vertex, list<adjacency_face> &planar_triangular_graph) :
          _vertex(vertex)
{
   get_hop_counts(vertex, _hop_counts);
}


som_task::run_representation_stage()
{
   double min_distance = numeric_limits<double>::max();
   unsigned int winning_vertex = 0;

   // Calculate the distances of all nearby items to the generalized mean vector at this node
   for (map<unsigned int, unsigned int>::const_iterator itr = hop_counts.begin(), end = hop_counts.end(); itr != end; ++itr)
   {
      double distance = 0;
      map<unsigned int, double> vertex_vector = vertex_vectors[itr->first];
      if (vertex_vector.size() * 10 > _generalized_mean_vector.size())
      {
         // Merge join
         for (map<unsigned int, double>::const_iterator property_itr = vertex_vector.begin(),
                                                        property_end = vertex_vector.end(),
                                                        weights_itr = _generalized_mean_vector.begin(),
                                                        weights_end = _generalized_mean_vector.end(),
                                                        property_itr != property_end &&
                                                        weights_itr != property_end &&
                                                        distance < min_distance;)
         {
            if (property_itr->first == weights_itr->second)
            {
               double vertex_property_value = property_itr->second;
               double weights_property_value = weights_itr->second;
 
               distance += vertex_property_value * (vertex_property_value - 2 * weights_property_value);
               ++property_itr;
               ++weights_itr;
            }
            else if (property_itr->first < weights_itr->first) 
            {
               distance += vertex_property_value * vertex_property_value;
               ++property_itr;
            }
            else ++weights_itr;
         }
      }
      else
      {
         // Lookup join
         for (map<unsigned int, double>::const_iterator property_itr = vertex_vector.begin(), property_end = vertex_vector.end();
                                                        property_itr != property_end &&
                                                        distance < min_distance;
                                                        ++property_itr)
         {
            double vertex_property_value = property_itr->second;
            double weights_property_value = 0;
 
            map<unsigned int, double>::const_iterator weights_property_itr = _generalized_mean_vector.find(property_itr->first);
            if (weights_property_itr != _generalized_mean_vector.end()) weights_property_value = weights_property_value.second;
 
            distance += vertex_property_value * (vertex_property_value - 2 * weights_property_value);
         }
      }
 
      // Compute the winning item that has the minimum distance
      if (distance < min_distance)
      {
         distance = min_distance;
         winning_vertex = itr->first;
      }
   }
 
   // Calculate the winning item contribution to the weighted mean of nearby items
   for (map<unsigned int, unsigned int>::const_iterator itr = hop_counts.begin(), end = hop_counts.end(); itr != end; ++itr)
   {
      unsigned int hop_count = itr->second;
      double weight = exp(-double(hop_count * hop_count) / double(_hop_count_radius * _hop_count_radius));
      
      if (weight > _min_weight)
      {
         send_weight_contribution(itr->first, weight);
         send_weighted_vector_contribution(itr->first, weight * winning_vertex_vector);
      }
   }
}

som_task::send_weight_contribution(unsigned int vertex, double weight)
{
}

som_task::receive_weight_contribution(double weight)
{
   _weight_contribution_sum += weight;
}

som_task::send_weighted_vector_contribution(unsigned int vertex, map<unsigned int, double> &weighted_vector)
{
}

som_task::receive_weight_contribution(map<unsigned int, double> &weighted_vector)
{
   _weighted_vector_sum += weighted_vector;
}
 
som_task::run_affection_stage()
{
   // Calculate the generalized mean vector at this node, given the accumulations
   _generalized_mean_vector = _weighted_vector_sum / _weight_contribution_sum;
}

