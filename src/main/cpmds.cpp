#include "cpmds.h"
#include <set>
#include <map>
#include <list>
#include <stack>
#include <queue>
#include <fstream>
#include <numeric>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/chrobak_payne_drawing.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/is_straight_line_drawing.hpp>
#include <boost/graph/planar_canonical_ordering.hpp>
#include <boost/graph/make_biconnected_planar.hpp>
#include <boost/graph/make_maximal_planar.hpp>

using namespace std;
using namespace boost;

// A circular packing MDS algorithm
cpmds::cpmds() : _flatten_item(0), _graph_iteration(0), _focus(125), _flattened_curvature_sum(0)
{
}

void cpmds::process_relation(unsigned int first_itm_id, unsigned int second_itm_id, double first_item_value, double second_item_value)
{
   cout << "Processing relation " << first_itm_id << "(" << first_item_value << ") and " <<
                                    second_itm_id << "(" << second_item_value << ")" << endl;
   if (first_itm_id == second_itm_id) return;
   if (first_itm_id > second_itm_id) swap(first_itm_id, second_itm_id);

   // Process one side at a time
   update_relation(first_itm_id, second_itm_id, first_item_value, second_item_value);
   update_item(first_itm_id, first_item_value, second_item_value);

   update_relation(second_itm_id, first_itm_id, second_item_value, first_item_value);
   update_item(second_itm_id, second_item_value, first_item_value);
}

// Utility for helping testing
const item_table_type &cpmds::items() const
{
   return item_table;
}

bool cpmds::has_neighbour(unsigned int id, item& neighbour) const
{
   
   list<unsigned int> perimeter;
   get_perimeter(id, perimeter);
   bool found = false;
   if (!perimeter.empty())
   {
      const item_pk& itm_pk = item_table.get<item_pk_type>();
      list<unsigned int>::const_iterator neighbour_itr = perimeter.begin();
      if (*neighbour_itr == 0) ++neighbour_itr;
      item_pk::const_iterator itm_itr = itm_pk.find(*neighbour_itr);
      neighbour = *itm_itr;
      found = true;
   }

   return found;
}

void cpmds::update_item(unsigned int id, double value, double compare_value)
{
   item i;

   item_pk& itm_pk = item_table.get<item_pk_type>();
   item_pk::const_iterator itm_itr = itm_pk.find(id);
   if (itm_itr == itm_pk.end())
   {
      i.id = id;
      i.location = point(10e23, 10e23);
      i.radius = 0;
      i.sample_count = 0;
      i.value_mean = 0;
      i.value_variance = 0;
      i.curvature = 2 * M_PI;
      i.is_located = false;
      i.curvature_delta = 0;
   }
   else
   { 
      // Update the item
      i = *itm_itr;
   }
   
   i.sample_count++;

   double value_mean = i.value_mean;
   i.value_mean = i.value_mean + (value - i.value_mean) / i.sample_count;

   double qsum_delta = (value - value_mean) * (value - i.value_mean);
   i.value_variance = i.value_variance + (qsum_delta - i.value_variance) / i.sample_count;

   if (itm_itr == itm_pk.end())
   {
      itm_pk.insert(i);
      map<unsigned int, item_adjustment> adjustments;
      insert_item(id, adjustments, true);
   }
   else
   {
      itm_pk.replace(itm_itr, i);
      on_radius_change(id, sqrt(itm_itr->sample_count));
   }
}

void cpmds::insert_item(unsigned int id, map<unsigned int, item_adjustment> &adjustments, bool update)
{
   if (update) adjustments.clear();

   // cout << "Inserting item " << id << endl;
   adjacency_face_pk& adf_pk = adjacency_face_table.get<adjacency_face_pk_type>();

   if (adjacency_face_table.empty())
   {
      // The very first faces 
      if (item_table.size() == 2)
      {
         //cout << "Adjacency face table size: " << adjacency_face_table.size() << endl;
         item_pk& itm_pk = item_table.get<item_pk_type>();
         unsigned int a_id = id;
         unsigned int b_id = item_table[0].id;
         item_pk::iterator a_itm_itr = itm_pk.find(a_id);
         item_pk::iterator b_itm_itr = itm_pk.find(b_id);
         item i_a = *a_itm_itr; 
         item i_b = *b_itm_itr; 

         i_a.connectivity += 1;
         i_b.connectivity += 1;
         i_a.curvature -= 2 * subtended_half_angle(i_a.radius, i_b.radius);
         i_b.curvature -= 2 * subtended_half_angle(i_b.radius, i_a.radius);

         pair<double, double> distance_weight = get_weighted_measure(a_id, b_id, WEIGHTED_DISTANCE);
         i_a.weight_sum += distance_weight.first;
         i_a.distance_weight_sum += distance_weight.second * distance_weight.first;
         i_b.weight_sum += distance_weight.first;
         i_b.distance_weight_sum += distance_weight.second * distance_weight.first;

         if (update)
         { 
           adf_pk.insert(make_face(item_table[0].id, id, 0));
           adf_pk.insert(make_face(id, item_table[0].id, 0));

           itm_pk.replace(a_itm_itr, i_a);
           itm_pk.replace(b_itm_itr, i_b);
         }
      }
   }
   else
   {
      adjacency_face_pk::iterator f_nearest_itr = find_k_nearest_face(id, 0);
      insert_item(id, f_nearest_itr, adjustments, update);
   }
}
      
void cpmds::insert_item(unsigned int id, adjacency_face_pk::iterator f_nearest_itr, map<unsigned int, item_adjustment> &adjustments, bool update)
{
   if (update) adjustments.clear();

   adjacency_face_pk& adf_pk = adjacency_face_table.get<adjacency_face_pk_type>();

   unsigned int a_id = f_nearest_itr->first_itm_id;
   unsigned int b_id = f_nearest_itr->second_itm_id;
   unsigned int c_id = f_nearest_itr->third_itm_id;

   // Look up the items
   item_pk& itm_pk = item_table.get<item_pk_type>();
   item_pk::iterator itm_itr = itm_pk.find(id);
   item_pk::iterator a_itr = itm_pk.find(a_id);
   item_pk::iterator b_itr = itm_pk.find(b_id);

   if (update)
   {
      // Create three new faces
      adf_pk.erase(f_nearest_itr);
      adf_pk.insert(make_face(a_id, b_id, id));
      adf_pk.insert(make_face(a_id, id, c_id));
      adf_pk.insert(make_face(id, b_id, c_id));
   }

   // Update the the angle sum for every item in the divided face
   double r = itm_itr->radius;
   double r_a = a_itr->radius;
   double r_b = b_itr->radius;

   item i = *itm_itr;
   item i_a = *a_itr;
   item i_b = *b_itr;

   pair<double, double> a_distance_weight = get_weighted_measure(id, a_id, WEIGHTED_DISTANCE);
   adjustments[id].curvature_delta -= subtended_angle(r, r_a, r_b);
   adjustments[id].connectivity_delta += 1;
   adjustments[id].weight_sum_delta += a_distance_weight.first;
   adjustments[id].distance_weight_sum_delta += a_distance_weight.second * a_distance_weight.first;

   adjustments[a_id].curvature_delta -= subtended_angle(r_a, r, r_b);
   adjustments[a_id].connectivity_delta += 1;
   adjustments[a_id].weight_sum_delta += a_distance_weight.first;
   adjustments[a_id].distance_weight_sum_delta += a_distance_weight.second * a_distance_weight.first;

   pair<double, double> b_distance_weight = get_weighted_measure(id, b_id, WEIGHTED_DISTANCE);
   adjustments[id].connectivity_delta += 1;
   adjustments[id].weight_sum_delta += b_distance_weight.first;
   adjustments[id].distance_weight_sum_delta += b_distance_weight.second * b_distance_weight.first;

   adjustments[b_id].curvature_delta -= subtended_angle(r_b, r, r_a);
   adjustments[b_id].connectivity_delta += 1;
   adjustments[b_id].weight_sum_delta += b_distance_weight.first;
   adjustments[b_id].distance_weight_sum_delta += b_distance_weight.second * b_distance_weight.first;

   if (c_id != 0)
   {

      item_pk::iterator c_itr = itm_pk.find(c_id);
      item i_c = *c_itr;
      double r_c = c_itr->radius;

      adjustments[id].curvature_delta -= subtended_angle(r, r_a, r_c)
                                       + subtended_angle(r, r_b, r_c);

      adjustments[a_id].curvature_delta -= subtended_angle(r_a, r, r_c)
                                         - subtended_angle(r_a, r_b, r_c);

      adjustments[b_id].curvature_delta -= subtended_angle(r_b, r, r_c)
                                         - subtended_angle(r_b, r_a, r_c);

      adjustments[c_id].curvature_delta -= subtended_angle(r_c, r, r_a)
                                         + subtended_angle(r_c, r, r_b)
                                         - subtended_angle(r_c, r_b, r_a);

      pair<double, double> c_distance_weight = get_weighted_measure(id, c_id, WEIGHTED_DISTANCE);
      adjustments[id].connectivity_delta += 1;
      adjustments[id].weight_sum_delta += c_distance_weight.first;
      adjustments[id].distance_weight_sum_delta += c_distance_weight.second * c_distance_weight.first;

      adjustments[c_id].connectivity_delta += 1;
      adjustments[c_id].weight_sum_delta += c_distance_weight.first;
      adjustments[c_id].distance_weight_sum_delta += c_distance_weight.second * c_distance_weight.first;

      if (update) itm_pk.replace(c_itr, i_c);
   }
   else
   {
      adjustments[id].curvature_delta -= subtended_half_angle(r, r_a)
                                       + subtended_half_angle(r, r_b);

      adjustments[a_id].curvature_delta -= subtended_half_angle(r_a, r)
                                         - subtended_half_angle(r_a, r_b);

      adjustments[b_id].curvature_delta -= subtended_half_angle(r_b, r)
                                         - subtended_half_angle(r_b, r_a);

   }

   if (update)
   {
      itm_pk.replace(itm_itr, i);
      itm_pk.replace(a_itr, i_a);
      itm_pk.replace(b_itr, i_b);

      // Update the perimeter item's angle sum
      for (map<unsigned int, item_adjustment>::const_iterator itr = adjustments.begin(), end = adjustments.end(); itr != end; ++itr)
      {
         itm_itr = itm_pk.find(itr->first);
         i = *itm_itr;
         i.curvature += itr->second.curvature_delta;
         i.curvature_delta += itr->second.curvature_delta;
         i.connectivity += itr->second.connectivity_delta;
         i.distance_weight_sum += itr->second.distance_weight_sum_delta;
         i.weight_sum += itr->second.weight_sum_delta;
         itm_pk.replace(itm_itr, i);
      }
   }   
}

// Find the nearest face to an item, given a k-hop region average
adjacency_face_pk::const_iterator cpmds::find_k_nearest_face(unsigned int id, unsigned int k) const
{
   // cout << "Inserting item " << id << endl;
   const adjacency_face_pk& adf_pk = adjacency_face_table.get<adjacency_face_pk_type>();

   adjacency_face_pk::const_iterator f_nearest_itr = adf_pk.end();

   unsigned int nearest_id = 0;
   if (k == 0)
   {
      // Find the nearest vertex
      const relation_distance_dist_idx& rel_distance_dist_idx = relation_distance_table.get<relation_distance_dist_idx_type>();
      relation_distance_dist_idx::const_iterator rel_nearest_itr = rel_distance_dist_idx.lower_bound(boost::make_tuple(id));
      nearest_id = rel_nearest_itr->second_itm_id;
   }
   else
   {
      // Find the nearest vertex, given a k-hop region average
      double min_distance = numeric_limits<double>::max();
      const item_pk& itm_pk = item_table.get<item_pk_type>();
      for (item_pk::const_iterator itr = itm_pk.begin(), end = itm_pk.end(); itr != end; ++itr)
      {
         pair<double, double> sub_measure = get_weighted_measure(itr->id, WEIGHTED_DISTANCE, AVG, id, k);
         if (sub_measure.second < min_distance)
         {
            min_distance = sub_measure.second;
            nearest_id = itr->id;
         }
      }
   }

   list<unsigned int> perimeter;
   get_perimeter(nearest_id, perimeter);

   double lowest_stress_delta = numeric_limits<double>::max();
   for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
   {
      list<unsigned int>::const_iterator next_itr = itr;
      ++next_itr;
      if (next_itr == perimeter.end()) next_itr = perimeter.begin();
     
      adjacency_face f = make_face(nearest_id, *itr, *next_itr);
      adjacency_face_pk::const_iterator f_itr = adf_pk.find(boost::make_tuple(f.first_itm_id, f.second_itm_id));
      
      map<unsigned int, item_adjustment> adjustments;
      const_cast<cpmds*>(this)->insert_item(id, f_itr, adjustments, false);

      double stress_delta = get_stress_delta(adjustments);

      if (stress_delta <= lowest_stress_delta)
      {
         lowest_stress_delta = stress_delta;
         f_nearest_itr = f_itr;
      }
     
   }

   return f_nearest_itr;
}

void cpmds::on_radius_change(unsigned int id, double r)
{
   list<unsigned int> perimeter;
   get_perimeter(id, perimeter);

   // Check it the first and last item in the perimeter are connected
   list<unsigned int>::const_iterator perimeter_edge_end_itr = perimeter.begin();

   if (!perimeter.empty())
   {

      item_pk& itm_pk = item_table.get<item_pk_type>();
      item_pk::iterator itm_itr = itm_pk.find(id);

      double r_old = itm_itr->radius;

      // Calculate the angle adjustments
      double central_curvature = 2 * M_PI;
      map<unsigned int, item_adjustment> adjustments;
      for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
      {
         if (*itr != 0)
         {
            item_adjustment ia;
            ia.id = *itr;
            ia.curvature_delta = 0;
            ia.weight_sum_delta = 0;
            ia.distance_weight_sum_delta = 0;
            ia.connectivity_delta = 0;
            adjustments[*itr] = ia;
         }
      }

      for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
      {
         list<unsigned int>::const_iterator prev_itr = itr;
         if (prev_itr == perimeter.begin())
            prev_itr = perimeter.end();
         --prev_itr;

         if (*itr != 0)
         {
            item_pk::iterator a_itm_itr = itm_pk.find(*itr);
            double r_a = a_itm_itr->radius;

            if (*prev_itr != 0)
            {
               item_pk::iterator b_itm_itr = itm_pk.find(*prev_itr);
               double r_b = b_itm_itr->radius;

               central_curvature -= subtended_angle(r, r_a, r_b); 
               adjustments[*itr].curvature_delta -= subtended_angle(r_a, r, r_b)
                                                  - subtended_angle(r_a, r_old, r_b);
               adjustments[*prev_itr].curvature_delta -= subtended_angle(r_b, r, r_a)
                                                       - subtended_angle(r_b, r_old, r_a);
            }
            else
            {
               central_curvature -= subtended_half_angle(r, r_a); 
               adjustments[*itr].curvature_delta -= subtended_half_angle(r_a, r)
                                                  - subtended_half_angle(r_a, r_old);
            }
         }
         else if (*prev_itr != 0)
         {
            item_pk::iterator b_itm_itr = itm_pk.find(*prev_itr);
            double r_b = b_itm_itr->radius;

            central_curvature -= subtended_half_angle(r, r_b); 
            adjustments[*prev_itr].curvature_delta -= subtended_half_angle(r_b, r)
                                                    - subtended_half_angle(r_b, r_old);
         }
      }

      // Update the items with the adjustments
      //cout << "Angle adjustments for size change: " << endl;
      item i = *itm_itr;
      i.curvature = central_curvature;
      i.radius = r;
      itm_pk.replace(itm_itr, i);

      for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
      {
         if (*itr != 0)
         {
            itm_itr = itm_pk.find(*itr);
            i = *itm_itr;
            i.curvature += adjustments[*itr].curvature_delta;
            i.curvature_delta += adjustments[*itr].curvature_delta;
            itm_pk.replace(itm_itr, i);
         }
      }
   }
}

void cpmds::get_increase_curvature_transform(unsigned int id, list<unsigned int>& fan_trfm_perimeter) const
{
   const adjacency_face_pk& adf_pk = adjacency_face_table.get<adjacency_face_pk_type>();
   double lowest_stress_delta = numeric_limits<double>::max();

   list<unsigned int> perimeter;
   get_perimeter(id, perimeter);

   // Push out the curvature, find the least similar adjacent node
   for (list<unsigned int>::const_iterator pivot_itr = perimeter.begin(), end = perimeter.end(); pivot_itr != end; ++pivot_itr)
   {
      // Seek to disconnect [id] from [pivot]
      bool potential_fan_trfm_clockwise = true;
      do
      {
         list<unsigned int>::const_iterator other_itr = pivot_itr;
         if (potential_fan_trfm_clockwise) 
         {
            if (other_itr == perimeter.begin()) other_itr = perimeter.end();
            --other_itr;
         }
         else
         {
            ++other_itr;
            if (other_itr == perimeter.end()) other_itr = perimeter.begin();
         }

         // The curvature of [pivot] and [other] will be significantly changed
         if (_flattened_items.find(*other_itr) == _flattened_items.end() &&
             _flattened_items.find(*pivot_itr) == _flattened_items.end())
         {

            list<unsigned int> pivot_perimeter;
            get_perimeter(*pivot_itr, pivot_perimeter);
            if (potential_fan_trfm_clockwise) reverse(pivot_perimeter.begin(), pivot_perimeter.end());

            list<unsigned int>::iterator pivot_id_itr = find(pivot_perimeter.begin(), pivot_perimeter.end(), id);
            rotate(pivot_perimeter.begin(), pivot_id_itr, pivot_perimeter.end());

            list<unsigned int> potential_fan_trfm_perimeter;
            potential_fan_trfm_perimeter.push_back(*pivot_itr);
            if (potential_fan_trfm_clockwise)
            {  
               potential_fan_trfm_perimeter.push_back(id);
               potential_fan_trfm_perimeter.push_back(*other_itr);
            }
            else
            {
               potential_fan_trfm_perimeter.push_back(*other_itr);
               potential_fan_trfm_perimeter.push_back(id);
            }

            adjacency_face f_tmp;
            for (list<unsigned int>::iterator pivot_last_itr = ++pivot_perimeter.begin(), pivot_other_itr = --pivot_perimeter.end(); pivot_last_itr != pivot_other_itr && potential_fan_trfm_perimeter.size() < 5; ++pivot_last_itr)
            {
               // Ensure that the joined neighbours of [id] are not already connected (separating triangle)
               // and that they are not both on the flattened boundary
               if (find_adjacency_face(*other_itr, *pivot_last_itr, f_tmp) != adf_pk.end()
                   || (_flattened_boundary_items.find(*other_itr) != _flattened_boundary_items.end() 
                       && _flattened_boundary_items.find(*pivot_last_itr) != _flattened_boundary_items.end()))
               {
                  break;
               }
            
               if (potential_fan_trfm_clockwise)
               {
                  potential_fan_trfm_perimeter.insert(++potential_fan_trfm_perimeter.begin(), *pivot_last_itr);
               }
               else
               {
                  potential_fan_trfm_perimeter.push_back(*pivot_last_itr);
               }

               if (_flattened_items.find(*pivot_last_itr) == _flattened_items.end())
               {
                  // TODO Tidy this up if possible
                  if (potential_fan_trfm_clockwise)
                     rotate(potential_fan_trfm_perimeter.begin(),
                            --potential_fan_trfm_perimeter.end(),
                            potential_fan_trfm_perimeter.end());
                  else
                     rotate(potential_fan_trfm_perimeter.begin(),
                            ++potential_fan_trfm_perimeter.begin(),
                            potential_fan_trfm_perimeter.end());

                  bool can_adjust = true;
                  map<unsigned int, item_adjustment> adjustments;
                  fan_transform(potential_fan_trfm_perimeter, adjustments);
                  if (potential_fan_trfm_perimeter.size() > 4)
                  {
                     // There will be a small adjustment to the flattened item, check to see if it is acceptable
                     list<unsigned int>::const_iterator minor_adjust_itr = pivot_last_itr;
                     --minor_adjust_itr;

                     if (_flattened_items.find(*minor_adjust_itr) != _flattened_items.end())
                     {

                        if (!is_beneficial_adjustment(*minor_adjust_itr, adjustments[*minor_adjust_itr].curvature_delta))
                        {
                           can_adjust = false;
                        }
                     }
                  }

                  if (can_adjust)
                  {
                     double stress_delta = get_stress_delta(adjustments, CURVATURE_STRESS);

                     if (stress_delta <= lowest_stress_delta)
                     {
                        lowest_stress_delta = stress_delta;
                        //cout << "Lowest stress delta - " << lowest_stress_delta << endl;
                        fan_trfm_perimeter.clear();
                        copy(potential_fan_trfm_perimeter.begin(), potential_fan_trfm_perimeter.end(), back_inserter(fan_trfm_perimeter));
                     }
                  }

                  if (potential_fan_trfm_clockwise)
                     rotate(potential_fan_trfm_perimeter.begin(),
                            ++potential_fan_trfm_perimeter.begin(),
                            potential_fan_trfm_perimeter.end());
                  else
                     rotate(potential_fan_trfm_perimeter.begin(),
                            --potential_fan_trfm_perimeter.end(),
                            potential_fan_trfm_perimeter.end());


               }
            }
         }

         potential_fan_trfm_clockwise = !potential_fan_trfm_clockwise;
      }
      while (potential_fan_trfm_clockwise == false);
   }
}

void cpmds::get_decrease_curvature_transform(unsigned int id, list<unsigned int>& fan_trfm_perimeter) const
{
   const adjacency_face_pk& adf_pk = adjacency_face_table.get<adjacency_face_pk_type>();

   // Store a list of all processed faces in a list
   vector<adjacency_face> processed_faces;
   list<unsigned int> candidate_faces;
   unsigned int last_iteration_face = 0;
   map<unsigned int, unsigned int> child_parent_path_map;
   set<unsigned int> discovered_vertices;

   // Start by including all the perimeter faces of the item that are modifiable
   list<unsigned int> perimeter;
   get_perimeter(id, perimeter);

   unsigned int prior_discovered_vertices_size = 0;
   for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
   {
      list<unsigned int>::const_iterator next_itr = itr;
      ++next_itr;
      if (next_itr == perimeter.end()) next_itr = perimeter.begin();

      if (_flattened_items.find(*itr) == _flattened_items.end() &&
          _flattened_items.find(*next_itr) == _flattened_items.end())
      {
         adjacency_face f_tmp;
         find_adjacency_face(*next_itr, *itr, f_tmp);

         processed_faces.push_back(f_tmp);
         discovered_vertices.insert(f_tmp.third_itm_id);
      }
   }

   if (discovered_vertices.empty())
   {
      // Handle the case where the two boundary neighbours are flattened exterior items
      map<unsigned int, list<unsigned int>::iterator>::const_iterator boundary_itr = _flattened_boundary_items.find(id);
      if (boundary_itr != _flattened_boundary_items.end())
      {
         list<unsigned int>::const_iterator itr = boundary_itr->second;

         list<unsigned int>::const_iterator prev_itr = itr;
         if (prev_itr == _flattened_boundary.begin()) prev_itr = _flattened_boundary.end();
         --prev_itr;

         list<unsigned int>::const_iterator next_itr = itr;
         ++next_itr;
         if (next_itr == _flattened_boundary.end()) next_itr = _flattened_boundary.begin();

         if (!is_interior_item(*prev_itr) && !is_interior_item(*next_itr))
         {
            // We can convert this item to an exterior item, it will just incrase the curvature of the
            //flattened exterior items which is fine.
            fan_trfm_perimeter.push_back(0);
            fan_trfm_perimeter.push_back(*next_itr);
            fan_trfm_perimeter.push_back(id);
            fan_trfm_perimeter.push_back(*prev_itr);
         }
      }
   }
   else
   {
      // Process each generation of faces until a non-boundary item or edge is found
      while (candidate_faces.size() < 3 && discovered_vertices.size() > prior_discovered_vertices_size )
      {
         // Check to see if any of the discovered faces present a usable item
         prior_discovered_vertices_size = discovered_vertices.size();
         for (unsigned int i = last_iteration_face, end = processed_faces.size(); i != end; ++i)
         {
            adjacency_face f_current = processed_faces[i];
            adjacency_face f_tmp; 
         
            if (f_current.third_itm_id == 0 ||
                _flattened_items.find(f_current.third_itm_id) == _flattened_items.end())
            {
               if ((f_current.third_itm_id == 0 ||
                   _flattened_boundary_items.find(f_current.third_itm_id) == _flattened_boundary_items.end())
                   //&& is_interior_item(f_current.third_itm_id)
                   && find_adjacency_face(f_current.third_itm_id, id, f_tmp) == adf_pk.end())
               {
                  candidate_faces.push_back(i);
               }

               find_adjacency_face(f_current.third_itm_id, f_current.second_itm_id, f_tmp);
               if (discovered_vertices.find(f_tmp.third_itm_id) == discovered_vertices.end())
               {
                  processed_faces.push_back(f_tmp);
                  child_parent_path_map.insert(make_pair(processed_faces.size() - 1, i));
                  discovered_vertices.insert(f_tmp.third_itm_id);
               }
               
               find_adjacency_face(f_current.first_itm_id, f_current.third_itm_id, f_tmp);
               if (discovered_vertices.find(f_tmp.third_itm_id) == discovered_vertices.end())
               {
                  processed_faces.push_back(f_tmp);
                  child_parent_path_map.insert(make_pair(processed_faces.size() - 1, i));
                  discovered_vertices.insert(f_tmp.third_itm_id);
               }
            }

            ++last_iteration_face;
         }
      }

      // Find the path that if transformed produces the lowest stress
      double lowest_stress_delta = numeric_limits<double>::max();
      for (list<unsigned int>::const_iterator itr = candidate_faces.begin(), end = candidate_faces.end(); itr != end; ++itr)
      {
         // Transform each path into a transform
         list<unsigned int> potential_fan_trfm_perimeter;
         list<unsigned int> potential_fan_trfm_clockwise;
         unsigned int child_id = *itr;

         adjacency_face f_first = processed_faces[child_id];
         adjacency_face f_current = f_first;
         potential_fan_trfm_perimeter.push_back(f_current.third_itm_id);
         potential_fan_trfm_perimeter.push_back(f_current.first_itm_id);
         potential_fan_trfm_clockwise.push_back(f_current.second_itm_id);

         // Build the transform perimeter from the face hierarchy
         map<unsigned int, unsigned int>::const_iterator path_itr = child_parent_path_map.find(child_id);
         while (path_itr != child_parent_path_map.end())
         {
            adjacency_face f_last = f_current;
            child_id = path_itr->second;
            f_current = processed_faces[child_id];

            if (f_current.third_itm_id == f_last.first_itm_id)
               potential_fan_trfm_perimeter.push_back(f_current.first_itm_id);
            else
               potential_fan_trfm_clockwise.push_back(f_current.second_itm_id);

            path_itr = child_parent_path_map.find(child_id);
         }

         potential_fan_trfm_clockwise.push_back(id);
         copy(potential_fan_trfm_clockwise.rbegin(),
              potential_fan_trfm_clockwise.rend(),
              back_inserter(potential_fan_trfm_perimeter));

         double stress_delta = 0;
         if (candidate_faces.size() > 1)
         {
            if (f_first.third_itm_id == 0)
            {
               // Always convert to an exterior item as the last resort 
               stress_delta = numeric_limits<double>::max();
            }
            else
            {
               // Find the best transform by calculating the induced stress
               map<unsigned int, item_adjustment> adjustments;
               fan_transform(potential_fan_trfm_perimeter, adjustments);

               stress_delta = get_stress_delta(adjustments, CURVATURE_STRESS);
            }
         }

         if (stress_delta <= lowest_stress_delta)
         {
            lowest_stress_delta = stress_delta;
            cout << "Lowest stress delta - " << lowest_stress_delta << endl;
            fan_trfm_perimeter.clear();
            copy(potential_fan_trfm_perimeter.begin(), potential_fan_trfm_perimeter.end(), back_inserter(fan_trfm_perimeter));
         }
      }
   }
}

void cpmds::get_perimeter(unsigned int id, list<unsigned int>& perimeter, unsigned int start_id, unsigned int end_id) const
{
   // First find all faces connected to this item
   const adjacency_face_pk& adf_pk = adjacency_face_table.get<adjacency_face_pk_type>();
   const adjacency_face_second_itm_idx& adf_second_itm_idx = adjacency_face_table.get<adjacency_face_second_itm_idx_type>();
   const adjacency_face_third_itm_idx& adf_third_itm_idx = adjacency_face_table.get<adjacency_face_third_itm_idx_type>();

   map<unsigned int, unsigned int> perimeter_edges;
   pair<adjacency_face_pk::const_iterator, adjacency_face_pk::const_iterator> first_range =
      adf_pk.equal_range(boost::make_tuple(id));
   for (; first_range.first != first_range.second; ++first_range.first)
   {
      perimeter_edges.insert(make_pair(first_range.first->second_itm_id, first_range.first->third_itm_id));
   }

   pair<adjacency_face_second_itm_idx::const_iterator, adjacency_face_second_itm_idx::const_iterator> second_range =
      adf_second_itm_idx.equal_range(boost::make_tuple(id));
   for (; second_range.first != second_range.second; ++second_range.first)
   {
      perimeter_edges.insert(make_pair(second_range.first->third_itm_id, second_range.first->first_itm_id));
   }

   pair<adjacency_face_third_itm_idx::const_iterator, adjacency_face_third_itm_idx::const_iterator> third_range =
      adf_third_itm_idx.equal_range(boost::make_tuple(id));
   for (; third_range.first != third_range.second; ++third_range.first)
   {
      perimeter_edges.insert(make_pair(third_range.first->first_itm_id, third_range.first->second_itm_id));
   }

   // Order the perimeter items cyclically
   if (!perimeter_edges.empty())
   {
      map<unsigned int, unsigned int>::iterator boundary_itm_itr = perimeter_edges.find(start_id);
      if (boundary_itm_itr != perimeter_edges.end())
      {
         // Start with the item before the boundary
         perimeter.push_back(start_id);
      }
      else 
      {
         perimeter.push_back(perimeter_edges.begin()->first);
      }
   }

   while (perimeter.size() < perimeter_edges.size())
   {
      perimeter.push_back(perimeter_edges[perimeter.back()]);
      if (perimeter.back() == end_id) break;
   }
}

void cpmds::remove_item(unsigned int id, map<unsigned int, item_adjustment> &adjustments, bool update)
{
   if (update) adjustments.clear();

   list<unsigned int> perimeter;
   get_perimeter(id, perimeter);

   adjacency_face_pk& adf_pk = adjacency_face_table.get<adjacency_face_pk_type>();

   if (update)
   {
      // Erase all the faces connected to this id
      for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
      {
         adjacency_face f_tmp;
         adf_pk.erase(find_adjacency_face(id, *itr, f_tmp));
      }
   }

   // Update the angle sum for all the perimeter items
   item_pk& itm_pk = item_table.get<item_pk_type>();
   item_pk::iterator itm_itr = itm_pk.find(id);
   double r = itm_itr->radius;
  
   {
      item_adjustment ia;
      ia.id = id;
      ia.curvature_delta = -itm_itr->curvature + 2 * M_PI;
      ia.weight_sum_delta = -itm_itr->weight_sum;
      ia.distance_weight_sum_delta = -itm_itr->distance_weight_sum;
      ia.connectivity_delta = -itm_itr->connectivity;
      adjustments[id] = ia;
   }

   if (perimeter.size() > 1)
   {
      for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
      {
         if (*itr != 0) 
         {
            item_adjustment ia;
            ia.id = *itr;
            ia.curvature_delta = 0;
            ia.weight_sum_delta = 0;
            ia.distance_weight_sum_delta = 0;
            ia.connectivity_delta = 0;
            adjustments[*itr] = ia;
         }
      }

      for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
      {
         list<unsigned int>::const_iterator next_itr = itr;
         ++next_itr;
         if (next_itr == perimeter.end()) next_itr = perimeter.begin();

         if (*itr != 0)
         {
            item_pk::iterator a_itm_itr = itm_pk.find(*itr);
            double r_a = a_itm_itr->radius;

            if (*next_itr != 0)
            {
               item_pk::iterator b_itm_itr = itm_pk.find(*next_itr);
               double r_b = b_itm_itr->radius;

               adjustments[*itr].curvature_delta      += subtended_angle(r_a, r, r_b);
               adjustments[*next_itr].curvature_delta += subtended_angle(r_b, r, r_a);
            }
            else
            {
               adjustments[*itr].curvature_delta      += subtended_half_angle(r_a, r);
            }

            pair<double, double> distance_weight = get_weighted_measure(id, *itr, WEIGHTED_DISTANCE);
            adjustments[*itr].connectivity_delta -= 1;
            adjustments[*itr].weight_sum_delta -= distance_weight.first;
            adjustments[*itr].distance_weight_sum_delta -= distance_weight.second * distance_weight.first;
         }
         else if (*next_itr != 0)
         {
            item_pk::iterator b_itm_itr = itm_pk.find(*next_itr);
            double r_b = b_itm_itr->radius;

            adjustments[*next_itr].curvature_delta += subtended_half_angle(r_b, r);
         }
      }
   }

   // Triangulate the left over void
   if (!perimeter.empty())
   {
      // Ensure the boundary edge (item id = 0) is the second item in the list, so it is *usually covered in the first pass
      // * this is not a guarantee, so care still has to be taken
      if (perimeter.front() == 0) rotate(perimeter.begin(), --perimeter.end(), perimeter.end());

      // Now triangulate the leftover hole evenly
      while (perimeter.size() > 2)
      {
         // Cover the middle items (odd indexes) with faces
         set<unsigned int> triangulated_items;
         for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end && boost::next(itr) != end; ++itr)
         {
            unsigned int a_id = *itr;
            unsigned int b_id = *boost::next(itr);
            unsigned int c_id = (boost::next(boost::next(itr)) != perimeter.end()) ? *boost::next(boost::next(itr)): *perimeter.begin();

            // Check that [i] and [i + 2] are not already connected
            adjacency_face f_tmp;
            if ((find_adjacency_face(c_id, a_id, f_tmp) == adf_pk.end() 
                || f_tmp.third_itm_id == id) &&
               (find_adjacency_face(a_id, c_id, f_tmp) == adf_pk.end()
                || f_tmp.third_itm_id != b_id))
            {
               if (update) adf_pk.insert(make_face(a_id, b_id, c_id));

               triangulated_items.insert(b_id);

               if (a_id != 0 && c_id != 0)
               {
                  item_pk::iterator a_itm_itr = itm_pk.find(a_id);
                  item_pk::iterator c_itm_itr = itm_pk.find(c_id);
                  double r_a = a_itm_itr->radius;
                  double r_c = c_itm_itr->radius;

                  // The final face does not create a new connection
                  if (perimeter.size() - triangulated_items.size() > 2)
                  {
                     pair<double, double> distance_weight = get_weighted_measure(a_id, c_id, WEIGHTED_DISTANCE);
                     adjustments[a_id].connectivity_delta += 1;
                     adjustments[a_id].weight_sum_delta += distance_weight.first;
                     adjustments[a_id].distance_weight_sum_delta += distance_weight.second * distance_weight.first;

                     adjustments[c_id].connectivity_delta += 1;
                     adjustments[c_id].weight_sum_delta += distance_weight.first;
                     adjustments[c_id].distance_weight_sum_delta += distance_weight.second * distance_weight.first;
                  }

                  // Update the angle sums
                  if (b_id != 0)
                  {
                     item_pk::iterator b_itm_itr = itm_pk.find(b_id);
                     double r_b = b_itm_itr->radius;
                     adjustments[a_itm_itr->id].curvature_delta -= subtended_angle(r_a, r_b, r_c);
                     adjustments[b_itm_itr->id].curvature_delta -= subtended_angle(r_b, r_c, r_a);
                     adjustments[c_itm_itr->id].curvature_delta -= subtended_angle(r_c, r_a, r_b);
                  }
                  else
                  {
                     adjustments[a_itm_itr->id].curvature_delta -= subtended_half_angle(r_a, r_c);
                     adjustments[c_itm_itr->id].curvature_delta -= subtended_half_angle(r_c, r_a);
                  }
               }
               else if (a_id != 0)
               {
                  item_pk::iterator a_itm_itr = itm_pk.find(a_id);
                  item_pk::iterator b_itm_itr = itm_pk.find(b_id);
                  double r_a = a_itm_itr->radius;
                  double r_b = b_itm_itr->radius;
                  adjustments[a_itm_itr->id].curvature_delta -= subtended_half_angle(r_a, r_b);
                  adjustments[b_itm_itr->id].curvature_delta -= subtended_half_angle(r_b, r_a);
               }
               else if (c_id != 0)
               {
                  item_pk::iterator b_itm_itr = itm_pk.find(b_id);
                  item_pk::iterator c_itm_itr = itm_pk.find(c_id);
                  double r_b = b_itm_itr->radius;
                  double r_c = c_itm_itr->radius;
                  adjustments[b_itm_itr->id].curvature_delta -= subtended_half_angle(r_b, r_c);
                  adjustments[c_itm_itr->id].curvature_delta -= subtended_half_angle(r_c, r_b);
               }

               ++itr;
            }
         }

         // The middle items (odd indexes) were covered, so remove them from the hole perimeter
         for (list<unsigned int>::iterator itr = perimeter.begin(); itr != perimeter.end();)
         {
            if (triangulated_items.find(*itr) != triangulated_items.end()) 
               itr = perimeter.erase(itr);
            else ++itr;
         }
      }
   } 

   if (update)
   {
      // Finally update the item's position
      item i = *itm_itr;
      i.location = point(10e23, 10e23);
      i.radius = 0;
      i.is_located = false;
      itm_pk.replace(itm_itr, i);

      // Update the perimeter item's angle sum
      for (map<unsigned int, item_adjustment>::const_iterator itr = adjustments.begin(), end = adjustments.end(); itr != end; ++itr)
      {
         itm_itr = itm_pk.find(itr->first);
         i = *itm_itr;
         i.curvature += itr->second.curvature_delta;
         i.curvature_delta += itr->second.curvature_delta;
         i.connectivity += itr->second.connectivity_delta;
         i.weight_sum += itr->second.weight_sum_delta;
         i.distance_weight_sum += itr->second.distance_weight_sum_delta;
         itm_pk.replace(itm_itr, i);
      }
   }
}


// Rotate the fanonals in a polygon from the first vertex to the last
void cpmds::fan_transform(const list<unsigned int> &perimeter, map<unsigned int, item_adjustment> &adjustments) const
{
   const_cast<cpmds*>(this)->fan_transform(perimeter, adjustments, false);
}

void cpmds::fan_transform(const list<unsigned int> &perimeter, map<unsigned int, item_adjustment> &adjustments, bool update)
{
   if (update) adjustments.clear();

   // Initialise the adjustments if necessary
   for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
   {
      if (*itr != 0 && adjustments.find(*itr) == adjustments.end()) 
      {
         item_adjustment ia;
         ia.id = *itr;
         ia.curvature_delta = 0;
         ia.weight_sum_delta = 0;
         ia.distance_weight_sum_delta = 0;
         ia.connectivity_delta = 0;
         adjustments[*itr] = ia;
      }
   }

   // First remove all the interior faces in the polygon
   adjacency_face_pk& adf_pk = adjacency_face_table.get<adjacency_face_pk_type>();
   item_pk& itm_pk = item_table.get<item_pk_type>();
   for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
   {
      list<unsigned int>::const_iterator next_itr = itr;
      ++next_itr;
      if (next_itr == perimeter.end()) next_itr = perimeter.begin();

      list<unsigned int>::const_iterator prev_itr = itr;
      if (prev_itr == perimeter.begin()) prev_itr = perimeter.end();
      --prev_itr;

      list<unsigned int>::const_iterator next_next_itr = next_itr;
      ++next_next_itr;
      if (next_next_itr == perimeter.end()) next_next_itr = perimeter.begin();
      
      adjacency_face f_inner;
      adjacency_face_pk::iterator f_inner_itr = find_adjacency_face(*itr, *next_itr, f_inner);

      if (f_inner_itr != adf_pk.end() && f_inner.third_itm_id != *prev_itr)
      {
         if (update)
         {
            adf_pk.erase(f_inner_itr);
    
            // Remove any boundary self joins that were deleted
            remove_boundary_self_join(f_inner.third_itm_id, *itr);

            if (f_inner.third_itm_id != *next_next_itr)
               remove_boundary_self_join(f_inner.third_itm_id, *next_itr);
         }

         unsigned int a_id = *itr;
         unsigned int b_id = *next_itr;
         unsigned int c_id = f_inner.third_itm_id;
         if (a_id != 0 && c_id != 0)
         {
            item_pk::iterator a_itm_itr = itm_pk.find(a_id);
            double r_a = a_itm_itr->radius;

            item_pk::iterator c_itm_itr = itm_pk.find(c_id);
            double r_c = c_itm_itr->radius;
          
            if (b_id != 0)
            {
               item_pk::iterator b_itm_itr = itm_pk.find(b_id);
               double r_b = b_itm_itr->radius;

               adjustments[a_id].curvature_delta += subtended_angle(r_a, r_c, r_b);
               adjustments[b_id].curvature_delta += subtended_angle(r_b, r_c, r_a);
               adjustments[c_id].curvature_delta += subtended_angle(r_c, r_a, r_b);
            }
            else
            {
               adjustments[a_id].curvature_delta += subtended_half_angle(r_a, r_c);
               adjustments[c_id].curvature_delta += subtended_half_angle(r_c, r_a);
            }

            // We discover this connection twice (from both directions), only process it once
            if (c_id > a_id)
            {
               //if (update) cout << "Disconnecting " << c_id << " and " << a_id << endl;
               pair<double, double> distance_weight = get_weighted_measure(c_id, a_id, WEIGHTED_DISTANCE);
               adjustments[c_id].connectivity_delta -= 1;
               adjustments[c_id].weight_sum_delta -= distance_weight.first;
               adjustments[c_id].distance_weight_sum_delta -= distance_weight.second * distance_weight.first;

               adjustments[a_id].connectivity_delta -= 1;
               adjustments[a_id].weight_sum_delta -= distance_weight.first;
               adjustments[a_id].distance_weight_sum_delta -= distance_weight.second * distance_weight.first;
            }
         }
         else if (a_id != 0 && b_id != 0)
         {
            item_pk::iterator a_itm_itr = itm_pk.find(a_id);
            double r_a = a_itm_itr->radius;

            item_pk::iterator b_itm_itr = itm_pk.find(b_id);
            double r_b = b_itm_itr->radius;
          
            adjustments[a_id].curvature_delta += subtended_half_angle(r_a, r_b);
            adjustments[b_id].curvature_delta += subtended_half_angle(r_b, r_a);
         }
         else if (b_id != 0 && c_id != 0)
         {
            item_pk::iterator b_itm_itr = itm_pk.find(b_id);
            double r_b = b_itm_itr->radius;

            item_pk::iterator c_itm_itr = itm_pk.find(c_id);
            double r_c = c_itm_itr->radius;
          
            adjustments[b_id].curvature_delta += subtended_half_angle(r_b, r_c);
            adjustments[c_id].curvature_delta += subtended_half_angle(r_c, r_b);
         }

         if (b_id != 0 && c_id != 0)
         {
            if (b_id > c_id && f_inner.third_itm_id != *next_next_itr)
            {
               //if (update) cout << "Disconnecting " << b_id << " and " << c_id << endl;
               pair<double, double> distance_weight = get_weighted_measure(b_id, c_id, WEIGHTED_DISTANCE);
               adjustments[b_id].connectivity_delta -= 1;
               adjustments[b_id].weight_sum_delta -= distance_weight.first;
               adjustments[b_id].distance_weight_sum_delta -= distance_weight.second * distance_weight.first;

               adjustments[c_id].connectivity_delta -= 1;
               adjustments[c_id].weight_sum_delta -= distance_weight.first;
               adjustments[c_id].distance_weight_sum_delta -= distance_weight.second * distance_weight.first;
            } 
         }
      }
   }

   // Retriangulate with the apex connecting to every other item on the boundary
   unsigned int apex_id = perimeter.front();
   for (list<unsigned int>::const_iterator itr = ++++perimeter.begin(),
                                             end = perimeter.end(),
                                             last_unconnected_itr = ++perimeter.begin(); itr != end; ++itr)
   {
      unsigned int first_itm_id, second_itm_id, third_itm_id;
      // Check that the apex is not already connected to the perimeter
      adjacency_face f_tmp;
      if (itr != --perimeter.end() &&
          find_adjacency_face(apex_id, *itr, f_tmp) != adf_pk.end())
      {
         list<unsigned int>::const_iterator next_itr = itr;
         ++next_itr;

         first_itm_id = *next_itr;
         second_itm_id = *last_unconnected_itr;
         third_itm_id = *itr;
      }
      else
      {
         first_itm_id = *itr;
         second_itm_id = apex_id;
         third_itm_id = *last_unconnected_itr;

         last_unconnected_itr = itr;
      }

      if (update) adf_pk.insert(make_face(first_itm_id, second_itm_id, third_itm_id));

      if (first_itm_id != 0 &&
          second_itm_id != 0)
      {
         item_pk::iterator a_itm_itr = itm_pk.find(first_itm_id);
         double r_a = a_itm_itr->radius;
         item_pk::iterator b_itm_itr = itm_pk.find(second_itm_id);
         double r_b = b_itm_itr->radius;

         if (third_itm_id != 0)
         {
            item_pk::iterator c_itm_itr = itm_pk.find(third_itm_id);
            double r_c = c_itm_itr->radius;
            adjustments[a_itm_itr->id].curvature_delta -= subtended_angle(r_a, r_b, r_c);
            adjustments[b_itm_itr->id].curvature_delta -= subtended_angle(r_b, r_c, r_a);
            adjustments[c_itm_itr->id].curvature_delta -= subtended_angle(r_c, r_a, r_b);
         }
         else
         {
            adjustments[a_itm_itr->id].curvature_delta -= subtended_half_angle(r_a, r_b);
            adjustments[b_itm_itr->id].curvature_delta -= subtended_half_angle(r_b, r_a);
         }

         if (itr != --perimeter.end())
         {
            pair<double, double> distance_weight = get_weighted_measure(first_itm_id, second_itm_id, WEIGHTED_DISTANCE);
            adjustments[first_itm_id].connectivity_delta += 1;
            adjustments[first_itm_id].weight_sum_delta += distance_weight.first;
            adjustments[first_itm_id].distance_weight_sum_delta += distance_weight.second * distance_weight.first;

            adjustments[second_itm_id].connectivity_delta += 1;
            adjustments[second_itm_id].weight_sum_delta += distance_weight.first;
            adjustments[second_itm_id].distance_weight_sum_delta += distance_weight.second * distance_weight.first;
         }
      }
      else if (first_itm_id != 0 && third_itm_id != 0)
      {
         item_pk::iterator a_itm_itr = itm_pk.find(first_itm_id);
         double r_a = a_itm_itr->radius;
         item_pk::iterator c_itm_itr = itm_pk.find(third_itm_id);
         double r_c = c_itm_itr->radius;

         adjustments[a_itm_itr->id].curvature_delta -= subtended_half_angle(r_a, r_c);
         adjustments[c_itm_itr->id].curvature_delta -= subtended_half_angle(r_c, r_a);
      }
      else if (second_itm_id != 0 && third_itm_id != 0)
      {
         item_pk::iterator b_itm_itr = itm_pk.find(second_itm_id);
         double r_b = b_itm_itr->radius;
         item_pk::iterator c_itm_itr = itm_pk.find(third_itm_id);
         double r_c = c_itm_itr->radius;

         adjustments[b_itm_itr->id].curvature_delta -= subtended_half_angle(r_b, r_c);
         adjustments[c_itm_itr->id].curvature_delta -= subtended_half_angle(r_c, r_b);
      }
   }

   if (update)
   {
      // Update the perimeter item's angle sum
      for (map<unsigned int, item_adjustment>::const_iterator itr = adjustments.begin(), end = adjustments.end(); itr != end; ++itr)
      {
         item_pk::iterator itm_itr = itm_pk.find(itr->first);
         item i = *itm_itr;
         i.curvature += itr->second.curvature_delta;
         i.curvature_delta += itr->second.curvature_delta;
         i.connectivity += itr->second.connectivity_delta;
         i.weight_sum += itr->second.weight_sum_delta;
         i.distance_weight_sum += itr->second.distance_weight_sum_delta;
         itm_pk.replace(itm_itr, i);
      }

      // Check for new boundary self joins
      // TODO - Is this necessary?
      for (list<unsigned int>::const_iterator itr = ++++perimeter.begin(), end = --perimeter.end(); itr != end; ++itr)
      {
         check_boundary_self_join(perimeter.front(), *itr);
      }
   }
}

void cpmds::update_relation(unsigned int first_itm_id, unsigned int second_itm_id, double first_item_value, double second_item_value)
{
   double distance = abs(first_item_value - second_item_value);

   relation_pk& relation_pk = relation_table.get<relation_pk_type>(); 
   relation_distance_pk& relation_distance_pk = relation_distance_table.get<relation_distance_pk_type>(); 

   // Create the relation record
   if (first_itm_id < second_itm_id) 
   {
      relation_pk::const_iterator r_itr = relation_pk.find(boost::make_tuple(first_itm_id, second_itm_id));
      if (r_itr == relation_table.end())
      {
         // Create the relation
         relation rel;
         rel.first_itm_id = first_itm_id;
         rel.second_itm_id = second_itm_id;
         rel.sample_count = 1;
         rel.distance = distance;
         rel.distance_variance = 0;
         rel.value_covariance = 0;
         relation_pk.insert(rel);
      }
      else
      {
         // Update the relation
         relation rel = *r_itr;
         rel.sample_count++;
         double value_mean = rel.distance;
         rel.distance = rel.distance + (distance - rel.distance) / rel.sample_count;

         double qsum_delta = (distance - value_mean) * (distance - rel.distance);
         rel.distance_variance = rel.distance_variance + (qsum_delta - rel.distance_variance) / rel.sample_count;
         relation_pk.replace(r_itr, rel);

         // Recalculate the covariance, note the items have the "old, n - 1" means
         item_pk& itm_pk = item_table.get<item_pk_type>();
         item_pk::const_iterator first_itm_itr = itm_pk.find(first_itm_id);
         item_pk::const_iterator second_itm_itr = itm_pk.find(second_itm_id);
         double value_product_old_mean =  rel.value_covariance + first_itm_itr->value_mean * second_itm_itr->value_mean;
         double value_product_mean = value_product_old_mean + ((first_item_value * second_item_value) - value_product_old_mean) / rel.sample_count;
         double first_value_mean = first_itm_itr->value_mean + (first_item_value - first_itm_itr->value_mean) / rel.sample_count;
         double second_value_mean = second_itm_itr->value_mean + (second_item_value - second_itm_itr->value_mean) / rel.sample_count;
 
         rel.value_covariance = value_product_mean - first_value_mean * second_value_mean;
      } 
   }

   // Update relation distance 
   relation_pk::const_iterator r_itr = relation_pk.find(boost::make_tuple(min(first_itm_id, second_itm_id), 
                                                                              max(first_itm_id, second_itm_id)));
   relation_distance_pk::const_iterator rd_itr = relation_distance_pk.find(boost::make_tuple(first_itm_id, second_itm_id));
   if (rd_itr != relation_distance_pk.end())
   {
      relation_distance rd = *rd_itr;
      rd.distance = r_itr->distance;
      relation_distance_pk.replace(rd_itr, rd);
   }
   else 
   {
      relation_distance rd;
      rd.first_itm_id = first_itm_id; 
      rd.second_itm_id = second_itm_id; 
      rd.distance = r_itr->distance;
      relation_distance_pk.insert(rd);
   }
}

void cpmds::get_inward_ordering(list<adjacency_face> &ordering)
{
   // First, define a vertex ordering
   set<unsigned int> processed_vertices;

   set<pair<unsigned int, unsigned int> > adjacency_edges;
   list<pair<unsigned int, unsigned int> > border_adjacency_edges;

   adjacency_face_third_itm_idx& adf_third_itm_idx = adjacency_face_table.get<adjacency_face_third_itm_idx_type>();

   pair<adjacency_face_third_itm_idx::const_iterator, adjacency_face_third_itm_idx::const_iterator> f_edge_range = adf_third_itm_idx.equal_range(boost::make_tuple(0));

   // Start with all the edge items 
   for (; f_edge_range.first != f_edge_range.second; ++f_edge_range.first)
   {
      pair<unsigned int, unsigned int> edge = make_pair(f_edge_range.first->first_itm_id, f_edge_range.first->second_itm_id);
      adjacency_edges.insert(edge);
      adjacency_edges.insert(make_pair(edge.second, edge.first));
      border_adjacency_edges.push_back(edge);
      
      if (processed_vertices.find(edge.second) == processed_vertices.end())
      {
         ordering.push_back(*f_edge_range.first);
         processed_vertices.insert(edge.second);
      } 
   }

   // For each edge, find adjacent edges that have not yet been iterated
   while (!border_adjacency_edges.empty())
   {
      list<pair<unsigned int, unsigned int> > previous_border_adjacency_edges = border_adjacency_edges;
      border_adjacency_edges.clear();

      for (list<pair<unsigned int, unsigned int> >::const_iterator edge_itr = previous_border_adjacency_edges.begin(), end = previous_border_adjacency_edges.end(); edge_itr != end; ++edge_itr)
      {
         adjacency_face f_next;
         find_adjacency_face(edge_itr->second, edge_itr->first, f_next);   

         pair<unsigned int, unsigned int> next_left_edge = make_pair(f_next.third_itm_id, f_next.first_itm_id);
         pair<unsigned int, unsigned int> next_right_edge = make_pair(f_next.second_itm_id, f_next.third_itm_id);

         if (adjacency_edges.find(next_left_edge) == adjacency_edges.end())
         {
            adjacency_edges.insert(next_left_edge);
            adjacency_edges.insert(make_pair(next_left_edge.second, next_left_edge.first));
            border_adjacency_edges.push_back(next_left_edge);
         }

         if (adjacency_edges.find(next_right_edge) == adjacency_edges.end())
         {
            adjacency_edges.insert(next_right_edge);
            adjacency_edges.insert(make_pair(next_right_edge.second, next_right_edge.first));
            border_adjacency_edges.push_back(next_right_edge);
         }

         if (processed_vertices.find(f_next.third_itm_id) == processed_vertices.end())
         {
            ordering.push_back(f_next);
            processed_vertices.insert(f_next.third_itm_id);
         } 
      }
   }
}

void cpmds::flatten_graph(graph_step_type step)
{
   revert_ideal_radii();

   bool all_items_processed = false;
   switch (step)
   {
      case GRAPH_STEP:
         while (!all_items_processed) all_items_processed = flatten_next_item(step);
      break;
      case ITEM_STEP:
      case MOVE_STEP:
         flatten_next_item(step);
      break;
      default:
      break;
   }

   output_statistics();

   layout_adjacency_graph();
}

void cpmds::output_statistics()
{
   //cout << "Stress (DIST, AVG): " << get_weighted_measure(DISTANCE, AVG).second << endl;
   cout << "Stress (RANK, AVG): " << get_weighted_measure(RANK, AVG).second << endl;
   //cout << "Stress (RANK, MIN): " << get_weighted_measure(RANK, MIN).second << endl;
   //cout << "Stress (RANK, MAX): " << get_weighted_measure(RANK, MAX).second << endl;
   cout << "Planar ? " << (is_planar() ? "TRUE" : "FALSE") << endl;
   cout << "Correct curvature ? " << (validate_curvature() ? "TRUE" : "FALSE") << endl;
   cout << "Correct stress ? " << (validate_distance() ? "TRUE" : "FALSE") << endl;

   // Output a histogram rank graph
   ofstream fout("/home/tgee/projects/c++/algo/data/cpmds-rankhist.dat"); 
   vector<unsigned int> rank_frequency(item_table.size(), 0);
  
   item_pk& itm_pk = item_table.get<item_pk_type>();
   for (item_pk::const_iterator itm_itr = itm_pk.begin(), itm_end = itm_pk.end(); itm_itr != itm_end; ++itm_itr)
   {
      list<unsigned int> perimeter; 
      get_perimeter(itm_itr->id, perimeter);
      for (list<unsigned int>::const_iterator per_itr = perimeter.begin(), per_end = perimeter.end(); per_itr != per_end; ++per_itr)
      {
         if (*per_itr != 0)
         {
            pair<double, double> rank = get_weighted_measure(itm_itr->id, *per_itr, RANK);
            ++rank_frequency[rank.second];
         }
      }
   }
   
   for (unsigned int i = 0; i < rank_frequency.size(); ++i)
   {
      fout << i << ", " << rank_frequency[i] << endl;
   }
}

bool cpmds::flatten_next_item(graph_step_type step)
{
   bool all_items_processed = false;

   if (_flattened_items.size() == item_table.size())
   {
      reset_flatten();
   }

   if (_flatten_item == 0)
   {
      if (_flattened_boundary.empty())
      {
         list<adjacency_face> face_order;
         get_inward_ordering(face_order);

         if (!face_order.empty())
         {
            unsigned int innermost_vertex = face_order.back().third_itm_id;
            if (innermost_vertex != 0)
            {
               _flattened_boundary.push_back(innermost_vertex);
               _flattened_boundary_items.insert(make_pair(innermost_vertex, _flattened_boundary.begin()));
            }
         }
      }

      // Get the connectivity of every vertex on the flattened boundary
      _boundary_connectivity.clear();
      for (list<unsigned int>::const_iterator itr = _flattened_boundary.begin(), end = _flattened_boundary.end(); itr != end; ++itr)
      {
         if (*itr != 0 && _flattened_items.find(*itr) == _flattened_items.end())
         {
            // Eliminate any items with a border self join connectivity greater than two
            // TODO Maintain a sorted list when calculating the self joins?
            unsigned int boundary_self_connectivity = 0;
            boundary_self_join_pk& bsj_pk = boundary_self_join_table.get<boundary_self_join_pk_type>();
            boundary_self_join_pk::const_iterator self_join_pk_itr, self_join_pk_end;
            tie(self_join_pk_itr, self_join_pk_end) = bsj_pk.equal_range(*itr); 
            boundary_self_connectivity += distance(self_join_pk_itr, self_join_pk_end);
            //for (;self_join_pk_itr != self_join_pk_end; ++self_join_pk_itr)
            //{
            //   if (!is_flattened_item(self_join_pk_itr->second_itm_id)) ++boundary_self_connectivity;
            //}

            boundary_self_join_reverse_idx& bsj_reverse_idx = boundary_self_join_table.get<boundary_self_join_reverse_idx_type>();
            boundary_self_join_reverse_idx::const_iterator self_join_r_idx_itr, self_join_r_idx_end;
            tie(self_join_r_idx_itr, self_join_r_idx_end) = bsj_reverse_idx.equal_range(*itr); 
            boundary_self_connectivity += distance(self_join_r_idx_itr, self_join_r_idx_end);
            //for (;self_join_r_idx_itr != self_join_r_idx_end; ++self_join_r_idx_itr)
            //{
            //   if (!is_flattened_item(self_join_r_idx_itr->second_itm_id)) ++boundary_self_connectivity;
            //}

            if (boundary_self_connectivity < 3)
            {
               // Prioritise the next item based on the constrained connectivity due to processed nodes
               list<unsigned int> boundary_perimeter;
               get_perimeter(*itr, boundary_perimeter); 
               unsigned int free_connectivity = 0;
               unsigned int constrained_connectivity = 0;

               // Subtract from this any perimeter items that have been processed and are fixed
               for (list<unsigned int>::const_iterator perimeter_itr = boundary_perimeter.begin(), perimeter_end = boundary_perimeter.end(); perimeter_itr != perimeter_end; ++perimeter_itr)
               {
                  if (*perimeter_itr != 0)
                  {
                     if (_flattened_items.find(*perimeter_itr) != _flattened_items.end())
                     {
                        constrained_connectivity++;
                     }
                     else
                     {
                        free_connectivity++;
                     }
                  }
               }

               _boundary_connectivity.insert(make_pair(*itr, make_pair(free_connectivity, constrained_connectivity)));
            }
         }
      }

      // Process the next item according to the free connectivity
      // TODO Avoid traversing the entire boundary
      unsigned int highest_constrained_connectivity = 0;
      unsigned int lowest_free_connectivity = numeric_limits<unsigned int>::max();
      for (map<unsigned int, pair<unsigned int, unsigned int> >::const_iterator itr = _boundary_connectivity.begin(), end = _boundary_connectivity.end(); itr != end; ++itr)
      {
         // Process the lowest connectivity first
         // TODO: Process highly constrained connectivity early as well
         unsigned int free_connectivity = itr->second.first;
         unsigned int constrained_connectivity = itr->second.second;
         if (constrained_connectivity > highest_constrained_connectivity)
         {
            highest_constrained_connectivity = constrained_connectivity;
            lowest_free_connectivity = free_connectivity;
            _flatten_item = itr->first;

            //cout << _flatten_item << " has highest constraints - " << highest_constrained_connectivity << 
            //                         " and free connections - " << lowest_free_connectivity << endl;
         }
         else if (constrained_connectivity == highest_constrained_connectivity &&
                  free_connectivity < lowest_free_connectivity)
         {
            lowest_free_connectivity = free_connectivity;
            _flatten_item = itr->first;

            //cout << _flatten_item << " has highest constraints - " << highest_constrained_connectivity << 
            //                         " and free connections - " << lowest_free_connectivity << endl;
         }
      }

      // Reset the curvature delta of the flattened_item
      item_pk& itm_pk = item_table.get<item_pk_type>();
      item_pk::iterator itm_itr = itm_pk.find(_flatten_item);
      item i = *itm_itr;

      i.curvature_delta = 0;
      itm_pk.replace(itm_itr, i);
   }

   if (_flatten_item != 0)
   {
      bool item_flattened = false;
      switch (step)
      {
         case GRAPH_STEP:
         case ITEM_STEP:
            while (!item_flattened) item_flattened = flatten_item(_flatten_item);
         break;
         case MOVE_STEP:
            item_flattened = flatten_item(_flatten_item);
         break;
         default:
         break;
      }

      if (item_flattened)
      {
         _flattened_items.insert(_flatten_item);
         update_flattened_boundary(_flatten_item);
         item_pk& itm_pk = item_table.get<item_pk_type>();
         item_pk::const_iterator itm_itr = itm_pk.find(_flatten_item);
         _flattened_curvature_sum += itm_itr->curvature;

         // Distribute the curvature surplus into the free surrounding items
         list<unsigned int> perimeter;
         get_perimeter(itm_itr->id, perimeter);
         unsigned int free_connectivity = 0;
         for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
         {
            if (*itr != 0 && _flattened_items.find(*itr) == _flattened_items.end()) 
            {
               ++free_connectivity;
            }
         }

         if (free_connectivity > 0)
         {
            double curvature_adjustment = itm_itr->curvature / free_connectivity;
            for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
            {
               if (*itr != 0 && _flattened_items.find(*itr) == _flattened_items.end())
               {
                  map<unsigned int, double>::iterator curvature_adjustment_itr = _curvature_adjustments.find(*itr);
                  if (curvature_adjustment_itr == _curvature_adjustments.end())
            
                  {
                     _curvature_adjustments.insert(make_pair(*itr, curvature_adjustment));
                  }
                  else
                  {
                     curvature_adjustment_itr->second += curvature_adjustment; 
                  }
               }
            }
         }

         cout << "Flattened Curvature Sum - " << _flattened_curvature_sum * 2 * M_PI  << endl;

         double curvature_adjustment_sum = 0;
         for (map<unsigned int, double>::const_iterator itr = _curvature_adjustments.begin(), end = _curvature_adjustments.end(); itr != end; ++itr)
         {
            curvature_adjustment_sum += itr->second;
         }

         cout << "Adjustment Curvature Sum - " << curvature_adjustment_sum * 2 * M_PI << endl;

         _flatten_item = 0;
      }
   }

   if (_flattened_items.size() == item_table.size())
   {
      all_items_processed = true;
      ++_graph_iteration;
   }

   return all_items_processed;
}

bool cpmds::flatten_item(unsigned int id)
{
   adjacency_face_pk& adf_pk = adjacency_face_table.get<adjacency_face_pk_type>();
   item_pk& itm_pk = item_table.get<item_pk_type>();
   item_pk::const_iterator itm_itr = itm_pk.find(id);

   cout << "flatten_item: # " << _graph_iteration << "-" << _flattened_items.size() << " Processing item " << id << " (" << itm_itr->curvature * 180 / M_PI << ")" << endl;

   bool is_acceptable_curvature = true;

   if (!is_interior_item(id))
   {
      // First try to convert to interior
      list<unsigned int> perimeter;
      get_perimeter(id, perimeter);
      const adjacency_face_third_itm_idx& adf_third_itm_idx = adjacency_face_table.get<adjacency_face_third_itm_idx_type>();

      adjacency_face_third_itm_idx::const_iterator f_outer_itr, f_outer_end_itr;
      tie(f_outer_itr, f_outer_end_itr) = adf_third_itm_idx.equal_range(boost::make_tuple(0)); 

      // There must be at least three exterior items to ensure planarity
      // Also ensure we didn't previously convert his item to an exterior item
      if (distance(f_outer_itr, f_outer_end_itr) > 3 &&
          _converted_exterior_items.find(id) == _converted_exterior_items.end())
      {
         // Ensure the items to be connected are not already adjacent
         // TODO: Let the perimeter boundary items be flattened as long as they still will have valid curvature
         adjacency_face f_tmp;
         if (//_flattened_items.find(*(++perimeter.begin())) == _flattened_items.end() &&
             //_flattened_items.find(perimeter.back()) == _flattened_items.end() &&
             find_adjacency_face(*(++perimeter.begin()), perimeter.back(), f_tmp) == adf_pk.end())
         {
            // Convert to an interior item by joining the two edge items
            cout << "Converting " << id << " to an interior item " << endl;
            list<unsigned int> fan_trfm_perimeter;
            map<unsigned int, item_adjustment> adjustments;
            fan_trfm_perimeter.push_back(*(++perimeter.begin()));
            fan_trfm_perimeter.push_back(id);
            fan_trfm_perimeter.push_back(perimeter.back());
            fan_trfm_perimeter.push_back(0);

            fan_transform(fan_trfm_perimeter, adjustments, true);
            is_acceptable_curvature = false;
         }
         else if (itm_itr->curvature < 0)
         {
            // Unable to convert to interior and the curvature is still unacceptable
            // Disconnect from the least similar adjacent perimeter item
            double lowest_stress_delta = numeric_limits<double>::max();
            list<unsigned int> fan_trfm_perimeter;
            bool clockwise_perimeter = false;

            unsigned int constrained_connectivity = 0;
            list<unsigned int> perimeter;
            get_perimeter(id, perimeter); 
            for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
            { 
               if (*itr != 0 && is_flattened_item(*itr)) ++constrained_connectivity;
            }

            do
            {
               list<unsigned int> potential_fan_trfm_perimeter;
               potential_fan_trfm_perimeter.push_back(0);

               adjacency_face f_outer, f_inner;
               if (clockwise_perimeter)
               {
                  find_adjacency_face(id, 0, f_outer);
                  potential_fan_trfm_perimeter.push_back(f_outer.third_itm_id);

                  find_adjacency_face(id, f_outer.third_itm_id, f_inner);
                  potential_fan_trfm_perimeter.push_back(f_inner.third_itm_id);
                  potential_fan_trfm_perimeter.push_back(id);
               }
               else
               {
                  find_adjacency_face(0, id, f_outer);
                  potential_fan_trfm_perimeter.push_back(id);

                  find_adjacency_face(f_outer.third_itm_id, id, f_inner);
                  potential_fan_trfm_perimeter.push_back(f_inner.third_itm_id);
                  potential_fan_trfm_perimeter.push_back(f_outer.third_itm_id);
               }

               item_pk::const_iterator disconnected_itm_itr = itm_pk.find(f_outer.third_itm_id);

               // Make sure we never disconnect this item from the flattened area
               if (disconnected_itm_itr->connectivity > 2 && 
                   (!is_flattened_item(f_outer.third_itm_id) || constrained_connectivity > 2))
               {
                  map<unsigned int, item_adjustment> adjustments;
                  fan_transform(potential_fan_trfm_perimeter, adjustments, false);

                  double stress_delta = get_stress_delta(adjustments, CURVATURE_STRESS);

                  if (stress_delta <= lowest_stress_delta)
                  {
                     lowest_stress_delta = stress_delta;
                     fan_trfm_perimeter.clear();
                     copy(potential_fan_trfm_perimeter.begin(), potential_fan_trfm_perimeter.end(), back_inserter(fan_trfm_perimeter));
                  }
               }

               clockwise_perimeter = !clockwise_perimeter;
            } 
            while (clockwise_perimeter != false);

            if (!fan_trfm_perimeter.empty())
            {
               map<unsigned int, item_adjustment> adjustments;
               fan_transform(fan_trfm_perimeter, adjustments, true);
               is_acceptable_curvature = false;

               // Update the flattened boundary in the case where the new boundary is flattened
               clockwise_perimeter = (*++fan_trfm_perimeter.begin() != id);
               list<unsigned int> new_boundary;
               new_boundary.push_back(*++++fan_trfm_perimeter.begin());
               unsigned int old_boundary = (clockwise_perimeter ? *++fan_trfm_perimeter.begin() : *++++++fan_trfm_perimeter.begin());
               
               if (is_flattened_item(new_boundary.front()))
               {
                  list<unsigned int>::iterator insert_itr = _flattened_boundary_items[old_boundary];
                  if (!clockwise_perimeter)
                  {
                     insert_itr++;
                     if (insert_itr == _flattened_boundary.end()) insert_itr = _flattened_boundary.begin();
                  }

                  for (list<unsigned int>::const_iterator new_boundary_itr = new_boundary.begin(), new_boundary_end = new_boundary.end(); new_boundary_itr != new_boundary_end; ++new_boundary_itr)
                  {
                     insert_itr = _flattened_boundary.insert(insert_itr, *new_boundary_itr);
                     _flattened_boundary_items.insert(make_pair(*new_boundary_itr, insert_itr));
                     ++insert_itr;
                  }

                  find_joins(new_boundary.begin(), new_boundary.end(),
                             _flattened_boundary.begin(), _flattened_boundary.end());
               }

               cout << "Disconnecting perimeter item " << id << " from " << old_boundary << endl;
            }
         }
      }
   }
   else 
   {
      // An interior item
      double curvature_tolerance = M_PI / 6;
      double adjusted_curvature = itm_itr->curvature;
      map<unsigned int, double>::const_iterator curvature_adjustment_itr = _curvature_adjustments.find(id);
      if (curvature_adjustment_itr != _curvature_adjustments.end())
      {
         adjusted_curvature += curvature_adjustment_itr->second;
      }

      list<unsigned int> fan_trfm_perimeter;

      if (adjusted_curvature + curvature_tolerance < 0) 
      {
         fan_trfm_perimeter.clear();
         get_increase_curvature_transform(id, fan_trfm_perimeter);

         // Push the best candidate away
         if (!fan_trfm_perimeter.empty())
         {
            cout << "Increasing " << id << " curvature (was " << itm_itr->curvature * 180 / M_PI << ", by moving " << fan_trfm_perimeter.back() << " away from " << *++fan_trfm_perimeter.begin() << ")." << endl;
            map<unsigned int, item_adjustment> adjustments;
            fan_transform(fan_trfm_perimeter, adjustments, true);
            is_acceptable_curvature = false;

            //double stress_delta = get_stress_delta(adjustments);
            //cout << "Actual stress delta - " << stress_delta << endl;

            // Handle an update to the boundary if a flattened item was altered
            if (fan_trfm_perimeter.size() > 4)
            {
               unsigned int adjusted_item = 0;
               unsigned int old_boundary = 0;

               list<unsigned int> new_boundary;
               new_boundary.push_back(*fan_trfm_perimeter.begin());

               if (id == *++fan_trfm_perimeter.begin())
               {
                  old_boundary = fan_trfm_perimeter.back();
                  adjusted_item = *++++fan_trfm_perimeter.begin();
               }
               else
               {
                  old_boundary = *++fan_trfm_perimeter.begin();
                  adjusted_item = *++++++fan_trfm_perimeter.begin();
               }

               if (_flattened_items.find(adjusted_item) != _flattened_items.end())
               {
                  if (_flattened_boundary_items.find(old_boundary) != _flattened_boundary_items.end())
                  {
                     list<unsigned int>::iterator insert_itr = _flattened_boundary_items[old_boundary];

                     // Check if the old boundary item is still on the boundary
                     bool is_old_boundary_on_boundary = false;
                     list<unsigned int> old_boundary_perimeter;
                     get_perimeter(old_boundary, old_boundary_perimeter);
                     for (list<unsigned int>::const_iterator itr = old_boundary_perimeter.begin(), end = old_boundary_perimeter.end(); itr != end; ++itr)
                     {
                        if (_flattened_items.find(*itr) != _flattened_items.end())
                        {
                           is_old_boundary_on_boundary = true;
                           break;
                        }
                     }

                     if (is_old_boundary_on_boundary)
                     {
                        if (*++fan_trfm_perimeter.begin() == id)
                        {
                           ++insert_itr;
                           if (insert_itr == _flattened_boundary.end()) insert_itr = _flattened_boundary.begin();
                        }
                     }
                     else
                     {
                        _flattened_boundary_items.erase(*insert_itr);
                        insert_itr = _flattened_boundary.erase(insert_itr);

                        // Update the self joins on the boundary
                        remove_boundary_self_joins(old_boundary);
                     }

                     for (list<unsigned int>::const_iterator new_boundary_itr = new_boundary.begin(), new_boundary_end = new_boundary.end(); new_boundary_itr != new_boundary_end; ++new_boundary_itr)
                     {
                        insert_itr = _flattened_boundary.insert(insert_itr, *new_boundary_itr);
                        _flattened_boundary_items.insert(make_pair(*new_boundary_itr, insert_itr));
                        ++insert_itr;
                     }

                     find_joins(new_boundary.begin(), new_boundary.end(),
                                _flattened_boundary.begin(), _flattened_boundary.end());
                  }
               }
            }
         }
      }
      else if (adjusted_curvature > curvature_tolerance) 
      {
         fan_trfm_perimeter.clear();
         get_decrease_curvature_transform(id, fan_trfm_perimeter);

         if (!fan_trfm_perimeter.empty())
         {

            // Calculate if this vertex can fit into the perimeter of [id]
            list<unsigned int>::const_iterator itr = find(fan_trfm_perimeter.begin(), fan_trfm_perimeter.end(), id);
            list<unsigned int>::const_iterator next_itr = itr;
            ++next_itr; 
            if (next_itr == fan_trfm_perimeter.end()) next_itr = fan_trfm_perimeter.begin();

            list<unsigned int>::const_iterator prev_itr = itr;
            if (prev_itr == fan_trfm_perimeter.begin()) prev_itr = fan_trfm_perimeter.end();
            --prev_itr; 

            item_pk::const_iterator a_itr = itm_pk.find(*next_itr);
            item_pk::const_iterator b_itr = itm_pk.find(*prev_itr);

            if (fan_trfm_perimeter.front() != 0)
            {
               item_pk::const_iterator c_itr = itm_pk.find(fan_trfm_perimeter.front());
               double r = itm_itr->radius;
               double r_a = a_itr->radius;
               double r_b = b_itr->radius;
               double r_c = c_itr->radius;

               double required_angle = subtended_angle(r, r_a, r_c)
                                    + subtended_angle(r, r_b, r_c)
                                    - subtended_angle(r, r_a, r_b);
               
               if (required_angle < adjusted_curvature + curvature_tolerance)
               {
                  cout << "Decreasing " << id << " curvature (was " << itm_itr->curvature * 180 / M_PI << ", by moving " << id << " towards " << fan_trfm_perimeter.front() << ")." << endl;

                  //cout << "Boundary:" << endl;
                  //for (list<unsigned int>::const_iterator itr = fan_trfm_perimeter.begin(), end = fan_trfm_perimeter.end(); itr != end; ++itr)
                  //{
                  //   cout << " --> " << *itr << endl;
                  //}
                  //cout << endl;

                  map<unsigned int, item_adjustment> adjustments;
                  fan_transform(fan_trfm_perimeter, adjustments, true);

                  //double stress_delta = 0;
                  //for (map<unsigned int, item_adjustment>::const_iterator itr = adjustments.begin(), end = adjustments.end(); itr != end; ++itr)
                  //{
                  //   item_pk::const_iterator adj_itm_itr = itm_pk.find(itr->first);
                  //   stress_delta += (itr->second.distance_weight_sum_delta - itr->second.weight_sum_delta * adj_itm_itr->distance_weight_sum / adj_itm_itr->weight_sum)
                  //                   / (adj_itm_itr->weight_sum + itr->second.weight_sum_delta);
                  //}
                  //stress_delta /= item_table.size();
                  //cout << "Actual stress delta - " << stress_delta << endl;

                  is_acceptable_curvature = false;

               }
            }
            else
            {
               cout << "Converting " << id << " to exterior item, curvature was " << itm_itr->curvature * 180 / M_PI << endl;

               map<unsigned int, item_adjustment> adjustments;
               fan_transform(fan_trfm_perimeter, adjustments, true);

               _converted_exterior_items.insert(id);
            }
         }
      }
   }

   return is_acceptable_curvature;
}

void cpmds::update_flattened_boundary(unsigned int id)
{
   if (_flattened_boundary_items.find(id) == _flattened_boundary_items.end()) return;

   // Update the boundary
   list<unsigned int>::iterator itr = _flattened_boundary_items[id];
         
   list<unsigned int>::iterator prev_itr = itr;
   if (prev_itr == _flattened_boundary.begin()) prev_itr = _flattened_boundary.end();
   --prev_itr;

   list<unsigned int>::iterator next_itr = itr;
   ++next_itr;
   if (next_itr == _flattened_boundary.end()) next_itr = _flattened_boundary.begin();

   // Find the new perimeter
   list<unsigned int> perimeter;
   get_perimeter(id, perimeter, *prev_itr, *next_itr);

   // If an exterior vertex, do not delete this vertex from the boundary
   list<unsigned int>::iterator edge_itr = find(perimeter.begin(), perimeter.end(), 0);
   if (edge_itr != perimeter.end()) *edge_itr = id;

   // Update the boundary
   _flattened_boundary_items.erase(*itr);
   list<unsigned int>::iterator insert_itr = _flattened_boundary.erase(itr);
   list<unsigned int>::iterator new_boundary_start_itr = perimeter.begin(); 
   list<unsigned int>::iterator new_boundary_end_itr = perimeter.end(); 

   if (_flattened_boundary.size() != 0 && perimeter.size() >= 2)
   {
      ++new_boundary_start_itr;
      --new_boundary_end_itr;
   }

   for (list<unsigned int>::const_iterator new_boundary_itr = new_boundary_start_itr, new_boundary_end = new_boundary_end_itr; new_boundary_itr != new_boundary_end; ++new_boundary_itr)
   {
      insert_itr = _flattened_boundary.insert(insert_itr, *new_boundary_itr);
      _flattened_boundary_items.insert(make_pair(*new_boundary_itr, insert_itr));
      ++insert_itr;
   }

   find_joins(new_boundary_start_itr, new_boundary_end_itr, 
              _flattened_boundary.begin(), _flattened_boundary.end());

   remove_boundary_self_joins(id);

   if (!is_interior_item(id))
   {
      // Add back in the outward joins to the boundary
      list<unsigned int>::iterator edge_prev_itr = edge_itr;
      if (edge_prev_itr == perimeter.begin()) edge_prev_itr = perimeter.end();
      --edge_prev_itr;

      list<unsigned int>::iterator edge_next_itr = edge_itr;
      ++edge_next_itr;
      if (edge_next_itr == perimeter.end()) edge_next_itr = perimeter.begin();

      check_boundary_self_join(id, *edge_prev_itr);
      check_boundary_self_join(id, *edge_next_itr);
   }

/*
   cout << "Flattened boundary: " << endl;
   for (list<unsigned int>::const_iterator itr = _flattened_boundary.begin(), end = _flattened_boundary.end(); itr != end; ++itr)
   {
      cout << *itr << endl;
   }

   cout << "Self joins: " << endl;
   boundary_self_join_pk& bsj_pk = boundary_self_join_table.get<boundary_self_join_pk_type>();
   for (boundary_self_join_pk::const_iterator itr = bsj_pk.begin(), end = bsj_pk.end(); itr != end; ++itr)
   {
      cout << itr->first_itm_id << " - " << itr->second_itm_id << endl;
   }
*/
}

bool cpmds::is_boundary_item(unsigned int id) const
{
   return _flattened_boundary_items.find(id) != _flattened_boundary_items.end();
}

bool cpmds::is_boundary_edge(unsigned int first_itm_id, unsigned int second_itm_id) const
{
   bool is_edge = false;
   map<unsigned int, list<unsigned int>::iterator>::const_iterator boundary_itr = _flattened_boundary_items.find(first_itm_id);
   if (boundary_itr != _flattened_boundary_items.end())
   {
      list<unsigned int>::const_iterator first_itr = boundary_itr->second;

      list<unsigned int>::const_iterator prev_itr = first_itr;
      if (prev_itr == _flattened_boundary.begin()) prev_itr = _flattened_boundary.end();
      --prev_itr;

      if (second_itm_id == *prev_itr) is_edge = true;
      else
      {
         list<unsigned int>::const_iterator next_itr = first_itr;
         ++next_itr;
         if (next_itr == _flattened_boundary.end()) next_itr = _flattened_boundary.begin();

         if (second_itm_id == *next_itr) is_edge = true;
      }

   }

   return is_edge;
}

bool cpmds::is_boundary_self_join(unsigned int first_itm_id, unsigned int second_itm_id) const
{
   const boundary_self_join_pk& bsj_pk = boundary_self_join_table.get<boundary_self_join_pk_type>();
   pair<unsigned int, unsigned int> lower_higher_id = first_itm_id > second_itm_id
                                                     ? make_pair(second_itm_id, first_itm_id) 
                                                     : make_pair(first_itm_id, second_itm_id);
   return bsj_pk.find(boost::make_tuple(lower_higher_id.first, lower_higher_id.second)) != bsj_pk.end();
}

void cpmds::remove_boundary_self_join(unsigned int first_itm_id, unsigned int second_itm_id)
{
   boundary_self_join_pk& bsj_pk = boundary_self_join_table.get<boundary_self_join_pk_type>();
   pair<unsigned int, unsigned int> lower_higher_id = first_itm_id > second_itm_id
                                                     ? make_pair(second_itm_id, first_itm_id) 
                                                     : make_pair(first_itm_id, second_itm_id);
   boundary_self_join_pk::iterator itr = bsj_pk.find(boost::make_tuple(lower_higher_id.first, lower_higher_id.second));
   if (itr != bsj_pk.end())
   {
      bsj_pk.erase(itr);
   }
}

void cpmds::remove_boundary_self_joins(unsigned int id)
{
   // Update the boundary self connectivity graph
   boundary_self_join_pk& bsj_pk = boundary_self_join_table.get<boundary_self_join_pk_type>();
   boundary_self_join_pk::const_iterator self_join_pk_itr, self_join_pk_end;
   tie(self_join_pk_itr, self_join_pk_end) = bsj_pk.equal_range(id); 
   bsj_pk.erase(self_join_pk_itr, self_join_pk_end);

   boundary_self_join_reverse_idx& bsj_reverse_idx = boundary_self_join_table.get<boundary_self_join_reverse_idx_type>();
   boundary_self_join_reverse_idx::const_iterator self_join_r_idx_itr, self_join_r_idx_end;
   tie(self_join_r_idx_itr, self_join_r_idx_end) = bsj_reverse_idx.equal_range(id); 
   bsj_reverse_idx.erase(self_join_r_idx_itr, self_join_r_idx_end);
}

void cpmds::check_boundary_self_join(unsigned int first_itm_id, unsigned int second_itm_id)
{
   boundary_self_join_pk& bsj_pk = boundary_self_join_table.get<boundary_self_join_pk_type>();
   pair<unsigned int, unsigned int> lower_higher_id = first_itm_id > second_itm_id
                                                     ? make_pair(second_itm_id, first_itm_id) 
                                                     : make_pair(first_itm_id, second_itm_id);
   if (bsj_pk.find(boost::make_tuple(lower_higher_id.first, lower_higher_id.second)) == bsj_pk.end())
   {
      // TODO, use a tree structure for the boundary to speed up this check
      if (_flattened_boundary_items.find(first_itm_id) != _flattened_boundary_items.end() &&
          _flattened_boundary_items.find(second_itm_id) != _flattened_boundary_items.end())
      {
         boundary_self_join bsj;
         bsj.first_itm_id = lower_higher_id.first;
         bsj.second_itm_id = lower_higher_id.second;
         bsj_pk.insert(bsj);
      }
   }
}

// Add any self joins between item sets
// TODO: Perform this in a more efficient way that brute force
void cpmds::find_joins(list<unsigned int>::const_iterator left_begin_itr,
                       list<unsigned int>::const_iterator left_end_itr,
                       list<unsigned int>::const_iterator right_begin_itr,
                       list<unsigned int>::const_iterator right_end_itr)
{
   adjacency_face_pk& adf_pk = adjacency_face_table.get<adjacency_face_pk_type>();
   boundary_self_join_pk& bsj_pk = boundary_self_join_table.get<boundary_self_join_pk_type>();

   for (list<unsigned int>::const_iterator left_itr = left_begin_itr; left_itr != left_end_itr; ++left_itr)
   {
      for (list<unsigned int>::const_iterator right_itr = right_begin_itr; right_itr != right_end_itr; ++right_itr)
      {
         adjacency_face f_tmp;
         if (find_adjacency_face(*left_itr, *right_itr, f_tmp) != adf_pk.end())
         {
            pair<unsigned int, unsigned int> lower_higher_id = *left_itr > *right_itr ? make_pair(*right_itr, *left_itr) 
                                                                                      : make_pair(*left_itr, *right_itr);

            if (lower_higher_id.first != 0)
            {
               if (bsj_pk.find(boost::make_tuple(lower_higher_id.first, lower_higher_id.second)) == bsj_pk.end())
               {
                  boundary_self_join bsj;
                  bsj.first_itm_id = lower_higher_id.first;
                  bsj.second_itm_id = lower_higher_id.second;
                  bsj_pk.insert(bsj);
               }
            }
         }
      }
   }
}

// Traverse the adjacency graph and recalculate the item positions
void cpmds::layout_adjacency_graph()
{
   item_pk& itm_pk = item_table.get<item_pk_type>();

   const unsigned int max_iterations = 1000;
   const double zero_curvature_angle_sum = 2 * M_PI;
   const double convergence_tolerance = 0.05;
   const double curvature_error_tolerance = 1e-4;
   const double radial_tolerance = 0.01;

   bool prevent_acceleration = true;
   double curvature_error_estimate = numeric_limits<double>::max();
   double error_reduction_factor = -1;

   // Start by solving for the required radii so the curvature is zero everywhere
   unsigned int iteration = 0;
   for (; curvature_error_estimate > curvature_error_tolerance && iteration < max_iterations; ++iteration)
   {
      double previous_error_estimate = curvature_error_estimate;
      double previous_error_reduction_factor = error_reduction_factor;
      bool   acceleration_possible = !prevent_acceleration;
      map<unsigned int, double> radius_delta;

      for (set<unsigned int>::const_iterator itr = _flattened_items.begin(), end = _flattened_items.end(); itr != end; ++itr)
      {
         if (is_interior_item(*itr))
         {
            item_pk::const_iterator itm_itr = itm_pk.find(*itr);
            double angle_sum = zero_curvature_angle_sum - itm_itr->curvature;
            double beta = sin(angle_sum / (2 * itm_itr->connectivity));

            // Calculate the Uniform Neighbour Model distance for this angle sum
            double unm_distance = itm_itr->radius * beta / (1 - beta);

            // Calculate analytically the radius required to have a zero curvature
            double gamma = sin(zero_curvature_angle_sum / (2 * itm_itr->connectivity));
            double r = unm_distance * (1 - gamma) / gamma;
            radius_delta[*itr] = r - itm_itr->radius;

            on_radius_change(*itr, r);
         }
      }

      // Calculate the error estimate as a mean square root of curvature
      curvature_error_estimate = 0;
      for (set<unsigned int>::const_iterator itr = _flattened_items.begin(), end = _flattened_items.end(); itr != end; ++itr)
      {
         if (is_interior_item(*itr))
         {
            item_pk::const_iterator itm_itr = itm_pk.find(*itr);
            curvature_error_estimate += itm_itr->curvature * itm_itr->curvature;                
         }
      }
      curvature_error_estimate = sqrt(curvature_error_estimate);
      error_reduction_factor = curvature_error_estimate / previous_error_estimate;

      prevent_acceleration = false;
      if (acceleration_possible && error_reduction_factor < 1)
      {
         curvature_error_estimate = error_reduction_factor * curvature_error_estimate;

         if (abs(error_reduction_factor - previous_error_reduction_factor) < convergence_tolerance)
         {
            // Perform super acceleration
            error_reduction_factor = error_reduction_factor / (1 - error_reduction_factor);
         }

         // Determine the largest error reduction factor that keeps all radii positive
         double largest_error_reduction_factor = 0;
         for (set<unsigned int>::const_iterator itr = _flattened_items.begin(), end = _flattened_items.end(); itr != end; ++itr)
         {
            if (is_interior_item(*itr))
            {
               item_pk::const_iterator itm_itr = itm_pk.find(*itr);
               double min_error_reduction_factor = -itm_itr->radius / radius_delta[*itr];

               if (min_error_reduction_factor > largest_error_reduction_factor)
               {
                  largest_error_reduction_factor = min_error_reduction_factor;
               }
            }
         }

         error_reduction_factor = min(error_reduction_factor, 0.5 * largest_error_reduction_factor);

         // Perform acceleration
         for (set<unsigned int>::const_iterator itr = _flattened_items.begin(), end = _flattened_items.end(); itr != end; ++itr)
         {
            if (is_interior_item(*itr))
            {
               item_pk::const_iterator itm_itr = itm_pk.find(*itr);
               double r = itm_itr->radius + error_reduction_factor * radius_delta[*itr];
               on_radius_change(*itr, r);
            }
         }

         prevent_acceleration = true;
      }
   }

   cout << "Solving curvature by perturbing radii, iterations - " << iteration << ", curvature_error_estimate - " << curvature_error_estimate << endl;

   // Layout entire perimeter of any flattened item
   map<unsigned int, pair<point, bool> > itm_locations; // id, location, located
   queue<unsigned int> pending_items;
   set<unsigned int> processed_parent_items;
   set<unsigned int> processed_child_items;
   map<unsigned int, unsigned int> child_parent_map;

   if (!_flattened_items.empty()) pending_items.push(*_flattened_items.begin());

   while (!pending_items.empty())
   {
      unsigned int id = pending_items.front();
      pending_items.pop();

      if (processed_parent_items.find(id) == processed_parent_items.end())
      {

         list<unsigned int> perimeter;
         get_perimeter(id, perimeter);

         item_pk::const_iterator located_first_itm_itr = itm_pk.find(id);
         double r_a = located_first_itm_itr->radius;

         // Find the parent that was located prior to this item
         map<unsigned int, unsigned int>::const_iterator child_parent_itr = child_parent_map.find(id);
         if (child_parent_itr == child_parent_map.end())
         {
            // The very first item
            itm_locations.insert(make_pair(id, make_pair(point(0, 0), true)));
            processed_child_items.insert(id);

            if (perimeter.front() == 0) perimeter.erase(perimeter.begin()); 
         
            if (!perimeter.empty())
            {
               unsigned int second_itm_id = perimeter.front();
               item_pk::const_iterator second_itm_itr = itm_pk.find(second_itm_id);

               double r_b = second_itm_itr->radius;
               itm_locations.insert(make_pair(second_itm_id, make_pair(point(r_a + r_b, 0), true)));
               processed_child_items.insert(second_itm_id);

               if (_flattened_items.find(second_itm_id) != _flattened_items.end())
               {
                  child_parent_map.insert(make_pair(second_itm_id, id));
                  pending_items.push(second_itm_id);
               }
            }
         }
         else
         {
            list<unsigned int>::iterator parent_itr = find(perimeter.begin(), perimeter.end(), child_parent_itr->second);
            
            // Put the parent at the beginning of the list
            rotate(perimeter.begin(), parent_itr, perimeter.end());
         }

         for(list<unsigned int>::const_iterator itr = ++perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
         {
            list<unsigned int>::const_iterator prior_itr = itr;
            --prior_itr;

            if (*itr != 0 && *prior_itr != 0 &&
                processed_child_items.find(*itr) == processed_child_items.end() &&
                (processed_child_items.find(*prior_itr) != processed_child_items.end() || 
                 processed_parent_items.find(*prior_itr) != processed_parent_items.end()))
            {
               // Locate the perimeter item 
               item_pk::const_iterator prior_itm_itr = itm_pk.find(*prior_itr);
               item_pk::const_iterator itm_itr = itm_pk.find(*itr);

               double r_b = prior_itm_itr->radius;
               double r_c = itm_itr->radius;

               point a = itm_locations[id].first;
               point b = itm_locations[*prior_itr].first;

               // Find the position of c assuming it is tangent to both a and b
               double ab_length = sqrt(pow(b.x - a.x, 2) + pow(b.y - a.y, 2));
               //if (ab_length > (r_a + 2 * r_c + r_b)) 
               if (ab_length - (r_a + r_b) > ((r_a + r_b) * radial_tolerance)) 
               {
                  // TODO place beside either A or B
                  cout << *itr << " cannot be placed. AB length = " << ab_length << ", r_a + r_b = " << r_a + r_b <<
                                  ", error % = " << (ab_length - (r_a + r_b)) / (r_a + r_b) <<  " tolerance = " << radial_tolerance << endl;
                  cout << "A = " << id << ", B = " << *prior_itr << endl;
                  while (!pending_items.empty()) pending_items.pop();
                  break;
               }
               else
               {
                  // cout << *itr << ", positional error % = " << (ab_length - (r_a + r_b)) / (r_a + r_b) <<  " tolerance = " << radial_tolerance << endl;
               }

               double cos_alpha = ((r_a + r_c) * (r_a + r_c) + ab_length * ab_length - (r_b + r_c) * (r_b + r_c))
                                 / (2. * ab_length * (r_a + r_c));
               double sin_alpha = sqrt(1 - pow(cos_alpha, 2));

               // To find the new point, we move along AB then along the normal
               point c(a.x + ((b.x - a.x) / ab_length) * (r_a + r_c) * cos_alpha 
                           - ((b.y - a.y) / ab_length) * (r_a + r_c) * sin_alpha,
                       a.y + ((b.y - a.y) / ab_length) * (r_a + r_c) * cos_alpha
                           + ((b.x - a.x) / ab_length) * (r_a + r_c) * sin_alpha);
                        
               // cout << "Laying out #" << processed_child_items.size() << " "  << *itr << " from " << id << " and " << *prior_itr << " at " << c << endl;
               itm_locations.insert(make_pair(*itr, make_pair(c, true)));

               if (_flattened_items.find(*itr) != _flattened_items.end())
               {
                  child_parent_map.insert(make_pair(*itr, id));
                  pending_items.push(*itr);
               }

               processed_child_items.insert(*itr);
            }
         }

         processed_parent_items.insert(id);
      }
   }

   // Reset previously calculated locations
   for (item_pk::iterator itr = itm_pk.begin(), end = itm_pk.end(); itr != end; ++itr)
   {
      if (itm_locations.find(itr->id) == itm_locations.end())
      {
         itm_locations.insert(make_pair(itr->id, make_pair(point(10e23, 10e23), false)));
      }
   }
 
   // Update the item table with the calculated locations
   for (map<unsigned int, pair<point, bool> >::const_iterator itr = itm_locations.begin(), end = itm_locations.end(); itr != end; ++itr)
   {
      item_pk::const_iterator itm_itr = itm_pk.find(itr->first);
      item i = *itm_itr;
      i.location = itr->second.first;
      i.is_located = itr->second.second;
      itm_pk.replace(itm_itr, i); 
   }
}

pair<double, double> cpmds::get_weighted_measure(unsigned int id, const adjacency_face& f, measure_type s_type, aggregate_type a_type) const
{
   list<unsigned int> vertices;
   vertices.push_back(f.first_itm_id);
   vertices.push_back(f.second_itm_id);
   if (f.third_itm_id != 0)
      vertices.push_back(f.third_itm_id);

   list<pair<double, double> > sub_measures;
   for (list<unsigned int>::const_iterator itr = vertices.begin(), end = vertices.end(); itr != end; ++itr)
   {
      sub_measures.push_back(get_weighted_measure(id, *itr, s_type));
   }

   pair<double, double> aggregate_measure = make_pair(0., 0.);
   if (a_type == AVG)
   {
      for (list<pair<double, double> >::const_iterator itr = sub_measures.begin(), end = sub_measures.end(); itr != end; ++itr)
      {
         aggregate_measure.first += itr->first;
         aggregate_measure.second += itr->first * itr->second;
      }

      aggregate_measure.second /= aggregate_measure.first;
   }
   else if (a_type == MAX)
   {
      for (list<pair<double, double> >::const_iterator itr = sub_measures.begin(), end = sub_measures.end(); itr != end; ++itr)
      {
         if (itr->second > aggregate_measure.second)
         {
            aggregate_measure.second = itr->second;
         }
      }
      aggregate_measure.first = 1.;
   }
   else if (a_type == MIN)
   {
      aggregate_measure.second = numeric_limits<double>::max(); 
      for (list<pair<double, double> >::const_iterator itr = sub_measures.begin(), end = sub_measures.end(); itr != end; ++itr)
      {
         if (itr->second < aggregate_measure.second)
         {
            aggregate_measure.second = itr->second;
         }
      }
      aggregate_measure.first = 1.;
   }
  
   return aggregate_measure;
}

adjacency_face cpmds::make_face(unsigned int first_itm_id, unsigned int second_itm_id, unsigned int third_itm_id)
{
   // By convention, faces have the lowest id last
   adjacency_face f;
   if (f.first_itm_id < second_itm_id)
   {
      if (f.first_itm_id < third_itm_id)
      {
         f.first_itm_id = second_itm_id;
         f.second_itm_id = third_itm_id;
         f.third_itm_id = first_itm_id;
      }
      else
      {
         f.first_itm_id = first_itm_id;
         f.second_itm_id = second_itm_id;
         f.third_itm_id = third_itm_id;
      }
   }
   else
   {
      if (second_itm_id < third_itm_id)
      {
         f.first_itm_id = third_itm_id;
         f.second_itm_id = first_itm_id;
         f.third_itm_id = second_itm_id;
      }
      else
      {
         f.first_itm_id = first_itm_id;
         f.second_itm_id = second_itm_id;
         f.third_itm_id = third_itm_id;
      }
   }

   return f;
}

adjacency_face_pk::const_iterator cpmds::find_adjacency_face(unsigned int first_itm_id, unsigned int second_itm_id, adjacency_face& f) const
{
   const adjacency_face_pk& adf_pk = adjacency_face_table.get<adjacency_face_pk_type>();
   const adjacency_face_second_itm_idx& adf_second_itm_idx = adjacency_face_table.get<adjacency_face_second_itm_idx_type>();
   const adjacency_face_third_itm_idx& adf_third_itm_idx = adjacency_face_table.get<adjacency_face_third_itm_idx_type>();

   adjacency_face_pk::const_iterator f_pk_itr = adf_pk.find(boost::make_tuple(first_itm_id,
                                                                                   second_itm_id));
   if (f_pk_itr != adf_pk.end())
   {
      f = *f_pk_itr;
   }
   else
   {
      adjacency_face_second_itm_idx::const_iterator f_second_itr = adf_second_itm_idx.find(boost::make_tuple(first_itm_id,
                                                                                                                  second_itm_id));
      if (f_second_itr != adf_second_itm_idx.end())
      {
         // Note that the cyclic item order is changed to maintain the correspondence with the request
         f.first_itm_id = first_itm_id;
         f.second_itm_id = second_itm_id;
         f.third_itm_id = f_second_itr->first_itm_id;
         f_pk_itr = adf_pk.find(boost::make_tuple(f_second_itr->first_itm_id, f_second_itr->second_itm_id));
      }
      else
      {
         adjacency_face_third_itm_idx::const_iterator f_third_itr = adf_third_itm_idx.find(boost::make_tuple(first_itm_id,
                                                                                                                  second_itm_id));
         if (f_third_itr != adf_third_itm_idx.end())
         {
            f.first_itm_id = first_itm_id;
            f.second_itm_id = second_itm_id;
            f.third_itm_id = f_third_itr->second_itm_id;
            f_pk_itr = adf_pk.find(boost::make_tuple(f_third_itr->first_itm_id, f_third_itr->second_itm_id));
         }
      }
   }

   return f_pk_itr;
}

// Find the angle subtended between AB and AC where A, B, C are the centres of three circles
double cpmds::subtended_angle(double r_a, double r_b, double r_c)
{
   return acos(1 - (2 * r_b * r_c) / ((r_a + r_c) * (r_a + r_b)));
}

// Find the angle subtended between the centre of circle A, the centre of circle B and it's tangent
double cpmds::subtended_half_angle(double r_a, double r_b)
{
   return asin(r_b / (r_a + r_b));
}

bool cpmds::is_interior_item(unsigned int id) const
{
   const adjacency_face_pk& adf_pk = adjacency_face_table.get<adjacency_face_pk_type>();
   adjacency_face edge_face;
   return (find_adjacency_face(id, 0, edge_face) == adf_pk.end() &&
           find_adjacency_face(0, id, edge_face) == adf_pk.end());
}

bool cpmds::is_flattened_item(unsigned int id) const
{
   return _flattened_items.find(id) != _flattened_items.end();
}

item cpmds::get_item_by_id(unsigned int id) const
{
   const item_pk& itm_pk = item_table.get<item_pk_type>();
   return *itm_pk.find(id);
}

void cpmds::get_graph_afl_drawing(map<unsigned int, point> &vertex_locations)
{
   const item_pk& itm_pk = item_table.get<item_pk_type>();
   if (itm_pk.empty()) return;

   // Start with an arbitary node, let's take the node with the smallest id
   unsigned int n0 = itm_pk.begin()->id;

   // Select [n1] to maximize h[0, 1]
   map<unsigned int, unsigned int> h0, h1, h2, h3, h4, h5;
   get_hop_counts(n0, h0);

   unsigned int n1 = n0;
   unsigned int h01 = 0;
   for (map<unsigned int, unsigned int>::const_iterator itr = h0.begin(), end = h0.end(); itr != end; ++itr)
   {
      if (itr->second > h01)
      {
         h01 = itr->second;
         n1 = itr->first;
      }
   }

   // Select [n2] to maximize h[1, 2]
   get_hop_counts(n1, h1);

   unsigned int n2 = n1;
   unsigned int h12 = 0;
   for (map<unsigned int, unsigned int>::const_iterator itr = h1.begin(), end = h1.end(); itr != end; ++itr)
   {
      if (itr->second > h12)
      {
         h12 = itr->second;
         n2 = itr->first;
      }
   }

   // Select n3 to minimuze h[1,3] - h[2,3], while maximizing h[1,3] + h[2, 3] in a tie
   get_hop_counts(n2, h2);

   unsigned int n3 = n1;
   unsigned int h13_23_min_diff = numeric_limits<unsigned int>::max();
   unsigned int h13_23_max_sum  = 0;
   for (map<unsigned int, unsigned int>::const_iterator h23_itr = h2.begin(), end = h2.end(); h23_itr != end; ++h23_itr)
   {
      map<unsigned int, unsigned int>::const_iterator h13_itr = h1.find(h23_itr->first);

      unsigned int h13_23_diff = h13_itr->second > h23_itr->second ?
                                 h13_itr->second - h23_itr->second :
                                 h23_itr->second - h13_itr->second;
      unsigned int h13_23_sum = h13_itr->second + h23_itr->second;

      if (h13_23_diff < h13_23_min_diff ||
          (h13_23_diff == h13_23_min_diff && h13_23_sum > h13_23_max_sum))
      {
         h13_23_min_diff = h13_23_diff;
         h13_23_max_sum = h13_23_sum;
         n3 = h23_itr->first;
      }
   }
   
   // Select n4 to minimuze h[1,4] - h[2,4], while maximizing h[3,4] in a tie
   get_hop_counts(n3, h3);

   unsigned int n4 = n1;
   unsigned int h14_24_min_diff = numeric_limits<unsigned int>::max();
   unsigned int h34_max = 0;
   for (map<unsigned int, unsigned int>::const_iterator h24_itr = h2.begin(), end = h2.end(); h24_itr != end; ++h24_itr)
   {
      map<unsigned int, unsigned int>::const_iterator h14_itr = h1.find(h24_itr->first);

      unsigned int h14_24_diff = h14_itr->second > h24_itr->second ?
                                 h14_itr->second - h24_itr->second :
                                 h24_itr->second - h14_itr->second;
      unsigned int h34 = h3[h24_itr->first];

      if (h14_24_diff < h14_24_min_diff ||
          (h14_24_diff == h14_24_min_diff && h34 > h34_max))
      {
         h14_24_min_diff = h14_24_diff;
         h34_max = h34;
         n4 = h24_itr->first;
      }
   }

   // Select n5 to minimize h[1,5] - h[2,5], while minimizing h[3,5] - h[4,5] in a tie
   get_hop_counts(n4, h4);

   unsigned int n5 = n1;
   unsigned int h15_25_min_diff = numeric_limits<unsigned int>::max();
   unsigned int h35_45_min_diff = numeric_limits<unsigned int>::max();
   for (map<unsigned int, unsigned int>::const_iterator h25_itr = h2.begin(), end = h2.end(); h25_itr != end; ++h25_itr)
   {
      map<unsigned int, unsigned int>::const_iterator h15_itr = h1.find(h25_itr->first);
      map<unsigned int, unsigned int>::const_iterator h35_itr = h3.find(h25_itr->first);
      map<unsigned int, unsigned int>::const_iterator h45_itr = h4.find(h25_itr->first);

      unsigned int h15_25_diff = h15_itr->second > h25_itr->second ?
                                 h15_itr->second - h25_itr->second :
                                 h25_itr->second - h15_itr->second;
      unsigned int h35_45_diff = h35_itr->second > h45_itr->second ?
                                 h35_itr->second - h45_itr->second :
                                 h45_itr->second - h35_itr->second;

      if (h15_25_diff < h15_25_min_diff ||
          (h15_25_diff == h15_25_min_diff && h35_45_diff < h35_45_min_diff))
      {
         h15_25_min_diff = h15_25_diff;
         h35_45_min_diff = h35_45_diff;
         n5 = h25_itr->first;
      }
   }

   // Get the hop_counts from n5
   get_hop_counts(n5, h5);
/*
   cout << "n0: " << n0 << endl;
   cout << "n1: " << n1 << endl;
   cout << "n2: " << n2 << endl;
   cout << "n3: " << n3 << endl;
   cout << "n4: " << n4 << endl;
   cout << "n5: " << n5 << endl;
*/

   // Finally layout the graph in polar co-ordinates using the hop counts from these reference nodes
   for (item_pk::const_iterator itr = itm_pk.begin(), end = itm_pk.end(); itr != end; ++itr)
   {
      float radius = h5[itr->id] * 100;
      /*
      float theta = h3[itr->id] != h4[itr->id] ?
                    atan((float(h1[itr->id]) - float(h2[itr->id])) / (float(h3[itr->id]) - float(h4[itr->id]))) :
                       h1[itr->id] > h2[itr->id] ?  
                       M_PI / 2 :
                       -M_PI / 2;
      */
      float theta = atan2((float(h1[itr->id]) - float(h2[itr->id])), (float(h3[itr->id]) - float(h4[itr->id])));
      vertex_locations.insert(make_pair(itr->id, point(radius * cos(theta), radius * sin(theta))));
   }
}

void cpmds::get_graph_straight_line_drawing(map<unsigned int, point> &vertex_locations)
{
   // Create a lookup between item ids and indices
   if (item_table.size() < 3) return;

   // Ensure the first three graph vertices are adjacent and exterior
   adjacency_face_third_itm_idx& adf_third_itm_idx = adjacency_face_table.get<adjacency_face_third_itm_idx_type>();
   adjacency_face_third_itm_idx::const_iterator f_outer_itr, f_outer_end_itr;
   for (tie(f_outer_itr, f_outer_end_itr) = adf_third_itm_idx.equal_range(boost::make_tuple(0)); f_outer_itr != f_outer_end_itr; ++f_outer_itr)
   {
      const adjacency_face_pk& adf_pk = adjacency_face_table.get<adjacency_face_pk_type>();
      map<unsigned int, unsigned int> id_to_index;
      vector<unsigned int> index_to_id;

      unsigned int first_itm_id = f_outer_itr->first_itm_id;
      unsigned int second_itm_id = f_outer_itr->second_itm_id;
      id_to_index.insert(make_pair(first_itm_id, 0));
      index_to_id.push_back(first_itm_id);
      id_to_index.insert(make_pair(second_itm_id, 1));
      index_to_id.push_back(second_itm_id);

      // Now enumerate all the other items
      const item_pk& itm_pk = item_table.get<item_pk_type>();
      unsigned int index = 2;
      for (item_pk::const_iterator itr = itm_pk.begin(), end = itm_pk.end(); itr != end; ++itr)
      {
         if (itr->id != first_itm_id && itr->id != second_itm_id)
         {
            id_to_index.insert(make_pair(itr->id, index));
            index_to_id.push_back(itr->id);
            index++;
         }
      }
      
      // Use boost to get a planar embedding
      typedef adjacency_list
         < setS,
           vecS,
           undirectedS,
           property<vertex_index_t, unsigned int>,
           property<edge_index_t, unsigned int>
         > graph;

      // Create the graph The functions
      // planar_canonical_ordering and chrobak_payne_straight_line_drawing both
      // require a maximal planar graph. 
      graph g(item_table.size());

      //cout << "Adjacencies:" << endl;
      for (adjacency_face_pk::const_iterator itr = adf_pk.begin(), end = adf_pk.end(); itr != end; ++itr)
      {
         //cout << itr->first_itm_id << ", " << itr->second_itm_id << ", " << itr->third_itm_id << endl;
         add_edge(id_to_index[itr->first_itm_id], id_to_index[itr->second_itm_id], g);
         if (itr->third_itm_id != 0)
         { 
            add_edge(id_to_index[itr->second_itm_id], id_to_index[itr->third_itm_id], g);
            add_edge(id_to_index[itr->third_itm_id], id_to_index[itr->first_itm_id], g);
         }
      }

      // Initialize the interior edge index
      property_map<graph, edge_index_t>::type e_index = get(edge_index, g);
      graph_traits<graph>::edges_size_type edge_count = 0;
      graph_traits<graph>::edge_iterator ei, ei_end;
      for(tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
         put(e_index, *ei, edge_count++);

      //Test for planarity; compute the planar embedding as a side-effect
      typedef list< graph_traits<graph>::edge_descriptor > vec_t;
      vector<vec_t> embedding(num_vertices(g));
      if (!boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                       boyer_myrvold_params::embedding = 
                                       &embedding[0]))
         cout << "Input graph is not planar" << endl;
  
      // Should already be bi-connected, but it can't hurt to be sure
      make_biconnected_planar(g, &embedding[0]);

      // Re-initialize the edge index, since we just added a few edges
      edge_count = 0;
      for(tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
         put(e_index, *ei, edge_count++);

      //Test for planarity again; compute the planar embedding as a side-effect
      if (!boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                        boyer_myrvold_params::embedding = 
                                       &embedding[0]))
         cout << "After calling make_biconnected, the graph is not planar" << endl;

      make_maximal_planar(g, &embedding[0]);

      // Re-initialize the edge index, since we just added a few edges
      edge_count = 0;
      for(tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
         put(e_index, *ei, edge_count++);

      // Test for planarity one final time; compute the planar embedding as a 
      // side-effect
      if (!boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                        boyer_myrvold_params::embedding = 
                                        &embedding[0]))
         cout << "After calling make_maximal_planar, the graph is not planar." << endl;
  
      // Find a canonical ordering
      list<graph_traits<graph>::vertex_descriptor> ordering;
      planar_canonical_ordering(g, &embedding[0], back_inserter(ordering));

      // Compute the straight line drawing
      if (ordering.size() != item_table.size())
      { 
         if (distance(f_outer_itr, f_outer_end_itr) == 1)
         {
            cout << "Failed to produce sensible ordering." << endl;
            cout << "Adjacency graph: " << endl;
            for (adjacency_face_pk::const_iterator itr = adf_pk.begin(), end = adf_pk.end(); itr != end; ++itr)
            {
               cout << itr->first_itm_id << ", " << itr->second_itm_id << ", " << itr->third_itm_id << endl;
            }
         }

         continue;
      }

      // Set up a property map to hold the mapping from vertices to points
      typedef iterator_property_map
         < vector<point>::iterator, 
           property_map<graph, vertex_index_t>::type 
         >
         straight_line_drawing_t;
      typedef property_map<graph, vertex_index_t>::type index_map_t;

      vector<point> straight_line_drawing_storage(num_vertices(g));

      straight_line_drawing_t straight_line_drawing
        (straight_line_drawing_storage.begin(), 
         get(vertex_index,g)
         );
      index_map_t index_map = get(vertex_index, g);

      chrobak_payne_straight_line_drawing(g, 
                                          embedding, 
                                          ordering.begin(),
                                          ordering.end(),
                                          straight_line_drawing
                                          );

      graph_traits<graph>::vertex_iterator vi, vi_end;
      for(tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi)
      {
         point p(get(straight_line_drawing, *vi));
         unsigned int id = index_to_id[index_map[*vi]];

         vertex_locations.insert(make_pair(id, p));
      }

      return;
   }
}

void cpmds::output_graph()
{
   // Output to vizgraph dot format
   ofstream fout("/home/tgee/projects/c++/algo/data/cpmds.dot"); 

   fout << "digraph \"cpmds\"" << endl;
   fout << "{" << endl;
   fout << "\tgraph [ model = \"subset\", overlap = false ];" << endl;
   fout << "\tnode [shape = circle ];" << endl;
   fout << endl;

   // Output the item positions from the planar embedding
   map<unsigned int, point> vertex_locations;
   get_graph_straight_line_drawing(vertex_locations);
   
   item_pk& itm_pk = item_table.get<item_pk_type>();
   for (item_pk::const_iterator itr = itm_pk.begin(), end = itm_pk.end(); itr != end; ++itr)
   {
      point p = vertex_locations[itr->id];
      fout << "\t_" << itr->id << " [pos=\"" << p.x * 100 << "," << p.y * 100 << "\"];" << endl;
   }

   fout << endl;
    
   // Output the edges according to the adjacency ;ist
   const adjacency_face_pk& adf_pk = adjacency_face_table.get<adjacency_face_pk_type>();
   for (adjacency_face_pk::const_iterator itr = adf_pk.begin(), end = adf_pk.end(); itr != end; ++itr)
   {
      fout << "\t_" << itr->first_itm_id;
      fout << "-> _" << itr->second_itm_id;

      if (itr->third_itm_id != 0)
      { 
         fout << "-> _" << itr->third_itm_id;
         fout << "-> _" << itr->first_itm_id;
      }

      fout << ";" << endl; 
   }
  
   fout << "}" << endl;
}

void cpmds::load_item(const item& itm)
{
   item_pk& itm_pk = item_table.get<item_pk_type>();
   itm_pk.insert(itm);
}

void cpmds::load_relation(const relation& rel)
{
   relation_table.insert(rel);
}

void cpmds::load_relation_distance(const relation_distance& rdt)
{
   relation_distance_table.insert(rdt);
}

void cpmds::load_adjacency_face(const adjacency_face& adf)
{
   adjacency_face_table.insert(adf);
}

void cpmds::reset_flatten()
{
   _flatten_item = 0;
   _flattened_boundary.clear();
   _flattened_boundary_items.clear();
   _flattened_items.clear();
   _converted_exterior_items.clear();
   boundary_self_join_table.clear();
   _flattened_curvature_sum = 0;
   _curvature_adjustments.clear();
}

pair<double, double> cpmds::get_weighted_measure(measure_type s_type, aggregate_type a_type) const
{
   list<pair<double, double> > sub_measures;
   const item_pk& itm_pk = item_table.get<item_pk_type>();
   for (item_pk::const_iterator itr = itm_pk.begin(), end = itm_pk.end(); itr != end; ++itr)
   {
      sub_measures.push_back(get_weighted_measure(itr->id, s_type, a_type));
   }

   pair<double, double> aggregate_measure = make_pair(0., 0.);
   if (a_type == AVG)
   {
      for (list<pair<double, double> >::const_iterator itr = sub_measures.begin(), end = sub_measures.end(); itr != end; ++itr)
      {
         aggregate_measure.first += itr->first;
         aggregate_measure.second += itr->first * itr->second;
      }

      aggregate_measure.second /= aggregate_measure.first;
   }
   else if (a_type == MAX)
   {
      for (list<pair<double, double> >::const_iterator itr = sub_measures.begin(), end = sub_measures.end(); itr != end; ++itr)
      {
         if (itr->second > aggregate_measure.second)
         {
            aggregate_measure.second = itr->second;
         }
      }
      aggregate_measure.first = 1.;
   }
   else if (a_type == MIN)
   {
      aggregate_measure.second = numeric_limits<double>::max(); 
      for (list<pair<double, double> >::const_iterator itr = sub_measures.begin(), end = sub_measures.end(); itr != end; ++itr)
      {
         if (itr->second < aggregate_measure.second)
         {
            aggregate_measure.second = itr->second;
         }
      }
      aggregate_measure.first = 1.;
   }
  
   return aggregate_measure;
}

pair<double, double> cpmds::get_weighted_measure(unsigned int first_itm_id, measure_type s_type, aggregate_type a_type, unsigned int second_itm_id, unsigned int max_hop_count) const
{
   list<pair<double, double> > sub_measures;

   if (second_itm_id == 0) second_itm_id = first_itm_id;

   if (second_itm_id == first_itm_id && max_hop_count == 1)
   {
      list<unsigned int> perimeter;
      get_perimeter(first_itm_id, perimeter);

      if (!perimeter.empty() &&
           perimeter.front() == 0) perimeter.erase(perimeter.begin());

      for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
      {
         sub_measures.push_back(get_weighted_measure(second_itm_id, *itr, s_type));
      }
   }
   else
   {
      map<unsigned int, unsigned int> hop_counts;
      get_hop_counts(first_itm_id, hop_counts, max_hop_count);

      for (map<unsigned int, unsigned int>::const_iterator itr = hop_counts.begin(), end = hop_counts.end(); itr != end; ++itr)
      {
         if (second_itm_id != itr->first)
         {
            sub_measures.push_back(get_weighted_measure(second_itm_id, itr->first, s_type));
         }
      }
   }

   pair<double, double> aggregate_measure = make_pair(0., 0.);
   if (a_type == AVG)
   {
      for (list<pair<double, double> >::const_iterator itr = sub_measures.begin(), end = sub_measures.end(); itr != end; ++itr)
      {
         aggregate_measure.first += itr->first;
         aggregate_measure.second += itr->first * itr->second;
      }

      aggregate_measure.second /= aggregate_measure.first;
   }
   else if (a_type == MAX)
   {
      for (list<pair<double, double> >::const_iterator itr = sub_measures.begin(), end = sub_measures.end(); itr != end; ++itr)
      {
         if (itr->second > aggregate_measure.second)
         {
            aggregate_measure.second = itr->second;
         }
      }
      aggregate_measure.first = 1.;
   }
   else if (a_type == MIN)
   {
      aggregate_measure.second = numeric_limits<double>::max(); 
      for (list<pair<double, double> >::const_iterator itr = sub_measures.begin(), end = sub_measures.end(); itr != end; ++itr)
      {
         if (itr->second < aggregate_measure.second)
         {
            aggregate_measure.second = itr->second;
         }
      }
      aggregate_measure.first = 1.;
   }
  
   return aggregate_measure;
}

pair<double, double> cpmds::get_weighted_measure(unsigned int first_itm_id, unsigned int second_itm_id, measure_type type) const
{
   pair<double, double> weighted_measure = make_pair(0., 0.);
   if (second_itm_id != 0)
   {

      if (type == RANK)
      {
         const relation_distance_dist_idx& rel_distance_dist_idx = relation_distance_table.get<relation_distance_dist_idx_type>();
         relation_distance_dist_idx::const_iterator rel_itr, rel_end_itr;
         tie(rel_itr, rel_end_itr) = rel_distance_dist_idx.equal_range(boost::make_tuple(first_itm_id));

         for (; rel_itr != rel_end_itr && rel_itr->second_itm_id != second_itm_id; ++rel_itr)
         {
            ++weighted_measure.second;
         }

         if (rel_itr == rel_end_itr) weighted_measure.first = 0;
         else weighted_measure.first = 1.; 
      }
      else if (type == DISTANCE || type == WEIGHTED_DISTANCE)
      {
         const relation_distance_pk& rel_distance_pk = relation_distance_table.get<relation_distance_pk_type>();
         relation_distance_pk::const_iterator rel_itr = rel_distance_pk.find(boost::make_tuple(first_itm_id, second_itm_id));

         if (rel_itr != rel_distance_pk.end())
         {
            weighted_measure.second = rel_itr->distance;
            if (type == WEIGHTED_DISTANCE) weighted_measure.first = rel_itr->weight;
            else weighted_measure.first = 1.;
         }
      }
   }

   return weighted_measure;
}

void cpmds::update_graph(graph_step_type step)
{
   // Get the current item list (note that updating items invalidates the index)
   if (_update_pending_items.empty())
   {
      item_pk& itm_pk = item_table.get<item_pk_type>();
      for (item_pk::const_iterator itr = itm_pk.begin(), end = itm_pk.end(); itr != end; ++itr)
      {
         _update_pending_items.insert(itr->id);
      }
   }

   bool all_items_processed;
   switch (step)
   {
      case GRAPH_STEP:
         while (!all_items_processed) all_items_processed = update_next_item();
         if (_focus > 0) --_focus;
         cout << "Focus : " << _focus << endl;
      break;
      case ITEM_STEP:
      case MOVE_STEP:
         update_next_item();
      break;
      default:
      break;
   }

   output_statistics();
}

bool cpmds::update_next_item()
{
   // Update items on the item list
   if (!_update_pending_items.empty())
   {
      unsigned int id = *_update_pending_items.begin();
      update_item(id);
      _update_pending_items.erase(id);
   }

   return _update_pending_items.empty();
}

void cpmds::update_item(unsigned int id)
{
   // This will potentially destroy the flatten process, so reset it
   reset_flatten();

   // See if moving the item will decrease the system stress
   adjacency_face_pk& adf_pk = adjacency_face_table.get<adjacency_face_pk_type>();
   adjacency_face_pk::iterator f_nearest_itr = find_k_nearest_face(id, _focus);

   bool move_item = false;
   if (f_nearest_itr != adf_pk.end() &&
       f_nearest_itr->first_itm_id != id &&
       f_nearest_itr->second_itm_id != id &&
       f_nearest_itr->third_itm_id != id)
   {
      map<unsigned int, item_adjustment> adjustments;
      double stress_delta = 0;
      if (_focus == 0)
      {
         remove_item(id, adjustments, false);
         insert_item(id, f_nearest_itr, adjustments, false);

         stress_delta = get_stress_delta(adjustments);
      }

      if (stress_delta <= 0)
      {
         if (!is_interior_item(id))
         {
            adjacency_face_third_itm_idx& adf_third_itm_idx = adjacency_face_table.get<adjacency_face_third_itm_idx_type>();
            adjacency_face_third_itm_idx::const_iterator f_outer_itr, f_outer_end_itr;
            tie(f_outer_itr, f_outer_end_itr) = adf_third_itm_idx.equal_range(boost::make_tuple(0)); 

            // There must be at least three exterior items to ensure planarity
            // Also ensure we didn't previously convert his item to an exterior item
            if (distance(f_outer_itr, f_outer_end_itr) > 3 &&
                !is_separating_item(id))
            {
               move_item = true;
            }
         }
         else move_item = true;
      }

      if (move_item)
      {
         cout << "#" << _update_pending_items.size() << " Stress delta for moving item " << id << " is " << stress_delta << endl; 
         // Move the item
         adjacency_face f_nearest = *f_nearest_itr;
         remove_item(id, adjustments, true);

         f_nearest_itr = adf_pk.find(boost::make_tuple(f_nearest.first_itm_id, f_nearest.second_itm_id));
         insert_item(id, f_nearest_itr, adjustments, true);
      }
   }
}

// A separating vertex is one that if remove will produce two graphs connected by a single vertex
bool cpmds::is_separating_item(unsigned int id) const
{
   bool is_separating = false;
   list<unsigned int> perimeter;
   get_perimeter(id, perimeter);
   
   unsigned int broken_perimeter_count = 0;
   unsigned int last_id = perimeter.back();
   bool is_last_interior = is_interior_item(perimeter.back());
   for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
   {
      // Ensure that removing this item does not produce two weakly connected graphs
      bool is_current_interior = false;
      unsigned int current_id = *itr;
      if (current_id == 0) current_id = id;
      if (*itr != 0) is_current_interior = is_interior_item(*itr);

      if (is_last_interior != is_current_interior) ++broken_perimeter_count; 
      else if (!is_current_interior)
      {
         adjacency_face tmp;
         find_adjacency_face(0, current_id, tmp);
         if (tmp.third_itm_id != last_id) broken_perimeter_count += 2;  
      }

      is_last_interior = is_current_interior;
      last_id = current_id;
   }

   // Check that the item does not connect two distinct boundaries
   if (broken_perimeter_count > 2) is_separating = true;

   // Also check that removing this item does not leave a perimeter vertex only connected to a single vertex
   if (get_connectivity(perimeter.back()) == 2) is_separating = true;
   if (get_connectivity(*(++perimeter.begin())) == 2) is_separating = true;

   return is_separating;
}

unsigned int cpmds::get_connectivity(unsigned int id) const
{
   unsigned int connectivity = 0;
   list<unsigned int> perimeter;
   get_perimeter(id, perimeter);

   connectivity = perimeter.size();
   if (perimeter.front() == 0) --connectivity;
   return connectivity; 
}

bool cpmds::is_planar() const
{
   // A very simple check, v - e + f = 1 ( where f does not include the outer face)
   // A necessary (but not sufficient) check for planarity
   const adjacency_face_third_itm_idx& adf_third_itm_idx = adjacency_face_table.get<adjacency_face_third_itm_idx_type>();
   adjacency_face_third_itm_idx::const_iterator f_inner_itr = adf_third_itm_idx.upper_bound(boost::make_tuple(0)); 
   adjacency_face_third_itm_idx::const_iterator f_end_itr = adf_third_itm_idx.end();
   
   unsigned int face_count = distance(f_inner_itr, f_end_itr);
   unsigned int edge_count = (face_count * 3 + adjacency_face_table.size() - face_count) / 2;
   unsigned int vertex_count = item_table.size();

   return (vertex_count - edge_count + face_count == 1); 
}

bool cpmds::validate_curvature() const
{
   bool valid_curvature = true;
   const item_pk& itm_pk = item_table.get<item_pk_type>();
   for (item_pk::const_iterator itr = itm_pk.begin(), end = itm_pk.end(); itr != end; ++itr)
   {
      valid_curvature &= validate_curvature(itr->id);
   }

   return valid_curvature;
}

bool cpmds::validate_distance() const
{
   bool valid_stress = true;
   const item_pk& itm_pk = item_table.get<item_pk_type>();
   for (item_pk::const_iterator itr = itm_pk.begin(), end = itm_pk.end(); itr != end; ++itr)
   {
      valid_stress &= validate_distance(itr->id);
   }

   return valid_stress;
}

bool cpmds::validate_curvature(unsigned int id) const
{
   bool valid_curvature = true; 
   double curvature = 2 * M_PI;

   list<unsigned int> perimeter;
   get_perimeter(id, perimeter);

   const item_pk& itm_pk = item_table.get<item_pk_type>();
   item_pk::const_iterator itm_itr = itm_pk.find(id);

   double r = itm_itr->radius;
   for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
   {
      list<unsigned int>::const_iterator next_itr = itr;
      ++next_itr;
      if (next_itr == perimeter.end()) next_itr = perimeter.begin();

      if (*itr != 0)
      {
         item_pk::const_iterator itm_a_itr = itm_pk.find(*itr);
         double r_a = itm_a_itr->radius;

         if (*next_itr != 0)
         {
            item_pk::const_iterator itm_b_itr = itm_pk.find(*next_itr);
            double r_b = itm_b_itr->radius;

            curvature -= subtended_angle(r, r_a, r_b);
         }
         else
         {
            curvature -= subtended_half_angle(r, r_a);
         }
      }
      else if (*next_itr != 0)
      {
         item_pk::const_iterator itm_b_itr = itm_pk.find(*next_itr);
         double r_b = itm_b_itr->radius;
         curvature -= subtended_half_angle(r, r_b);
      }
   }

   // Check for a difference bigger than a tenth of a degree
   if (abs(itm_itr->curvature - curvature) > M_PI / 1800)
   {
      valid_curvature = false;
      cout << "Item " << id << " has curvature " << itm_itr->curvature << " compared to theoretical " << curvature << endl;
   }

   return valid_curvature;
}

bool cpmds::validate_distance(unsigned int id) const
{
   bool valid_distance = true; 

   list<unsigned int> perimeter;
   get_perimeter(id, perimeter);

   const item_pk& itm_pk = item_table.get<item_pk_type>();
   item_pk::const_iterator itm_itr = itm_pk.find(id);

   double connectivity = 0;
   double weight_sum = 0;
   double distance_weight_sum = 0;
   for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
   {
      if (*itr != 0)
      {
         ++connectivity;
         pair<double, double> distance_weight = get_weighted_measure(id, *itr, WEIGHTED_DISTANCE);
         weight_sum += distance_weight.first;
         distance_weight_sum += distance_weight.first * distance_weight.second;
      }
   }

   // Check for a stress difference bigger than a tenth of a degree
   if (itm_itr->connectivity != connectivity)
   {
      valid_distance = false;
      cout << "Item " << id << " has connectivity " << itm_itr->connectivity << " compared to theoretical " << connectivity << endl;
   }
   if (abs(itm_itr->weight_sum - weight_sum) > 10e-5)
   {
      valid_distance = false;
      cout << "Item " << id << " has weight sum " << itm_itr->weight_sum << " compared to theoretical " << weight_sum << endl;
   }
   if (abs(itm_itr->distance_weight_sum - distance_weight_sum) > 10e-5)
   {
      valid_distance = false;
      cout << "Item " << id << " has distance weight sum " << itm_itr->distance_weight_sum << " compared to theoretical " << distance_weight_sum << endl;
   }

   return valid_distance;
}

bool cpmds::is_beneficial_adjustment(unsigned int id, double curvature_adjustment) const
{
   bool is_beneficial = false;
   const item_pk& itm_pk = item_table.get<item_pk_type>();
   item_pk::const_iterator itm_itr = itm_pk.find(id);
   double adjusted_curvature = itm_itr->curvature;

   map<unsigned int, double>::const_iterator curvature_adjustment_itr = _curvature_adjustments.find(id);
   if (curvature_adjustment_itr != _curvature_adjustments.end())
   {
      adjusted_curvature += curvature_adjustment_itr->second;
   }

   if ((adjusted_curvature > 0 && curvature_adjustment < 0) ||
       (adjusted_curvature < 0 && curvature_adjustment > 0)) is_beneficial = true;

   if (is_beneficial) cout << "Adjusting flattened item " << id << " is OK !" << endl;

   return is_beneficial;
}

void cpmds::revert_ideal_radii()
{
   // Change all the item radii to the theoretical ideal
   const item_pk& itm_pk = item_table.get<item_pk_type>();
   for (item_pk::const_iterator itr = itm_pk.begin(), end = itm_pk.end(); itr != end; ++itr)
   {
      on_radius_change(itr->id, sqrt(itr->sample_count));
   }
}

// Calculate the hop counts from a vertex
void cpmds::get_hop_counts(unsigned int id, map<unsigned int, unsigned int> &hop_counts, unsigned int max_hop_count) const
{
   set<unsigned int> processed_vertices;
   
   queue<unsigned int> unprocessed_vertices;

   hop_counts[id] = 0;
   unprocessed_vertices.push(id);

   while (!unprocessed_vertices.empty())   
   {
      unsigned int current_id = unprocessed_vertices.front();
      unprocessed_vertices.pop();

      if (processed_vertices.find(current_id) == processed_vertices.end())
      {
         list<unsigned int> perimeter;
         get_perimeter(current_id, perimeter);  

         unsigned int current_hop_count = hop_counts[current_id];
        
         if (max_hop_count != 0 && current_hop_count < max_hop_count)
         {
            for (list<unsigned int>::const_iterator itr = perimeter.begin(), end = perimeter.end(); itr != end; ++itr)
            {
               if (*itr != 0 && hop_counts.find(*itr) == hop_counts.end())
               {
                  hop_counts[*itr] = current_hop_count + 1;
                  unprocessed_vertices.push(*itr);
               }
            }
         }

         processed_vertices.insert(current_id);
      }
   }
}

double cpmds::get_stress_delta(const map<unsigned int, item_adjustment> &adjustments, stress_type s_type) const
{
   double stress_delta = 0;
   const item_pk& itm_pk = item_table.get<item_pk_type>();
   for (map<unsigned int, item_adjustment>::const_iterator itr = adjustments.begin(), end = adjustments.end(); itr != end; ++itr)
   {
      item_pk::const_iterator adj_itm_itr = itm_pk.find(itr->first);
      switch (s_type)
      {
         case DISTANCE_STRESS:
            stress_delta += (itr->second.distance_weight_sum_delta - itr->second.weight_sum_delta * adj_itm_itr->distance_weight_sum / adj_itm_itr->weight_sum) / (adj_itm_itr->weight_sum + itr->second.weight_sum_delta);
         break;
         case CURVATURE_STRESS:
            stress_delta += pow(adj_itm_itr->curvature + itr->second.curvature_delta, 2.0) - pow(adj_itm_itr->curvature, 2);
         break;
         default:
         break;
      }
   }

   return stress_delta;
}


