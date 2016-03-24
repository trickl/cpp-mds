#include "pmds.h"
#include <stdexcept>
#include <set>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

using namespace std;
using namespace boost;
using namespace boost::numeric::ublas;

pmds::pmds(unsigned int v_max,
           unsigned int h_max,
           double spring_constant,
           unsigned int max_polar_search_depth)
   : _v_max(v_max), _h_max(h_max), _spring_constant(spring_constant),
     _max_polar_search_depth(max_polar_search_depth)
{
}

void pmds::process_relation(unsigned int first_item_id, unsigned int second_item_id, float first_item_value, float second_item_value)
{
   if (first_item_id == second_item_id) return;
   if (first_item_id > second_item_id) swap(first_item_id, second_item_id);

   // Process one side at a time
   update_relation(first_item_id, second_item_id, first_item_value, second_item_value);
   update_item(first_item_id, first_item_value, second_item_value);
   adjust_item(first_item_id);

   update_relation(second_item_id, first_item_id, second_item_value, first_item_value);
   update_item(second_item_id, second_item_value, first_item_value);
   adjust_item(second_item_id);
}

void pmds::adjust_item(unsigned int id)
{
   item_pk_idx&     item_pk     = item_table.get<item_pk_idx_type>();
   item_pk_idx::iterator item_itr = item_pk.find(id);
   if (isnan(item_itr->location[0]) || item_table.size() == 1) return;

   relation_pk_idx& relation_pk = relation_table.get<relation_pk_idx_type>(); 
   set<unsigned int> item_comparison_list;

   // Find the nearest v_max locations from the relations
   double max_distance = numeric_limits<double>::infinity();
   relation_distance_dist_idx& rel_distance_dist_idx = relation_distance_table.get<relation_distance_dist_idx_type>();

   unsigned int v_count = 0;
   relation_distance_dist_idx::const_iterator rel_distance_itr = rel_distance_dist_idx.lower_bound(boost::make_tuple(id, 0));
   relation_distance_dist_idx::const_iterator rel_distance_upper_itr = rel_distance_dist_idx.upper_bound(boost::make_tuple(id, max_distance));
   for (; rel_distance_itr != rel_distance_upper_itr && v_count < _v_max; ++rel_distance_itr, ++v_count)
   {
      item_comparison_list.insert(rel_distance_itr->second_item_id);
   }

   // Sample from the furthest locations using the relations (to flatten out the mds).
   unsigned int s_count = 0;
   relation_distance_dist_idx::const_reverse_iterator rel_distance_ritr(rel_distance_dist_idx.upper_bound(boost::make_tuple(id, max_distance)));
   relation_distance_dist_idx::const_reverse_iterator rel_distance_lower_ritr(rel_distance_dist_idx.lower_bound(boost::make_tuple(id, 0)));
   for (; rel_distance_ritr != rel_distance_lower_ritr && s_count < _v_max; ++rel_distance_ritr, ++s_count)
   {
      item_comparison_list.insert(rel_distance_ritr->second_item_id);
   }

   // Find the nearest h_max locations from the spatial index
   if (!item_comparison_list.empty())
   {
      // Define minimum bounding rectangle for a local search 
      point mbr_bottom_left, mbr_top_right; 
      double mbr_half_size = rel_distance_lower_ritr->distance_mean * 2.;
      mbr_bottom_left[0] = item_itr->location[0] - mbr_half_size;
      mbr_bottom_left[1] = item_itr->location[1] - mbr_half_size;
      mbr_top_right[0] = item_itr->location[0] + mbr_half_size;
      mbr_top_right[1] = item_itr->location[1] + mbr_half_size;
      unsigned int h_count = 0;
      for(item_spatial_idx_type::const_iterator itr = item_spatial_idx.begin(mbr_bottom_left, mbr_top_right), end = item_spatial_idx.end(); itr != end && h_count < _h_max; ++itr, ++h_count)
      {
         if (itr->second != id) item_comparison_list.insert(itr->second);
      }
   }

   // Calculate the force on each particle due to all comparison items
   float net_force[2];
   fill(net_force, net_force + 2, 0);

   // Comparisons will be made between the locations in V and S
   for (set<unsigned int>::const_iterator itr = item_comparison_list.begin(), end = item_comparison_list.end(); itr != end; ++itr)
   {
      unsigned int second_item_id = *itr;

      // Force is proportional to difference between higher and lower dimensional distances
      item_pk_idx::const_iterator second_item_itr = item_pk.find(second_item_id);

      float loc_dist[2];
      loc_dist[0] = item_itr->location[0] - second_item_itr->location[0];
      loc_dist[1] = item_itr->location[1] - second_item_itr->location[1];

      double distance = sqrt((loc_dist[0] * loc_dist[0]) +
                             (loc_dist[1] * loc_dist[1]));

      // Default the direction if co-located
      if (distance == 0)
      {
         loc_dist[0] = 1.;
         loc_dist[1] = 0.;
         distance = 1.;
      }
      
      // Find the relation distance, initialise to the item average
      relation_pk_idx::const_iterator r_itr = relation_pk.find(boost::make_tuple(min(id, second_item_id), 
                                                                                        max(id, second_item_id)));
      double relation_distance;
      if (r_itr == relation_table.end())
      {
         relation_distance = second_item_itr->distance_mean;
      }
      else 
      {
         // Use a small sample weighted version of this distance
         const unsigned int small_sample_size = 1; // TODO: Figure this out from heuristics
         double linear_interp = double(small_sample_size) / (double(small_sample_size + r_itr->sample_count));
         relation_distance = linear_interp * second_item_itr->distance_mean + (1 - linear_interp) * r_itr->distance_mean; 
      }

      double force_size = _spring_constant * (relation_distance - distance);

      double scaled_force = force_size / (distance * double(item_comparison_list.size()));
      net_force[0] += loc_dist[0] * scaled_force;
      net_force[1] += loc_dist[1] * scaled_force;
   }

   item updated_item = *item_itr;
   updated_item.location[0] += net_force[0];
   updated_item.location[1] += net_force[1];
   
   if (!isnan(item_itr->location[0]))
   {
     item_spatial_idx.erase(item_itr->location);
     item_spatial_idx.insert(make_pair(updated_item.location, updated_item.id));
   }

   item_pk.replace(item_itr, updated_item);
}

// Utility for helping testing
matrix<double, column_major> pmds::item_locations_as_matrix() const
{
   const item_pk_idx& item_pk = item_table.get<item_pk_idx_type>();

   // Get the maximum id
   matrix<double, column_major> X(0, 2);
   if (item_table.size() > 0)
   {
      unsigned max_id = 0;
      for (unsigned int i = 0; i < item_table.size(); ++i)
      {
         max_id = max(max_id, item_table[i].id);
      }

      X.resize(max_id + 1, 2);
      for (unsigned int i = 0; i <= max_id; ++i)
      {
         item_pk_idx::const_iterator item_itr = item_pk.find(i);
         X(i, 0) = item_itr->location[0];
         X(i, 1) = item_itr->location[1];
      } 
   }

   return X;
}

void pmds::update_item(unsigned int id, float value, float compare_value)
{
   float distance = abs(value - compare_value);

   item i;

   item_pk_idx& item_pk = item_table.get<item_pk_idx_type>();
   item_pk_idx::const_iterator item_itr = item_pk.find(id);
   if (item_itr == item_pk.end())
   {
      i.id = id;
      i.location[0] = numeric_limits<float>::quiet_NaN();
      i.location[1] = numeric_limits<float>::quiet_NaN();
      i.sample_count = 0;
      i.distance_mean = 0;
      i.value_mean = 0;
      i.value_variance = 0;
   }
   else
   { 
      // Update the item
      i = *item_itr;
   }

   // Locate the item if possible
   if (isnan(i.location[0])) 
   {
      locate_item(id, i.location);

      if (!isnan(i.location[0]))
      {
         item_spatial_idx.insert(make_pair(i.location, id));

         // Create the relationship distance records
         relation_pk_idx& relation_pk = relation_table.get<relation_pk_idx_type>(); 
         relation_distance_pk_idx& relation_distance_pk = relation_distance_table.get<relation_distance_pk_idx_type>(); 
         pair<relation_pk_idx::const_iterator, relation_pk_idx::const_iterator>  range_itr = relation_pk.equal_range(boost::make_tuple(id));
         for (; range_itr.first != range_itr.second; ++range_itr.first)
         {
            relation_distance rd;
            rd.first_item_id = range_itr.first->second_item_id; 
            rd.second_item_id = id; 
            rd.distance_mean = range_itr.first->distance_mean;
            relation_distance_pk.insert(rd);
         }

         relation_reverse_idx& relation_rpk_idx = relation_table.get<relation_reverse_idx_type>(); 
         pair<relation_reverse_idx::const_iterator, relation_reverse_idx::const_iterator>  reverse_range_itr = relation_rpk_idx.equal_range(boost::make_tuple(id));
         for (; reverse_range_itr.first != reverse_range_itr.second; ++reverse_range_itr.first)
         {
            relation_distance rd;
            rd.first_item_id = reverse_range_itr.first->first_item_id; 
            rd.second_item_id = id; 
            rd.distance_mean = reverse_range_itr.first->distance_mean;
            relation_distance_pk.insert(rd);
         }
      }
   }
   
   i.sample_count++;
   i.distance_mean = i.distance_mean + (distance - i.distance_mean) / i.sample_count;

   float value_mean = i.value_mean;
   i.value_mean = i.value_mean + (value - i.value_mean) / i.sample_count;

   float qsum_delta = (value - value_mean) * (value - i.value_mean);
   i.value_variance = i.value_variance + (qsum_delta - i.value_variance) / i.sample_count;

   if (item_itr == item_pk.end())
   {
      item_pk.insert(i);
   }
   else
   {
      item_pk.replace(item_itr, i);
   }
}

void pmds::locate_item(unsigned int id, point& location)
{
   item_pk_idx& item_pk = item_table.get<item_pk_idx_type>();

   if (item_table.empty())
   {
      // The very first item
      location[0] = 0.;
      location[1] = 0.;
   }
   else
   {
      location[0] = numeric_limits<float>::quiet_NaN();
      location[1] = numeric_limits<float>::quiet_NaN();

      // Find the nearest neighbour and get all the neighbours
      relation_distance_dist_idx& rel_distance_dist_idx = relation_distance_table.get<relation_distance_dist_idx_type>();

      set<unsigned int> item_comparison_list;
      relation_distance_dist_idx::const_iterator rel_distance_itr = rel_distance_dist_idx.lower_bound(boost::make_tuple(id, 0));
      relation_distance_dist_idx::const_iterator rel_distance_upper_itr = rel_distance_dist_idx.upper_bound(boost::make_tuple(id, numeric_limits<double>::infinity()));
      relation_distance_dist_idx::const_iterator rel_nearest_itr = rel_distance_itr;

      if (rel_distance_itr != rel_distance_upper_itr)
      { 
         for (; rel_distance_itr != rel_distance_upper_itr; ++rel_distance_itr)
         {
            item_comparison_list.insert(rel_distance_itr->second_item_id);
         }

         item_pk_idx::const_iterator nearest_item_itr = item_pk.find(rel_nearest_itr->second_item_id);
         item_comparison_list.erase(nearest_item_itr->id);

         // Calculate initial position from relations using a polar search
         relation_pk_idx& relation_pk = relation_table.get<relation_pk_idx_type>(); 
         double angle = 0;
         double search_angle = M_PI / 2.;
         double radius = rel_nearest_itr->distance_mean;

         // Only locate the item once a minimum number of related items are available
         size_t min_relations_for_location = min(size_t(sqrt(relation_distance_table.size())), size_t(log10(sqrt(relation_distance_table.size())) + 1));
         if (item_comparison_list.size() + 1 >= min_relations_for_location) 
         {
            if (!item_comparison_list.empty())
            { 
               for (unsigned int search_iteration = 0; search_iteration < _max_polar_search_depth; ++search_iteration)
               {
                  // Find the stress of both configurations
                  double stress = 0;
                  for (unsigned int configuration = 0; configuration < 2; ++configuration)
                  {
                     for (set<unsigned int>::const_iterator itr = item_comparison_list.begin(), end = item_comparison_list.end(); itr != end; ++itr)
                     {
                        // Position according according to the nearest and lowest stress configuration
                        float direction_vector[2];
                        float query_angle = configuration ? angle - search_angle : angle + search_angle;
                        direction_vector[0] = radius * cos(query_angle);
                        direction_vector[1] = radius * sin(query_angle);
                  
                        location[0] = nearest_item_itr->location[0] + direction_vector[0];
                        location[1] = nearest_item_itr->location[1] + direction_vector[1];
   
                           unsigned int second_item_id = *itr;
                        item_pk_idx::const_iterator second_item_itr = item_pk.find(second_item_id);
   
                        float loc_dist[2];
                        loc_dist[0] = location[0] - second_item_itr->location[0];
                        loc_dist[1] = location[1] - second_item_itr->location[1];
   
                        double distance = sqrt((loc_dist[0] * loc_dist[0]) +
                                               (loc_dist[1] * loc_dist[1]));

                        // Find the relation distance, initialise to the item average
                        relation_pk_idx::const_iterator r_itr = relation_pk.find(boost::make_tuple(min(id, second_item_id), 
                                                                                           max(id, second_item_id)));
                           
                        double relation_distance;
                        if (r_itr == relation_table.end())
                        {
                           relation_distance = second_item_itr->distance_mean;
                        }
                        else 
                        {
                           // Use a small sample weighted version of this distance
                           const unsigned int small_sample_size = 1; // TODO: Figure this out from heuristics
                           double linear_interp = double(small_sample_size) / (double(small_sample_size + r_itr->sample_count));
                           relation_distance = linear_interp * second_item_itr->distance_mean + (1 - linear_interp) * r_itr->distance_mean; 
                        }

                        stress += pow(relation_distance - distance, 2.);
                     }

                     stress = -stress;
                  }

                  // Change the angle to the lower stress configuration 
                  angle = stress < 0 ? angle + search_angle : angle - search_angle;
                  search_angle /= 2;
               }
            }

            float direction_vector[2];
            direction_vector[0] = radius * cos(angle);
            direction_vector[1] = radius * sin(angle);

            location[0] = nearest_item_itr->location[0] + direction_vector[0];
            location[1] = nearest_item_itr->location[1] + direction_vector[1];
         }
      }
   }
}

void pmds::update_relation(unsigned int first_item_id, unsigned int second_item_id, float first_item_value, float second_item_value)
{
   float distance = abs(first_item_value - second_item_value);

   relation_pk_idx& relation_pk = relation_table.get<relation_pk_idx_type>(); 
   relation_distance_pk_idx& relation_distance_pk = relation_distance_table.get<relation_distance_pk_idx_type>(); 

   // Create the relation record
   if (first_item_id < second_item_id) 
   {
      relation_pk_idx::const_iterator r_itr = relation_pk.find(boost::make_tuple(first_item_id, second_item_id));
      if (r_itr == relation_table.end())
      {
         // Create the relation
         relation rel;
         rel.first_item_id = first_item_id;
         rel.second_item_id = second_item_id;
         rel.sample_count = 1;
         rel.distance_mean = distance;
         rel.distance_variance = 0;
         rel.value_covariance = 0;
         relation_pk.insert(rel);
      }
      else
      {
         // Update the relation
         relation rel = *r_itr;
         rel.sample_count++;
         float value_mean = rel.distance_mean;
         rel.distance_mean = rel.distance_mean + (distance - rel.distance_mean) / rel.sample_count;

         float qsum_delta = (distance - value_mean) * (distance - rel.distance_mean);
         rel.distance_variance = rel.distance_variance + (qsum_delta - rel.distance_variance) / rel.sample_count;
         relation_pk.replace(r_itr, rel);

         // Recalculate the covariance, note the items have the "old, n - 1" means
         item_pk_idx& item_pk = item_table.get<item_pk_idx_type>();
         item_pk_idx::const_iterator first_item_itr = item_pk.find(first_item_id);
         item_pk_idx::const_iterator second_item_itr = item_pk.find(second_item_id);
         double value_product_old_mean =  rel.value_covariance + first_item_itr->value_mean * second_item_itr->value_mean;
         double value_product_mean = value_product_old_mean + ((first_item_value * second_item_value) - value_product_old_mean) / rel.sample_count;
         double first_value_mean = first_item_itr->value_mean + (first_item_value - first_item_itr->value_mean) / rel.sample_count;
         double second_value_mean = second_item_itr->value_mean + (second_item_value - second_item_itr->value_mean) / rel.sample_count;
 
         rel.value_covariance = value_product_mean - first_value_mean * second_value_mean;
      } 
   }

   // Update relation distance 
   relation_pk_idx::const_iterator r_itr = relation_pk.find(boost::make_tuple(min(first_item_id, second_item_id), 
                                                                              max(first_item_id, second_item_id)));
   relation_distance_pk_idx::const_iterator rd_itr = relation_distance_pk.find(boost::make_tuple(first_item_id, second_item_id));
   if (rd_itr != relation_distance_pk.end())
   {
      relation_distance rd = *rd_itr;
      rd.distance_mean = r_itr->distance_mean;
      relation_distance_pk.replace(rd_itr, rd);
   }
   else 
   {
      item_pk_idx& item_pk = item_table.get<item_pk_idx_type>();
      item_pk_idx::iterator second_item_itr = item_pk.find(second_item_id);
      if (second_item_itr != item_pk.end() && !isnan(second_item_itr->location[0]))
      {
         relation_distance rd;
         rd.first_item_id = first_item_id; 
         rd.second_item_id = second_item_id; 
         rd.distance_mean = r_itr->distance_mean;
         relation_distance_pk.insert(rd);
      }
   }
}

