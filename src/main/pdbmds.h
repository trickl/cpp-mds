#ifndef _PDBMDS_H
#define _PDBMDS_H

#include <boost/numeric/ublas/fwd.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>
#include <geometry/point/point.h>
#include <vector>
#include <map>

// Progressive Database MDS Algorithm
// Using a database backend for persistent storage and low memory footprint
// Designed for incremental relation based updates and very scalable
// Only for reducing data into a two-dimensional space (visual representation).
// Based on a force-directed iterative model
class pdbmds
{
public:
   struct rating
   {
      unsigned int user_id;
      unsigned int item_id;
      double item_value;

      rating(unsigned int user_id, unsigned int item_id, double item_value)
         : user_id(user_id), item_id(item_id), item_value(item_value) {};
   };

   struct item_detail
   {
      unsigned int id;
      point location;
      unsigned int sample_count;
      float value_mean;
      std::string title;

      item_detail(unsigned int id, point location, unsigned int sample_count, float value_mean, std::string title)
         : id(id), location(location), sample_count(sample_count), value_mean(value_mean), title(title) {};
   };

   enum distance_model_type
   {
      abs_difference_distance_model = 0,
      correlation_distance_model = 1
   };

   pdbmds(unsigned int max_closest_relations = 10,
        unsigned int max_furthest_relations = 5,
        unsigned int max_local_relations = 5,
        double mds_boundary_factor = 2,
        double velocity_dampening = 0.5,
        unsigned int max_polar_search_depth = 5,
        unsigned int small_sample_size = 1,
        distance_model_type distance_model = abs_difference_distance_model,
        bool output_sql = false);

   void queue_ratings(std::vector<rating> ratings);
   void queue_adjustments(std::vector<unsigned int> ids);

   void process_rating_queue(unsigned int batch_size, bool recursive = true, bool validate = false);
   void process_adjustment_queue(unsigned int batch_size);

   void get_items(std::map<unsigned int, item_detail> &items) const;
   void get_items(double min_x, double min_y, double max_x, double max_y, double resolution, std::map<unsigned int, item_detail> &items) const;
   void get_hull(std::vector<std::pair<unsigned int, point> > &ids) const;

private:

   const unsigned int _max_closest_relations;
   const unsigned int _max_furthest_relations;
   const unsigned int _max_local_relations;
   const double _mds_boundary_factor;
   const double _velocity_dampening;
   const unsigned int _max_polar_search_depth;
   const unsigned int _small_sample_size;
   const distance_model_type _distance_model;
   const bool _output_sql;
};
          

#endif
