#ifndef _PMDS_H
#define _PMDS_H

#include <boost/numeric/ublas/fwd.hpp>
#include <gsl/gsl_rng.h>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>
#include <ssrc/spatial/kd_tree.h> 
#include <geometry/point/point.h>

// Symmetric, first item id < second item id
struct relation
{
   unsigned int first_item_id;
   unsigned int second_item_id;
   unsigned int sample_count;
   float        distance_mean;
   float        distance_variance;
   float        value_covariance;
};

// Cartesian, contains duplicates [id, second_id], [second_id, id]
// For the purpose od distance indexes
// An entry only exists in the table if the item is located, i.e.
// has a valid location
struct relation_distance
{
   unsigned int first_item_id;
   unsigned int second_item_id;
   float        distance_mean;
};

struct relation_pk_idx_type{};
struct relation_reverse_idx_type{};


// First item id < Second item id, Not Range Searchable
typedef boost::multi_index::multi_index_container<
   relation,
   boost::multi_index::indexed_by<
      boost::multi_index::ordered_unique<boost::multi_index::tag<relation_pk_idx_type>,
         boost::multi_index::composite_key<
            relation,
            boost::multi_index::member<relation, unsigned int, &relation::first_item_id>,
            boost::multi_index::member<relation, unsigned int, &relation::second_item_id>
         >
      >,
      boost::multi_index::ordered_unique<boost::multi_index::tag<relation_reverse_idx_type>,
         boost::multi_index::composite_key<
            relation,
            boost::multi_index::member<relation, unsigned int, &relation::second_item_id>,
            boost::multi_index::member<relation, unsigned int, &relation::first_item_id>
         >
      >
   >
> relation_table_type;

typedef relation_table_type::index<relation_pk_idx_type>::type relation_pk_idx;
typedef relation_table_type::index<relation_reverse_idx_type>::type relation_reverse_idx;

struct relation_distance_pk_idx_type{};
struct relation_distance_dist_idx_type{};

// Range Searchable
typedef boost::multi_index::multi_index_container<
   relation_distance,
   boost::multi_index::indexed_by<
      boost::multi_index::ordered_unique<boost::multi_index::tag<relation_distance_pk_idx_type>,
         boost::multi_index::composite_key<
            relation_distance,
            boost::multi_index::member<relation_distance, unsigned int, &relation_distance::first_item_id>,
            boost::multi_index::member<relation_distance, unsigned int, &relation_distance::second_item_id>
         >
      >,
      boost::multi_index::ordered_non_unique<boost::multi_index::tag<relation_distance_dist_idx_type>,
         boost::multi_index::composite_key<
            relation_distance,
            boost::multi_index::member<relation_distance, unsigned int, &relation_distance::first_item_id>,
            boost::multi_index::member<relation_distance, float, &relation_distance::distance_mean>
         >
      >
   >
> relation_distance_table_type;

typedef relation_distance_table_type::index<relation_distance_pk_idx_type>::type relation_distance_pk_idx;
typedef relation_distance_table_type::index<relation_distance_dist_idx_type>::type relation_distance_dist_idx;

struct item
{
   unsigned int id;
   point location;
   unsigned int sample_count;
   float value_mean;
   float value_variance;
   float distance_mean;
};

struct item_pk_idx_type{};

typedef boost::multi_index::multi_index_container<
   item,
   boost::multi_index::indexed_by<
      boost::multi_index::random_access<>,
      boost::multi_index::hashed_unique<boost::multi_index::tag<item_pk_idx_type>,
         boost::multi_index::member<item, unsigned int, &item::id>
      >
   >
> item_table_type;

typedef item_table_type::index<item_pk_idx_type>::type item_pk_idx;

typedef ssrc::spatial::kd_tree<point, unsigned int> item_spatial_idx_type; 

// Hybrid Proprietary Progressive MDS Algorithm
// Designed for incremental relation based updates and very scalable
// Only for reducing data into a two-dimensional space (visual representation).
// Based on a force-directed iterative model
// R is a n x n relational matrix of dissimilarities (assumed to be Euclidean relations)
// p is the dimensionality of the target space 
// X are the projected points
class pmds
{
public:
   pmds(unsigned int v_max = 10,
        unsigned int h_max = 10,
        double spring_constant = 1.,
        unsigned int max_polar_search_depth = 5);

   void process_relation(unsigned int first_item_id, unsigned int second_item_id, float first_item_value, float second_item_value);
   void adjust_item(unsigned int id);

   boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> item_locations_as_matrix() const;

private:
   void update_item(unsigned int id, float value, float compare_value);
   void update_relation(unsigned int first_item_id, unsigned int second_item_id, float value, float compare_value);
   void locate_item(unsigned int id, point& location);

   gsl_rng* _rng;
   const unsigned int _v_max;
   const unsigned int _h_max;
   const double _spring_constant;
   const unsigned int _max_polar_search_depth;

   relation_table_type          relation_table;
   relation_distance_table_type relation_distance_table;

   item_table_type       item_table;
   item_spatial_idx_type item_spatial_idx; 
};
          

#endif
