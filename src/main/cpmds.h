#ifndef _CPMDS_H
#define _CPMDS_H

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>
#include <geometry/point/point.h>
#include <geometry/line/line.h>
#include <set>
#include <queue>
#include <map>
#include <list>

// Symmetric, first item id < second item id
struct relation
{
   unsigned int first_itm_id;
   unsigned int second_itm_id;
   unsigned int sample_count;
   double        distance;
   double        distance_variance;
   double        value_covariance;
};

// Cartesian, contains duplicates [id, second_itm_id], [second_itm_id, id]
// For the purpose od distance indexes
// An entry only exists in the table if the item is located, i.e.
// has a valid location
struct relation_distance
{
   unsigned int first_itm_id;
   unsigned int second_itm_id;
   double        distance;
   double        weight;
};

struct item
{
   unsigned int id;
   point location;
   unsigned int sample_count;
   unsigned int connectivity;
   double distance_weight_sum;
   double weight_sum;
   double value_mean;
   double value_variance;
   double radius;
   double curvature;
   double curvature_delta;
   bool   is_located;
};

struct item_adjustment
{
   unsigned int id;
   int connectivity_delta;
   double distance_weight_sum_delta;
   double weight_sum_delta;
   double curvature_delta;
};

struct adjacency_face
{
   unsigned int first_itm_id;
   unsigned int second_itm_id;
   unsigned int third_itm_id;
};

struct boundary_self_join
{
   unsigned int first_itm_id;
   unsigned int second_itm_id;
};

struct relation_pk_type{};
struct relation_reverse_idx_type{};

// First item id < Second item id, Not Range Searchable
typedef boost::multi_index::multi_index_container<
   relation,
   boost::multi_index::indexed_by<
      boost::multi_index::ordered_unique<boost::multi_index::tag<relation_pk_type>,
         boost::multi_index::composite_key<
            relation,
            boost::multi_index::member<relation, unsigned int, &relation::first_itm_id>,
            boost::multi_index::member<relation, unsigned int, &relation::second_itm_id>
         >
      >,
      boost::multi_index::ordered_unique<boost::multi_index::tag<relation_reverse_idx_type>,
         boost::multi_index::composite_key<
            relation,
            boost::multi_index::member<relation, unsigned int, &relation::second_itm_id>,
            boost::multi_index::member<relation, unsigned int, &relation::first_itm_id>
         >
      >
   >
> relation_table_type;

typedef relation_table_type::index<relation_pk_type>::type relation_pk;
typedef relation_table_type::index<relation_reverse_idx_type>::type relation_reverse_idx;

struct relation_distance_pk_type{};
struct relation_distance_dist_idx_type{};

// Range Searchable
typedef boost::multi_index::multi_index_container<
   relation_distance,
   boost::multi_index::indexed_by<
      boost::multi_index::ordered_unique<boost::multi_index::tag<relation_distance_pk_type>,
         boost::multi_index::composite_key<
            relation_distance,
            boost::multi_index::member<relation_distance, unsigned int, &relation_distance::first_itm_id>,
            boost::multi_index::member<relation_distance, unsigned int, &relation_distance::second_itm_id>
         >
      >,
      boost::multi_index::ordered_non_unique<boost::multi_index::tag<relation_distance_dist_idx_type>,
         boost::multi_index::composite_key<
            relation_distance,
            boost::multi_index::member<relation_distance, unsigned int, &relation_distance::first_itm_id>,
            boost::multi_index::member<relation_distance, double, &relation_distance::distance>
         >
      >
   >
> relation_distance_table_type;

typedef relation_distance_table_type::index<relation_distance_pk_type>::type relation_distance_pk;
typedef relation_distance_table_type::index<relation_distance_dist_idx_type>::type relation_distance_dist_idx;

struct item_pk_type{};
struct item_curvature_idx_type{};

typedef boost::multi_index::multi_index_container<
   item,
   boost::multi_index::indexed_by<
      boost::multi_index::random_access<>,
      boost::multi_index::hashed_unique<boost::multi_index::tag<item_pk_type>,
         boost::multi_index::member<item, unsigned int, &item::id>
      >,
      boost::multi_index::ordered_non_unique<boost::multi_index::tag<item_curvature_idx_type>,
         boost::multi_index::member<item, double, &item::curvature>
      >
   >
> item_table_type;

typedef item_table_type::index<item_pk_type>::type item_pk;
typedef item_table_type::index<item_curvature_idx_type>::type item_curvature_idx;

struct adjacency_face_pk_type{};
struct adjacency_face_second_itm_idx_type{};
struct adjacency_face_third_itm_idx_type{};

// First item id < Second item id, Not Range Searchable
typedef boost::multi_index::multi_index_container<
   adjacency_face,
   boost::multi_index::indexed_by<
      boost::multi_index::ordered_unique<boost::multi_index::tag<adjacency_face_pk_type>,
         boost::multi_index::composite_key<
            adjacency_face,
            boost::multi_index::member<adjacency_face, unsigned int, &adjacency_face::first_itm_id>,
            boost::multi_index::member<adjacency_face, unsigned int, &adjacency_face::second_itm_id>
         >
      >,
      boost::multi_index::ordered_unique<boost::multi_index::tag<adjacency_face_second_itm_idx_type>,
         boost::multi_index::composite_key<
            adjacency_face,
            boost::multi_index::member<adjacency_face, unsigned int, &adjacency_face::second_itm_id>,
            boost::multi_index::member<adjacency_face, unsigned int, &adjacency_face::third_itm_id>
         >
      >,
      boost::multi_index::ordered_unique<boost::multi_index::tag<adjacency_face_third_itm_idx_type>,
         boost::multi_index::composite_key<
            adjacency_face,
            boost::multi_index::member<adjacency_face, unsigned int, &adjacency_face::third_itm_id>,
            boost::multi_index::member<adjacency_face, unsigned int, &adjacency_face::first_itm_id>
         >
      >
   >
> adjacency_face_table_type;


typedef adjacency_face_table_type::index<adjacency_face_pk_type>::type adjacency_face_pk;
typedef adjacency_face_table_type::index<adjacency_face_second_itm_idx_type>::type adjacency_face_second_itm_idx;
typedef adjacency_face_table_type::index<adjacency_face_third_itm_idx_type>::type adjacency_face_third_itm_idx;

struct boundary_self_join_pk_type{};
struct boundary_self_join_reverse_idx_type{};

typedef boost::multi_index::multi_index_container<
   boundary_self_join,
   boost::multi_index::indexed_by<
      boost::multi_index::ordered_unique<boost::multi_index::tag<boundary_self_join_pk_type>,
         boost::multi_index::composite_key<
            boundary_self_join,
            boost::multi_index::member<boundary_self_join, unsigned int, &boundary_self_join::first_itm_id>,
            boost::multi_index::member<boundary_self_join, unsigned int, &boundary_self_join::second_itm_id>
         >
      >,
      boost::multi_index::ordered_unique<boost::multi_index::tag<boundary_self_join_reverse_idx_type>,
         boost::multi_index::composite_key<
            boundary_self_join,
            boost::multi_index::member<boundary_self_join, unsigned int, &boundary_self_join::second_itm_id>,
            boost::multi_index::member<boundary_self_join, unsigned int, &boundary_self_join::first_itm_id>
         >
      >
   >
> boundary_self_join_table_type;

typedef boundary_self_join_table_type::index<boundary_self_join_pk_type>::type boundary_self_join_pk;
typedef boundary_self_join_table_type::index<boundary_self_join_reverse_idx_type>::type boundary_self_join_reverse_idx;

// Circle Packing MDS Algorithm
class cpmds
{
public:
   enum graph_step_type
   {
      GRAPH_STEP,
      ITEM_STEP,
      MOVE_STEP
   };

   enum measure_type
   {
      WEIGHTED_DISTANCE,
      DISTANCE,
      RANK
   };

   enum stress_type
   {
      CURVATURE_STRESS,
      DISTANCE_STRESS
   };

   enum aggregate_type
   {
      AVG,
      MAX,
      MIN
   };

public:
   cpmds();
   void process_relation(unsigned int first_itm_id, unsigned int second_itm_id, double first_item_value, double second_item_value);

   const item_table_type &items() const;
   item get_item_by_id(unsigned int id) const;
   bool is_interior_item(unsigned int id) const;
   bool is_flattened_item(unsigned int id) const;
   bool is_boundary_item(unsigned int id) const;
   bool is_boundary_self_join(unsigned int first_itm_id, unsigned int second_itm_id) const;
   bool is_boundary_edge(unsigned int first_itm_id, unsigned int second_itm_id) const;
   bool has_neighbour(unsigned int id, item& neighbour) const;
   void get_perimeter(unsigned int id, std::list<unsigned int>& perimeter, unsigned int start_id = 0, unsigned int end_id = 0) const;
   void get_graph_straight_line_drawing(std::map<unsigned int, point> &vertex_locations);
   void get_graph_afl_drawing(std::map<unsigned int, point> &vertex_locations);
   void get_hop_counts(unsigned int id, std::map<unsigned int, unsigned int> &hop_counts, unsigned int max_hop_count = 0) const;

   void load_item(const item& itm);
   void load_relation(const relation& rel);
   void load_relation_distance(const relation_distance& rdt);
   void load_adjacency_face(const adjacency_face& adf);

   void flatten_graph(graph_step_type step = GRAPH_STEP);
   void update_graph(graph_step_type step = GRAPH_STEP);
   void output_graph();
   void output_statistics();
   void reset_flatten();

   std::pair<double, double> get_weighted_measure(measure_type s_type = WEIGHTED_DISTANCE, aggregate_type a_type = AVG) const;
   std::pair<double, double> get_weighted_measure(unsigned int first_itm_id, measure_type s_type = WEIGHTED_DISTANCE, aggregate_type a_type = AVG, unsigned int second_itm_id = 0, unsigned int max_hop_count = 1) const;
   std::pair<double, double> get_weighted_measure(unsigned int first_itm_id, unsigned int second_itm_id, measure_type type = WEIGHTED_DISTANCE) const;
   std::pair<double, double> get_weighted_measure(unsigned int id, const adjacency_face& f, measure_type type = WEIGHTED_DISTANCE, aggregate_type a_type = AVG) const;
   double get_stress_delta(const std::map<unsigned int, item_adjustment> &adjustments, stress_type = DISTANCE_STRESS) const;

private:

   void get_inward_ordering(std::list<adjacency_face> &ordering);

   bool flatten_next_item(graph_step_type step);
   bool flatten_item(unsigned int id);

   void insert_item(unsigned int id, std::map<unsigned int, item_adjustment> &adjustments, bool update);
   void insert_item(unsigned int id, adjacency_face_pk::iterator f_nearest_itr, std::map<unsigned int, item_adjustment> &adjustments, bool update);
   void remove_item(unsigned int id, std::map<unsigned int, item_adjustment> &adjustments, bool update);

   bool update_next_item();
   void update_item(unsigned int id);
   void update_item(unsigned int id, double value, double compare_value);

   void get_decrease_curvature_transform(unsigned int id, std::list<unsigned int>& fan_trfm_perimeter) const;
   void get_increase_curvature_transform(unsigned int id, std::list<unsigned int>& fan_trfm_perimeter) const;
   bool is_beneficial_adjustment(unsigned int id, double curvature_adjustment) const;
   void fan_transform(const std::list<unsigned int> &perimeter, std::map<unsigned int, item_adjustment> &adjustments, bool update);
   void fan_transform(const std::list<unsigned int> &perimeter, std::map<unsigned int, item_adjustment> &adjustments) const;

   void update_flattened_boundary(unsigned int id);

   void update_relation(unsigned int first_itm_id, unsigned int second_itm_id, double value, double compare_value);
   void on_radius_change(unsigned int id, double r);
   void layout_adjacency_graph();

   adjacency_face_pk::const_iterator find_adjacency_face(unsigned int first_itm_id, unsigned int second_itm_id, adjacency_face& f) const;
   adjacency_face_pk::const_iterator find_k_nearest_face(unsigned int id, unsigned int k) const;

   void find_joins(std::list<unsigned int>::const_iterator left_itr,
                   std::list<unsigned int>::const_iterator left_end_itr,
                   std::list<unsigned int>::const_iterator right_itr,
                   std::list<unsigned int>::const_iterator right_end_itr);
   void remove_boundary_self_join(unsigned int first_itm_id, unsigned int second_itm_id);
   void remove_boundary_self_joins(unsigned int id);
   void check_boundary_self_join(unsigned int first_itm_id, unsigned int second_itm_id);
   bool is_planar() const;
   bool is_separating_item(unsigned int id) const;
   unsigned int get_connectivity(unsigned int id) const;
   bool validate_curvature() const;
   bool validate_curvature(unsigned int id) const;
   bool validate_distance() const;
   bool validate_distance(unsigned int id) const;
   void revert_ideal_radii();

   static double subtended_angle(double r_a, double r_b, double r_c);
   static double subtended_half_angle(double r_a, double r_b);
   static adjacency_face make_face(unsigned int first_itm_id, unsigned int second_itm_id, unsigned int third_itm_id);

   // System tables
   adjacency_face_table_type adjacency_face_table;
   relation_table_type relation_table;
   relation_distance_table_type relation_distance_table;
   item_table_type item_table;
   boundary_self_join_table_type boundary_self_join_table;

   // For analysis
   std::list<unsigned int> _flattened_boundary;
   std::map<unsigned int, std::list<unsigned int>::iterator > _flattened_boundary_items;
   std::set<unsigned int> _flattened_items;
   std::set<unsigned int> _update_pending_items;
   std::set<unsigned int> _converted_exterior_items;
   unsigned int _flatten_item;
   unsigned int _graph_iteration;
   std::map<unsigned int, std::pair<unsigned int, unsigned int> > _boundary_connectivity;
   unsigned int _focus;
   double _flattened_curvature_sum;
   std::map<unsigned int, double> _curvature_adjustments;
};

#endif
