#include "isomap.h"
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_cdf.h>

using namespace std;
using namespace boost;
using namespace boost::numeric::ublas;

// Flow property for edges
struct edge_flow_t {
  typedef edge_property_tag kind;
};

void rkisomap(const symmetric_matrix<double, lower, column_major> &R,
              symmetric_matrix<double, lower, column_major> &S,
              unsigned int k,
              double m,
              double flow_intolerance)
{
   if (R.size1() != R.size2()) throw invalid_argument("Relation matrix R must be square.");
   const unsigned int n = R.size1(); // Number of points

   // Create a neighbourhood graph using all the known dissimilarities
   // Directed, as we create both an in flow and out flow for each edge (for robustness)
   typedef adjacency_list< listS, vecS, directedS, no_property, property< edge_weight_t, double, property< boost::edge_flow_t, int > >, disallow_parallel_edge_tag> graph_t;
   typedef graph_t::vertex_descriptor vertex_t;
   typedef graph_t::edge_descriptor edge_t;
   graph_t graph(n);

   // Add an edge for each dissimilarity, weighted by the distance
   property_map<graph_t, edge_weight_t>::type weight_map = get(edge_weight, graph);
   property_map<graph_t, boost::edge_flow_t>::type   flow_map   = get(boost::edge_flow_t(), graph);
   for (unsigned int i = 0; i < n; ++i)
   {
      // Sort the nearest neighbours to this point
      multimap<double, unsigned int> sorted_distances;
      for (unsigned int j = 0; j < n; ++j)
      {
         if (i != j) sorted_distances.insert(std::pair<double, unsigned int>(R(i, j), j));
      }

      unsigned int pos = 0; // Add the k nearest neighbours to the graph
      for (map<double, unsigned int>::const_iterator itr = sorted_distances.begin(), end = sorted_distances.end();
           pos < k && itr != end; ++itr, ++pos)
      {
         // Add an edge weighted with the distance
         edge_t e;
         bool exists;

         tie(e, exists) = add_edge(i, itr->second, graph);
         weight_map[e] = pow(R(i, itr->second), m);
         flow_map[e]   = 0;
      }
   }

   // Create a distance matrix of the shortest path between two nodes
   // using Dijkstra's shortest paths
   S.resize(n);
   bool is_robust = false;
   int max_robust_iterations = 3;
   for (int robust_itr = 0; is_robust == false && robust_itr <= max_robust_iterations; ++robust_itr)
   {
      property_map<graph_t, vertex_index_t>::type vertex_index_map = get(vertex_index, graph);

      for (unsigned int i = 0; i < n; ++i)
      {
         std::vector<vertex_t> predecessors(num_vertices(graph));
         std::vector<double> shortest_paths(num_vertices(graph));

         // TODO: Given symmetrical dissimilarities, only necessary to find shortest paths
         // to nodes > j + 1
         dijkstra_shortest_paths(graph, vertex(i, graph), &predecessors[0], &shortest_paths[0],
                                 weight_map, vertex_index_map, 
                                 less<double>(), // Distance compare
                                 closed_plus<double>(), // Distance combine
                                 numeric_limits<double>::infinity(), // Max distance
                                 0.0, // Zero distance
                                 default_dijkstra_visitor() // Vistor Pattern 
                                 );

         graph_traits<graph_t>::vertex_iterator vi, vend;
         for (tie(vi, vend) = vertices(graph); vi != vend; ++vi)
         {
            // Transformed distance matrix is the shortest path from vertex i to every other vertex
            S(i, vertex_index_map[*vi]) = pow(shortest_paths[*vi], 1.0 / m);
         }

         // Now calculate the total flows in the graph in order to remove any outliers
         // that are connecting manifolds. This adds the robustness.
         for (unsigned int j = 0; j < n; ++j)
         {
            if (i != j)
            {
               edge_t e;
               bool exists;
               tie(e, exists)  = edge(vertex_index_map[predecessors[j]], j, graph);

               // Number of shortest paths through this edge
               if (exists)
               { 
                 flow_map[e]++;
               }
            }
         }
      }

      // Calculate the distribution of flow across edges
      int flow_array[num_edges(graph)];
      long edge_index = 0;
      graph_traits<graph_t>::edge_iterator e_i, e_end;
      for (tie(e_i, e_end) = edges(graph); e_i != e_end; ++e_i)
      {
         int edge_flow = flow_map[*e_i];
         flow_array[edge_index++] = edge_flow;
      }

      double flow_sigma = gsl_stats_int_sd(flow_array, 1, edge_index);
      double flow_mean  = gsl_stats_int_mean(flow_array, 1, edge_index);

      list<edge_t> edges_pending_removal;
      for (tie(e_i, e_end) = edges(graph); e_i != e_end; ++e_i)
      {
         // Cumulative distribution function P
         double cumulative_density = gsl_cdf_gaussian_P(flow_map[*e_i] - flow_mean, flow_sigma);

         if (flow_intolerance > 1 - cumulative_density)
         {
            // This edge has too much flow going through it, check that we don't orphan either vertex by removal
            if (out_degree(source(*e_i, graph), graph) > 1 && out_degree(target(*e_i, graph), graph) > 1)
            {
               edges_pending_removal.push_back(*e_i);
            }
         }
      }  
   
      if (edges_pending_removal.size() > 0)
      {
         cout << "Iteration: " << robust_itr << " failed robust test." << endl;
         // Safe removal without invalidating the iterator
         for (list<edge_t>::const_iterator itr = edges_pending_removal.begin(), end = edges_pending_removal.end(); itr != end; ++itr)
         {
            remove_edge(*itr, graph);
         }

         // Reset the flow for all remaining edges
         for (tie(e_i, e_end) = edges(graph); e_i != e_end; ++e_i)
         {
            flow_map[*e_i] = 0;
         }
      } 
      else
      {
         cout << "Iteration: " << robust_itr << " passed robust test." << endl;
         is_robust = true;
      }
   }
}
