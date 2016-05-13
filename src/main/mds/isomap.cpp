#include "isomap.h"
#include <stdexcept>
#include <vector>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include "cmds.h"

using namespace std;
using namespace boost;
using namespace boost::numeric::ublas;

void isomap(const symmetric_matrix<double, lower, column_major> &R,
            symmetric_matrix<double, lower, column_major> &S,
            unsigned int k)
{
   if (R.size1() != R.size2()) throw invalid_argument("Relation matrix R must be square.");
   const unsigned int n = R.size1(); // Number of points

   // Create a neighbourhood graph using all the known dissimilarities
   // Undirected as we assume symmetrical dissimlarities and we need a weight for each edge
   // No parallel edges allowed
   typedef adjacency_list< listS, vecS, undirectedS, no_property, property< edge_weight_t, double>, disallow_parallel_edge_tag> graph_t;
   typedef graph_t::vertex_descriptor vertex_t;
   typedef graph_t::edge_descriptor edge_t;
   graph_t graph(n);

   // Add an edge for each dissimilarity, weighted by the distance
   property_map<graph_t, edge_weight_t>::type weight_map = get(edge_weight, graph);
   for (unsigned int i = 0; i < n; ++i)
   {
      // Sort the nearest neighbours to this point
      multimap<double, unsigned int> distance_index_map;
      for (unsigned int j = 0; j < n; ++j)
      {
         if (i != j) distance_index_map.insert(std::pair<double, unsigned int>(R(i, j), j));
      }

      unsigned int pos = 0; // Add the k nearest neighbours to the graph
      for (map<double, unsigned int>::const_iterator itr = distance_index_map.begin(), end = distance_index_map.end();
           pos < k && itr != end; ++itr, ++pos)
      {
         // Add an edge weighted with the distance
         edge_t e;
         bool exists;

         // Note undirected and no parallel edges, so this allows flow in and out
         tie(e, exists) = add_edge(i, itr->second, graph);
         weight_map[e] = R(i, itr->second);
      }
   }

   // Create a distance matrix of the shortest path between two nodes
   // using Dijkstra's shortest paths
   S.resize(n);
   for (unsigned int i = 0; i < n; ++i)
   {
      std::vector<graph_t::vertex_descriptor> adjacent_nodes(num_vertices(graph));
      std::vector<double> shortest_paths(num_vertices(graph));
      property_map<graph_t, vertex_index_t>::type index_map = get(vertex_index, graph);

      // TODO: Given symmetrical dissimilarities, only necessary to find shortest paths
      // to nodes > j + 1

      dijkstra_shortest_paths(graph, vertex(i, graph), &adjacent_nodes[0], &shortest_paths[0],
                              weight_map, index_map, 
                              less<double>(), // Distance compare
                              closed_plus<double>(), // Distance combine
                              numeric_limits<double>::infinity(), // Max distance
                              0.0, // Zero distance
                              default_dijkstra_visitor() // Vistor Pattern
                              );


      graph_traits<graph_t>::vertex_iterator vi, vend;
      for (tie(vi, vend) = vertices(graph); vi != vend; ++vi)
      {
         S(i, index_map[*vi]) = shortest_paths[*vi];
      }
   } 
}
