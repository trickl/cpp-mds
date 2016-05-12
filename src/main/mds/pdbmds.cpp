#include "pdbmds.h"
#include <mysql++.h>
#include <stdexcept>
#include <set>

using namespace std;
using namespace boost;
using namespace mysqlpp;

// Establishes a connection to a MySQL database server
void connect_to_db(Connection& con)
{
   // TODO: Hide password
   con.connect("mds", "127.0.0.1", "Sleeper", "kifaru", 3306);
   if (!con) throw runtime_error(string("Database connection failed: ") + con.error());
};

pdbmds::pdbmds(unsigned int max_closest_relations,
               unsigned int max_furthest_relations,
               unsigned int max_local_relations,
               double mds_boundary_factor,
               double velocity_dampening,
               unsigned int max_polar_search_depth,
               unsigned int small_sample_size,
               distance_model_type distance_model,
               bool output_sql)
   : _max_closest_relations(max_closest_relations),
     _max_furthest_relations(max_furthest_relations),
     _max_local_relations(max_local_relations),
     _mds_boundary_factor(mds_boundary_factor),
     _velocity_dampening(velocity_dampening),
     _max_polar_search_depth(max_polar_search_depth),
     _small_sample_size(small_sample_size),
     _distance_model(distance_model),
     _output_sql(output_sql)
{
}

void pdbmds::queue_ratings(std::vector<rating> ratings)
{
   try
   {
      // Establish the connection to the database server.
      Connection con(use_exceptions);
      connect_to_db(con);

      for (std::vector<rating>::const_iterator itr = ratings.begin(), end = ratings.end(); itr != end; ++itr)
      {           
         // Put the relations onto the queue
         Query query = con.query();
         query << "CALL queue_rating(" << itr->user_id << ", " << itr->item_id << ", " << itr->item_value << ")";
         if (_output_sql) cout << query.str() << " \\\\" << endl;

         query.execute();

         if (con.errnum())
         {
         }
      }
   }
   catch (const BadQuery& er)
   {
      // Handle any query errors
      cerr << "Query error: " << er.what() << endl;
      throw;
   }
   catch (const Exception& er)
   {
      // Catch-all for any other MySQL++ exceptions
      cerr << "Error: " << er.what() << endl;
      throw;
   }
}

void pdbmds::process_rating_queue(unsigned int batch_size, bool recursive, bool validate)
{
   try
   {
      // Establish the connection to the database server.
      Connection con(use_exceptions);
      connect_to_db(con);

      bool queue_empty = false;
      do 
      {
         // Process the relation queue
         Query query = con.query();
         query << "CALL process_rating_queue(" << batch_size << ", " << _max_polar_search_depth << ", " << _small_sample_size << ", " << _distance_model << ")";
         if (_output_sql) cout << query.str() << " \\\\" << endl;
 
         query.execute();

         if (con.errnum())
         {
         }

         if (recursive)
         {
            query = con.query();
            query << "SELECT COUNT(1) FROM rtq_rating_queue" << endl;
            StoreQueryResult res = query.store();
            if (!res || !res.num_rows() == 1 || int(res[0][0]) == 0)
            {
               queue_empty = true;
            } 
         }

         // Perform test
         if (validate)
         {
            query = con.query();
            query <<

               "select truth.id, \n"
               "       truth.sample_count as true_sample_count, \n"
               "       itm.sample_count as itm_sample_count, \n"
               "       truth.value_mean as true_value_mean, \n"
               "       itm.value_mean as itm_value_mean, \n"
               "       truth.value_variance as true_value_variance, \n"
               "       itm.value_variance as itm_value_variance \n"
               "from itm_item itm \n"
               "right outer join ( \n"
               "select item_id as id, count(1) as sample_count, avg(value) as value_mean, variance(value) as value_variance \n"
               "from rat_rating group by item_id ) truth on (truth.id = itm.id) \n"
               "where abs(truth.sample_count - itm.sample_count) > 0 OR \n"
               "abs(truth.value_mean - itm.value_mean) / truth.value_mean > 0.001 OR \n"
               "abs(truth.value_variance - itm.value_variance) / truth.value_variance > 0.001";

            UseQueryResult res = query.use();
            if (con.errnum()) throw runtime_error(con.error());
            else if (res.fetch_row())
            {
               throw runtime_error("Validate Item Failed.");
            }

            query = con.query();
            query <<

               "select truth.first_item_id, \n"
               "       truth.second_item_id, \n"
               "       truth.sample_count as true_sample_count, \n"
               "       rel.sample_count as rel_sample_count, \n"
               "       truth.value_abs_diff_mean as true_value_abs_diff_mean, \n"
               "       rel.value_abs_diff_mean as rel_value_abs_diff_mean, \n"
               "       truth.value_abs_diff_variance as true_value_abs_diff_variance, \n"
               "       rel.value_abs_diff_variance as rel_value_abs_diff_variance, \n"
               "       truth.first_value_mean as true_first_value_mean,\n"
               "       rel.first_value_mean as rel_first_value_mean,\n"
               "       truth.second_value_mean as true_second_value_mean,\n"
               "       rel.second_value_mean as rel_second_value_mean,\n"
               "       truth.first_value_variance as true_first_value_variance,\n"
               "       rel.first_value_variance as rel_first_value_variance,\n"
               "       truth.second_value_variance as true_second_value_variance,\n"
               "       rel.second_value_variance as rel_second_value_variance,\n"
               "       truth.value_covariance as true_value_covariance,\n"
               "       rel.value_covariance as rel_value_covariance\n"
               "       from rel_relation rel \n"
               "RIGHT OUTER JOIN \n"
               "(SELECT a.item_id as first_item_id, \n"
               "        b.item_id as second_item_id, \n"
               "        COUNT(1) as sample_count, \n"
               "        avg(abs(a.value - b.value)) as value_abs_diff_mean, \n"
               "        variance(abs(a.value - b.value)) as value_abs_diff_variance, \n"
               "        avg(a.value) as first_value_mean, \n"
               "        avg(b.value) as second_value_mean, \n" 
               "        variance(a.value) as first_value_variance, \n"
               "        variance(b.value) as second_value_variance, \n" 
               "        avg(a.value * b.value) - avg(a.value) * avg(b.value) as value_covariance \n"
               "FROM rat_rating a \n"
               "INNER JOIN rat_rating b ON (a.user_id = b.user_id and a.item_id < b.item_id) \n"
               "GROUP BY a.item_id, b.item_id) as truth ON (truth.first_item_id = rel.first_item_id AND \n"
               "                                            truth.second_item_id = rel.second_item_id) \n"
               "where abs(truth.sample_count - rel.sample_count) > 0 OR \n"
               "abs(truth.value_abs_diff_mean - rel.value_abs_diff_mean) / truth.value_abs_diff_mean > 0.001 OR \n"
               "abs(truth.value_abs_diff_variance - rel.value_abs_diff_variance) / truth.value_abs_diff_variance > 0.001 OR \n"
               "abs(truth.first_value_mean - rel.first_value_mean) / truth.first_value_mean > 0.001 OR \n"
               "abs(truth.second_value_mean - rel.second_value_mean) / truth.second_value_mean > 0.001 OR \n"
               "abs(truth.first_value_variance - rel.first_value_variance) / truth.first_value_variance > 0.001 OR \n"
               "abs(truth.second_value_variance - rel.second_value_variance) / truth.second_value_variance > 0.001 OR \n"
               "abs(truth.value_covariance - rel.value_covariance) / truth.value_covariance > 0.001";


            res = query.use();
            if (con.errnum()) throw runtime_error(con.error());
            else if (res.fetch_row())
            {
               throw runtime_error("Validate Relation Failed.");
            }
         }
      } while (recursive && !queue_empty);
   }
   catch (const BadQuery& er)
   {
      // Handle any query errors
      cerr << "Query error: " << er.what() << endl;
      throw;
   }
   catch (const Exception& er)
   {
      // Catch-all for any other MySQL++ exceptions
      cerr << "Error: " << er.what() << endl;
      throw;
   }
}

void pdbmds::queue_adjustments(std::vector<unsigned int> ids)
{
   if (ids.empty()) return;

   try
   {
      // Establish the connection to the database server.
      Connection con(use_exceptions);
      connect_to_db(con);

      // Put the item on the adjustment queue
      for (std::vector<unsigned int>::const_iterator itr = ids.begin(), end = ids.end(); itr != end; ++itr)
      {
         Query query = con.query();
         query << "CALL queue_adjustment( " << *itr << ")";
         if (_output_sql) cout << query.str() << " \\\\" << endl;

         query.execute();

         if (con.errnum())
         {
         }
      }
   }
   catch (const BadQuery& er)
   {
      // Handle any query errors
      cerr << "Query error: " << er.what() << endl;
      throw;
   }
   catch (const Exception& er)
   {
      // Catch-all for any other MySQL++ exceptions
      cerr << "Error: " << er.what() << endl;
      throw;
   }
}

void pdbmds::process_adjustment_queue(unsigned int batch_size)
{
   try
   {
      // Establish the connection to the database server.
      Connection con(use_exceptions);
      connect_to_db(con);

      // Process the item adjustment queue
      Query query = con.query();
      query << "CALL process_adjustment_queue(" << batch_size << ", " << _max_closest_relations << ", " << _max_furthest_relations << ", " <<_max_local_relations << ", " << _mds_boundary_factor << ", " << _velocity_dampening << ")";
      if (_output_sql) cout << query.str() << " \\\\" << endl;
 
      query.execute();

      if (con.errnum())
      {
      }
   }
   catch (const BadQuery& er)
   {
      // Handle any query errors
      cerr << "Query error: " << er.what() << endl;
      throw;
   }
   catch (const Exception& er)
   {
      // Catch-all for any other MySQL++ exceptions
      cerr << "Error: " << er.what() << endl;
      throw;
   }
}


// Utility for helping testing
void pdbmds::get_items(std::map<unsigned int, item_detail> &items) const
{
   // Get the item locations from the item table
   try
   {
      // Establish the connection to the database server.
      Connection con(use_exceptions);
      connect_to_db(con);

      // Put the item on the adjustment queue
      Query query = con.query();
      query << "SELECT itm.id, X(itm.location), Y(itm.location), itm.sample_count, itm.value_mean, itd.title "
               "FROM itm_item itm "
               "INNER JOIN itd_item_detail itd ON (itm.id = itd.id)";

      if (UseQueryResult res = query.use())
      {
         while (Row row = res.fetch_row())
         {
            items.insert(make_pair(row.at(0),
                                       (item_detail(int(row.at(0)),
                                        point(row.at(1), row.at(2)),
                                        int(row.at(3)),
                                        float(row.at(4)),
                                        string(row.at(5))))));
         }
         if (con.errnum())
         {
         }
      }
   }
   catch (const BadQuery& er)
   {
      // Handle any query errors
      cerr << "Query error: " << er.what() << endl;
      throw;
   }
   catch (const Exception& er)
   {
      // Catch-all for any other MySQL++ exceptions
      cerr << "Error: " << er.what() << endl;
      throw;
   }
}

void pdbmds::get_items(double min_x, double min_y, double max_x, double max_y, double resolution, std::map<unsigned int, item_detail> &items) const
{
   // Get the item locations from the item table
   try
   {
      // Establish the connection to the database server.
      Connection con(use_exceptions);
      connect_to_db(con);

      // Put the item on the adjustment queue
      Query query = con.query();
      query << "SELECT itm.id, X(itm.location), Y(itm.location), itm.sample_count, itm.value_mean, itd.title "
               "FROM itm_item itm "
               "INNER JOIN itd_item_detail itd ON (itm.id = itd.id) "
               "WHERE X(itm.location) >= " << min_x << " "
               "AND Y(itm.location) >=  " << min_y << " "
               "AND X(itm.location) <=  " << max_x << " "
               "AND Y(itm.location) <=  " << max_y << " ";

      if (UseQueryResult res = query.use())
      {
         while (Row row = res.fetch_row())
         {
            items.insert(make_pair(row.at(0),
                                       (item_detail(int(row.at(0)),
                                        point(row.at(1), row.at(2)),
                                        int(row.at(3)),
                                        float(row.at(4)),
                                        string(row.at(5))))));
         }
         if (con.errnum())
         {
         }
      }
   }
   catch (const BadQuery& er)
   {
      // Handle any query errors
      cerr << "Query error: " << er.what() << endl;
      throw;
   }
   catch (const Exception& er)
   {
      // Catch-all for any other MySQL++ exceptions
      cerr << "Error: " << er.what() << endl;
      throw;
   }
}

void pdbmds::get_hull(vector<pair<unsigned int, point> > &hull_ids) const
{
   // Get the item locations from the item table
   try
   {
      // Establish the connection to the database server.
      Connection con(use_exceptions);
      connect_to_db(con);

      // Put the item on the adjustment queue
      Query query = con.query();
      query << "SELECT itm.id, X(itm.location), Y(itm.location) "
               "FROM chp_convex_hull_points chp "
               "INNER JOIN itm_item itm ON (chp.point_id = itm.id) "
               "ORDER BY chp.cyclic_order";

      if (UseQueryResult res = query.use())
      {
         while (Row row = res.fetch_row())
         {
            hull_ids.push_back(make_pair((unsigned int) row[0], point(row[1], row[2])));
         }
         if (con.errnum())
         {
         }
      }
   }
   catch (const BadQuery& er)
   {
      // Handle any query errors
      cerr << "Query error: " << er.what() << endl;
      throw;
   }
   catch (const Exception& er)
   {
      // Catch-all for any other MySQL++ exceptions
      cerr << "Error: " << er.what() << endl;
      throw;
   }
}

