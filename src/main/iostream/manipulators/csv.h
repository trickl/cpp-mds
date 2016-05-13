#ifndef _CSV_H
#define _CSV_H

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <algorithm>

class csv_ostream_manipulator
{
public:
   csv_ostream_manipulator(std::ostream &os = std::cout) : _os(os) {};

   // Copy manipulator params
   csv_ostream_manipulator &copy(const csv_ostream_manipulator& rhs) {return *this;};

   template<class T>
   std::ostream& operator <<(const std::vector<T> &data)
   {
      for (unsigned int i = 0; i < data.size(); i++)
      {
         unsigned int j = 0;
         for (; j < T::size() - 1; j++)
         {
            _os << data[i](j) << ", ";
         }

         _os << data[i](j) << std::endl;
      }

      return this->_os;
   }

   template<class T>
   std::ostream& operator <<(const T &data)
   {
      for (unsigned int i = 0; i < data.size1(); i++)
      {
         unsigned int j = 0;
         for (;j < data.size2() - 1; j++)
         {
            _os << data(i, j) << ", ";
         }

         _os << data(i, j) << std::endl;
      }

      return this->_os;
   }


private:
   std::ostream &_os;
};

class csv_istream_manipulator
{
public:
   csv_istream_manipulator(std::istream &is = std::cin) : _is(is) {};

   // Copy manipulator params
   csv_istream_manipulator &copy(const csv_istream_manipulator& rhs) {return *this;}; 

   void read_cells(std::vector<std::vector<std::string> >& cells, long &max_col) const
   {
      std::string line;
      while(getline(_is, line))
      {
         max_col = 1;
         cells.push_back(std::vector<std::string>());

         while (line.find(",") != std::string::npos)
         {
            cells.back().push_back(line.substr(0, line.find(",")));
            line = line.substr(line.find(",") + 1);
            max_col++;
         }

         cells.back().push_back(line);
      }
   }

   template <class T>
   const std::istream &operator >>(T &data) const
   {
      std::vector<std::vector<std::string> > cells;
      long max_col = 0;

      read_cells(cells, max_col);

      // Try to convert all cells to doubles
      data.resize(cells.size(), max_col);

      for (unsigned long i = 0; i < cells.size(); i++)
      {
         for (unsigned long j = 0; j < cells[i].size(); j++)
         {
            double d = 0;
            std::istringstream token;
            token.str(cells[i][j]);
            token >> d;
            data(i, j) = d;
         }
      }

      return this->_is;
   }

   template <class T>
   const std::istream &operator >>(std::vector<T> &data) const
   {
      std::vector<std::vector<std::string> > cells;

      long max_col = 0;
      read_cells(cells, max_col);

      // Try to convert all cells to doubles
      data.resize(cells.size());

      for (unsigned long i = 0; i < cells.size(); i++)
      {
         for (unsigned long j = 0; j < cells[i].size(); j++)
         {
            double d = 0;
            std::istringstream token;
            token.str(cells[i][j]);
            token >> d;
            data[i](j) = d;
         }
      }

      return this->_is;
   }

private:
   std::istream &_is;
};

// Parameterless construction, e.g. cout << csv << M
// N.B. We need a separate function to call the constructor as you cannot
// get the pointer to a constructor directly.
csv_ostream_manipulator operator<<(std::ostream &os, csv_ostream_manipulator (*f_ptr)(std::ostream&))
{
   return f_ptr(os);
}

csv_ostream_manipulator csv(std::ostream& os)
{
   return csv_ostream_manipulator(os);
}

// Explicit construction allows parameters, e.g. cout << csv_ostream_manipulator(params) << M
csv_ostream_manipulator operator<<(std::ostream &os, csv_ostream_manipulator rhs)
{
   // Must construct with ostream reference and copy manipulator parameters
   return csv_ostream_manipulator(os).copy(rhs);
}

// Parameterless construction, e.g. cin >> csv >> M
csv_istream_manipulator operator>>(std::istream &is, csv_istream_manipulator (*f_ptr)(std::istream&))
{
   return f_ptr(is);  
}

csv_istream_manipulator csv(std::istream& is)
{
   return csv_istream_manipulator(is);
}

// Explicit construction allows parameters, e.g. cin >> csv_istream_manipulator(params) >> M
csv_istream_manipulator operator>>(std::istream &is, csv_istream_manipulator rhs)
{
   // Must construct with istream reference and copy manipulator parameters
   return csv_istream_manipulator(is).copy(rhs);
}

#endif
