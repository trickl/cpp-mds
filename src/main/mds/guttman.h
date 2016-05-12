#ifndef _GUTTMAN_H
#define _GUTTMAN_H

#include <boost/numeric/ublas/fwd.hpp>

void guttman_transform(const boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> &X,
                       const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &R,
                       const boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower, boost::numeric::ublas::column_major> &W,
                       boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> &G);

#endif
