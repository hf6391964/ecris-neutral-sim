// Copyright (c) 2015 GeometryFactory
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_IO_STL_READER_H
#define CGAL_IO_STL_READER_H

#include <CGAL/array.h>
#include <boost/cstdint.hpp> 
#include <vector>
#include <map>
#include <iostream>
#include <string>

namespace CGAL{

  bool
  read_STL( std::istream& input,
            std::vector< cpp11::array<double,3> >& points,
            std::vector< cpp11::array<int,3> >& facets,
            bool verbose = false);
} // namespace CGAL

#endif // CGAL_IO_STL_READER_H
