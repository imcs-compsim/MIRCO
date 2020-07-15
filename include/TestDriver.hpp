/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author  H. Carter Edwards  <hcedwar@sandia.gov>
 */

#ifndef main_TestDriver_hpp
#define main_TestDriver_hpp

#include <map>
#include <string>
#include <iosfwd>

#include <util/Parallel.hpp>

namespace phdmesh {

typedef
void (*TestSubprogram)( ParallelMachine , std::istream & );

typedef std::map< std::string , TestSubprogram > TestDriverMap ;

int test_driver( ParallelMachine , std::istream & , const TestDriverMap & );
int test_driver( ParallelMachine , const TestDriverMap & ,
                 int argc , const char * const * argv );

}

#endif

