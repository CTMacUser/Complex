//  Boost Complex Numbers, general ops header file  --------------------------//

//  Copyright 2013 Daryle Walker.
//  Distributed under the Boost Software License, Version 1.0.  (See the
//  accompanying file LICENSE_1_0.txt or a copy at
//  <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

/** \file   boost/math/complex.hpp
    \brief  Collection of support functions for complex numbers.

    \author  Daryle Walker

    \version  0.1

    \copyright  Boost Software License, version 1.0

    Contains the declarations (and definitions) of function templates, including
    operators, that work with the `complex_it` and `complex_rt` class templates,
    that are not core to either template, or interact with both templates.

    \warning  This library requires C++2011 features.
 */

#ifndef BOOST_MATH_COMPLEX_HPP
#define BOOST_MATH_COMPLEX_HPP

#include <complex>

#include "boost/math/complex_it.hpp"
#include "boost/math/complex_rt.hpp"
#include <boost/math/octonion.hpp>
#include <boost/math/quaternion.hpp>


// Flag to check if complex-number implmentations are packed
#if 0
#define BOOST_MATH_COMPLEX_IS_PACKED  1
#else
#define BOOST_MATH_COMPLEX_IS_PACKED  0
#endif
/** \def  BOOST_MATH_COMPLEX_IS_PACKED
    \brief  Flag for packed implementations on complex numbers.

    If this pre-processor flag is set to non-zero, then the implmentations of
    #boost::math::complex_rt and #boost::math::complex_it don't have any padding
    at the beginning, middle, or end of their structures.  This means that:

    - `sizeof( complex_?t<T, R> ) == ( sizeof(T) << R )`
    - Given a `complex_?t<T, R>` object `x`, you can find a given component *i*
      with `reintepret_cast< T(&)[1 << R] >( x )[ i ]`.  (This property is always
      true with `complex_it<T, ?>` when `T` is standard-layout.)
    - Given an array segment of `complex_?t<T, R>` objects starting at *p*, you
      can find component *i* of element *n* with `reinterpret_cast<T *>( p )[ (n
      << R) + i ]`.
 */


//! Name-space for all Boost (non-macro) items
namespace boost
{
//! Name-space for various mathematical functions, macros, and types in Boost.
namespace math
{


// Put stuff here.


}  // namespace math
}  // namespace boost


#endif // BOOST_MATH_COMPLEX_HPP
