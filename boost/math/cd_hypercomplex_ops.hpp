//  Boost.Complex Library, cd_hypercomplex_ops.hpp header file  --------------//

//  Copyright 2012 Daryle Walker.
//  Distributed under the Boost Software License, Version 1.0.  (See the
//  accompanying file LICENSE_1_0.txt or a copy at
//  <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

/** \file
    \brief  Core header for modeling Cayley-Dickson hypercomplex numbers.

    \author  Daryle Walker

    \version  0.5

    \copyright  Boost Software License, version 1.0

    Contains the standard operations, including operators, for the class
    templates modelling Cayley-Dickson hypercomplex numbers.  The algorithms
    used slightly differ depending if a given class stores its components
    iteratively or recursively.
 */

#ifndef BOOST_MATH_CD_HYPERCOMPLEX_OPS_HPP
#define BOOST_MATH_CD_HYPERCOMPLEX_OPS_HPP

#include "boost/math/cd_hypercomplex_core.hpp"

#include <algorithm>
#include <cstddef>
#include <utility>


namespace boost
{
namespace math
{


//  Hypercomplex number equality operator definitions  -----------------------//

/** \brief Equality comparison

    Compares two hypercomplex values for equality.  Two values are equal when
    every pair of corresponding components are equal.  When the values differ
    in length, the longer value's trailing components are compared against
    zero, acting as if the shorter value was zero-extended.

    The result can be expressed as:
    - Single: `R == Q`
    - Iterative: `(c0 == d0) && (c1 == d1) && ... && (cN == dN)`
    - Recursive: `(B0 == A0) && (B1 == A1)`

    \pre  `operator ==` must be supported between objects of types `T` and `U`.
    \pre  `T` and `U` must support Boolean conversion.  (It may be `explicit`.)

    \param[in] l  The left operand.
    \param[in] r  The right operand.

    \throws  Anything equality-comparison may throw.  If `T::dimensions` and
             `U::dimensions` differ; then any anything that `operator bool` of
             the longer value may throw too.

    \returns  <code>!(<var>l</var> - <var>r</var>)</code>.
 */
template < typename T, std::size_t R, typename U, std::size_t S >
auto operator ==( cdh_complex_ai<T, R> const &l, cdh_complex_ai<U, S> const &r )
 -> decltype( std::declval<T const &>() == std::declval<U const &>() )
{
    using std::find_if;

    auto const  common_length = std::min( l.size(), r.size() );

    return std::equal( l.begin(), l.begin() + common_length, r.begin() )
     && ( find_if(l.begin() + common_length, l.end(), [] ( T const &ll )
        { return static_cast<bool>(ll); }) == l.end() )
     && ( find_if(r.begin() + common_length, r.end(), [] ( U const &rr )
        { return static_cast<bool>(rr); }) == r.end() );
}

//! \overload
template < typename T, typename U >
inline constexpr
auto operator ==( cdh_complex_ar<T, 0> const &l, cdh_complex_ar<U, 0> const &r )
 -> decltype( std::declval<T const &>() == std::declval<U const &>() )
{ return l.r[0] == r.r[0]; }

//! \overload
template < typename T, typename U, std::size_t R >
inline constexpr
auto operator ==( cdh_complex_ar<T, R> const &l, cdh_complex_ar<U, R> const &r )
 -> decltype( std::declval<T const &>() == std::declval<U const &>() )
{ return (l.b[ 0 ] == r.b[ 0 ]) && (l.b[ 1 ] == r.b[ 1 ]); }

//! \overload
template < typename T, std::size_t R, typename U, std::size_t S,
 typename std::enable_if<(R < S)>::type... >
inline constexpr
auto operator ==( cdh_complex_ar<T, R> const &l, cdh_complex_ar<U, S> const &r )
 -> decltype( std::declval<T const &>() == std::declval<U const &>() )
{ return (l == r.b[ 0 ]) && !r.b[1]; }

//! \overload
template < typename T, std::size_t R, typename U, std::size_t S,
 typename std::enable_if<(R > S)>::type... >
inline constexpr
auto operator ==( cdh_complex_ar<T, R> const &l, cdh_complex_ar<U, S> const &r )
 -> decltype( std::declval<T const &>() == std::declval<U const &>() )
{ return (l.b[ 0 ] == r) && !l.b[1]; }

/** \brief Inequality comparison

    Compares two hypercomplex values for inequality.  Two values are unequal
    when at least one pair of corresponding components is unequal.  When the
    values differ in length, the longer value's trailing components are compared
    against zero, acting as if the shorter value was zero-extended.

    The result can be expressed as:
    - Single: `R != Q`
    - Iterative: `(c0 != d0) || (c1 != d1) || ... || (cN != dN)`
    - Recursive: `(B0 != A0) || (B1 != A1)`

    \pre  `operator ==` (and `operator !=`) must be supported between objects of
          types `T` and `U`.
    \pre  `T` and `U` must support Boolean conversion.  (It may be `explicit`.)

    \param[in] l  The left operand.
    \param[in] r  The right operand.

    \throws  Anything (in)equality-comparison may throw.  If `T::dimensions` and
             `U::dimensions` differ; then any anything that `operator bool` of
             the longer value may throw too.

    \returns  <code>!(<var>l</var> == <var>r</var>)</code>.
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline
auto operator !=( cdh_complex_ai<T, R> const &l, cdh_complex_ai<U, S> const &r )
 -> decltype( std::declval<T const &>() == std::declval<U const &>() )
{ return !(l == r); }

//! \overload
template < typename T, typename U >
inline constexpr
auto operator !=( cdh_complex_ar<T, 0> const &l, cdh_complex_ar<U, 0> const &r )
 -> decltype( std::declval<T const &>() != std::declval<U const &>() )
{ return l.r[0] != r.r[0]; }

//! \overload
template < typename T, typename U, std::size_t R >
inline constexpr
auto operator !=( cdh_complex_ar<T, R> const &l, cdh_complex_ar<U, R> const &r )
 -> decltype( std::declval<T const &>() != std::declval<U const &>() )
{ return (l.b[ 0 ] != r.b[ 0 ]) || (l.b[ 1 ] != r.b[ 1 ]); }

//! \overload
template < typename T, std::size_t R, typename U, std::size_t S,
 typename std::enable_if<(R < S)>::type... >
inline constexpr
auto operator !=( cdh_complex_ar<T, R> const &l, cdh_complex_ar<U, S> const &r )
 -> decltype( std::declval<T const &>() != std::declval<U const &>() )
{ return (l != r.b[ 0 ]) || r.b[1]; }

//! \overload
template < typename T, std::size_t R, typename U, std::size_t S,
 typename std::enable_if<(R > S)>::type... >
inline constexpr
auto operator !=( cdh_complex_ar<T, R> const &l, cdh_complex_ar<U, S> const &r )
 -> decltype( std::declval<T const &>() != std::declval<U const &>() )
{ return (l.b[ 0 ] != r) || l.b[1]; }


}  // namespace math
}  // namespace boost


#endif  // BOOST_MATH_CD_HYPERCOMPLEX_OPS_HPP
