//  Implementation Detail, index_tuple11.hpp header file  --------------------//

//  Copyright 2012 Daryle Walker.
//  Distributed under the Boost Software License, Version 1.0.  (See the
//  accompanying file LICENSE_1_0.txt or a copy at
//  <http://www.boost.org/LICENSE_1_0.txt>.)

//  There's no corresponding Boost library (yet); it's an implementation detail.

/** \file
    \brief  Common header for the index list tuple technique.

    \author  Daryle Walker
    \author  Several C++ programmers on the web that I based this code on.
             (I need suggestions on who to thank here.)

    \version  0.1

    \copyright  Boost Software License, version 1.0

    Contains the declarations and definitions of the class and function
    templates needed for the list-of-numbers variadic parameter technique.

    This construct needs C++2011 enabled for your compiler.
 */

#ifndef BOOST_DETAIL_INDEX_TUPLE11_HPP
#define BOOST_DETAIL_INDEX_TUPLE11_HPP

#include <cstddef>


namespace boost
{
namespace detail
{


/** \brief  Meta-programming provider of orderly numeric variadic list.

    You provide a properly-build object of an instantiation of this class
    template as a parameter to a (member) function template, and said function
    has a formal template parameter that is a variadic list of `size_t`, then
    function has a list of ordered numbers available for the template-based
    expansions in its implementation.

    As the numeric variadic list will almost always be "0, 1, ..., X", you
    should use following meta-function class templates and/or function template
    to create objects from this template.

    \tparam Indices  The list of numbers (type `std::size_t`) that a function
                     template will use in its definition after receiving an
                     object of an instantiation of this class template.
 */
template < std::size_t ...Indices >
struct index_tuple
{
    //! Meta-programming helper for building other `index_tuple` instantiations.
    typedef index_tuple< Indices..., sizeof...(Indices) >  next;
};

/** \brief  Meta-programming function to create conventional `index_tuple`
            instantiations.

    Provides a type-alias that creates an instantiation of `index_tuple` where
    its numeric variadic list is "0, 1, ..., `Size` - 1".  (If `Size` is 1, then
    the list is just "0".)  It's recursively defined, so there is a
    specialization for the base case.

    \tparam Size  The number of elements in the template parameter list of the
                  `index_tuple` instantiation created here.
 */
template < std::size_t Size >
struct build_indices
{
    //! A type-alias to `index_tuple\<0, 1, ..., Size - 1\>`.
    typedef typename build_indices<Size - 1u>::type::next  type;
};

/** \brief  Base case for `build_indices`.

    See the notes from the main version for more.  When `Size` is 0, then the
    list is empty.
 */
template < >
struct build_indices< 0u >
{
    //! A type-alias to the empty `index_tuple` instantiation.
    typedef index_tuple<>  type;
};

/** \brief  Type-based variant of `build_indices`.

    When the list length of `index_tuple` has to match a list of types, then
    use this variant.

    \tparam T  The list of types to be enumerated.
 */
template < typename ...T >
using build_indices_from_types = build_indices< sizeof...(T) >;

/** \brief  Creates an `index_tuple` object for later (meta-)programming.

    Given a numeric limit, create an object of an `index_tuple` instantiation
    with the numbers 0 through the limit (including 0, but not the limit) in
    order as the actual template parameters of the type.  Passing this object to
    a (member) function template lets said function use the numeric list in its
    implementation.

    \tparam Size  The number of elements in the template parameter list of the
                  `index_tuple` instantiation created here.

    \throws  No-throw.

    \returns  An object of type `index_tuple\<0, ..., Size - 1\>`.  The object
              is empty of data.
 */
template < std::size_t Size >
inline constexpr
auto  make_indices() noexcept -> typename build_indices<Size>::type
{ return {}; }


}  // namespace detail
}  // namespace boost


#endif  // BOOST_DETAIL_INDEX_TUPLE11_HPP
