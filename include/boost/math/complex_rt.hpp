//  Boost Complex Numbers, recursive-mode header file  -----------------------//

//  Copyright 2013 Daryle Walker.
//  Distributed under the Boost Software License, Version 1.0.  (See the
//  accompanying file LICENSE_1_0.txt or a copy at
//  <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

/** \file   boost/math/complex_rt.hpp
    \brief  A complex-number class template with recursive composition.

    \author  Daryle Walker

    \version  0.1

    \copyright  Boost Software License, version 1.0

    Contains the declaration (and definition) of class template `complex_rt`,
    that models Cayley-Dickson hypercomplex numbers (complex, quaternions,
    octonions, etc.) by storing the components and synthesizing the related
    operations recursively.

    \warning  This library requires C++2011 features.
 */

#ifndef BOOST_MATH_COMPLEX_RT_HPP
#define BOOST_MATH_COMPLEX_RT_HPP

#include <cstddef>
#include <ostream>
#include <sstream>
#include <tuple>
#include <type_traits>

// Put includes from Boost here.


namespace boost
{
namespace math
{


//  Hypercomplex number (recursive) class template definitions  --------------//

// Forward declaration since the base case is the (partial) specialization.
template < typename Number, std::size_t Rank >
struct complex_rt;

/** \brief  A complex number class nesting its components, base case.

This class template models hypercomplex numbers built with Cayley-Dickson
construction, with the storage and operations defined recursively.  This means
that the storage for the base case is a single component and the algorithms for
the support operations are trivial.

(The `Rank` is fixed as 0 in the base case.)

See the description for #boost::math::complex_it for more.

    \pre  Same as #boost::math::complex_it

    \tparam Number  The component type
 */
template < typename Number >
struct complex_rt< Number, 0u >
{
    // Core types
    //! \copydoc  #boost::math::complex_it::size_type
    typedef std::size_t  size_type;
    //! \copydoc  #boost::math::complex_it::value_type
    typedef Number      value_type;

    // Sizing parameters
    //! \copybrief  #boost::math::complex_it::rank
    static constexpr  size_type  rank = 0u;
    //! \copybrief  #boost::math::complex_it::static_size
    static constexpr  size_type  static_size = 1u;

    // Support types
    //! A type-alias for compatibility with other instantiations.
    typedef complex_rt  barrage_type;

    // Component(s) access
    //! \copydoc  #boost::math::complex_it::operator[](size_type)const
    constexpr
    auto  operator []( size_type ) const noexcept -> value_type const &
    { return r; }
    //! \overload
    auto  operator []( size_type ) noexcept -> value_type &  { return r; }

    //! \copydoc  #boost::math::complex_it::lower_barrage
    constexpr
    auto  lower_barrage() const noexcept -> barrage_type const &
    { return *this; }
    //! \overload
    auto  lower_barrage() noexcept -> barrage_type &  { return *this; }
    //! \copydoc  #boost::math::complex_it::upper_barrage
    constexpr
    auto  upper_barrage() const noexcept -> barrage_type const &
    { return *this; }
    //! \overload
    auto  upper_barrage() noexcept -> barrage_type &  { return *this; }

    // Conditions
    //! \copydoc  #boost::math::complex_it::operator bool()const
    explicit constexpr  operator bool() const  { return static_cast<bool>(r); }

private:
    // Member data
    value_type  r;
};

/** \brief  A complex number class nesting its components, recursive case.

This class template models hypercomplex numbers built with Cayley-Dickson
construction, with the storage and operations defined recursively.  This means
that the storage for the general case is a pair of `complex_rt` objects of
one-rank lower and the algorithms for the support operations are recursive.

See the description for #boost::math::complex_it for more.

    \pre  Same as #boost::math::complex_it

    \tparam Number  The component type
    \tparam Rank    The Cayley-Dickson construction level
 */
template < typename Number, std::size_t Rank >
struct complex_rt
{
    // Core types
    //! \copydoc  #boost::math::complex_it::size_type
    typedef std::size_t  size_type;
    //! \copydoc  #boost::math::complex_it::value_type
    typedef Number      value_type;

    // Sizing parameters
    //! \copybrief  #boost::math::complex_it::rank
    static constexpr  size_type  rank = Rank;
    //! \copybrief  #boost::math::complex_it::static_size
    static constexpr  size_type  static_size = 1ULL << rank;

    // Support types
    //! The type immediately lower in Cayley-Dickson construction.
    typedef complex_rt<value_type, rank - 1u>  barrage_type;

    // Component(s) access
    //! \copydoc  #boost::math::complex_it::operator[](size_type)const
    constexpr
    auto  operator []( size_type i ) const noexcept -> value_type const &
    { return (i >= static_size / 2u) ? b[1][i - static_size / 2u] : b[0][i]; }
    //! \overload
    auto  operator []( size_type i ) noexcept -> value_type &
    { return (i >= static_size / 2u) ? b[1][i - static_size / 2u] : b[0][i]; }

    //! \copydoc  #boost::math::complex_it::lower_barrage
    constexpr
    auto  lower_barrage() const noexcept -> barrage_type const &
    { return b[0]; }
    //! \overload
    auto  lower_barrage() noexcept -> barrage_type &  { return b[0]; }
    //! \copydoc  #boost::math::complex_it::upper_barrage
    constexpr
    auto  upper_barrage() const noexcept -> barrage_type const &
    { return b[1]; }
    //! \overload
    auto  upper_barrage() noexcept -> barrage_type &  { return b[1]; }

    // Conditions
    //! \copydoc  #boost::math::complex_it::operator bool()const
    explicit constexpr  operator bool() const  { return b[0] || b[1]; }

private:
    // Member data
    barrage_type  b[ 2 ];
};


//  Class-static data member definitions  ------------------------------------//

/** Gives access to an implied template parameter.  The base level of 0
    represents real numbers.
 */
template < typename Number >
constexpr
typename complex_rt<Number, 0u>::size_type  complex_rt<Number, 0u>::rank;

/** There is only one component at the base level.
 */
template < typename Number >
constexpr
typename complex_rt<Number, 0u>::size_type  complex_rt<Number, 0u>::static_size;

/** Gives access to a template parameter.  Past the base level, 1 is the
    (traditional) complex number level, 2 for quaternions, 3 for octonions, etc.
 */
template < typename Number, std::size_t Rank >
constexpr
typename complex_rt<Number, Rank>::size_type  complex_rt<Number, Rank>::rank;

/** The component count doubles when going to the next rank.
 */
template < typename Number, std::size_t Rank >
constexpr
typename complex_rt<Number, Rank>::size_type
  complex_rt<Number, Rank>::static_size;


//  Equality operators  ------------------------------------------------------//

/** \brief  Equality comparison for `complex_rt`.

Compares two `complex_rt` objects for equality.  If the two objects differ in
length, the upper barrage(s) of the longer object will be checked against the
additive identity.

    \relates  #boost::math::complex_rt

    \pre  `std::declval<T>() == std::declval<U>()` is well-formed and the result
          type is (contextually) convertible to `bool`.

    \param l  The left-side argument.
    \param r  The right-side argument.

    \retval true   If the corresponding barrages between *l* and *r* compare
                   as equal; and when their lengths don't match, all the excess
                   components of the longer object are zero.
    \retval false  Otherwise.
 */
template < typename T, typename U >
inline constexpr
bool  operator ==( complex_rt<T, 0u> const &l, complex_rt<U, 0u> const &r )
{ return static_cast<bool>(l[ 0 ] == r[ 0 ]); }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, typename U, std::size_t R >
inline constexpr
bool  operator ==( complex_rt<T, R> const &l, complex_rt<U, R> const &r )
{
    return ( l.lower_barrage() == r.lower_barrage() ) && ( l.upper_barrage() ==
     r.upper_barrage() );
}

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline constexpr
auto  operator ==( complex_rt<T, R> const &l, complex_rt<U, S> const &r )
 -> typename std::enable_if<(R < S), bool>::type
{ return (l == r.lower_barrage()) && !r.upper_barrage(); }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline constexpr
auto  operator ==( complex_rt<T, R> const &l, complex_rt<U, S> const &r )
 -> typename std::enable_if<(R > S), bool>::type
{ return (l.lower_barrage() == r) && !l.upper_barrage(); }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T >
inline constexpr
bool  operator ==( complex_rt<T, 0u> const &l, T const &r )
{ return l[0] == r; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T >
inline constexpr
bool  operator ==( T const &l, complex_rt<T, 0u> const &r )
{ return l == r[0]; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline constexpr
bool  operator ==( complex_rt<T, R> const &l, T const &r )
{ return (l.lower_barrage() == r) && !l.upper_barrage(); }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline constexpr
bool  operator ==( T const &l, complex_rt<T, R> const &r )
{ return (l == r.lower_barrage()) && !r.upper_barrage(); }

/** \brief  Inequality comparison for `complex_rt`.

Compares two `complex_rt` objects for inequality.  If the two objects differ in
length, the upper barrage(s) of the longer object will be checked against the
additive identity.

    \relates  #boost::math::complex_rt

    \pre  `std::declval<T>() == std::declval<U>()` is well-formed and the result
          type is (contextually) convertible to `bool`.

    \param l  The left-side argument.
    \param r  The right-side argument.

    \retval true   If there is a mismatch between corresponding barrages
                   between *l* and *r*; or when their lengths differ, there is
                   at least one non-zero value amoung the excess components.
    \retval false  Otherwise.
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline constexpr
bool  operator !=( complex_rt<T, R> const &l, complex_rt<U, S> const &r )
{ return not operator ==(l, r); }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline constexpr
bool  operator !=( complex_rt<T, R> const &l, T const &r )
{ return not operator ==(l, r); }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline constexpr
bool  operator !=( T const &l, complex_rt<T, R> const &r )
{ return not operator ==(l, r); }


//  Input/output operators  --------------------------------------------------//

/** \brief  Output-streaming for `complex_rt`.

Outputs a string representation of the given number to the given stream.  When
the rank of the number is zero, the sole (real) component is written; otherwise
the lower- and upper-barrages, surrounded by parentheses, are written.

    \relatesalso  #boost::math::complex_rt

    \pre  `std::declval<basic_ostream<Ch,Tr> &>() << std::declval<T>()` is
          well-formed and the result type is the same as the first operand.

    \param[in,out] o  The stream to send the output
    \param[in]     x  The complex number to be written

    \returns  `o`
 */
template < typename Ch, class Tr, typename T >
inline
std::basic_ostream<Ch, Tr> &
operator <<( std::basic_ostream<Ch, Tr> &o, complex_rt<T, 0u> const &x )
{ return o << x[0]; }

/** \overload
    \relatesalso  #boost::math::complex_rt
 */
template < typename Ch, class Tr, typename T, std::size_t R >
std::basic_ostream<Ch, Tr> &
operator <<( std::basic_ostream<Ch, Tr> &o, complex_rt<T, R> const &x )
{
    std::basic_ostringstream<Ch, Tr>  s;

    s.flags( o.flags() );
    s.imbue( o.getloc() );
    s.precision( o.precision() );
    return o << '(' << x.lower_barrage() << ',' << x.upper_barrage() << ')';
}


}  // namespace math
}  // namespace boost


//  Specializations from the STD namespace  ----------------------------------//

//! The standard namespace
namespace std
{
    //! Provides the number of components in a `complex_rt` object.
    template < typename T, size_t R >
    class tuple_size< boost::math::complex_rt<T, R> >
        : public integral_constant< size_t, boost::math::complex_rt<T,
           R>::static_size >
    { };

    //! Provides the type of each component in a `complex_rt` object.
    template < size_t I, typename T, size_t R >
    class tuple_element< I, boost::math::complex_rt<T, R> >
    {
        static_assert( I < boost::math::complex_rt<T, R>::static_size,
         "Index too large" );

    public:
        //! The component type, index-independent for complex numbers.
        typedef typename boost::math::complex_rt<T, R>::value_type  type;
    };

}  // namespace std


#endif // BOOST_MATH_COMPLEX_RT_HPP
