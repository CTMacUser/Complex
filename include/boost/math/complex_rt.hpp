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
#include <utility>

#include "boost/math/complex_it.hpp"


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

    // Constructors
    /** \brief  Default-construction

    The compiler-generated definition is used to keep `complex_rt` trivial (when
    #value_type is trivial).

        \post  When used for value-initialization, `(*this)[0] == value_type{}`.
        \post  When used for default-initialization, the component is left
               uninitialized if #value_type uses trivial default-construction,
               otherwise it acts the same as value-initialization.
     */
    complex_rt() = default;
    /** \brief  Single-real constructor

    Constructs a `complex_rt` object from a single real number.  Can act as a
    conversion.

        \param[in] r  The real number to convert.

        \post  `(*this)[0] == r`.
     */
    constexpr  complex_rt( value_type const &r )  : r{ r }  {}

    /** \brief  Convert from a `complex_rt`, same #rank, different #value_type

    Constructs a `complex_rt` object from one with a different #value_type, but
    the same #rank.  Acts as a conversion.

        \pre  The component type for `s` has to be implicitly convertible to
              #value_type.

        \param[in] s  The real number copy source.

        \post  `(*this)[0] == s[0]`.

     */
    template <
        typename  T,
        typename    = typename std::enable_if<std::is_convertible<T,
         value_type>::value>::type
    >
    constexpr  complex_rt( complex_rt<T, rank> const &s )  : r{ s[0] }  {}
    /** \brief  Convert from a longer `complex_rt` object.

    Constructs a `complex_rt` from the first (i.e. real) component of an
    object with a higher #rank.  The conversion is `explicit` to prevent
    ambiguities with the same-length conversion constructor.

        \pre  The component type for `senior` has to implicitly convertible to
              #value_type.

        \param[in] senior  The source to copy.

        \post  `(*this)[0] == senior[0]`.
     */
    template <
        typename   T,
        size_type  R,
        typename     = typename std::enable_if<R>::type,
        typename     = typename std::enable_if<std::is_convertible<T,
         value_type>::value>::type
    >
    explicit constexpr  complex_rt( complex_rt<T, R> const &senior )
        : r{ senior[0] }
    {}
    /** \brief  Convert from a `complex_rt` object with explicitly-convertible
                component types.

    Constructs a `complex_rt` from the components of another object, but the
    respective component types are *not* implicitly convertible.  Obviously,
    the composite-level constructor is also `explicit`.

        \pre  The component type for `cc` has to be explicitly convertible to
              #value_type.

        \param[in] cc  The source to copy.

        \post  `(*this)[0] == static_cast<value_type>(cc[0])`.
     */
    template <
        typename   T,
        size_type  R,
        typename     = typename std::enable_if<not std::is_convertible<T,
         value_type>::value>::type
    >
    explicit  complex_rt( complex_rt<T, R> const &cc )
        : r{ static_cast<value_type>(cc[ 0 ]) }
    {}

    /** \brief  Convert from a `complex_it`.

    Constructs a `complex_rt` object with one from the other philosophy.  To
    prevent ambiguities, conversion is marked `explicit`.  The source object can
    be of any size; excess components will be truncated.

        \param[in] c  The `complex_it` object to convert.

        \post  `(*this)[0] == static_cast<value_type>(c[0])`.
     */
    template < typename T, size_type R >
    explicit constexpr  complex_rt( complex_it<T, R> const &c )
        : r{ static_cast<value_type>(c[ 0 ]) }
    {}

    // More operators
    /** \brief  Convert to a `complex_it`.

    Creates a `complex_it` object with one from this philosophy.  To prevent
    ambiguities, conversion is marked `explicit`.  The destination object can be
    of any size; excess components will be value-initialized.

        \returns  A value `x` such that (given destination component type `T`):
                  - `x[k] == static_cast<T>( (*this)[k] )`.
                    for all indices `k` that are valid for both `x` and `*this`.
                  - `x[k] == T{}` for all valid indices `k >=` #static_size.
     */
    template < typename T, size_type R >
    explicit constexpr  operator complex_it<T, R>() const
    { return complex_it<T, R>{static_cast<T>( r )}; }

private:
    // Authorize use of hidden constructor below for complex_rt<Number, 1u>
    friend class complex_rt<value_type, rank + 1u>;

    // Hidden constructor to take a specific component from a complex_it object.
    template < size_type I, typename T, size_type R >
    explicit constexpr
    complex_rt( std::integral_constant<size_type, I>, complex_it<T, R> const
     &c )
        : r{ (I < complex_it<T, R>::static_size) ? static_cast<value_type>(c[I])
          : value_type{} }
    {}

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

    // Constructors
    /** \brief  Default-construction

    The compiler-generated definition is used to keep `complex_rt` trivial (when
    #value_type is trivial).

        \post  When used for value-initialization, `(*this)[i] == value_type{}`
               for every valid `i`.
        \post  When used for default-initialization, each component is left
               uninitialized if #value_type uses trivial default-construction,
               otherwise it acts the same as value-initialization.
     */
    complex_rt() = default;
    /** \brief  Single-real conversion

    Constructs a `complex_rt` object from a single real number.  Can act as a
    conversion.

        \param[in] r  The real number to convert.

        \post  `(*this)[0] == r`.
        \post  For `0 < i <` #static_size, `(*this)[i] == value_type{}`.
     */
    constexpr  complex_rt( value_type const &r )  : b{ {r} }  {}
    /** \brief  List-of-reals constructor

    Constructs a `complex_rt` object from a list of real numbers.

        \pre  `sizeof...(u)` \<= #static_size - 2.
        \pre  Each parameter in `u` must have a (non-narrowing) implicit
              conversion to #value_type.

        \param[in] r  The real component.
        \param[in] i  The classic imaginary component.
        \param[in] u  The higher-order imaginary components.  May be empty.

        \post  `(*this)[0] == r`.
        \post  `(*this)[1] == i`.
        \post  For `1 < k <= sizeof...(u)+1`, `(*this)[k] == Explode(u; k-2)`.
        \post  For `sizeof...(u) + 1 < k <` #static_size, `(*this)[k] ==
               value_type{}`.
     */
    template < typename ...Args >
    constexpr  complex_rt( value_type const &r, value_type const &i, Args const
     &...u )
        : complex_rt{ complex_rt<value_type, 0u>{r}, complex_rt<value_type,
          0u>{i}, complex_rt<Args, 0u>{u}... }
    {}

    /** \brief  Convert from a `complex_rt`, same #rank, different #value_type

    Constructs a `complex_rt` object from one with a different #value_type, but
    the same #rank.  Acts as a conversion.

        \pre  The component type for `s` has to be implicitly convertible to
              #value_type.

        \param[in] s  The complex number copy source.

        \post  `(*this)[k] == s[k]` for all valid indices `k`.
     */
    template <
        typename  T,
        typename    = typename std::enable_if<std::is_convertible<T,
         value_type>::value>::type
    >
    constexpr  complex_rt( complex_rt<T, rank> const &s )
        : b{ s.lower_barrage(), s.upper_barrage() }
    {}
    /** \brief  Convert from one or two barrages.

    Constructs a `complex_rt` from one or two barrages.  It is not required that
    one or both barrages share #value_type.  Can act as a conversion.

        \pre  The component types for `l` and `u` each have to be implicitly
              convertible to #value_type.

        \param[in] l  The complex number copy source for the lower barrage.
        \param[in] u  The complex number copy source for the upper barrage.  If
                      not given, it defaults to `barrage_type{}`.

        \post  `this->lower_barrage() == l`.
        \post  `this->upper_barrage() == u`.
     */
    template <
        typename  T,
        typename  U = value_type,
        typename    = typename std::enable_if<std::is_convertible<T,
         value_type>::value>::type,
        typename    = typename std::enable_if<std::is_convertible<U,
         value_type>::value>::type
    >
    constexpr  complex_rt( complex_rt<T, rank - 1u> const &l,
     complex_rt<U, rank - 1u> const &u = {} )  : b{ l, u }  {}
    /** \brief  Convert from list of sub-barrages.

    Constructs a `complex_rt` from one or more `complex_rt` objects that are
    smaller than `barrage_type`.  It is not required that the sub-barrages share
    #value_type, nor use the same component type among themselves.  Can act as a
    conversion.

        \pre  The component types for `first` and each `rest` all have to be
              implicitly convertible to #value_type.
        \pre  The total number of components from `first` and each `rest`
              combined cannot exceed #static_size.

        \param[in] first  The first set of components to copy.
        \param[in] rest   The remaining sets of components to copy.  May be
                          empty.

        \post  Given `auto const x = std::tuple_cat(first, rest...);` and `using
               TT = decltype(x);`:
               - For `0 <= k < Min(static_size, std::tuple_size<TT>::value)`,
                 `get<k>(*this) == get<k>(x)`.
               - For `std::tuple_size<TT>::value <= k < static_size`,
                 `get<k>(*this) == value_type{}`.
     */
    template <
        typename    T,
        size_type   R,
        typename ...U,
        typename      = typename std::enable_if<(R < rank - 1u)>::type,
        typename      = typename std::enable_if<((1u + sizeof...( U )) <= (1ULL
         << ( rank - R )))>::type,
        typename      = typename std::enable_if<std::is_convertible<T,
         value_type>::value>::type
    >
    constexpr
    complex_rt( complex_rt<T, R> const &first, complex_rt<U, R> const &...rest )
        : complex_rt{ std::integral_constant<size_type, 0u>{}, first, rest... }
    {}
    /** \brief  Convert from a longer `complex_rt` object.

    Constructs a `complex_rt` from the first (i.e. real-ward) components of an
    object with a higher #rank.  The conversion is `explicit` to prevent
    ambiguities with the same-length and (sub-)barrage conversion constructors.

        \pre  The component type for `senior` has to implicitly convertible to
              #value_type.

        \param[in] senior  The source to copy.

        \post  For `0 <= k <` #static_size, `(*this)[k] == senior[k]`.
     */
    template <
        typename   T,
        size_type  R,
        typename     = typename std::enable_if<(R > rank)>::type,
        typename     = typename std::enable_if<std::is_convertible<T,
         value_type>::value>::type
    >
    explicit constexpr  complex_rt( complex_rt<T, R> const &senior )
        : complex_rt{ senior.lower_barrage() }
    {}
    /** \brief  Convert from `complex_rt` objects with explicitly-convertible
                component types.

    Constructs a `complex_rt` from the components of another object, but the
    respective component types are *not* implicitly convertible.  Obviously,
    the composite-level constructor is also `explicit`.

        \pre  The component type for `cc` has to be explicitly convertible to
              #value_type.

        \param[in] cc  The source to copy.

        \post  For `0 <= k < Min(static_size, decltype(cc)::static_size)`,
               `(*this)[k] == static_cast<value_type>(cc[k])`.
        \post  For `decltype(cc)::static_size <= k <` #static_size,
               `(*this)[k] == value_type{}`.
     */
    template <
        typename   T,
        size_type  R,
        typename     = typename std::enable_if<not std::is_convertible<T,
         value_type>::value>::type
    >
    explicit  complex_rt( complex_rt<T, R> const &cc )
        : b{}
    {
        auto  i = 1ULL << std::min( rank, R );

        while ( i-- )
            ( *this )[ i ] = static_cast<value_type>( cc[i] );
    }

    /** \brief  Convert from a `complex_it`.

    Constructs a `complex_rt` object with one from the other philosophy.  To
    prevent ambiguities, conversion is marked `explicit`.  The source object can
    be of any size; either excess components will be truncated or missing
    components will be value-initialized.

        \pre  There is some sort of conversion from `decltype(c)::value_type` to
              #value_type.

        \param[in] c  The `complex_it` object to convert.

        \post  For `0 <= k < Min( decltype(c)::static_size,` #static_size `)`,
               `(*this)[k] == static_cast<value_type>(c[k])`.
        \post  For `decltype(c)::static_size <= k <` #static_size,
               `(*this)[k] == value_type{}`.
     */
    template < typename T, size_type R >
    explicit constexpr  complex_rt( complex_it<T, R> const &c )
        : complex_rt{ std::integral_constant<size_type, 0u>{}, c }
    {}

    // More operators
    /** \brief  Convert to a `complex_it`.

    Creates a `complex_it` object with one from this philosophy.  To prevent
    ambiguities, conversion is marked `explicit`.  The destination object can be
    of any size; excess components will be value-initialized.

        \returns  A value `x` such that:
                  - `x[k] == static_cast<decltype(x)::value_type>((*this)[k])`
                    for all indices `k` that are valid for both `x` and `*this`.
                  - `x[k] == decltype(x)::value_type{}` for `k >=` #static_size.
     */
    template < typename T, size_type R >
    explicit  operator complex_it<T, R>() const
    {
        complex_it<T, R>  result{};
        auto              b = begin( result );
        auto const        e = end( result );
        size_type         i = 0u;

        while ( (static_size > i) && (e != b) )
            *b++ = static_cast<T>( (*this)[i++] );
        return result;
    }

private:
    // Authorize use of hidden constructors below for complex_rt<Number,Rank+1u>
    friend class complex_rt<value_type, rank + 1u>;

    // Hidden constructor to take specific components from a complex_it object.
    template < size_type Start, typename T, size_type R >
    explicit constexpr
    complex_rt( std::integral_constant<size_type, Start>, complex_it<T, R> const
     &c )
        : b{ barrage_type{std::integral_constant<size_type, Start>{}, c},
          barrage_type{std::integral_constant<size_type, Start + static_size /
          2u>{}, c} }
    {}

    // Hidden constructor to value-initialize when sources run dry early.
    template < size_type Start >
    explicit constexpr
    complex_rt( std::integral_constant<size_type, Start> )  : b{}  {}

    // Hidden constructor to take the first barrage coming.
    template < typename T >
    explicit //constexpr
    complex_rt( std::integral_constant<size_type, 0u>, complex_rt<T, rank - 1u>
     const &first )  : b{ first }  {}

    // Hidden constructor to take the first two barrages coming.
    template < typename T, typename U, typename ...V >
    explicit //constexpr
    complex_rt( std::integral_constant<size_type, 0u>,
     complex_rt<T, rank - 1u> const &first,
     complex_rt<U, rank - 1u> const &second,
     complex_rt<V, rank - 1u> const &... )
        : b{ first, second }
    {}

    // Hidden constructor to take the later barrages.
    template < size_type Start, typename T, typename ...U >
    explicit constexpr
    complex_rt( std::integral_constant<size_type, Start>,
     complex_rt<T, rank - 1u> const & ,
     complex_rt<U, rank - 1u> const &...rest )
        : complex_rt{ std::integral_constant<size_type, Start - static_size /
          2u>{}, rest... }
    {}

    // Hidden constructor to handle sub-barrages.
    template <
        size_type    Start,
        typename     T,
        size_type    R,
        typename  ...U,
        typename       = typename std::enable_if<(R < rank - 1u)>::type
    >
    explicit constexpr
    complex_rt( std::integral_constant<size_type, Start> marker, complex_rt<T,
     R> const &first, complex_rt<U, R> const &...rest )
        : b{ barrage_type{marker, first, rest...}, barrage_type{
          std::integral_constant<size_type, Start + static_size / 2u>{}, first,
          rest...} }
    {}

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


//  Object support functions  ------------------------------------------------//

/** \brief  Swap routine for `complex_rt`.

Exchanges the state of the two given objects.

    \relatesalso  #boost::math::complex_rt

    \pre  There is a swapping routine, called `swap`, either in namespace `std`
          for built-ins, or found via ADL for other types.

    \param a  The first object to have its state exchanged.
    \param b  The second object to have its state exchanged.

    \throws Whatever  the element-level swap does.

    \post  `a` is equivalent to the old state of `b`, while `b` is equivalent to
           the old state of `a`.
 */
template < typename T >
inline
void  swap( complex_rt<T, 0u> &a, complex_rt<T, 0u> &b )
 noexcept( detail::is_swap_nothrow<T>() )
{
    using std::swap;

    swap( a[0], b[0] );
}

/** \overload
    \relatesalso  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline
void  swap( complex_rt<T, R> &a, complex_rt<T, R> &b )
 noexcept( detail::is_swap_nothrow<T>() )
{
    swap( a.lower_barrage(), b.lower_barrage() );
    swap( a.upper_barrage(), b.upper_barrage() );
}


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
                   at least one non-zero value among the excess components.
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
