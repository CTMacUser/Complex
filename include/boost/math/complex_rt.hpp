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

#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
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
template < typename Number, std::size_t Rank = 1u >
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
    //! \copybrief  #boost::math::complex_it::has_padding
    static constexpr  bool       has_padding = sizeof( complex_rt ) >
     static_size * sizeof( value_type );

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

    //! \copydoc  #boost::math::complex_it::real()const
    constexpr  auto  real() const -> value_type   { return r; }
    //! \copydoc  #boost::math::complex_it::real(value_type const&)
               void  real( value_type const &r )  { this->r = r; }
    //! \copydoc  #boost::math::complex_it::imag()const
    constexpr  auto  imag() const -> value_type   { return {}; }

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

    //! \copydoc  #boost::math::complex_it::unreal()const
    constexpr  auto  unreal() const -> complex_rt  { return {}; }
    //! \copydoc  #boost::math::complex_it::unreal(complex_it const&)
               void  unreal( complex_rt const & )  {}

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
    \tparam Rank    The Cayley-Dickson construction level.  If not given, it
                    defaults to 1, in order to model regular complex numbers.
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
    //! \copybrief  #boost::math::complex_it::has_padding
    static constexpr  bool       has_padding = sizeof( complex_rt ) >
     static_size * sizeof( value_type );

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

    //! \copydoc  #boost::math::complex_it::real()const
    constexpr  auto  real() const -> value_type   { return b[0].real(); }
    //! \copydoc  #boost::math::complex_it::real(value_type const&)
               void  real( value_type const &r )  { b[0].real(r); }
    //! \copydoc  #boost::math::complex_it::imag()const
    constexpr  auto  imag() const -> value_type   { return operator[](1); }
    //! \copydoc  #boost::math::complex_it::imag(value_type const&)
               void  imag( value_type const &i )  { operator[](1) = i; }

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

    //! \copydoc  #boost::math::complex_it::unreal()const
    constexpr
    auto  unreal() const -> complex_rt
    { return {b[ 0 ].unreal(), b[ 1 ]}; }
    //! \copydoc  #boost::math::complex_it::unreal(complex_it const&)
    void  unreal( complex_rt const &u )
    {
        b[ 0 ].unreal( u.lower_barrage() );
        b[ 1 ] = u.upper_barrage();
    }

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

//! \copydetails  #boost::math::complex_it::has_padding
template < typename Number >
constexpr
bool  complex_rt<Number, 0u>::has_padding;

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

//! \copydetails  #boost::math::complex_it::has_padding
template < typename Number, std::size_t Rank >
constexpr
bool  complex_rt<Number, Rank>::has_padding;


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


//  Tuple interface functions  -----------------------------------------------//

/** \brief  Accesses a component of a `complex_rt`.

Extracts the given component of the given complex number object.

    \pre  0 \<= *I* \< #boost::math::complex_rt::static_size.

    \tparam I  The index of the desired component.

    \param c  The complex number object containing the component.

    \returns  A reference to the desired component.
 */
template < std::size_t I, typename T, std::size_t R >
inline constexpr
auto  get( complex_rt<T, R> const &c ) noexcept -> T const &
{
    static_assert( I < complex_rt<T, R>::static_size, "Index too large" );

    return c[ I ];
}

//! \overload
template < std::size_t I, typename T, std::size_t R >
inline
auto  get( complex_rt<T, R> &c ) noexcept -> T &
{ return const_cast<T &>(get<I>( const_cast<complex_rt<T,R> const &>(c) )); }

//! \overload
template < std::size_t I, typename T, std::size_t R >
inline
auto  get( complex_rt<T, R> &&c ) noexcept -> T &&
{ return std::forward<T>(get<I>( c )); }


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


//  Addition operators  ------------------------------------------------------//

/** \brief  Identity operator

Returns the given value (possibly normalized).

The definition can be broken down as:
- Real: `+r`
- Component-wise: `{ +c[0], +c[1], ..., +c[2^Rank - 1] }`
- Barrage-wise: `{ +Lower; +Upper }`

    \relates  #boost::math::complex_rt

    \pre  `+declval<T>()` is well-formed.

    \param[in] x  The input value.

    \returns  `x`.
 */
template < typename T >
inline constexpr
auto  operator +( complex_rt<T, 0u> const &x )
 -> complex_rt<decltype( +std::declval<T>() ), 0u>
{ return {+x[ 0 ]}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline constexpr
auto  operator +( complex_rt<T, R> const &x )
 -> complex_rt<decltype( +std::declval<T>() ), R>
{ return {+x.lower_barrage(), +x.upper_barrage()}; }

/** \brief  Addition

Calculates the sum of the given values.

The definition can be broken down as:
- Real: `rA + rB`
- Component-wise: `{ cA[0] + cB[0], ..., cA[ 2^Rank - 1 ] + cB[ 2^Rank - 1 ] }`
- Barrage-wise: `{ lowerA + lowerB; upperA + upperB }`

    \relates  #boost::math::complex_rt

    \pre  `declval<T>() + declval<U>()` is well-formed.  Let's call the result
          type `SS`.
    \pre  `+declval<T>()` and `+declval<U>()` are well-formed.
    \pre  Both `T` and `U` can implicitly be assigned to a `SS`.

    \param[in] augend  The first value to be added.
    \param[in] addend  The second value to be added.

    \returns  The sum of `augend` and `addend`.
 */
template < typename T, typename U >
inline constexpr
auto  operator +( complex_rt<T, 0u> const &augend, complex_rt<U, 0u> const
 &addend )
 -> complex_rt<decltype( std::declval<T>() + std::declval<U>() ), 0u>
{ return {augend[ 0 ] + addend[ 0 ]}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, typename U, std::size_t R >
inline constexpr
auto  operator +( complex_rt<T,R> const &augend, complex_rt<U,R> const &addend )
 -> complex_rt<decltype( std::declval<T>() + std::declval<U>() ), R>
{
    return { augend.lower_barrage() + addend.lower_barrage(),
     augend.upper_barrage() + addend.upper_barrage() };
}

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline constexpr
auto  operator +( complex_rt<T,R> const &augend, complex_rt<U,S> const &addend )
 -> typename std::enable_if< (R < S), complex_rt<decltype( std::declval<T>() +
 std::declval<U>() ), S> >::type
{ return {augend + addend.lower_barrage(), +addend.upper_barrage()}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline constexpr
auto  operator +( complex_rt<T,R> const &augend, complex_rt<U,S> const &addend )
 -> typename std::enable_if< (R > S), complex_rt<decltype( std::declval<T>() +
 std::declval<U>() ), R> >::type
{ return {augend.lower_barrage() + addend, +augend.upper_barrage()}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T >
inline constexpr
auto  operator +( T const &augend, complex_rt<T, 0u> const &addend )
 -> complex_rt<decltype( std::declval<T>() + std::declval<T>() ), 0u>
{ return {augend + addend[ 0 ]}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline constexpr
auto  operator +( T const &augend, complex_rt<T, R> const &addend )
 -> complex_rt<decltype( std::declval<T>() + std::declval<T>() ), R>
{ return {augend + addend.lower_barrage(), +addend.upper_barrage()}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T >
inline constexpr
auto  operator +( complex_rt<T, 0u> const &augend, T const &addend )
 -> complex_rt<decltype( std::declval<T>() + std::declval<T>() ), 0u>
{ return {augend[ 0 ] + addend}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline constexpr
auto  operator +( complex_rt<T, R> const &augend, T const &addend )
 -> complex_rt<decltype( std::declval<T>() + std::declval<T>() ), R>
{ return {augend.lower_barrage() + addend, +augend.upper_barrage()}; }

/** \brief  Add-and-assign

Calculates the sum of the given objects into the first.

    \relates  #boost::math::complex_rt

    \pre  `declval<T &>() += declval<U>()` is well-formed.
    \pre  The rank of `addend` doesn't exceed that of `augend_sum`.

    \param[in,out] augend_sum  The first value to be added, and the location of
                               the future sum.
    \param[in]     addend      The second value to be added.

    \returns  A reference to post-addition `augend_sum`.
 */
template < typename T, typename U >
inline
auto  operator +=(complex_rt<T, 0u> &augend_sum,complex_rt<U, 0u> const &addend)
 -> complex_rt<T, 0u> &
{ return augend_sum[0] += addend[0], augend_sum; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, typename U, std::size_t R >
inline
auto  operator +=( complex_rt<T, R> &augend_sum,complex_rt<U, R> const &addend )
 -> complex_rt<T, R> &
{
    augend_sum.lower_barrage() += addend.lower_barrage();
    augend_sum.upper_barrage() += addend.upper_barrage();
    return augend_sum;
}

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline
auto  operator +=( complex_rt<T, R> &augend_sum,complex_rt<U, S> const &addend )
 -> typename std::enable_if< (R > S), complex_rt<T, R> >::type &
{ return augend_sum.lower_barrage() += addend, augend_sum; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T >
inline
auto  operator +=( complex_rt<T, 0u> &augend_sum, T const &addend )
 -> complex_rt<T, 0u> &
{ return augend_sum[0] += addend, augend_sum; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline
auto  operator +=( complex_rt<T, R> &augend_sum, T const &addend )
 -> complex_rt<T, R> &
{ return augend_sum.lower_barrage() += addend, augend_sum; }

/** \brief  Pre-increment

Applies the successor function to the given object, turning it into the sum of
the original value and one.  (This is the same inappropriate application that
the built-in floating types get.)  Only the real part is affected.

    \relates  #boost::math::complex_rt

    \pre  `++declval<T &>()` is well-formed.

    \param[in,out] augend_sum  The object to be affected.

    \returns  A reference to post-increment `augend_sum`.
 */
template < typename T >
inline
auto  operator ++( complex_rt<T, 0u> &augend_sum ) -> complex_rt<T, 0u> &
{ return ++augend_sum[0], augend_sum; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline
auto  operator ++( complex_rt<T, R> &augend_sum ) -> complex_rt<T, R> &
{ return ++augend_sum.lower_barrage(), augend_sum; }

/** \brief  Post-increment

Applies the successor function to the given object, turning it into the sum of
the original value and one.  (This is the same inappropriate application that
the built-in floating types get.)  Only the real part is affected.

    \relates  #boost::math::complex_rt

    \pre  `declval<T &>()++` is well-formed.

    \param[in,out] augend_sum  The object to be affected.

    \returns  A copy of pre-increment `augend_sum`.
 */
template < typename T >
inline
auto  operator ++( complex_rt<T, 0u> &augend_sum, int ) -> complex_rt<T, 0u>
{ return {augend_sum[ 0 ]++}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline
auto  operator ++( complex_rt<T, R> &augend_sum, int ) -> complex_rt<T, R>
{ return {augend_sum.lower_barrage()++, augend_sum.upper_barrage()}; }


//  Subtraction operators  ---------------------------------------------------//

/** \brief  Negation operator

Returns the additive inverse of the given value.

The definition can be broken down as:
- Real: `-r`
- Component-wise: `{ -c[0], -c[1], ..., -c[2^Rank - 1] }`
- Barrage-wise: `{ -Lower; -Upper }`

    \relates  #boost::math::complex_rt

    \pre  `-declval<T>()` is well-formed.

    \param[in] x  The input value.

    \returns  A value `y` such that `x + y == decltype(x){}`.
 */
template < typename T >
inline constexpr
auto  operator -( complex_rt<T, 0u> const &x )
 -> complex_rt<decltype( -std::declval<T>() ), 0u>
{ return {-x[ 0 ]}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline constexpr
auto  operator -( complex_rt<T, R> const &x )
 -> complex_rt<decltype( -std::declval<T>() ), R>
{ return {-x.lower_barrage(), -x.upper_barrage()}; }

/** \brief  Subtraction

Calculates the difference of the given values.

The definition can be broken down as:
- Real: `rA - rB`
- Component-wise: `{ cA[0] - cB[0], ..., cA[ 2^Rank - 1 ] - cB[ 2^Rank - 1 ] }`
- Barrage-wise: `{ lowerA - lowerB; upperA - upperB }`

    \relates  #boost::math::complex_rt

    \pre  `declval<T>() - declval<U>()` is well-formed.  Let's call the result
          type `SS`.
    \pre  `+declval<T>()` and `-declval<U>()` are well-formed.
    \pre  Both `T` and `U` can implicitly be assigned to a `SS`.

    \param[in] minuend     The value to be subtracted from.
    \param[in] subtrahend  The value to be subtracted.

    \returns  The difference of the `subtrahend` from the `minuend`.
 */
template < typename T, typename U >
inline constexpr
auto  operator -( complex_rt<T, 0u> const &minuend, complex_rt<U, 0u> const
 &subtrahend )
 -> complex_rt<decltype( std::declval<T>() - std::declval<U>() ), 0u>
{ return {minuend[ 0 ] - subtrahend[ 0 ]}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, typename U, std::size_t R >
inline constexpr
auto  operator -( complex_rt<T, R> const &minuend, complex_rt<U, R> const
 &subtrahend )
 -> complex_rt<decltype( std::declval<T>() - std::declval<U>() ), R>
{
    return { minuend.lower_barrage() - subtrahend.lower_barrage(),
     minuend.upper_barrage() - subtrahend.upper_barrage() };
}

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline constexpr
auto  operator -( complex_rt<T, R> const &minuend, complex_rt<U, S> const
 &subtrahend )
 -> typename std::enable_if< (R < S), complex_rt<decltype( std::declval<T>() -
 std::declval<U>() ), S> >::type
{ return {minuend - subtrahend.lower_barrage(), -subtrahend.upper_barrage()}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline constexpr
auto  operator -( complex_rt<T, R> const &minuend, complex_rt<U, S> const
 &subtrahend )
 -> typename std::enable_if< (R > S), complex_rt<decltype( std::declval<T>() -
 std::declval<U>() ), R> >::type
{ return {minuend.lower_barrage() - subtrahend, +minuend.upper_barrage()}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T >
inline constexpr
auto  operator -( T const &minuend, complex_rt<T, 0u> const &subtrahend )
 -> complex_rt<decltype( std::declval<T>() - std::declval<T>() ), 0u>
{ return {minuend - subtrahend[ 0 ]}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline constexpr
auto  operator -( T const &minuend, complex_rt<T, R> const &subtrahend )
 -> complex_rt<decltype( std::declval<T>() - std::declval<T>() ), R>
{ return {minuend - subtrahend.lower_barrage(), -subtrahend.upper_barrage()}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T >
inline constexpr
auto  operator -( complex_rt<T, 0u> const &minuend, T const &subtrahend )
 -> complex_rt<decltype( std::declval<T>() - std::declval<T>() ), 0u>
{ return {minuend[ 0 ] - subtrahend}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline constexpr
auto  operator -( complex_rt<T, R> const &minuend, T const &subtrahend )
 -> complex_rt<decltype( std::declval<T>() - std::declval<T>() ), R>
{ return {minuend.lower_barrage() - subtrahend, +minuend.upper_barrage()}; }

/** \brief  Subtract-and-assign

Calculates the difference of the given objects into the first.

    \relates  #boost::math::complex_rt

    \pre  `declval<T &>() -= declval<U>()` is well-formed.
    \pre  The rank of `subtrahend` doesn't exceed that of `minuend_difference`.

    \param[in,out] minuend_difference  The value to be subtracted from, and the
                                       location of the future difference.
    \param[in]     subtrahend          The value to be subtracted.

    \returns  A reference to post-subtraction `minuend_difference`.
 */
template < typename T, typename U >
inline
auto  operator -=( complex_rt<T, 0u> &minuend_difference, complex_rt<U, 0u>
 const &subtrahend )
 -> complex_rt<T, 0u> &
{ return minuend_difference[0] -= subtrahend[0], minuend_difference; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, typename U, std::size_t R >
inline
auto  operator -=( complex_rt<T, R> &minuend_difference, complex_rt<U, R> const
 &subtrahend )
 -> complex_rt<T, R> &
{
    minuend_difference.lower_barrage() -= subtrahend.lower_barrage();
    minuend_difference.upper_barrage() -= subtrahend.upper_barrage();
    return minuend_difference;
}

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline
auto  operator -=( complex_rt<T, R> &minuend_difference, complex_rt<U, S> const
 &subtrahend )
 -> typename std::enable_if< (R > S), complex_rt<T, R> >::type &
{ return minuend_difference.lower_barrage() -= subtrahend, minuend_difference; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T >
inline
auto  operator -=( complex_rt<T, 0u> &minuend_difference, T const &subtrahend )
 -> complex_rt<T, 0u> &
{ return minuend_difference[0] -= subtrahend, minuend_difference; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline
auto  operator -=( complex_rt<T, R> &minuend_difference, T const &subtrahend )
 -> complex_rt<T, R> &
{ return minuend_difference.lower_barrage() -= subtrahend, minuend_difference; }

/** \brief  Pre-decrement

Applies the predecessor function to the given object, turning it into the
difference of the original value and one.  (This is the same inappropriate
application that the built-in floating types get.)  Only the real part is
affected.

    \relates  #boost::math::complex_rt

    \pre  `--declval<T &>()` is well-formed.

    \param[in,out] minuend_difference  The object to be affected.

    \returns  A reference to post-decrement `minuend_difference`.
 */
template < typename T >
inline
auto  operator --( complex_rt<T, 0u> &minuend_difference ) -> complex_rt<T,0u> &
{ return --minuend_difference[0], minuend_difference; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline
auto  operator --( complex_rt<T, R> &minuend_difference ) -> complex_rt<T, R> &
{ return --minuend_difference.lower_barrage(), minuend_difference; }

/** \brief  Post-decrement

Applies the predecessor function to the given object, turning it into the
difference of the original value and one.  (This is the same inappropriate
application that the built-in floating types get.)  Only the real part is
affected.

    \relates  #boost::math::complex_rt

    \pre  `declval<T &>()--` is well-formed.

    \param[in,out] minuend_difference  The object to be affected.

    \returns  A copy of pre-decrement `minuend_difference`.
 */
template < typename T >
inline
auto  operator --(complex_rt<T, 0u> &minuend_difference,int) -> complex_rt<T,0u>
{ return {minuend_difference[ 0 ]--}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline
auto  operator --( complex_rt<T, R> &minuend_difference,int ) -> complex_rt<T,R>
{
    return { minuend_difference.lower_barrage()--,
     minuend_difference.upper_barrage() };
}


//  Hypercomplex condition functions  ----------------------------------------//

/** \brief  Complex conjugate, in operator form

Returns the complex conjugate of the given value.  This is commonly given as
`conj(x)` in computer code, but the `~x` notation is reminiscent of the compact
notations of this operation in prose.  (The `operator ~` would be otherwise
unused, anyway.)

    \relates  #boost::math::complex_rt

    \param[in] x  The input value.

    \returns  `Conj(x)`.
 */
template < typename T >
inline constexpr
auto  operator ~( complex_rt<T, 0u> const &x ) -> complex_rt<T, 0u>
{ return {+x[ 0 ]}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline constexpr
auto  operator ~( complex_rt<T, R> const &x ) -> complex_rt<T, R>
{ return {~x.lower_barrage(), -x.upper_barrage()}; }

/** \brief  Complex conjugate

Returns the complex conjugate of the given value.  This function is a core
operation for complex numbers.  The component type must *not* have a version of
this function (or at least one that differs from the identity function), since
that implementation will be ignored.

The definition can be broken down as:
- Real: `+r`
- Component-wise: `{ +c[0], -c[1], ..., -c[2^Rank - 1] }`
- Barrage-wise: `{ Conj(Lower); -Upper }`

This function has to act as an anti-involution (modulo real-life computing
issues), which means:
- Self-inverse: `Conj( Conj(x) ) == x`.
- Linear w/ Adding: `Conj( x + y ) == Conj( x ) + Conj( y )`.
- Linear w/ Scale: `Conj( scalar * x ) == scalar * Conj( x )`.
- Reversed w/ Multiplying: `Conj( x * y ) == Conj( y ) * Conj( x )`.

    \relatesalso  #boost::math::complex_rt

    \param[in] x  The input value.

    \returns  The reflection of `x` on the real axis.
 */
template < typename T, std::size_t R >
inline constexpr
auto  conj( complex_rt<T, R> const &x ) -> complex_rt<T, R>
{ return ~x; }

/** \brief  Cayley norm

Returns the Cayley norm of the given value.  This is formed by multiplying a
value by its conjugate, which results in the square of the conventional norm
(i.e. Euclidean norm, absolute value, or distance).

The Cayley norm is always a non-negative real number.  It is zero only when the
input value is zero.  Real-life computational arithmetic may result in negative
values or a zero-norm from a non-zero value due to overflow, wraparound, or
other condition events.

The definition can be broken down as:
- Real: `r^2`
- Component-wise: `c[0]^2 + c[1]^2 +...+ c[2^Rank - 1]^2`
- Barrage-wise: `Norm(Lower) + Norm(Upper)`

    \relatesalso  #boost::math::complex_rt

    \param[in] x  The input value.

    \returns  `Norm(x) := Conj(x) * x`.
 */
template < typename T >
inline constexpr
auto  norm( complex_rt<T, 0u> const &x )
 -> decltype( std::declval<T>() * std::declval<T>() )
{ return x[0] * x[0]; }

/** \overload
    \relatesalso  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline constexpr
auto  norm( complex_rt<T, R> const &x )
 -> decltype( std::declval<T>() * std::declval<T>() )
{ return norm(x.lower_barrage()) + norm(x.upper_barrage()); }


//  Multiplication operators  ------------------------------------------------//

/** \brief  Multiplication, scalar

Calculates the product of the given values.  One of the factors is a real
scalar; multiplication between a scalar and a hypercomplex number is easy to
define.  (It's like scalar-vector multiplication.)

The definition can be broken down as:
- Real: `scalar * r`
- Component-wise: `{ scalar * c[0], ..., scalar * c[ 2^Rank - 1 ] }`
- Barrage-wise: `{ scalar * Lower, scalar * Upper }`

    \relates  #boost::math::complex_rt

    \pre  `declval<T>() * declval<T>()` is well-formed.

    \param[in] multiplicand  The first factor to be multiplied.
    \param[in] multiplier    The second factor to be multiplied.

    \returns  The product of `multiplicand` and `multiplier`.
 */
template < typename T >
inline constexpr
auto  operator *( T const &multiplicand, complex_rt<T, 0u> const &multiplier )
 -> complex_rt<decltype( std::declval<T>() * std::declval<T>() ), 0u>
{ return {multiplicand * multiplier[ 0 ]}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline constexpr
auto  operator *( T const &multiplicand, complex_rt<T, R> const &multiplier )
 -> complex_rt<decltype( std::declval<T>() * std::declval<T>() ), R>
{
    return { multiplicand * multiplier.lower_barrage(),
     multiplicand * multiplier.upper_barrage() };
}

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T >
inline constexpr
auto  operator *( complex_rt<T, 0u> const &multiplicand, T const &multiplier )
 -> complex_rt<decltype( std::declval<T>() * std::declval<T>() ), 0u>
{ return {multiplicand[ 0 ] * multiplier}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline constexpr
auto  operator *( complex_rt<T, R> const &multiplicand, T const &multiplier )
 -> complex_rt<decltype( std::declval<T>() * std::declval<T>() ), R>
{
    return { multiplicand.lower_barrage() * multiplier,
     multiplicand.upper_barrage() * multiplier };
}

/** \brief  Multiplication, Cayley

Calculates the product of the given values, where both are hypercomplex.  This
is the key operation when going up the Cayley-Dickson construction ladder.

Note that this operation is non-commutative starting from quaternions, and
non-associative starting from octonions.  (It's always power-associative.)

The definition can be broken down as:
- Real: `rA * rB`
- Component-wise: polynomial multiplication, where the two unit elements in
  each product term combine.  The combined unit is dependent on the input units
  involved, and (usually) their order.
- Barrage-wise: `{ lowerA * lowerB - Conj(upperB) * upperA, upperB * lowerA +
  upperA * Conj(lowerB) }`

There are eight combinations of the barrage definition, 3 binary decisions, that
still result in a purely-real Cayley-norm.
- Swap the factors of the first term of the product's lower barrage.
- Swap the factors of the second term of the product's lower barrage.
- Switch which factor of the lower barrage's second term is conjugated.

(There are four combinations for the product's upper barrage's formula, but
they're fixed, determined from the four states that the lower barrage's second
term can take.)  The combination used in library is the one used by Boost's
Quaternion and Octonion libraries, and the Cayley-Dickson Construction page on
Wikipedia.  It was also used by R. Shafer in his 1954 paper "On the algebras
formed by the Cayley-Dickson process."

    \relates  #boost::math::complex_rt

    \pre  `declval<T>() * declval<U>()` is well-formed.
    \pre  The various additive/subtractive operators are well-formed.

    \param[in] multiplicand  The first factor to be multiplied.
    \param[in] multiplier    The second factor to be multiplied.

    \returns  The product of `multiplicand` and `multiplier`.
 */
template < typename T, typename U >
inline constexpr
auto  operator *( complex_rt<T, 0u> const &multiplicand, complex_rt<U, 0u> const
 &multiplier )
 -> complex_rt<decltype( std::declval<T>() * std::declval<U>() ), 0u>
{ return {multiplicand[ 0 ] * multiplier[ 0 ]}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, typename U, std::size_t R >
inline constexpr
auto  operator *( complex_rt<T, R> const &multiplicand, complex_rt<U, R> const
 &multiplier )
 -> complex_rt<decltype( std::declval<T>() * std::declval<U>() ), R>
{
    return { multiplicand.lower_barrage() * multiplier.lower_barrage() -
     ~multiplier.upper_barrage() * multiplicand.upper_barrage(),
     multiplier.upper_barrage() * multiplicand.lower_barrage() +
     multiplicand.upper_barrage() * ~multiplier.lower_barrage() };
}

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline constexpr
auto  operator *( complex_rt<T, R> const &multiplicand, complex_rt<U, S> const
 &multiplier )
 -> typename std::enable_if< (R < S), complex_rt<decltype( std::declval<T>() *
 std::declval<U>() ), S> >::type
{
    return { multiplicand * multiplier.lower_barrage(),
     multiplier.upper_barrage() * multiplicand };
}

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline constexpr
auto  operator *( complex_rt<T, R> const &multiplicand, complex_rt<U, S> const
 &multiplier )
 -> typename std::enable_if< (R > S), complex_rt<decltype( std::declval<T>() *
 std::declval<U>() ), R> >::type
{
    return { multiplicand.lower_barrage() * multiplier,
     multiplicand.upper_barrage() * ~multiplier };
}

/** \brief  Multiply-and-assign, scalar

Calculates the product of the given objects into the first.

    \relates  #boost::math::complex_rt

    \pre  `declval<T &>() *= declval<T>()` is well-formed.

    \param[in,out] multiplicand_product  The first factor to be multiplied, and
                                         the location of the future product.
    \param[in]     multiplier            The second factor to be multiplied.

    \returns  A reference to post-multiplication `multiplicand_product`.
 */
template < typename T >
inline
auto  operator *=( complex_rt<T,0u> &multiplicand_product, T const &multiplier )
 -> complex_rt<T, 0u> &
{ return multiplicand_product[0] *= multiplier, multiplicand_product; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline
auto  operator *=( complex_rt<T, R> &multiplicand_product, T const &multiplier )
 -> complex_rt<T, R> &
{
    multiplicand_product.lower_barrage() *= multiplier;
    multiplicand_product.upper_barrage() *= multiplier;
    return multiplicand_product;
}

/** \brief  Multiply-and-assign, Cayley

Calculates the product of the given objects into the first.

    \relates  #boost::math::complex_rt

    \pre  `declval<T &>() *= declval<U>()` is well-formed.
    \pre  The various additive/subtractive operators are well-formed.
    \pre  `multiplier`'s rank doesn't exceed that of `multiplicand_product`.

    \param[in,out] multiplicand_product  The first factor to be added, and the
                                         location of the future product.
    \param[in]     multiplier            The second factor to be multiplied.

    \returns  A reference to post-multiplication `multiplicand_product`.
 */
template < typename T, typename U >
inline
auto  operator *=( complex_rt<T, 0u> &multiplicand_product, complex_rt<U, 0u>
 const &multiplier )
 -> complex_rt<T, 0u> &
{ return multiplicand_product[0] *= multiplier[0], multiplicand_product; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, typename U, std::size_t R >
auto  operator *=( complex_rt<T, R> &multiplicand_product, complex_rt<U, R>
 const &multiplier )
 -> complex_rt<T, R> &
{
    auto const  mp_copy = multiplicand_product;

    multiplicand_product.lower_barrage() *=  multiplier.lower_barrage();
    multiplicand_product.upper_barrage() *= ~multiplier.lower_barrage();
    multiplicand_product.lower_barrage() -= ~multiplier.upper_barrage() *
     mp_copy.upper_barrage();
    multiplicand_product.upper_barrage() +=  multiplier.upper_barrage() *
     mp_copy.lower_barrage();
    return multiplicand_product;
}

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline
auto  operator *=( complex_rt<T, R> &multiplicand_product, complex_rt<U, S>
 const &multiplier )
 -> typename std::enable_if< (R > S), complex_rt<T, R> >::type &
{
    multiplicand_product.lower_barrage() *=  multiplier;
    multiplicand_product.upper_barrage() *= ~multiplier;
    return multiplicand_product;
}


//  Division operators  ------------------------------------------------------//

/** \brief  Division, scalar

Calculates the quotient of the given values.  The divisor is a real scalar;
division by a scalar upon a hypercomplex number is easy to define.  (It's like
element-wise division.)

A scalar dividend and hypercomplex divisor is covered under Cayley division.

The type of division, quotient-and-remainder (integer) or exact-quotient
(floating or rational), matches that of the type's `value_type`.

The definition can be broken down as:
- Real: `r / scalar`
- Component-wise: `{ c[0] / scalar, ..., c[ 2^Rank - 1 ] / scalar }`
- Barrage-wise: `{ Lower / scalar, Upper / scalar }`

    \relates  #boost::math::complex_rt

    \pre  `declval<T>() / declval<T>()` is well-formed.
    \pre  `divisor` is *not* zero.

    \param[in] dividend  The value to be divided.
    \param[in] divisor   The value to divide by.

    \returns  The quotient of `dividend` and `divisor`.
 */
template < typename T >
inline constexpr
auto  operator /( complex_rt<T, 0u> const &dividend, T const &divisor )
 -> complex_rt<decltype( std::declval<T>() / std::declval<T>() ), 0u>
{ return {dividend[ 0 ] / divisor}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline constexpr
auto  operator /( complex_rt<T, R> const &dividend, T const &divisor )
 -> complex_rt<decltype( std::declval<T>() / std::declval<T>() ), R>
{
    return { dividend.lower_barrage() / divisor, dividend.upper_barrage() /
     divisor };
}

/** \brief  Division, Cayley

Calculates the quotient of the given values, where both operands are
hypercomplex.  Division works via multiplication of a reciprocal (which makes it
an exact-quotient type of division).  The dividend is the left-factor (a.k.a.
the multiplicand) in this altered multiplication, and the reciprocal of the
divisor is the right factor (a.k.a. the multiplier).

The problem with just computing the divisor's reciprocal is integer types.  The
reciprocal of a hypercomplex number is its conjugate divided by its (scalar)
Cayley norm.  The norm is always larger than any of the components, so the
reciprocal always reduces to zero in integer arithmetic.  To prevent this, the
dividend should be multiplied by the divisor's conjugate before the scalar
division by the divisor's norm, to ensure that the first product is large enough
to give a useful answer after the division.  This method of applying the
operations can work in general, but it risks overflow, so it's safer to not use
it for types that don't need it.

This function uses traits from `std::numeric_limits` to determine if a
component type is integral before determining the order of operations.

    \relates  #boost::math::complex_rt

    \pre  `declval<T>() / declval<U>()` is well-formed.
    \pre  `divisor` is *not* zero.
    \pre  The various additive, subtractive, and multiplicative operators are
          well-formed.

    \param[in] dividend  The value to be divided.
    \param[in] divisor   The value to divide by.

    \returns  The quotient of `dividend` and `divisor`.
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline
auto  operator /( complex_rt<T, R> const &dividend, complex_rt<U, S> const
 &divisor )
 -> typename std::enable_if< std::numeric_limits<U>::is_integer,
 complex_rt<decltype( std::declval<T>() / std::declval<U>() ), ( R < S ) ? S :
 R> >::type  // for integers
{
    return (dividend * conj( divisor )) / static_cast<decltype(
     std::declval<T>() * std::declval<U>() )>(norm( divisor ));
}

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline
auto  operator /( complex_rt<T, R> const &dividend, complex_rt<U, S> const
 &divisor )
 -> typename std::enable_if< not std::numeric_limits<U>::is_integer,
 complex_rt<decltype( std::declval<T>() / std::declval<U>() ), ( R < S ) ? S :
 R> >::type
{ return dividend * (conj( divisor ) / norm( divisor )); }  // for non-integers

/** \brief  Modulo, scalar

Calculates the remainder of the given values.  The divisor is a real scalar;
division by a scalar upon a hypercomplex number is easy to define.  (It's like
element-wise division.)

A scalar dividend and hypercomplex divisor is covered under Cayley division.

The type's `value_type` should use quotient-and-remainder division.

The definition can be broken down as:
- Real: `r % scalar`
- Component-wise: `{ c[0] % scalar, ..., c[ 2^Rank - 1 ] % scalar }`
- Barrage-wise: `{ Lower % scalar, Upper % scalar }`

    \relates  #boost::math::complex_rt

    \pre  `declval<T>() % declval<T>()` is well-formed.
    \pre  `divisor` is *not* zero.

    \param[in] dividend  The value to be divided.
    \param[in] divisor   The value to divide by.

    \returns  The remainder from `divisor` dividing `dividend`.
 */
template < typename T >
inline constexpr
auto  operator %( complex_rt<T, 0u> const &dividend, T const &divisor )
 -> complex_rt<decltype( std::declval<T>() % std::declval<T>() ), 0u>
{ return {dividend[ 0 ] % divisor}; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline constexpr
auto  operator %( complex_rt<T, R> const &dividend, T const &divisor )
 -> complex_rt<decltype( std::declval<T>() % std::declval<T>() ), R>
{
    return { dividend.lower_barrage() % divisor, dividend.upper_barrage() %
     divisor };
}

/** \brief  Modulo, Cayley

Calculates the remainder of the given values, where both operands are
hypercomplex.  Division works via multiplication of a reciprocal (which makes it
an exact-quotient type of division).  The dividend is the left-factor (a.k.a.
the multiplicand) in this altered multiplication, and the reciprocal of the
divisor is the right factor (a.k.a. the multiplier).

To prevent underflow, hypercomplex numbers with integer-type components need to
multiply the dividend and the conjugate of the divisor before applying scalar
division with the divisor's norm.  (When the component type is not an integer
type, the divisor's conjugate and norm can be combined first to prevent
overflow.)  Since fractions are involved even when using integer components, the
result will usually be truncated from the true result.  This modulo operation
returns that trunctated part; it has no real meaning in Cayley division, since
it's an exact-quotient style.  (It can be used to recover the dividend given the
divisor, quotient, and this remainder.)

    \relates  #boost::math::complex_rt

    \pre  `declval<T>() % declval<U>()` is well-formed.
    \pre  `divisor` is *not* zero.
    \pre  The various additive, subtractive, and multiplicative operators are
          well-formed.

    \param[in] dividend  The value to be divided.
    \param[in] divisor   The value to divide by.

    \returns  The remainder from `divisor` dividing `dividend`.
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline
auto  operator %( complex_rt<T, R> const &dividend, complex_rt<U, S> const
 &divisor )
 -> complex_rt<decltype( std::declval<T>() % std::declval<U>() ), ( R < S ) ?
 S : R>
{ return dividend - ((dividend * conj(divisor)) / norm( divisor )) * divisor; }

/** \brief  Divide-and-assign

Calculates the quotient of the given objects into the first.

The type of division, quotient-and-remainder (integer) or exact-quotient
(floating or rational), matches that of the type's `value_type`.

    \relates  #boost::math::complex_rt

    \pre  `declval<T &>() /= declval<X>()` is well-formed, where `X` is `T` for
          scalar divisors, and `U` for hypercomplex ones.
    \pre  `divisor` is *not* zero.
    \pre  The various additive, subtractive, and multiplicative operators are
          well-formed when the divisor is hypercomplex.
    \pre  When `divisor` is a hypercomplex number, its rank doesn't exceed that
          of `dividend_quotient`.

    \param[in,out] dividend_quotient  The value to be divided, and the location
                                      of the future quotient.
    \param[in]     divisor            The value to divide by.

    \returns  A reference to post-division `dividend_quotient`.
 */
template < typename T >
inline
auto  operator /=( complex_rt<T, 0u> &dividend_quotient, T const &divisor )
 -> complex_rt<T, 0u> &
{ return dividend_quotient[0] /= divisor, dividend_quotient; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline
auto  operator /=( complex_rt<T, R> &dividend_quotient, T const &divisor )
 -> complex_rt<T, R> &
{
    dividend_quotient.lower_barrage() /= divisor;
    dividend_quotient.upper_barrage() /= divisor;
    return dividend_quotient;
}

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline
auto  operator /=( complex_rt<T, R> &dividend_quotient, complex_rt<U, S> const
 &divisor )
 -> typename std::enable_if< (R >= S) && std::numeric_limits<U>::is_integer,
 complex_rt<T, R> >::type &
{
    dividend_quotient *= conj( divisor );
    dividend_quotient /= norm( divisor );
    return dividend_quotient;
}

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline
auto  operator /=( complex_rt<T, R> &dividend_quotient, complex_rt<U, S> const
 &divisor )
 -> typename std::enable_if< (R >= S) && not std::numeric_limits<U>::is_integer,
 complex_rt<T, R> >::type &
{
    dividend_quotient *= conj( divisor ) / norm( divisor );
    return dividend_quotient;
}

/** \brief  Modulo-and-assign

Calculates the remainder of the given objects into the first.

The type's `value_type` should use quotient-and-remainder division.

    \relates  #boost::math::complex_rt

    \pre  `declval<T &>() %= declval<X>()` is well-formed, where `X` is `T` for
          scalar divisors, and `U` for hypercomplex ones.
    \pre  `divisor` is *not* zero.
    \pre  The various additive, subtractive, and multiplicative operators are
          well-formed when the divisor is hypercomplex.
    \pre  When `divisor` is a hypercomplex number, its rank doesn't exceed that
          of `dividend_remainder`.

    \param[in,out] dividend_remainder  The value to be divided, and the location
                                       of the future remainder.
    \param[in]     divisor             The value to divide by.

    \returns  A reference to post-division `dividend_remainder`.
 */
template < typename T >
inline
auto  operator %=( complex_rt<T, 0u> &dividend_remainder, T const &divisor )
 -> complex_rt<T, 0u> &
{ return dividend_remainder[0] %= divisor, dividend_remainder; }

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R >
inline
auto  operator %=( complex_rt<T, R> &dividend_remainder, T const &divisor )
 -> complex_rt<T, R> &
{
    dividend_remainder.lower_barrage() %= divisor;
    dividend_remainder.upper_barrage() %= divisor;
    return dividend_remainder;
}

/** \overload
    \relates  #boost::math::complex_rt
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline
auto  operator %=( complex_rt<T, R> &dividend_remainder, complex_rt<U, S> const
 &divisor )
 -> typename std::enable_if< (R >= S), complex_rt<T, R> >::type &
{ return dividend_remainder -= (dividend_remainder / divisor) * divisor; }


//  Component functions  -----------------------------------------------------//

/** \brief  Real part

Returns the real part of the given value.

The definition can be broken down as:
- Real: `r`
- Component-wise: `c[0]`
- Barrage-wise: `Re(Lower)`

    \param[in] x  The input value.

    \returns  `Re(x) := ( x + Conj(x) ) / 2`.
 */
template < typename T, std::size_t R >
inline
auto  real( complex_rt<T, R> const &x ) -> T
{ return x.real(); }

/** \brief  Imaginary part

Returns the imaginary part of the given value.  It doesn't have a natural
definition from Cayley-Dickson construction.

    \param[in] x  The input value.

    \returns  `Im(x) := Re( -i * x )`, where *i* is the classic imaginary unit.
 */
template < typename T, std::size_t R >
inline
auto  imag( complex_rt<T, R> const &x ) -> T
{ return x.imag(); }

/** \brief  Unreal part

Returns the unreal part of the given value.

The definition can be broken down as:
- Real: `0`
- Component-wise: `{ 0, c[1], ..., c[2^Rank - 1] }`
- Barrage-wise: `{ Ur(Lower); Upper }`

    \param[in] x  The input value.

    \returns  `Ur(x) := ( x - Conj(x) ) / 2`.
 */
template < typename T, std::size_t R >
inline
auto  unreal( complex_rt<T, R> const &x ) -> complex_rt<T, R>
{ return x.unreal(); }


//  Norm/distance functions  -------------------------------------------------//

/** \brief  Taxicab (L-1) norm / Manhattan distance

Returns the First order Lebesgue norm when treating the components of the given
number like a vector.  It is defined as:

- `Root<1>( |c[0]|^1 + ... + |c[2 ^ Rank - 1]|^1 )`

which equals the sum of the absolute value of the components.

Boost Quaternion and Octonion libraries call this function "l1," which is a
small letter "L" followed by the numeral "1."  This is hard to discern with some
fonts, since the characters look like each other (and the vertical line "|").
This is why this library uses a different name for the function.

    \pre  `declval<T>() < declval<T>()` is well-formed.
    \pre  `-declval<T>()` is well-formed.
    \pre  `+declval<T>()` is well-formed.
    \pre  `declval<T>() + declval<T>()` is well-formed.

    \param[in] x  The input value.

    \returns  `||x||_1 := Sum{i}( |c[i]| )`.
 */
template < typename T >
inline constexpr
T  taxi( complex_rt<T, 0u> const &x )
{ return (x[ 0 ] < T{}) ? -x[0] : +x[0]; }

//! \overload
template < typename T, std::size_t R >
inline constexpr
T  taxi( complex_rt<T, R> const &x )
{ return taxi(x.lower_barrage()) + taxi(x.upper_barrage()); }

/** \brief  Euclidean (L-2) norm / Euclidean distance

Returns the Second order Lebesgue norm when treating the components of the given
number like a vector.  It is also called the absolute value and the magnitude.
It is defined as:

- `Root<2>( |c[0]|^2 + ... + |c[2 ^ Rank - 1]|^2 )`

which equals the square root of the Cayley-norm.

    \pre  `sqrt( declval<T>() )` is well-formed, and `sqrt` can be found through
          ADL.
    \pre  For any given regular `T` value `y`, `y^2 == |y|^2`.

    \param[in] x  The input value.

    \returns  `||x||_2 := Sqrt( Sum{i}(|c[i]|^2) )`.
 */
template < typename T, std::size_t R >
inline
auto  abs( complex_rt<T, R> const &x )
 -> decltype( sqrt(std::declval<T>() * std::declval<T>()) )
{
    using std::sqrt;

    return sqrt( norm(x) );
}

/** \brief  Maximum (L-infinity) norm / Chebyshev distance

Returns the Infinite order Lebesgue norm when treating the components of the
given number like a vector.  It is also called the uniform norm, infinity norm,
and supremum (sup) norm.  It is defined as:

- `Limit{p -> Infinity}(Root<p>( |c[0]|^p + ... + |c[2 ^ Rank - 1]|^p ))`

which equals the maximum of the components' absolute values.

    \pre  `abs( declval<T>() )` is well-formed, and `abs` can be found through
          ADL.
    \pre  `declval<T>() < declval<T>()` is well-formed.

    \param[in] x  The input value.

    \returns  `||x||_infinity := Max{i}( |c[i]| )`.
 */
template < typename T >
inline
T  sup( complex_rt<T, 0u> const &x )
{
    using std::abs;

    return abs( x[0] );
}

//! \overload
template < typename T, std::size_t R >
inline
T  sup( complex_rt<T, R> const &x )
{ return std::max(sup( x.lower_barrage() ), sup( x.upper_barrage() )); }

/** \brief  Sign / Unit-vector

Returns the direction of a hypercomplex value from the origin, expressed as a
hypercomplex value with its magnitude nullified, i.e. scaled to 1.  When the
input value is zero, zero is the returned sign.

    \pre  `abs(x)` is well-formed.
    \pre  `decltype(x)` should support the exact-quotient style of division.

    \param[in] x  The input value.

    \returns  `Sgn( x ) := x / Abs( x )`, when `x` is not zero.
    \returns  0, when `x` is zero.
 */
template < typename T, std::size_t R >
inline
auto  sgn( complex_rt<T, R> const &x ) -> complex_rt<T, R>
{ return x ? x / abs(x) : x; }


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
