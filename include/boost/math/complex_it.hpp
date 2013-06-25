//  Boost Complex Numbers, iterative-mode header file  -----------------------//

//  Copyright 2013 Daryle Walker.
//  Distributed under the Boost Software License, Version 1.0.  (See the
//  accompanying file LICENSE_1_0.txt or a copy at
//  <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

/** \file   boost/math/complex_it.hpp
    \brief  A complex-number class template with iterative composition.

    \author  Daryle Walker

    \version  0.1

    \copyright  Boost Software License, version 1.0

    Contains the declaration (and definition) of class template `complex_it`,
    that models Cayley-Dickson hypercomplex numbers (complex, quaternions,
    octonions, etc.) by storing the components in an array and synthesizing the
    related operations iteratively.

    \warning  This library requires C++2011 features.
 */

#ifndef BOOST_MATH_COMPLEX_IT_HPP
#define BOOST_MATH_COMPLEX_IT_HPP

#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iterator>
#include <numeric>
#include <ostream>
#include <sstream>
#include <tuple>
#include <type_traits>
#include <utility>

// Put includes from Boost here.


namespace boost
{
namespace math
{


//  Hypercomplex number (iterative) class template definition  ---------------//

/** \brief  A complex number class storing its components in an array.

This class template models hypercomplex numbers built with Cayley-Dickson
construction, with the storage and operations defined iteratively.  This means
that the storage is a flat, one-dimensional C-level array and the algorithms
used for support operations use looping.

This class should satisfy the literal, standard-layout, trivially-copyable,
trivial, and/or P.O.D. properties iff the component type also satisfies the
corresponding property.  The type supports default-, copy-, and/or
move-construction, destruction, copy- and/or move-assignment if the component
type supports the corresponding operation.

The component type should be a regular type, but can generate condition events
(divide-by-zero, over-flow, etc.) or irregular values (NaN, infinities, etc.)
with certain combinations of operands.  Depending on the component type and/or
irregularity, encountering an irregularity may make further computations
useless.

Individual operations may have additional restrictions beyond those given below.

    \pre  `Number` has to meet the requirements of Section 26.2
          [numeric.requirements] of the C++2011 standard.
    \pre  The additive identity of `Number` has to be reachable via
          value-initialization (i.e. `Number{}`).
    \pre  The multiplicative identity of `Number` has to be reachable via a
          single application of `operator ++()` from the additive identity.
    \pre  `Number` has to support (contextual) Boolean conversion, with the
          additive identity mapping to `false` and all other regular values
          mapping to `true`.
    \pre  Either `Number` does not have a conjugate function (nominally
          called `conj`), or it has one but it's equivalent to the identity
          function.
    \pre  0 \<= `Rank` \< log2( maximum allowable size )

    \tparam Number  The component type
    \tparam Rank    The Cayley-Dickson construction level.  If not given, it
                    defaults to 1, in order to model regular complex numbers.
 */
template < typename Number, std::size_t Rank = 1u >
struct complex_it
{
    // Core types
    //! The type for size-based meta-data and access indices.
    typedef std::size_t  size_type;
    //! The component type.  Gives access to a template parameter.
    typedef Number      value_type;

    // Sizing parameters
    //! The rung of Cayley-Dickson construction.
    static constexpr  size_type  rank = Rank;
    //! The total number of components, of type #value_type.
    static constexpr  size_type  static_size = 1ULL << rank;
    //! Does this type have padding bytes?
    static constexpr  bool       has_padding = sizeof( complex_it ) >
     static_size * sizeof( value_type );

    // Support types
    /** \brief  The type immediately lower in Cayley-Dickson construction.

    Past the scalar base case, a new level of C.D. construction can be seen as a
    composite of two objects of the current level.

    When #rank is zero, this type-alias points to the current `complex_it`
    instantiation as a degenerate case.
     */
    typedef complex_it<value_type, rank - !!rank>  barrage_type;

    // Component(s) access
    /** \brief  Access to component data.

    The supplied index value is **not** bounds-checked.

        \pre  *i* \< #static_size

        \param[in] i  The index of the selected component.

        \throws  Nothing.

        \returns  A reference to the given component.
     */
    constexpr
    auto  operator []( size_type i ) const noexcept -> value_type const &
    { return c[i]; }
    //! \overload
    auto  operator []( size_type i ) noexcept -> value_type &  { return c[i]; }

    /** \brief    Real-component inspector
        \returns  `(*this)[0]`.
     */
    constexpr  auto  real() const -> value_type   { return c[0]; }
    /** \brief     Real-component mutator
        \param[in] r  The new real component value.
        \post      `(*this)[0] == r`.
     */
               void  real( value_type const &r )  { c[0] = r; }
    /** \brief    Imaginary-component inspector
        \returns  `(*this)[1]`.  (Or `value_type{}` when #rank is 0.)
     */
    constexpr  auto  imag() const -> value_type
    { return rank ? c[1] : value_type{}; }  // Hope the choice is optimized out.
    /** \brief     Imaginary-component mutator
        \pre       #rank \> 0.
        \param[in] i  The new (classic) imaginary component value.
        \post      `(*this)[1] == i`.
     */
               void  imag( value_type const &i )  { c[1] = i; }

    /** \brief  The lower (i.e. real-ward) decomposed half of this value

    When #rank is 0, `*this` is returned as a degenerate case.

    Cayley-Dickson construction builds a multi-valued scalar from composing two
    instances from the immediately-lower level.  This function decomposes a
    value to its two lower-ranked instances and returns the one made from the
    lower-numbered indices.

    The iteration-based class template version is `const` and returns a value,
    while the recursion-based class template versions return references, in
    `const` and mutable versions.

        \see  #upper_barrage()const

        \returns  `Lower( *this ) := { c[0]; c[1]; ...; c[2^(Rank - 1) - 1] }`
     */
    auto  lower_barrage() const -> barrage_type
    {
        barrage_type  result;

        std::copy( &c[0], &c[static_size / 2u + !rank], &result[0] );
        return result;
    }
    /** \brief     Mutate the lower decomposed half of this object
        \see       #upper_barrage(barrage_type const&)
        \see       #lower_barrage()const
        \param[in] b  The new value of the targeted half.
        \post      `this->lower_barrage() == b`.
     */
    void  lower_barrage( barrage_type const &b )
    { std::copy(&b[ 0 ], &b[ static_size / 2u + !rank ], &c[ 0 ]); }
    /** \brief  The upper (i.e. imaginary-ward) decomposed half of this value

    When #rank is 0, `*this` is returned as a degenerate case.

    Cayley-Dickson construction builds a multi-valued scalar from composing two
    instances from the immediately-lower level.  This function decomposes a
    value to its two lower-ranked instances and returns the one made from the
    higher-numbered indices.

    The iteration-based class template version is `const` and returns a value,
    while the recursion-based class template versions return references, in
    `const` and mutable versions.

        \see  #lower_barrage()const

        \returns  `Upper( *this ) := { c[2^(Rank - 1)]; ...; c[2^Rank - 1] }`
     */
    auto  upper_barrage() const -> barrage_type
    {
        barrage_type  result;

        std::copy( &c[static_size / 2u], &c[static_size], &result[0] );
        return result;
    }
    /** \brief     Mutate the upper decomposed half of this object
        \see       #lower_barrage(barrage_type const&)
        \see       #upper_barrage()const
        \param[in] b  The new value of the targeted half.
        \post      `this->upper_barrage() == b`.
     */
    void  upper_barrage( barrage_type const &b )
    { std::copy(&b[0 ], &b[static_size / 2u + !rank ], &c[static_size / 2u ]); }

    /** \brief  Unreal-components inspector

    The generalization of `imag` from regular-complex numbers to higher-level
    hypercomplex numbers is done by considering all the non-real components at
    once.  Since the unreal component vector has to stay a power-of-2 in length,
    the real-component spot is zero-filled.  (This also keeps each component
    associated with its corresponding hypercomplex unit.)

        \returns  A value `y` such that:
                  - `y[0] == value_type{}`.
                  - For 0 \< `k` \< #static_size, `(*this)[k] == y[k]`.
     */
    auto  unreal() const -> complex_it
    {
        complex_it  result = *this;

        result[ 0 ] = {};
        return result;
    }
    /** \brief     Unreal-components mutator
        \param[in] u  The new unreal component values.
        \post      `this->unreal() == u.unreal()`.
                   (The real parts don't change/get-copied!)
     */
    void  unreal( complex_it const &u )
    { std::copy(&u[ 1 ], &u[ static_size ], &c[ 1 ]); }

    // Conditions
    /** \brief  Boolean conversion

    Makes values of this type suitable for Boolean situations.  The member
    function is marked `explicit`, so it should be usable for contextual
    Boolean calls.  It uses contextual Boolean conversions of its components in
    turn.

    The definition can be broken down as:
    - Real: `static_cast<bool>(r)`
    - Component-wise: `c[0] || c[1] || ... || c[2^Rank - 1]`
    - Barrage-wise: `Lower || Upper`

        \returns  `Norm( *this ) != 0`
     */
    explicit
    operator bool() const
    { for (auto const &cc : c) if (cc) return true; return false; }

    // Constructors
    /** \brief  Default-construction

    The compiler-generated definition is used to keep `complex_it` trivial (when
    #value_type is trivial).

        \post  When used for value-initialization, `(*this)[i] == value_type{}`
               for every valid `i`.
        \post  When used for default-initialization, each component is left
               uninitialized if #value_type uses trivial default-construction,
               otherwise it acts the same as value-initialization.
     */
    complex_it() = default;
    /** \brief  List-of-reals constructor / Single-real conversion

    Constructs a `complex_it` object from a list of real numbers.  Can act as a
    conversion when given one real number.

        \pre  `sizeof...(i)` \< #static_size.
        \pre  Each parameter in `i` must have a (non-narrowing) implicit
              conversion to #value_type.

        \param[in] r  The real component.
        \param[in] i  The imaginary components.  May be empty.

        \post  `(*this)[0] == r`.
        \post  For `0 < k <= sizeof...(i)`, `(*this)[k] == Explode(i; k - 1)`.
        \post  For `sizeof...(i) < k <` #static_size, `(*this)[k] ==
               value_type{}`.
     */
    template < typename ...Args >
    constexpr  complex_it( value_type const &r, Args const &...i )
     : c{ r, i... }  {}

    /** \brief  Convert from one or more non-longer `complex_it` objects.

    Constructs a `complex_it` from a concatenation of components from a given
    list of `complex_it` objects.  The given sources do not need a common
    component type.  Can act as a conversion.

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
        typename      = typename std::enable_if<(R <= rank)>::type,
        typename      = typename std::enable_if<((1u + sizeof...( U )) <= (1ULL
         << ( rank - R )))>::type,
        typename      = typename std::enable_if<std::is_convertible<T,
         value_type>::value>::type
    >
    complex_it( complex_it<T, R> const &first, complex_it<U, R> const &...rest )
    { copy_barrages(0u, first, rest...); }
    /** \brief  Convert from a longer `complex_it` object.

    Constructs a `complex_it` from the first (i.e. real-ward) components of an
    object with a higher #rank.  The conversion is `explicit` to prevent
    ambiguities with the non-longer conversion constructor.

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
    explicit  complex_it( complex_it<T, R> const &senior )
    { std::copy(&senior[ 0 ], &senior[ static_size ], &c[ 0 ]); }
    /** \brief  Convert from `complex_it` objects with explicitly-convertible
                component types.

    Constructs a `complex_it` from the components of another object, but the
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
    explicit  complex_it( complex_it<T, R> const &cc )
        : c{}
    {
        std::transform( &cc[0], &cc[1ULL << std::min( rank, R )], &c[0],
         [](T const &t){return static_cast<value_type>( t );} );
    }

private:
    // Implements the cross-copy/barrage/sub-barrage constructor, base case
    void  copy_barrages( size_type start )
    { std::fill(&c[ start ], &c[ static_size ], value_type{}); }

    // Implements the cross-copy/barrage/sub-barrage constructor, main case
    // (Assumes `R <= rank` and `start + ((1 + sizeof...(U)) << (rank - R)) <=
    // static_size`.)
    template < typename T, size_type R, typename ...U >
    void  copy_barrages( size_type start, complex_it<T, R> const &first,
     complex_it<U, R> const &...rest )
    {
        constexpr auto  barrage_length = static_size >> ( rank - R );

        std::copy( &first[0], &first[barrage_length], &c[start] );
        copy_barrages( start + barrage_length, rest... );
    }

    // Member data
    value_type  c[ static_size ];
};


//  Class-static data member definitions  ------------------------------------//

/** Gives access to a template parameter.  The base level of 0 represents real
    numbers, 1 for complex numbers, 2 for quaternions, 3 for octonions, etc.
 */
template < typename Number, std::size_t Rank >
constexpr
typename complex_it<Number, Rank>::size_type  complex_it<Number, Rank>::rank;

/** The component count doubles as rank increases, starting at 1 for the base
    level.
 */
template < typename Number, std::size_t Rank >
constexpr
typename complex_it<Number, Rank>::size_type
  complex_it<Number, Rank>::static_size;

/** When this is `false`, then a pointer to an array segment of this type can be
    converted (with `reinterpret_cast`) to a pointer to an array segment of
    #value_type with #static_size times the elements.  (Or an array of this type
    can be converted to an array reference of `value_type`.)  This mirrors the
    same requirement on `std::complex`.

    You can also convert a pointer to one object of this type to a pointer to an
    array segment of #value_type, #static_size in length.  (Or a reference to
    this object to an array of `static_size` `value_type` objects.)  But even if
    `has_padding` is `true`, this conversion is still legal if this type is
    considered standard-layout (which requires `value_type` to be
    standard-layout).
 */
template < typename Number, std::size_t Rank >
constexpr
bool  complex_it<Number, Rank>::has_padding;


//  Implementation details  --------------------------------------------------//

//! \cond
namespace detail
{
    //! Detect if a type's swap (found via ADL for non-built-ins) throws.
    template < typename T, typename U = T >
    inline constexpr
    bool  is_swap_nothrow() noexcept
    {
        using std::swap;

        return noexcept( swap(std::declval<T &>(), std::declval<U &>()) );
    }

}  // namespace detail
//! \endcond


//  Object support functions  ------------------------------------------------//

/** \brief  Forward iteration over components, start point

Generates a start point for iterating over the given object's component data in
a forward direction.  Progression through the elements starts from the real
component and goes through a total of #boost::math::complex_it::static_size
elements.

This function is exclusive to iteration-based hypercomplex number types.

    \relatesalso  #boost::math::complex_it

    \see  #boost::math::end(boost::math::complex_it<T,R>const&)

    \param c  The object to have its begin-iterator accessed.

    \throws  Nothing.

    \returns  An iterator pointing to the first element (i.e. the real
              component).
 */
template < typename T, std::size_t R >
inline constexpr
auto  begin( complex_it<T, R> const &c ) noexcept -> T const *
{ return &c[0]; }

/** \overload
    \relatesalso  #boost::math::complex_it
 */
template < typename T, std::size_t R >
inline
auto  begin( complex_it<T, R> &c ) noexcept -> T *
{ return &c[0]; }

/** \brief  Forward iteration over components, end point

Generates an end point for iterating over the given object's component data in
a forward direction.

This function is exclusive to iteration-based hypercomplex number types.

    \relatesalso  #boost::math::complex_it

    \see  #boost::math::begin(boost::math::complex_it<T,R>const&)

    \param c  The object to have its end-iterator accessed.

    \throws  Nothing.

    \returns  An iterator pointing to past the last element.
 */
template < typename T, std::size_t R >
inline constexpr
auto  end( complex_it<T, R> const &c ) noexcept -> T const *
{ return begin(c) + complex_it<T, R>::static_size; }

/** \overload
    \relatesalso  #boost::math::complex_it
 */
template < typename T, std::size_t R >
inline
auto  end( complex_it<T, R> &c ) noexcept -> T *
{ return begin(c) + complex_it<T, R>::static_size; }

/** \brief  Swap routine for `complex_it`.

Exchanges the state of the two given objects.

    \relatesalso  #boost::math::complex_it

    \pre  There is a swapping routine, called `swap`, either in namespace `std`
          for built-ins, or found via ADL for other types.

    \param a  The first object to have its state exchanged.
    \param b  The second object to have its state exchanged.

    \throws Whatever  the element-level swap does.

    \post  `a` is equivalent to the old state of `b`, while `b` is equivalent to
           the old state of `a`.
 */
template < typename T, std::size_t R >
inline
void  swap( complex_it<T, R> &a, complex_it<T, R> &b )
 noexcept( detail::is_swap_nothrow<T>() )
{ std::swap_ranges(begin( a ), end( a ), begin( b )); }


//  Tuple interface functions  -----------------------------------------------//

/** \brief  Accesses a component of a `complex_it`.

Extracts the given component of the given complex number object.

    \pre  0 \<= *I* \< #boost::math::complex_it::static_size.

    \tparam I  The index of the desired component.

    \param c  The complex number object containing the component.

    \returns  A reference to the desired component.
 */
template < std::size_t I, typename T, std::size_t R >
inline constexpr
auto  get( complex_it<T, R> const &c ) noexcept -> T const &
{
    static_assert( I < complex_it<T, R>::static_size, "Index too large" );

    return c[ I ];
}

//! \overload
template < std::size_t I, typename T, std::size_t R >
inline
auto  get( complex_it<T, R> &c ) noexcept -> T &
{ return const_cast<T &>(get<I>( const_cast<complex_it<T,R> const &>(c) )); }

//! \overload
template < std::size_t I, typename T, std::size_t R >
inline
auto  get( complex_it<T, R> &&c ) noexcept -> T &&
{ return std::forward<T>(get<I>( c )); }


//  Equality operators  ------------------------------------------------------//

/** \brief  Inequality comparison for `complex_it`.

Compares two `complex_it` objects for inequality.  If the two objects differ in
length, the remaining components of the longer object will be checked against
the additive identity.

    \relates  #boost::math::complex_it

    \pre  `std::declval<T>() != std::declval<U>()` is well-formed and the result
          type is (contextually) convertible to `bool`.

    \param l  The left-side argument.
    \param r  The right-side argument.

    \retval false  If every component in *l* is equal to its corresponding
                   component in *r*; and when one object is longer, the excess
                   components in that object are all zero.
    \retval true   Otherwise.
 */
template < typename T, std::size_t R, typename U, std::size_t S >
bool  operator !=( complex_it<T, R> const &l, complex_it<U, S> const &r )
{
    auto        lb = begin( l );
    auto const  le = end( l );
    auto        rb = begin( r );
    auto const  re = end( r );

    while ( (le != lb) && (re != rb) )
        if ( *lb++ != *rb++ )
            return true;
    while ( le != lb )
        if ( *lb++ )
            return true;
    while ( re != rb )
        if ( *rb++ )
            return true;
    return false;
}

/** \overload
    \relates  #boost::math::complex_it
 */
template < typename T, std::size_t R >
bool  operator !=( complex_it<T, R> const &l, T const &r )
{
    auto        lb = begin( l );
    auto const  le = end( l );

    if ( *lb != r )
        return true;
    while ( le != ++lb )
        if ( *lb )
            return true;
    return false;
}

/** \overload
    \relates  #boost::math::complex_it
 */
template < typename T, std::size_t R >
inline
bool  operator !=( T const &l, complex_it<T, R> const &r )
{ return operator !=(r, l); }

/** \brief  Equality comparison for `complex_it`.

Compares two `complex_it` objects for equality.  If the two objects differ in
length, the remaining components of the longer object will be checked against
the additive identity.

    \relates  #boost::math::complex_it

    \pre  `std::declval<T>() != std::declval<U>()` is well-formed and the result
          type is (contextually) convertible to `bool`.

    \param l  The left-side argument.
    \param r  The right-side argument.

    \retval true   If the corresponding components between *l* and *r* compare
                   as equal; and when their lengths don't match, all the excess
                   components of the longer object are zero.
    \retval false  Otherwise.
 */
template < typename T, std::size_t R, typename U, std::size_t S >
inline
bool  operator ==( complex_it<T, R> const &l, complex_it<U, S> const &r )
{ return not operator !=(l, r); }

/** \overload
    \relates  #boost::math::complex_it
 */
template < typename T, std::size_t R >
inline
bool  operator ==( complex_it<T, R> const &l, T const &r )
{ return not operator !=(l, r); }

/** \overload
    \relates  #boost::math::complex_it
 */
template < typename T, std::size_t R >
inline
bool  operator ==( T const &l, complex_it<T, R> const &r )
{ return not operator !=(l, r); }


//  Input/output operators  --------------------------------------------------//

/** \brief  Output-streaming for `complex_it`.

Outputs a string representation of the given number to the given stream.  When
the rank of the number is zero, the sole (real) component is written; otherwise
a comma-separated list of the components, surrounded by parentheses, is written.

    \relatesalso  #boost::math::complex_it

    \pre  `std::declval<basic_ostream<Ch,Tr> &>() << std::declval<T>()` is
          well-formed and the result type is the same as the first operand.

    \param[in,out] o  The stream to send the output
    \param[in]     x  The complex number to be written

    \returns  `o`
 */
template < typename Ch, class Tr, typename T >
inline
std::basic_ostream<Ch, Tr> &
operator <<( std::basic_ostream<Ch, Tr> &o, complex_it<T, 0u> const &x )
{ return o << x[0]; }

/** \overload
    \relatesalso  #boost::math::complex_it
 */
template < typename Ch, class Tr, typename T, std::size_t R >
std::basic_ostream<Ch, Tr> &
operator <<( std::basic_ostream<Ch, Tr> &o, complex_it<T, R> const &x )
{
    std::basic_ostringstream<Ch, Tr>  s;
    auto                              b = begin( x );
    auto const                        e = end( x );

    s.flags( o.flags() );
    s.imbue( o.getloc() );
    s.precision( o.precision() );
    s << '(' << *b++;
    while ( e != b )
        s << ',' << *b++;
    s << ')';
    return o << s.str();
}


//  Addition operators  ------------------------------------------------------//

/** \brief  Identity operator

Returns the given value (possibly normalized).

The definition can be broken down as:
- Real: `+r`
- Component-wise: `{ +c[0], +c[1], ..., +c[2^Rank - 1] }`
- Barrage-wise: `{ +Lower; +Upper }`

    \relates  #boost::math::complex_it

    \pre  `+declval<T>()` is well-formed.

    \param[in] x  The input value.

    \returns  `x`.
 */
template < typename T, std::size_t R >
auto  operator +( complex_it<T, R> const &x )
 -> complex_it<decltype( +std::declval<T>() ), R>
{
    decltype( +x )  result;
    auto            rb = begin( result );

    for ( auto const &xx : x )
        *rb++ = +xx;
    return result;
}

/** \brief  Addition

Calculates the sum of the given values.

The definition can be broken down as:
- Real: `rA + rB`
- Component-wise: `{ cA[0] + cB[0], ..., cA[ 2^Rank - 1 ] + cB[ 2^Rank - 1 ] }`
- Barrage-wise: `{ lowerA + lowerB; upperA + upperB }`

    \relates  #boost::math::complex_it

    \pre  `declval<T>() + declval<U>()` is well-formed.  Let's call the result
          type `SS`.
    \pre  `+declval<T>()` and `+declval<U>()` are well-formed.
    \pre  Both `T` and `U` can implicitly be assigned to a `SS`.

    \param[in] augend  The first value to be added.
    \param[in] addend  The second value to be added.

    \returns  The sum of `augend` and `addend`.
 */
template < typename T, std::size_t R, typename U, std::size_t S >
auto  operator +( complex_it<T,R> const &augend, complex_it<U,S> const &addend )
 -> complex_it<decltype( std::declval<T>() + std::declval<U>() ), ( (R<S)?S:R )>
{
    decltype( augend + addend )  sum;
    auto                         sb = begin( sum );
    auto                        aub = begin( augend );
    auto const                  aue = end( augend );
    auto                        adb = begin( addend );
    auto const                  ade = end( addend );

    while ( (aue != aub) && (ade != adb) )
        *sb++ = *aub++ + *adb++;
    while ( aue != aub )
        *sb++ = +*aub++;
    while ( ade != adb )
        *sb++ = +*adb++;
    // assert(sb == end(sum)); assert(aub == aue); assert(adb == ade);
    return sum;
}

/** \overload
    \relates  #boost::math::complex_it
 */
template < typename T, std::size_t R >
auto  operator +( T const &augend, complex_it<T, R> const &addend )
 -> complex_it<decltype( std::declval<T>() + std::declval<T>() ), R>
{
    decltype( augend + addend )  sum = +addend;

    sum[ 0 ] = augend + addend[ 0 ];
    return sum;
}

/** \overload
    \relates  #boost::math::complex_it
 */
template < typename T, std::size_t R >
auto  operator +( complex_it<T, R> const &augend, T const &addend )
 -> complex_it<decltype( std::declval<T>() + std::declval<T>() ), R>
{
    decltype( augend + addend )  sum = augend;

    sum[ 0 ] = augend[ 0 ] + addend;
    return sum;
}

/** \brief  Add-and-assign

Calculates the sum of the given objects into the first.

    \relates  #boost::math::complex_it

    \pre  `declval<T &>() += declval<U>()` is well-formed.
    \pre  The rank of `addend` doesn't exceed that of `augend_sum`.

    \param[in,out] augend_sum  The first value to be added, and the location of
                               the future sum.
    \param[in]     addend      The second value to be added.

    \returns  A reference to post-addition `augend_sum`.
 */
template < typename T, std::size_t R, typename U, std::size_t S >
auto  operator +=( complex_it<T,R> &augend_sum, complex_it<U,S> const &addend )
 -> typename std::enable_if< (R >= S), complex_it<T, R> >::type &
{
    auto  asb = begin( augend_sum );

    for ( auto const &x : addend )
        *asb++ += x;
    // assert(asb <= end(augend_sum));
    return augend_sum;
}

/** \overload
    \relates  #boost::math::complex_it
 */
template < typename T, std::size_t R >
inline
auto  operator +=( complex_it<T, R> &augend_sum, T const &addend )
 -> complex_it<T, R> &
{ return augend_sum[0] += addend, augend_sum; }

/** \brief  Pre-increment

Applies the successor function to the given object, turning it into the sum of
the original value and one.  (This is the same inappropriate application that
the built-in floating types get.)  Only the real part is affected.

    \relates  #boost::math::complex_it

    \pre  `++declval<T &>()` is well-formed.

    \param[in,out] augend_sum  The object to be affected.

    \returns  A reference to post-increment `augend_sum`.
 */
template < typename T, std::size_t R >
inline
auto  operator ++( complex_it<T, R> &augend_sum ) -> complex_it<T, R> &
{ return ++augend_sum[0], augend_sum; }

/** \brief  Post-increment

Applies the successor function to the given object, turning it into the sum of
the original value and one.  (This is the same inappropriate application that
the built-in floating types get.)  Only the real part is affected.

    \relates  #boost::math::complex_it

    \pre  `declval<T &>()++` is well-formed.

    \param[in,out] augend_sum  The object to be affected.

    \returns  A copy of pre-increment `augend_sum`.
 */
template < typename T, std::size_t R >
auto  operator ++( complex_it<T, R> &augend_sum, int ) -> complex_it<T, R>
{
    auto const  copy = augend_sum;

    augend_sum[ 0 ]++;
    return copy;
}


//  Subtraction operators  ---------------------------------------------------//

/** \brief  Negation operator

Returns the additive inverse of the given value.

The definition can be broken down as:
- Real: `-r`
- Component-wise: `{ -c[0], -c[1], ..., -c[2^Rank - 1] }`
- Barrage-wise: `{ -Lower; -Upper }`

    \relates  #boost::math::complex_it

    \pre  `-declval<T>()` is well-formed.

    \param[in] x  The input value.

    \returns  A value `y` such that `x + y == decltype(x){}`.
 */
template < typename T, std::size_t R >
auto  operator -( complex_it<T, R> const &x )
 -> complex_it<decltype( -std::declval<T>() ), R>
{
    decltype( -x )  result;
    auto            rb = begin( result );

    for ( auto const &xx : x )
        *rb++ = -xx;
    return result;
}

/** \brief  Subtraction

Calculates the difference of the given values.

The definition can be broken down as:
- Real: `rA - rB`
- Component-wise: `{ cA[0] - cB[0], ..., cA[ 2^Rank - 1 ] - cB[ 2^Rank - 1 ] }`
- Barrage-wise: `{ lowerA - lowerB; upperA - upperB }`

    \relates  #boost::math::complex_it

    \pre  `declval<T>() - declval<U>()` is well-formed.  Let's call the result
          type `SS`.
    \pre  `+declval<T>()` and `-declval<U>()` are well-formed.
    \pre  Both `T` and `U` can implicitly be assigned to a `SS`.

    \param[in] minuend     The value to be subtracted from.
    \param[in] subtrahend  The value to be subtracted.

    \returns  The difference of the `subtrahend` from the `minuend`.
 */
template < typename T, std::size_t R, typename U, std::size_t S >
auto  operator -(complex_it<T,R>const &minuend,complex_it<U,S>const &subtrahend)
 -> complex_it<decltype( std::declval<T>() + std::declval<U>() ), ( (R<S)?S:R )>
{
    decltype( minuend - subtrahend )  difference;
    auto                              db = begin( difference );
    auto                              mb = begin( minuend );
    auto const                        me = end( minuend );
    auto                              sb = begin( subtrahend );
    auto const                        se = end( subtrahend );

    while ( (me != mb) && (se != sb) )
        *db++ = *mb++ - *sb++;
    while ( me != mb )
        *db++ = +*mb++;
    while ( se != sb )
        *db++ = -*sb++;
    // assert(db == end(difference)); assert(mb == me); assert(sb == se);
    return difference;
}

/** \overload
    \relates  #boost::math::complex_it
 */
template < typename T, std::size_t R >
auto  operator -( T const &minuend, complex_it<T, R> const &subtrahend )
 -> complex_it<decltype( std::declval<T>() - std::declval<T>() ), R>
{
    decltype( minuend - subtrahend )  difference = -subtrahend;

    difference[ 0 ] = minuend - subtrahend[ 0 ];
    return difference;
}

/** \overload
    \relates  #boost::math::complex_it
 */
template < typename T, std::size_t R >
auto  operator -( complex_it<T, R> const &minuend, T const &subtrahend )
 -> complex_it<decltype( std::declval<T>() - std::declval<T>() ), R>
{
    decltype( minuend - subtrahend )  difference = minuend;

    difference[ 0 ] = minuend[ 0 ] - subtrahend;
    return difference;
}

/** \brief  Subtract-and-assign

Calculates the difference of the given objects into the first.

    \relates  #boost::math::complex_it

    \pre  `declval<T &>() -= declval<U>()` is well-formed.
    \pre  The rank of `subtrahend` doesn't exceed that of `minuend_difference`.

    \param[in,out] minuend_difference  The value to be subtracted from, and the
                                       location of the future difference.
    \param[in]     subtrahend          The value to be subtracted.

    \returns  A reference to post-subtraction `minuend_difference`.
 */
template < typename T, std::size_t R, typename U, std::size_t S >
auto  operator -=( complex_it<T, R> &minuend_difference, complex_it<U, S> const
 &subtrahend )
 -> typename std::enable_if< (R >= S), complex_it<T, R> >::type &
{
    auto  mdb = begin( minuend_difference );

    for ( auto const &x : subtrahend )
        *mdb++ -= x;
    // assert(mdb <= end(minuend_difference));
    return minuend_difference;
}

/** \overload
    \relates  #boost::math::complex_it
 */
template < typename T, std::size_t R >
inline
auto  operator -=( complex_it<T, R> &minuend_difference, T const &subtrahend )
 -> complex_it<T, R> &
{ return minuend_difference[0] -= subtrahend, minuend_difference; }

/** \brief  Pre-decrement

Applies the predecessor function to the given object, turning it into the
difference of the original value and one.  (This is the same inappropriate
application that the built-in floating types get.)  Only the real part is
affected.

    \relates  #boost::math::complex_it

    \pre  `--declval<T &>()` is well-formed.

    \param[in,out] minuend_difference  The object to be affected.

    \returns  A reference to post-decrement `minuend_difference`.
 */
template < typename T, std::size_t R >
inline
auto  operator --( complex_it<T, R> &minuend_difference ) -> complex_it<T, R> &
{ return --minuend_difference[0], minuend_difference; }

/** \brief  Post-decrement

Applies the predecessor function to the given object, turning it into the
difference of the original value and one.  (This is the same inappropriate
application that the built-in floating types get.)  Only the real part is
affected.

    \relates  #boost::math::complex_it

    \pre  `declval<T &>()--` is well-formed.

    \param[in,out] minuend_difference  The object to be affected.

    \returns  A copy of pre-decrement `minuend_difference`.
 */
template < typename T, std::size_t R >
auto  operator --(complex_it<T, R> &minuend_difference, int) -> complex_it<T, R>
{
    auto const  copy = minuend_difference;

    minuend_difference[ 0 ]--;
    return copy;
}


//  Hypercomplex condition operator and functions  ---------------------------//

/** \brief  Complex conjugate, in operator form

Returns the complex conjugate of the given value.  This is commonly given as
`conj(x)` in computer code, but the `~x` notation is reminiscent of the compact
notations of this operation in prose.  (The `operator ~` would be otherwise
unused, anyway.)

    \relates  #boost::math::complex_it

    \param[in] x  The input value.

    \returns  `Conj(x)`.
 */
template < typename T, std::size_t R >
auto  operator ~( complex_it<T, R> const &x ) -> complex_it<T, R>
{
    decltype( operator~(x) )  result;  // "decltype(~x)" causes ICE
    auto                      rb = begin( result );
    auto const                re = end( result );
    auto                      xb = begin( x );

    *rb++ = +*xb++;
    while ( re != rb )
        *rb++ = -*xb++;
    return result;
}

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

    \relatesalso  #boost::math::complex_it

    \param[in] x  The input value.

    \returns  The reflection of `x` on the real axis.
 */
template < typename T, std::size_t R >
inline
auto  conj( complex_it<T, R> const &x ) -> complex_it<T, R>
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

    \relatesalso  #boost::math::complex_it

    \param[in] x  The input value.

    \returns  `Norm(x) := Conj(x) * x`.
 */
template < typename T, std::size_t R >
inline
auto  norm( complex_it<T, R> const &x )
 -> decltype( std::declval<T>() * std::declval<T>() )
{ return std::inner_product(begin(x), end(x), begin(x), decltype(norm(x)){}); }


//  Multiplication operators  ------------------------------------------------//

/** \brief  Multiplication, scalar

Calculates the product of the given values.  One of the factors is a real
scalar; multiplication between a scalar and a hypercomplex number is easy to
define.  (It's like scalar-vector multiplication.)

The definition can be broken down as:
- Real: `scalar * r`
- Component-wise: `{ scalar * c[0], ..., scalar * c[ 2^Rank - 1 ] }`
- Barrage-wise: `{ scalar * Lower, scalar * Upper }`

    \relates  #boost::math::complex_it

    \pre  `declval<T>() * declval<T>()` is well-formed.

    \param[in] multiplicand  The first factor to be multiplied.
    \param[in] multiplier    The second factor to be multiplied.

    \returns  The product of `multiplicand` and `multiplier`.
 */
template < typename T, std::size_t R >
auto  operator *( T const &multiplicand, complex_it<T, R> const &multiplier )
 -> complex_it<decltype( std::declval<T>() * std::declval<T>() ), R>
{
    decltype( multiplicand * multiplier )  product;
    auto                                   pb = begin( product );

    for ( auto const &mm : multiplier )
        *pb++ = multiplicand * mm;
    // assert( pb == end(product) );
    return product;
}

/** \overload
    \relates  #boost::math::complex_it
 */
template < typename T, std::size_t R >
auto  operator *( complex_it<T, R> const &multiplicand, T const &multiplier )
 -> complex_it<decltype( std::declval<T>() * std::declval<T>() ), R>
{
    decltype( multiplicand * multiplier )  product;
    auto                                   pb = begin( product );

    for ( auto const &mm : multiplicand )
        *pb++ = mm * multiplier;
    // assert( pb == end(product) );
    return product;
}

/** \brief  Multiply-and-assign, scalar

Calculates the product of the given objects into the first.

    \relates  #boost::math::complex_it

    \pre  `declval<T &>() *= declval<T>()` is well-formed.

    \param[in,out] multiplicand_product  The first factor to be multiplied, and
                                         the location of the future product.
    \param[in]     multiplier            The second factor to be multiplied.

    \returns  A reference to post-multiplication `multiplicand_product`.
 */
template < typename T, std::size_t R >
auto  operator *=( complex_it<T, R> &multiplicand_product, T const &multiplier )
 -> complex_it<T, R> &
{
    for ( auto &mp : multiplicand_product )
        mp *= multiplier;
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

    \relates  #boost::math::complex_it

    \pre  `declval<T>() / declval<T>()` is well-formed.
    \pre  `divisor` is *not* zero.

    \param[in] dividend  The value to be divided.
    \param[in] divisor   The value to divide by.

    \returns  The quotient of `dividend` and `divisor`.
 */
template < typename T, std::size_t R >
auto  operator /( complex_it<T, R> const &dividend, T const &divisor )
 -> complex_it<decltype( std::declval<T>() / std::declval<T>() ), R>
{
    decltype( dividend / divisor )  quotient;
    auto                            qb = begin( quotient );

    for ( auto const &dd : dividend )
        *qb++ = dd / divisor;
    // assert( qb == end(quotient) );
    return quotient;
}

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

    \relates  #boost::math::complex_it

    \pre  `declval<T>() % declval<T>()` is well-formed.
    \pre  `divisor` is *not* zero.

    \param[in] dividend  The value to be divided.
    \param[in] divisor   The value to divide by.

    \returns  The remainder from `divisor` dividing `dividend`.
 */
template < typename T, std::size_t R >
auto  operator %( complex_it<T, R> const &dividend, T const &divisor )
 -> complex_it<decltype( std::declval<T>() % std::declval<T>() ), R>
{
    decltype( dividend % divisor )  remainder;
    auto                            rb = begin( remainder );

    for ( auto const &dd : dividend )
        *rb++ = dd % divisor;
    // assert( rb == end(remainder) );
    return remainder;
}

/** \brief  Divide-and-assign, scalar

Calculates the quotient of the given objects into the first.

The type of division, quotient-and-remainder (integer) or exact-quotient
(floating or rational), matches that of the type's `value_type`.

    \relates  #boost::math::complex_it

    \pre  `declval<T &>() /= declval<T>()` is well-formed.
    \pre  `divisor` is *not* zero.

    \param[in,out] dividend_quotient  The value to be divided, and the location
                                      of the future quotient.
    \param[in]     divisor            The value to divide by.

    \returns  A reference to post-division `dividend_quotient`.
 */
template < typename T, std::size_t R >
auto  operator /=( complex_it<T, R> &dividend_quotient, T const &divisor )
 -> complex_it<T, R> &
{
    for ( auto &dq : dividend_quotient )
        dq /= divisor;
    return dividend_quotient;
}

/** \brief  Modulo-and-assign, scalar

Calculates the remainder of the given objects into the first.

The type's `value_type` should use quotient-and-remainder division.

    \relates  #boost::math::complex_it

    \pre  `declval<T &>() %= declval<T>()` is well-formed.
    \pre  `divisor` is *not* zero.

    \param[in,out] dividend_remainder  The value to be divided, and the location
                                       of the future remainder.
    \param[in]     divisor             The value to divide by.

    \returns  A reference to post-division `dividend_remainder`.
 */
template < typename T, std::size_t R >
auto  operator %=( complex_it<T, R> &dividend_remainder, T const &divisor )
 -> complex_it<T, R> &
{
    for ( auto &dr : dividend_remainder )
        dr %= divisor;
    return dividend_remainder;
}


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
auto  real( complex_it<T, R> const &x ) -> T
{ return x.real(); }

/** \brief  Imaginary part

Returns the imaginary part of the given value.  It doesn't have a natural
definition from Cayley-Dickson construction.

    \param[in] x  The input value.

    \returns  `Im(x) := Re( -i * x )`, where *i* is the classic imaginary unit.
 */
template < typename T, std::size_t R >
inline
auto  imag( complex_it<T, R> const &x ) -> T
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
auto  unreal( complex_it<T, R> const &x ) -> complex_it<T, R>
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
    \pre  `declval<T &>() -= declval<T>()` is well-formed.
    \pre  `declval<T &>() += declval<T>()` is well-formed.

    \param[in] x  The input value.

    \returns  `||x||_1 := Sum{i}( |c[i]| )`.
 */
template < typename T, std::size_t R >
T  taxi( complex_it<T, R> const &x )
{
    T  result{};

    for ( auto const &xx : x )
        ( xx < T{} ) ? ( result -= xx ) : ( result += xx );
    return result;
}

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
auto  abs( complex_it<T, R> const &x )
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
template < typename T, std::size_t R >
T  sup( complex_it<T, R> const &x )
{
    using std::abs;

    auto        xb = begin( x );
    T           result = abs( *xb++ );
    auto const  xe = end( x );

    while ( xe != xb )
        result = std::max( result, abs(*xb++) );
    return result;
}

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
auto  sgn( complex_it<T, R> const &x ) -> complex_it<T, R>
{ return x ? x / abs(x) : x; }


}  // namespace math
}  // namespace boost


//  Specializations from the STD namespace  ----------------------------------//

//! The standard namespace
namespace std
{
    //! Provides the number of components in a `complex_it` object.
    template < typename T, size_t R >
    class tuple_size< boost::math::complex_it<T, R> >
        : public integral_constant< size_t, boost::math::complex_it<T,
           R>::static_size >
    { };

    //! Provides the type of each component in a `complex_it` object.
    template < size_t I, typename T, size_t R >
    class tuple_element< I, boost::math::complex_it<T, R> >
    {
        static_assert( I < boost::math::complex_it<T, R>::static_size,
         "Index too large" );

    public:
        //! The component type, index-independent for complex numbers.
        typedef typename boost::math::complex_it<T, R>::value_type  type;
    };

}  // namespace std


#endif // BOOST_MATH_COMPLEX_IT_HPP
