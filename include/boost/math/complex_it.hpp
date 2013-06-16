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
#include <cstddef>
#include <iterator>
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
    \tparam Rank    The Cayley-Dickson construction level
 */
template < typename Number, std::size_t Rank >
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

    /** \brief  The lower (i.e. real-ward) decomposed half of this value

    When #rank is 0, `*this` is returned as a degenerate case.

    Cayley-Dickson construction builds a multi-valued scalar from composing two
    instances from the immediately-lower level.  This function decomposes a
    value to its two lower-ranked instances and returns the one made from the
    lower-numbered indices.

    The iteration-based class template version is `const` and returns a value,
    while the recursion-based class template versions return references, in
    `const` and mutable versions.

        \see  #upper_barrage

        \returns  `Lower( *this ) := { c[0]; c[1]; ...; c[2^(Rank - 1) - 1] }`
     */
    auto  lower_barrage() const -> barrage_type
    {
        barrage_type  result;

        std::copy( &c[0], &c[static_size / 2u + !rank], &result[0] );
        return result;
    }
    /** \brief  The upper (i.e. imaginary-ward) decomposed half of this value

    When #rank is 0, `*this` is returned as a degenerate case.

    Cayley-Dickson construction builds a multi-valued scalar from composing two
    instances from the immediately-lower level.  This function decomposes a
    value to its two lower-ranked instances and returns the one made from the
    higher-numbered indices.

    The iteration-based class template version is `const` and returns a value,
    while the recursion-based class template versions return references, in
    `const` and mutable versions.

        \see  #lower_barrage

        \returns  `Upper( *this ) := { c[2^(Rank - 1)]; ...; c[2^Rank - 1] }`
     */
    auto  upper_barrage() const -> barrage_type
    {
        barrage_type  result;

        std::copy( &c[static_size / 2u], &c[static_size], &result[0] );
        return result;
    }

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
