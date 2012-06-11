//  Boost.Complex Library, cd_hypercomplex_core.hpp header file  -------------//

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

    Contains the declarations and definitions of several class templates that
    model Cayley-Dickson hypercomplex numbers.  There are two pairs of opposing
    ways to model: iterative vs. recursive storage, and constructors vs.
    aggregate-initialization; the four combinations are modeled here.  All the
    applicable math operators are provided, plus some core condition function
    templates.
 */

#ifndef BOOST_MATH_CD_HYPERCOMPLEX_CORE_HPP
#define BOOST_MATH_CD_HYPERCOMPLEX_CORE_HPP

#include "boost/detail/index_tuple11.hpp"

#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>


namespace boost
{
namespace math
{


//  Forward declarations  ----------------------------------------------------//

//! Cayley-Dickson hypercomplex numbers, aggregate/iterative model.
template < typename CommutativeRing, std::size_t Rank >
struct cdh_complex_ai;

//! Cayley-Dickson hypercomplex numbers, aggregate/recursive model.
template < typename CommutativeRing, std::size_t Rank >
struct cdh_complex_ar;

//! Cayley-Dickson hypercomplex numbers, aggregate/recursive model, base case.
template < typename CommutativeRing >
struct cdh_complex_ar<CommutativeRing, 0>;

//template < typename CommutativeRing, std::size_t Rank >
//class cdh_complex_ci;

//template < typename CommutativeRing, std::size_t Rank >
//class cdh_complex_cr;

//template < typename CommutativeRing >
//class cdh_complex_cr<CommutativeRing, 0>;


//  Hypercomplex number class template definitions  --------------------------//

/** This model of C.D. hypercomplex numbers stores its components directly in a
    flat, top-level one-dimensional array.  There are no constructors; the type
    is an aggregate class.  The aggregate property may lead to time-savings when
    creating objects of this type, since there will be no initialization code
    run at declaration time (unless the component type has a non-trivial default
    constructor).

    This class will satisfy the standard-layout, trivially-copyable, trivial,
    and/or P.O.D. properties iff the component type also satisfies the
    corresponding property.  The type supports default-, copy-, and/or
    move-construction, destruction, copy- and/or move-assignment if the
    component type also supports the corresponding special operation.

    The component type should be a regular type, however its operations may
    cause condition events and/or irregular values on some combinations of
    operands.  Condition events are things like: divide-by-zero, overflow,
    underflow, and/or wrap-around.  Irregular values include: not-a-number,
    (positive) infinity, and negative infinity.  The component type should meet
    all expectations when these events and values are avoided, and may meet them
    even when not avoided (e.g. wrap-around, the infinities).  Operands with
    irregular values tend to be viral; continuing operations will result in
    answers with more irregular components.

    \pre  `CommutativeRing` has to meet the requirements given in section 26.2,
          *Numeric type requirements*, of the C++11 standard
          [numeric.requirements].
    \pre  `CommutativeRing` has to model a commutative ring as best as possible.
          - (Binary) `+` has to do the addition operation.
          - (Binary) `\*` has to do the multiplication operation.
          - (Unary) `-` has to generate the additive inverse.
          - Default/value-initialization has to return the additive identity
            (a.k.a. 0).
          - The multiplicative identity (a.k.a. 1) has to be reachable from the
            additive one via a single application of the `++` operator.
          - There should be a Boolean conversion operator, with 0 mapping to
            `false` and all other regular values to `true`.
          - The related operators; (unary) `+`, (binary) `-`, `==`, `!=`, `\<`,
            `\>`, `\<=`, `\>=`, `/`, `%` (optional), `++` (both), `--` (both),
            `+=`, `-=`, `\*=`, `%=` (optional), and `/=`; should be available
            and keep their conventional connections to each other and the
            primary operators.
    \pre  0 \<= `Rank` \< `std\::numeric_limits\<unsigned long\>\::digits`. 

    \tparam CommutativeRing  The component type.
    \tparam Rank             The rung of Cayley-Dickson construction used.  The
                             number of components is based on this.
 */
template < typename CommutativeRing, std::size_t Rank >
struct cdh_complex_ai
{
    // Template parameters
    /** \brief  Type for each component.

        This type is used to store each component of this multi-valued number.
        It is an alias for the first parameter of this class template.
     */
    typedef CommutativeRing       value_type;
    //! The level of the Cayley-Dickson construction modeled.
    static constexpr std::size_t  rank = Rank;

    // Support types & constants
    /** \brief  Type for the elements of the data array.

        This is an alias for the element type of the array storing the
        components of this multi-valued number.  It aliases the component type
        iff the array directly stores the data.
     */
    typedef value_type            element_type;
    //! The number of elements directly stored in the data array.
    static constexpr std::size_t  element_count = 1UL << rank;
    //! The number of components stored in the data array.
    static constexpr std::size_t  dimensions = element_count;

    // Member data
    /** \brief The component data.

        This member sub-object stores each component.  For iterative models,
        each array element is a component.  The first element is the real part
        and the later elements (if any) are the various **components** of the
        unreal part.
     */
    element_type  c[ element_count ];

    // Range-for/Iteration support
    //! The beginning of the immutable-access iteration range.
    constexpr
    auto  begin() const noexcept -> value_type const *;
    //! The beginning of the mutable-access iteration range.
    auto  begin() noexcept -> value_type *;
    //! The end of the mutable-access iteration range.
    auto  end() noexcept -> value_type *;
    //! The end of the immutable-access iteration range.
    constexpr
    auto  end() const noexcept -> value_type const *;

    //! The number of components contained
    constexpr
    auto  size() const noexcept -> std::size_t;

    // Operators
    //! Boolean conversion.
    explicit operator bool() const
     noexcept( std::is_nothrow_constructible<bool, value_type const &>::value );

    //! Cross-conversion: same-or-greater size, implicitly convertible types.
    template <
        typename U, std::size_t S,
        typename std::enable_if<
            std::is_convertible<CommutativeRing, U>::value && (S >= Rank)
        >::type*...
    >
    constexpr
    operator cdh_complex_ai<U, S>() const
     noexcept(
         std::is_nothrow_constructible<U, value_type const &>::value
         && (( S == rank ) || std::is_nothrow_default_constructible<U>::value)
         && std::is_nothrow_move_constructible<U>::value
     );
    //! Cross-conversion: lesser size and/or explicitly convertible types.
    template <
        typename U, std::size_t S,
        typename std::enable_if<
            !( std::is_convertible<CommutativeRing, U>::value && (S >= Rank) ),
            decltype(static_cast<U>( std::declval<CommutativeRing const &>() ))
        >::type*...
    >
    explicit constexpr
    operator cdh_complex_ai<U, S>() const
     noexcept(
         noexcept( static_cast<U>(std::declval<value_type const &>()) )
         && (( S <= rank ) || std::is_nothrow_default_constructible<U>::value)
         && std::is_nothrow_move_constructible<U>::value
     );
};

/** This partial specialization is the base case for the recursively-storing
    C.D. hypercomplex number class template.  It stores a single component
    (directly as a top-level array).  There are no constructors; the type is an
    aggregate class.  The aggregate property may lead to time-savings when
    creating objects of this type, since there will be no initialization code
    run at declaration time (unless the component type has a non-trivial default
    constructor).

    See the notes for the related `#cdh_complex_ar` class template for further
    information.

    \pre  `CommutativeRing` has the same pre-conditions as listed for the
          `#cdh_complex_ar` class template.  (The `Rank` parameter is fixed as
          `0` for this partial specialization.)

    \tparam CommutativeRing  The component type.
 */
template < typename CommutativeRing >
struct cdh_complex_ar<CommutativeRing, 0>
{
    // Template parameters
    //! \copydoc cdh_complex_ai::value_type
    typedef CommutativeRing       value_type;
    //! \copybrief cdh_complex_ai::rank
    static constexpr std::size_t  rank = 0;

    // Support types & constants
    //! \copydoc cdh_complex_ai::element_type
    typedef value_type            element_type;
    //! \copybrief cdh_complex_ai::element_count
    static constexpr std::size_t  element_count = 1;
    //! \copybrief cdh_complex_ai::dimensions
    static constexpr std::size_t  dimensions = element_count;

    // Member data
    /** \brief The component data.

        This member sub-object stores the component.  The recursive model's base
        case uses a single component as the sole array element.  That element is
        the **real** part.  The unreal part isn't supported.
     */
    element_type  r[ element_count ];

    // Range-for/Iteration support
    //! \copybrief cdh_complex_ai::begin()const
    constexpr
    auto  begin() const noexcept -> value_type const *;
    //! \copybrief cdh_complex_ai::begin()
    auto  begin() noexcept -> value_type *;
    //! \copybrief cdh_complex_ai::end()
    auto  end() noexcept -> value_type *;
    //! \copybrief cdh_complex_ai::end()const
    constexpr
    auto  end() const noexcept -> value_type const *;

    //! \copybrief cdh_complex_ai::size()const
    constexpr
    auto  size() const noexcept -> std::size_t;

    //! Iterate a function call over the component, mutable access.
    template < typename Function >
    void  iterate( Function &&f );
    //! Iterate a function call over the component, immutable access.
    template < typename Function >
    void  iterate( Function &&f ) const;

    // Operators
    //! \copybrief cdh_complex_ai::operator bool()const
    constexpr explicit
    operator bool() const
     noexcept( std::is_nothrow_constructible<bool, value_type const &>::value );

    //! Cross-conversion: same size, implicitly convertible types.
    template <
        typename U,
        typename std::enable_if<std::is_convertible<CommutativeRing, U>::value
        >::type*...
    >
    constexpr
    operator cdh_complex_ar<U, 0>() const
     noexcept(
         std::is_nothrow_constructible<U, value_type const &>::value
         && std::is_nothrow_move_constructible<U>::value
     );
    //! Cross-conversion: same size, explicitly convertible types.
    template <
        typename U,
        typename std::enable_if<
            not std::is_convertible<CommutativeRing, U>::value,
            decltype(static_cast<U>( std::declval<CommutativeRing const &>() ))
        >::type*...
    >
    explicit constexpr
    operator cdh_complex_ar<U, 0>() const
     noexcept(
         noexcept( static_cast<U>(std::declval<value_type const &>()) )
         && std::is_nothrow_move_constructible<U>::value
     );
    //! Cross-conversion: greater size, same or implicitly convertible types.
    template <
        typename U, std::size_t S,
        typename std::enable_if<std::is_convertible<CommutativeRing, U>::value
         && S>::type*...
    >
    constexpr
    operator cdh_complex_ar<U, S>() const
     noexcept(
         std::is_nothrow_constructible<U, value_type const &>::value
         && std::is_nothrow_default_constructible<U>::value
         && std::is_nothrow_move_constructible<U>::value
     );
    //! Cross-conversion: greater size, explicitly convertible types.
    template <
        typename U, std::size_t S,
        typename std::enable_if<
            !std::is_convertible<CommutativeRing, U>::value && S,
            decltype(static_cast<U>( std::declval<CommutativeRing const &>() ))
        >::type*...
    >
    explicit constexpr
    operator cdh_complex_ar<U, S>() const
     noexcept(
         noexcept( static_cast<U>(std::declval<value_type const &>()) )
         && std::is_nothrow_default_constructible<U>::value
         && std::is_nothrow_move_constructible<U>::value
     );
};

/** This model of C.D. hypercomplex numbers stores its components indirectly;
    the top-level array is two elements of the previous C.D. rung, resulting in
    the components being indirectly stored as a tree.  There are no constructors;
    the type is an aggregate class.  The aggregate property may lead to
    time-savings when creating objects of this type, since there will be no
    initialization code run at declaration time (unless the component type has a
    non-trivial default constructor).

    See the notes for the related `#cdh_complex_ar` class template for further
    information.

    \pre  `CommutativeRing` and `Rank` have the same pre-conditions as listed
          for the `#cdh_complex_ar` class template.

    \tparam CommutativeRing  The component type.
    \tparam Rank             The rung of Cayley-Dickson construction used.  The
                             number of components is based on this.
 */
template < typename CommutativeRing, std::size_t Rank >
struct cdh_complex_ar
{
    // Template parameters
    //! \copydoc cdh_complex_ai::value_type
    typedef CommutativeRing       value_type;
    //! \copybrief cdh_complex_ai::rank
    static constexpr std::size_t  rank = Rank;

    // Support types & constants
    //! \copydoc cdh_complex_ai::element_type
    typedef cdh_complex_ar<value_type, rank - 1u>  element_type;
    //! \copybrief cdh_complex_ai::element_count
    static constexpr std::size_t                   element_count = 2;
    //! \copybrief cdh_complex_ai::dimensions
    static constexpr std::size_t                   dimensions = element_count *
     element_type::dimensions;

    // Member data
    /** \brief The component data.

        This member sub-object indirectly stores each component.  For the
        recursive model's case, each of the two array elements stores half of
        the components, nicknamed a **barrage**.  The real part is the frontmost
        component of the first element, while the unreal part is (split between
        the back components of the first element and) the second element.
     */
    element_type  b[ element_count ];

    // Iteration support
    //! \copybrief cdh_complex_ai::size()const
    constexpr
    auto  size() const noexcept -> std::size_t;

    //! Iterate a function call over the components, mutable access.
    template < typename Function >
    void  iterate( Function &&f );
    //! Iterate a function call over the components, immutable access.
    template < typename Function >
    void  iterate( Function &&f ) const;

    // Operators
    //! \copybrief cdh_complex_ai::operator bool()const
    constexpr explicit
    operator bool() const
     noexcept( std::is_nothrow_constructible<bool, value_type const &>::value );

    //! Cross-conversion: same size, implicitly convertible types.
    template <
        typename U,
        typename std::enable_if<std::is_convertible<CommutativeRing, U>::value
        >::type*...
    >
    constexpr
    operator cdh_complex_ar<U, rank>() const
     noexcept(
         std::is_nothrow_constructible<U, value_type const &>::value
         && std::is_nothrow_move_constructible<U>::value
     );
    //! Cross-conversion: same size, explicitly convertible types.
    template <
        typename U,
        typename std::enable_if<
            not std::is_convertible<CommutativeRing, U>::value,
            decltype( static_cast<U>(std::declval<CommutativeRing const &>()) )
        >::type*...
    >
    explicit constexpr
    operator cdh_complex_ar<U, rank>() const
     noexcept(
         noexcept( static_cast<U>(std::declval<value_type const &>()) )
         && std::is_nothrow_move_constructible<U>::value
     );
    //! Cross-conversion: greater size, same or implicitly convertible types.
    template <
        typename U, std::size_t S,
        typename std::enable_if<std::is_convertible<CommutativeRing, U>::value
         && (S > Rank)>::type*...
    >
    constexpr
    operator cdh_complex_ar<U, S>() const
     noexcept(
         std::is_nothrow_constructible<U, value_type const &>::value
         && std::is_nothrow_default_constructible<U>::value
         && std::is_nothrow_move_constructible<U>::value
     );
    //! Cross-conversion: greater size, explicitly convertible types.
    template <
        typename U, std::size_t S,
        typename std::enable_if<
            !std::is_convertible<CommutativeRing, U>::value && (S > Rank),
            decltype( static_cast<U>(std::declval<CommutativeRing const &>()) )
        >::type*...
    >
    explicit constexpr
    operator cdh_complex_ar<U, S>() const
     noexcept(
         noexcept( static_cast<U>(std::declval<value_type const &>()) )
         && std::is_nothrow_default_constructible<U>::value
         && std::is_nothrow_move_constructible<U>::value
     );
    //! Cross-conversion: immediately-smaller size, same type.
    explicit constexpr
    operator cdh_complex_ar<value_type, rank - 1u>() const
     noexcept( std::is_nothrow_move_constructible<value_type>::value );
    //! Cross-conversion: smaller sizes, same type.
    template <std::size_t S, typename std::enable_if<(S < Rank - 1u)>::type*...>
    explicit constexpr
    operator cdh_complex_ar<value_type, S>() const
     noexcept( std::is_nothrow_move_constructible<value_type>::value );
    //! Cross-conversion: lesser size, implicitly convertible types.
    template <
        typename U, std::size_t S,
        typename std::enable_if<
            not std::is_same<CommutativeRing, U>::value
            && std::is_convertible<CommutativeRing, U>::value
            && (S < Rank)
        >::type*...
    >
    explicit constexpr
    operator cdh_complex_ar<U, S>() const
     noexcept(
         std::is_nothrow_constructible<U, value_type const &>::value
         && std::is_nothrow_move_constructible<U>::value
     );
    //! Cross-conversion: lesser size, explicitly convertible types.
    template <
        typename U, std::size_t S,
        typename std::enable_if<
            !std::is_convertible<CommutativeRing, U>::value && (S < Rank),
            decltype( static_cast<U>(std::declval<CommutativeRing const &>()) )
        >::type*...
    >
    explicit constexpr
    operator cdh_complex_ar<U, S>() const
     noexcept(
         noexcept( static_cast<U>(std::declval<value_type const &>()) )
         && std::is_nothrow_move_constructible<U>::value
     );
};


//  Hypercomplex number class-static member definitions  ---------------------//

/** This member indicates the number of times Cayley-Dickson construction has
    been applied to the component type.  A rank of zero just makes this type a
    wrapper around `#value_type`, and each subsequent rank doubles the number
    of components.  It is an alias for the second parameter of this class
    template.
 */
template < typename CommutativeRing, std::size_t Rank >
constexpr
std::size_t  cdh_complex_ai<CommutativeRing, Rank>::rank;

/** Since this type uses the iterative model for component storage, this value
    always equals `#dimensions`.
 */
template < typename CommutativeRing, std::size_t Rank >
constexpr
std::size_t  cdh_complex_ai<CommutativeRing, Rank>::element_count;

/** Cayley-Dickson construction doubles the number of components between levels,
    and we always start with one component at the lowest level (i.e. a simple
    wrapper class); therefore this value is two to the power of `#rank`.
 */
template < typename CommutativeRing, std::size_t Rank >
constexpr
std::size_t  cdh_complex_ai<CommutativeRing, Rank>::dimensions;

//! \copydetails cdh_complex_ai::rank
template < typename CommutativeRing >
constexpr
std::size_t  cdh_complex_ar<CommutativeRing, 0>::rank;

/** Since this type is a base case for the recursive model for component
    storage, this value always equals 1 (which is also the value of
    `#dimensions`).
 */
template < typename CommutativeRing >
constexpr
std::size_t  cdh_complex_ar<CommutativeRing, 0>::element_count;

//! \copydetails cdh_complex_ai::dimensions
template < typename CommutativeRing >
constexpr
std::size_t  cdh_complex_ar<CommutativeRing, 0>::dimensions;

//! \copydetails cdh_complex_ai::rank
template < typename CommutativeRing, std::size_t Rank >
constexpr
std::size_t  cdh_complex_ar<CommutativeRing, Rank>::rank;

/** Since this type uses the recursive model for component storage, this value
    always equals 2, representing the doubling of components at each step of
    Cayley-Dickson construction.  (Therefore, `#dimensions` is doubled between
    steps.)
 */
template < typename CommutativeRing, std::size_t Rank >
constexpr
std::size_t  cdh_complex_ar<CommutativeRing, Rank>::element_count;

//! \copydetails cdh_complex_ai::dimensions
template < typename CommutativeRing, std::size_t Rank >
constexpr
std::size_t  cdh_complex_ar<CommutativeRing, Rank>::dimensions;


//  Hypercomplex number class iteration methods definitions  -----------------//

/** Iteration proceeds from the real component to the last indexed component of
    the highest supported rank.  The iterator type is a pointer, so all the
    random-access iterator semantics apply.

    The various sequence/container type-aliases aren't supported, though, just
    the minimum needed to enable objects of this type for the ranged-based `for`
    statement.

    \returns  A pointer to the real component, immutable access.
 */
template < typename CommutativeRing, std::size_t Rank >
inline constexpr
auto
cdh_complex_ai<CommutativeRing, Rank>::begin() const noexcept
 -> value_type const *
{ return &c[0]; }

/** \see begin()const

    \returns  A pointer to the real component, mutable access.
 */
template < typename CommutativeRing, std::size_t Rank >
inline
auto
cdh_complex_ai<CommutativeRing, Rank>::begin() noexcept -> value_type *
{ return &c[0]; }

/** \see begin()const

    \returns  A pointer to one-past the highest-indexed component of the highest
              supported rank, mutable access.
 */
template < typename CommutativeRing, std::size_t Rank >
auto
cdh_complex_ai<CommutativeRing, Rank>::end() noexcept -> value_type *
{ return &c[dimensions]; }

/** \see begin()const

    \returns  A pointer to one-past the highest-indexed component of the highest
              supported rank, immutable access.
 */
template < typename CommutativeRing, std::size_t Rank >
inline constexpr
auto
cdh_complex_ai<CommutativeRing, Rank>::end() const noexcept
 -> value_type const *
{ return &c[dimensions]; }

/** This member function provides access to the number of components in a way
    that is compatible with a container interface, and you just need an object's
    name instead of the type's name.

    \returns  `#dimensions`.
 */
template < typename CommutativeRing, std::size_t Rank >
inline constexpr
auto
cdh_complex_ai<CommutativeRing, Rank>::size() const noexcept -> std::size_t
{ return dimensions; }

//! \copydetails cdh_complex_ai::begin()const
template < typename CommutativeRing >
inline constexpr
auto
cdh_complex_ar<CommutativeRing, 0>::begin() const noexcept -> value_type const *
{ return &r[0]; }

//! \copydetails cdh_complex_ai::begin()
template < typename CommutativeRing >
inline
auto
cdh_complex_ar<CommutativeRing, 0>::begin() noexcept -> value_type *
{ return &r[0]; }

//! \copydetails cdh_complex_ai::end()
template < typename CommutativeRing >
auto
cdh_complex_ar<CommutativeRing, 0>::end() noexcept -> value_type *
{ return &r[dimensions]; }

//! \copydetails cdh_complex_ai::end()const
template < typename CommutativeRing >
inline constexpr
auto
cdh_complex_ar<CommutativeRing, 0>::end() const noexcept -> value_type const *
{ return &r[dimensions]; }

//! \copydetails cdh_complex_ai::size()const
template < typename CommutativeRing >
inline constexpr
auto
cdh_complex_ar<CommutativeRing, 0>::size() const noexcept -> std::size_t
{ return dimensions; }

//! \copydetails cdh_complex_ai::size()const
template < typename CommutativeRing, std::size_t Rank >
inline constexpr
auto
cdh_complex_ar<CommutativeRing, Rank>::size() const noexcept -> std::size_t
{ return dimensions; }

/** Calls the given function, using a (mutable l-value) reference to the stored
    component as the sole parameter.

    \tparam Function  The type of `f`.

    \param[in] f  The function to be called.  It may be a function pointer, a
                  function object, or a lambda.

    \throws  Whatever calling `f` with the stored component can throw.
 */
template < typename CommutativeRing >
template < typename Function >
inline
void  cdh_complex_ar<CommutativeRing, 0>::iterate( Function &&f )
{ std::forward<Function>(f)(r[ 0 ]); }

/** Calls the given function, using a (`const` l-value) reference to the stored
    component as the sole parameter.

    \tparam Function  The type of `f`.

    \param[in] f  The function to be called.  It may be a function pointer, a
                  function object, or a lambda.

    \throws  Whatever calling `f` with the stored component can throw.
 */
template < typename CommutativeRing >
template < typename Function >
inline
void  cdh_complex_ar<CommutativeRing, 0>::iterate( Function &&f ) const
{ std::forward<Function>(f)(r[ 0 ]); }

/** Calls the given function, using (mutable l-value) references to each of the
    stored components, in order of increasing index, as the sole parameter.

    \tparam Function  The type of `f`.

    \param[in] f  The function to be called.  It may be a function pointer, a
                  function object, or a lambda.

    \throws  Whatever calling `f` with a component can throw.
 */
template < typename CommutativeRing, std::size_t Rank >
template < typename Function >
void  cdh_complex_ar<CommutativeRing, Rank>::iterate( Function &&f )
{
    b[ 0 ].iterate( std::forward<Function>(f) );
    b[ 1 ].iterate( std::forward<Function>(f) );
}

/** Calls the given function, using a (`const` l-value) reference to each of the
    stored components, in order of increasing index, as the sole parameter.

    \tparam Function  The type of `f`.

    \param[in] f  The function to be called.  It may be a function pointer, a
                  function object, or a lambda.

    \throws  Whatever calling `f` with a component can throw.
 */
template < typename CommutativeRing, std::size_t Rank >
template < typename Function >
void  cdh_complex_ar<CommutativeRing, Rank>::iterate( Function &&f ) const
{
    b[ 0 ].iterate( std::forward<Function>(f) );
    b[ 1 ].iterate( std::forward<Function>(f) );
}


//  Hypercomplex number class operator member definitions  -------------------//

/** This member operator provides a contextual conversion to `bool`, enabling
    this type to be used in language constructs that use objects/values in
    Boolean contexts.  This conversion is marked `explicit` so it can be used
    for language constructs without participating in general numeric contexts.

    Example language constructs needing contextual conversion to `bool` include:
    the built-in `!`, `&&`, and `||` operators; the first part of the `?:`
    operator (but not the other two parts); the `noexcept` operator; and the
    test-condition part of `if`, `while`, and regular-`for` statements.

    The result can be expressed as:
    - Single: `static_cast\<bool\>( R )`
    - Iterative: `c0 || c1 || ... || cN`
    - Recursive: `B0 || B1`

    \pre  `#value_type` has to support Boolean conversion.  The key conversion
          operation may be marked `explicit`.

    \throws  Anything that Boolean conversion for `value_type` may throw.

    \returns  *Norm*(`\*this`) != 0.  (Provided that `value_type` doesn't have
              zero divisors; the squaring of any regular value maps to the
              positive cone; and the actual computation doesn't involve
              irregular values, doesn't produce irregular values, or have
              condition events occur.)
 */
template < typename CommutativeRing, std::size_t Rank >
cdh_complex_ai<CommutativeRing, Rank>::operator bool() const
 noexcept( std::is_nothrow_constructible<bool, value_type const &>::value )
{
    for ( auto const &x : c )
        if ( x )
            return true;
    return false;
}

//! \copydetails cdh_complex_ai::operator bool()const
template < typename CommutativeRing >
inline constexpr
cdh_complex_ar<CommutativeRing, 0>::operator bool() const
 noexcept( std::is_nothrow_constructible<bool, value_type const &>::value )
{ return static_cast<bool>(r[ 0 ]); }

//! \copydetails cdh_complex_ai::operator bool()const
template < typename CommutativeRing, std::size_t Rank >
inline constexpr
cdh_complex_ar<CommutativeRing, Rank>::operator bool() const
 noexcept( std::is_nothrow_constructible<bool, value_type const &>::value )
{ return b[0] || b[1]; }


//  Hypercomplex number status function definitions  -------------------------//

/** \brief  Find the minimum Cayley-Dickson construction rung required for the
            given value.

    Inspects the given value for the highest index that corresponds to a nonzero
    component, then returns the smallest rung that supports that index.  Returns
    0 when the value (i.e. all components) is zero.

    The result can be expressed as:
    - Single: `0`
    - Iterative: adapt recursive version
    - Recursive: `!B1 ? dynamic_rank(B0) : static_rank([B0;B1])`

    \pre  `T` must support Boolean conversion.  (It may be `explicit`.)

    \param[in] x  The value to inspect.

    \throws  Anything that Boolean conversion for `value_type` may throw.

    \returns  0 if `x` is zero or otherwise purely real; else `n`, where it is
              the smallest integer that 2<sup>n</sup> - 1 is greater than or
              equal to the largest index with a nonzero component.
 */
template < typename T, std::size_t R >
auto  dynamic_rank( cdh_complex_ai<T, R> const &x )
 noexcept( std::is_nothrow_constructible<bool, T const &>::value )
 -> std::size_t
{
    auto  r = R;

    for ( auto  d = x.size() ; r ; --r, d /= 2 )
        for ( auto  i = d / 2 ; i < d ; ++i )
            if ( x.c[i] )
                return r;
    return r;
}

//! \overload
template < typename T >
inline constexpr
auto  dynamic_rank( cdh_complex_ar<T, 0> const & ) noexcept -> std::size_t
{ return 0; }

//! \overload
template < typename T, std::size_t R >
inline constexpr
auto  dynamic_rank( cdh_complex_ar<T, R> const &x )
 noexcept( std::is_nothrow_constructible<bool, T const &>::value )
 -> std::size_t
{ return x.b[1] ? R : dynamic_rank(x.b[ 0 ]); }


//  Hypercomplex number access function definitions  -------------------------//

/** \brief  Extract the component with the given index from the given value.

    This function and its overloads provides an interface compatible with
    `std::tuple` for accessing elements.

    \pre  0 \<= `I` \< `decltype(x)\::dimensions`.

    \tparam I  The index of the desired component.

    \param[in] x  A reference to the hypercomplex number object to be inspected.

    \throws  Nothing.  However, a compile-time error will occur if `I` exceeds
             its bounds.

    \returns  A reference to the desired component.  It will match the r-value
              vs. l-value vs. `const` l-value state as the source object.
 */
template < std::size_t I, typename T, std::size_t R >
inline constexpr
auto  get( cdh_complex_ai<T, R> const &x ) noexcept -> T const &
{
    static_assert( I < cdh_complex_ai<T,R>::dimensions, "index out of bounds" );

    return x.c[ I ];
}

//! \overload
template < std::size_t I, typename T, std::size_t R >
inline constexpr
auto  get( cdh_complex_ai<T, R> &x ) noexcept -> T &
{
    static_assert( I < cdh_complex_ai<T,R>::dimensions, "index out of bounds" );

    return x.c[ I ];
}

//! \overload
template < std::size_t I, typename T, std::size_t R >
inline constexpr
auto  get( cdh_complex_ai<T, R> &&x ) noexcept -> T &&
{
    static_assert( I < cdh_complex_ai<T,R>::dimensions, "index out of bounds" );

    return static_cast<T &&>( x.c[I] );
}

//! \overload
template < std::size_t I, typename T >
inline constexpr
auto  get( cdh_complex_ar<T, 0> const &x ) noexcept -> T const &
{
    static_assert( !I, "index out of bounds" );

    return x.r[ I ];
}

//! \overload
template < std::size_t I, typename T >
inline constexpr
auto  get( cdh_complex_ar<T, 0> &x ) noexcept -> T &
{
    static_assert( !I, "index out of bounds" );

    return x.r[ I ];
}

//! \overload
template < std::size_t I, typename T >
inline constexpr
auto  get( cdh_complex_ar<T, 0> &&x ) noexcept -> T &&
{
    static_assert( !I, "index out of bounds" );

    return static_cast<T &&>( x.r[I] );
}

//! \overload
template < std::size_t I, typename T, std::size_t R >
inline constexpr
auto  get( cdh_complex_ar<T, R> const &x ) noexcept -> T const &
{
    static_assert( I < cdh_complex_ar<T,R>::dimensions, "index out of bounds" );

    return get<I % (cdh_complex_ar<T,R>::dimensions / 2u)>( x.b[I >= ( x.size()
     / 2u )] );
}

//! \overload
template < std::size_t I, typename T, std::size_t R >
inline constexpr
auto  get( cdh_complex_ar<T, R> &x ) noexcept -> T &
{
    static_assert( I < cdh_complex_ar<T,R>::dimensions, "index out of bounds" );

    return get<I % (cdh_complex_ar<T,R>::dimensions / 2u)>( x.b[I >= ( x.size()
     / 2u )] );
}

//! \overload
template < std::size_t I, typename T, std::size_t R >
inline constexpr
auto  get( cdh_complex_ar<T, R> &&x ) noexcept -> T &&
{ return static_cast<T &&>(get<I>( x )); }


//  Implementation details  --------------------------------------------------//

//! \cond
namespace detail
{
    // Create a `cdh_complex_ai` object from an initializer.
    template < typename T, std::size_t R, class Source, std::size_t ...Indices >
    inline constexpr
    auto  implicitly_convert_to_ai( Source const &x,
     boost::detail::index_tuple<Indices...> )
     noexcept(
         std::is_nothrow_constructible<
             T,
             typename Source::value_type const &
         >::value
         && ( (sizeof...( Indices ) == ( 1UL << R ))
             || std::is_nothrow_default_constructible<T>::value
         )
         && std::is_nothrow_move_constructible<T>::value
     )
     -> boost::math::cdh_complex_ai<T, R>
    { return {{ boost::math::get<Indices>(x)... }}; }

    // Create a `cdh_complex_ai` object from an initializer, but using
    // explicit conversion.
    template < typename T, std::size_t R, class Source, std::size_t ...Indices >
    inline constexpr
    auto  explicitly_convert_to_ai( Source const &x,
     boost::detail::index_tuple<Indices...> )
     noexcept(
         noexcept(
             static_cast<T>(std::declval<typename Source::value_type const &>())
         )
         && ( (sizeof...( Indices ) >= ( 1UL << R ))
             || std::is_nothrow_default_constructible<T>::value
         )
         && std::is_nothrow_move_constructible<T>::value
     )
     -> boost::math::cdh_complex_ai<T, R>
    { return {{ static_cast<T>(boost::math::get<Indices>( x ))... }}; }

}  // namespace detail
//! \endcond


//  Hypercomplex number class, more operator member definitions  -------------//

/** Convert this object to an object of another instantiation of this class
    template.  This particular member function handles the cases of implicit
    conversion can be mechanically and logically done:  when the destination
    type has at least as many components as the source type, and the source
    component type can be implicitly converted to the destination component
    type.

    \pre  `U` can be implicitly converted to `#value_type`.
    \pre  `S` \>= `#rank`.

    \tparam U  The component type of the returned value.
    \tparam S  The rank (i.e. Cayley-Dickson rung) of the returned value.

    \throws  Whatever conversions from `value_type` to `U` may throw, plus any
             exceptions from returning `U` objects and/or creating default
             value-initialized `U` objects.

    \returns  An object `x` such that `get\<I\>(\*this) == get\<I\>(x)` for all
              `I` values the source and destination have in common.  If the
              destination value has more components, then all those later
              components have values of zero.
 */
template < typename CommutativeRing, std::size_t Rank >
template <
    typename U, std::size_t S,
    typename std::enable_if<
        std::is_convertible<CommutativeRing, U>::value && (S >= Rank)
    >::type*...
>
inline constexpr
cdh_complex_ai<CommutativeRing, Rank>::operator cdh_complex_ai<U, S>() const
 noexcept(
     std::is_nothrow_constructible<U, value_type const &>::value
     && (( S == rank ) || std::is_nothrow_default_constructible<U>::value)
     && std::is_nothrow_move_constructible<U>::value
 )
{
    return detail::implicitly_convert_to_ai<U, S>( *this,
     boost::detail::make_indices<dimensions>() );
}

/** Convert this object to an object of another instantiation of this class
    template.  This particular member function handles any destination type
    that is convertible from the source type, even if by explicit conversion,
    and can handle any difference in component length.  (This function blocks
    use by destination types that are covered by the implicit-conversion
    operator.)  The conversions done by this function are marked explicit.  If
    the destination type uses fewer components than the source, then only the
    ones with a lower index are kept. 

    \pre  `U` can be explicitly converted to `#value_type`.

    \tparam U  The component type of the returned value.
    \tparam S  The rank (i.e. Cayley-Dickson rung) of the returned value.

    \throws  Whatever conversions from `value_type` to `U` may throw, plus any
             exceptions from returning `U` objects and/or creating default
             value-initialized `U` objects.

    \returns  An object `x` such that `get\<I\>(\*this) == get\<I\>(x)` for all
              `I` values the source and destination have in common.  If the
              destination value has more components, then all those later
              components have values of zero.
 */
template < typename CommutativeRing, std::size_t Rank >
template <
    typename U, std::size_t S,
    typename std::enable_if<
        !( std::is_convertible<CommutativeRing, U>::value && (S >= Rank) ),
        decltype(static_cast<U>( std::declval<CommutativeRing const &>() ))
    >::type*...
>
inline constexpr
cdh_complex_ai<CommutativeRing, Rank>::operator cdh_complex_ai<U, S>() const
 noexcept(
     noexcept( static_cast<U>(std::declval<value_type const &>()) )
     && (( S <= rank ) || std::is_nothrow_default_constructible<U>::value)
     && std::is_nothrow_move_constructible<U>::value
 )
{
    return detail::explicitly_convert_to_ai<U, S>( *this,
     boost::detail::make_indices<(S < Rank) ? (1UL << S) : dimensions>() );
}

/** Convert this object to an object of another instantiation of this class
    template, where both the source and destination types support exactly one
    component.  This particular member function works as an implicit conversion,
    and therefore requires the source type to be able to implicitly convert to
    the destination type.

    \pre  `U` can be implicitly converted to `#value_type`.

    \tparam U  The component type of the returned value.

    \throws  Whatever conversions from `value_type` to `U` may throw, plus any
             exceptions from returning `U` objects.

    \returns  An object `x` such that `get\<0\>(\*this) == get\<0\>(x)`.
 */
template < typename CommutativeRing >
template <
    typename U,
    typename std::enable_if<std::is_convertible<CommutativeRing, U>::value
    >::type*...
>
inline constexpr
cdh_complex_ar<CommutativeRing, 0>::operator cdh_complex_ar<U, 0>() const
 noexcept(
     std::is_nothrow_constructible<U, value_type const &>::value
     && std::is_nothrow_move_constructible<U>::value
 )
{ return {{ r[0] }}; }

/** Convert this object to an object of another instantiation of this class
    template, where both the source and destination types support exactly one
    component.  This particular member function works as an explicit conversion,
    and therefore requires the source type to be able to explicitly convert to
    the destination type.

    \pre  `U` can be explicitly converted to `#value_type`.

    \tparam U  The component type of the returned value.

    \throws  Whatever conversions from `value_type` to `U` may throw, plus any
             exceptions from returning `U` objects.

    \returns  An object `x` such that `get\<0\>(\*this) == get\<0\>(x)`.
 */
template < typename CommutativeRing >
template <
    typename U,
    typename std::enable_if<
        not std::is_convertible<CommutativeRing, U>::value,
        decltype(static_cast<U>( std::declval<CommutativeRing const &>() ))
    >::type*...
>
inline constexpr
cdh_complex_ar<CommutativeRing, 0>::operator cdh_complex_ar<U, 0>() const
 noexcept(
     noexcept( static_cast<U>(std::declval<value_type const &>()) )
     && std::is_nothrow_move_constructible<U>::value
 )
{ return {{ static_cast<U>(r[ 0 ]) }}; }

/** Convert this object to an object of another instantiation of this class
    template, where the source type supports a single component and the
    destination type supports mutliple components.  This particular member
    function works as an implicit conversion, and therefore requires the source
    type to be able to implicitly convert to the destination type.

    \pre  `U` can be implicitly converted to `#value_type`.
    \pre  `S` \>= `#rank`.

    \tparam U  The component type of the returned value.
    \tparam S  The rank (i.e. Cayley-Dickson rung) of the returned value.

    \throws  Whatever conversions from `value_type` to `U` may throw, plus any
             exceptions from returning `U` objects and/or creating default
             value-initialized `U` objects.

    \returns  An object `x` such that `get\<0\>(\*this) == get\<0\>(x)`  and
              later components of the destination have values of zero.
 */
template < typename CommutativeRing >
template <
    typename U, std::size_t S,
    typename std::enable_if<std::is_convertible<CommutativeRing, U>::value &&
     S>::type*...
>
inline constexpr
cdh_complex_ar<CommutativeRing, 0>::operator cdh_complex_ar<U, S>() const
 noexcept(
     std::is_nothrow_constructible<U, value_type const &>::value
     && std::is_nothrow_default_constructible<U>::value
     && std::is_nothrow_move_constructible<U>::value
 )
{ return {{ *this, {} }}; }

/** Convert this object to an object of another instantiation of this class
    template, where the source type supports a single component and the
    destination type supports mutliple components.  This particular member
    function works as an explicit conversion, and therefore requires the source
    type to be able to explicitly convert to the destination type.

    \pre  `U` can be explicitly converted to `#value_type`.
    \pre  `S` \>= `#rank`.

    \tparam U  The component type of the returned value.
    \tparam S  The rank (i.e. Cayley-Dickson rung) of the returned value.

    \throws  Whatever conversions from `value_type` to `U` may throw, plus any
             exceptions from returning `U` objects and/or creating default
             value-initialized `U` objects.

    \returns  An object `x` such that `get\<0\>(\*this) == get\<0\>(x)` and
              later components of the destination have values of zero.
 */
template < typename CommutativeRing >
template <
    typename U, std::size_t S,
    typename std::enable_if<
        !std::is_convertible<CommutativeRing, U>::value && S,
        decltype(static_cast<U>( std::declval<CommutativeRing const &>() ))
    >::type*...
>
inline constexpr
cdh_complex_ar<CommutativeRing, 0>::operator cdh_complex_ar<U, S>() const
 noexcept(
     noexcept( static_cast<U>(std::declval<value_type const &>()) )
     && std::is_nothrow_default_constructible<U>::value
     && std::is_nothrow_move_constructible<U>::value
 )
{ return {{ this->operator cdh_complex_ar<U, S - 1u>(), {} }}; }

/** Convert this object to an object of another instantiation of this class
    template, where both the source and destination types support the same
    number of components.  This particular member function works as an implicit
    conversion, and therefore requires the source type to be able to implicitly
    convert to the destination type.

    \pre  `U` can be implicitly converted to `#value_type`.

    \tparam U  The component type of the returned value.

    \throws  Whatever conversions from `value_type` to `U` may throw, plus any
             exceptions from returning `U` objects.

    \returns  An object `x` such that `get\<I\>(\*this) == get\<I\>(x)` for all
              valid `I` values.
 */
template < typename CommutativeRing, std::size_t Rank >
template <
    typename U,
    typename std::enable_if<std::is_convertible<CommutativeRing, U>::value
    >::type*...
>
inline constexpr
cdh_complex_ar<CommutativeRing, Rank>::operator cdh_complex_ar<U, rank>() const
 noexcept(
     std::is_nothrow_constructible<U, value_type const &>::value
     && std::is_nothrow_move_constructible<U>::value
 )
{ return {{ b[0], b[1] }}; }

/** Convert this object to an object of another instantiation of this class
    template, where both the source and destination types support the same
    number of components.  This particular member function works as an explicit
    conversion, and therefore requires the source type to be able to explicitly
    convert to the destination type.

    \pre  `U` can be explicitly converted to `#value_type`.

    \tparam U  The component type of the returned value.

    \throws  Whatever conversions from `value_type` to `U` may throw, plus any
             exceptions from returning `U` objects.

    \returns  An object `x` such that `get\<I\>(\*this) == get\<I\>(x)` for all
              valid `I` values.
 */
template < typename CommutativeRing, std::size_t Rank >
template <
    typename U,
    typename std::enable_if<
        not std::is_convertible<CommutativeRing, U>::value,
        decltype( static_cast<U>(std::declval<CommutativeRing const &>()) )
    >::type*...
>
inline constexpr
cdh_complex_ar<CommutativeRing, Rank>::operator cdh_complex_ar<U, rank>() const
 noexcept(
     noexcept( static_cast<U>(std::declval<value_type const &>()) )
     && std::is_nothrow_move_constructible<U>::value
 )
{
    return { {b[0].operator cdh_complex_ar<U, rank - 1u>(),
     b[1].operator cdh_complex_ar<U, rank - 1u>()} };
}

/** Convert this object to an object of another instantiation of this class
    template, where the source type supports fewer (but still multiple)
    components than the destination type.  This particular member function works
    as an implicit conversion, and therefore requires the source type to be able
    to implicitly convert to the destination type.

    \pre  `U` can be implicitly converted to `#value_type`.
    \pre  `S` \>= `#rank`.

    \tparam U  The component type of the returned value.
    \tparam S  The rank (i.e. Cayley-Dickson rung) of the returned value.

    \throws  Whatever conversions from `value_type` to `U` may throw, plus any
             exceptions from returning `U` objects and/or creating default
             value-initialized `U` objects.

    \returns  An object `x` such that `get\<I\>(\*this) == get\<I\>(x)` for all
              `I` values the source and destination have in common.  Later
              components of the destination have values of zero.
 */
template < typename CommutativeRing, std::size_t Rank >
template <
    typename U, std::size_t S,
    typename std::enable_if<std::is_convertible<CommutativeRing, U>::value &&
     (S > Rank)>::type*...
>
inline constexpr
cdh_complex_ar<CommutativeRing, Rank>::operator cdh_complex_ar<U, S>() const
 noexcept(
     std::is_nothrow_constructible<U, value_type const &>::value
     && std::is_nothrow_default_constructible<U>::value
     && std::is_nothrow_move_constructible<U>::value
 )
{ return {{ *this, {} }}; }

/** Convert this object to an object of another instantiation of this class
    template, where the source type supports fewer (but still multiple)
    components than the destination type.  This particular member function works
    as an explicit conversion, and therefore requires the source type to be able
    to explicitly convert to the destination type.

    \pre  `U` can be explicitly converted to `#value_type`.
    \pre  `S` \>= `#rank`.

    \tparam U  The component type of the returned value.
    \tparam S  The rank (i.e. Cayley-Dickson rung) of the returned value.

    \throws  Whatever conversions from `value_type` to `U` may throw, plus any
             exceptions from returning `U` objects and/or creating default
             value-initialized `U` objects.

    \returns  An object `x` such that `get\<I\>(\*this) == get\<I\>(x)` for all
              `I` values the source and destination have in common.  Later
              components of the destination have values of zero.
 */
template < typename CommutativeRing, std::size_t Rank >
template <
    typename U, std::size_t S,
    typename std::enable_if<
        !std::is_convertible<CommutativeRing, U>::value && (S > Rank),
        decltype( static_cast<U>(std::declval<CommutativeRing const &>()) )
    >::type*...
>
inline constexpr
cdh_complex_ar<CommutativeRing, Rank>::operator cdh_complex_ar<U, S>() const
 noexcept(
     noexcept( static_cast<U>(std::declval<value_type const &>()) )
     && std::is_nothrow_default_constructible<U>::value
     && std::is_nothrow_move_constructible<U>::value
 )
{ return {{ this->operator cdh_complex_ar<U, S - 1u>(), {} }}; }

/** Convert this object to an object of another instantiation of this class
    template, with the source and destination using the same component type, but
    the destination is one rung lower in Cayley-Dickson construction.  Although
    same-type (the component type) conversion is implicit, this member function
    is marked explicit since there is a loss of information.  Only components of
    the lower half of valid indices are kept.

    \throws  Any exceptions from returning `#value_type` objects.

    \returns  The lower-barrage of `\*this`.
 */
template < typename CommutativeRing, std::size_t Rank >
inline constexpr
cdh_complex_ar<CommutativeRing, Rank>::operator cdh_complex_ar<value_type,
 rank - 1u>()
 const noexcept( std::is_nothrow_move_constructible<value_type>::value )
{ return b[0]; }

/** Convert this object to an object of another instantiation of this class
    template, with the source and destination using the same component type, but
    the destination is multiple rungs lower in Cayley-Dickson construction.
    Although same-type (the component type) conversion is implicit, this member
    function is marked explicit since there is a loss of information.  Only
    components of the lowest section of valid indices are kept.

    \throws  Any exceptions from returning `#value_type` objects.

    \returns  The lowest section of the lower-barrage of `\*this`.
 */
template < typename CommutativeRing, std::size_t Rank >
template < std::size_t S, typename std::enable_if<(S < Rank - 1u)>::type*... >
inline constexpr
cdh_complex_ar<CommutativeRing, Rank>::operator cdh_complex_ar<value_type, S>()
 const noexcept( std::is_nothrow_move_constructible<value_type>::value )
{ return b[0].operator cdh_complex_ar<value_type, S>(); }

/** Convert this object to an object of another instantiation of this class
    template, where the source type supports more components than the
    destination type.  Although the source type has to be able to implicitly
    convert to the destination type, this member function is marked explicit due
    to the loss of information.

    \pre  `U` can be implicitly converted to `#value_type`.
    \pre  `S` \>= `#rank`.

    \tparam U  The component type of the returned value.
    \tparam S  The rank (i.e. Cayley-Dickson rung) of the returned value.

    \throws  Whatever conversions from `value_type` to `U` may throw, plus any
             exceptions from returning `U` objects.

    \returns  An object `x` such that `get\<I\>(\*this) == get\<I\>(x)` for all
              `I` values the source and destination have in common.
 */
template < typename CommutativeRing, std::size_t Rank >
template <
    typename U, std::size_t S,
    typename std::enable_if<
        not std::is_same<CommutativeRing, U>::value &&
        std::is_convertible<CommutativeRing, U>::value && (S < Rank)
    >::type*...
>
inline constexpr
cdh_complex_ar<CommutativeRing, Rank>::operator cdh_complex_ar<U, S>() const
 noexcept(
     std::is_nothrow_constructible<U, value_type const &>::value
     && std::is_nothrow_move_constructible<U>::value
 )
{ return b[0].operator cdh_complex_ar<U, S>(); }

/** Convert this object to an object of another instantiation of this class
    template, where the source type supports more components than the
    destination type.  This particular member function works as an explicit
    conversion, and therefore requires the source type to be able to explicitly
    convert to the destination type.

    \pre  `U` can be explicitly converted to `#value_type`.
    \pre  `S` \>= `#rank`.

    \tparam U  The component type of the returned value.
    \tparam S  The rank (i.e. Cayley-Dickson rung) of the returned value.

    \throws  Whatever conversions from `value_type` to `U` may throw, plus any
             exceptions from returning `U` objects.

    \returns  An object `x` such that `get\<I\>(\*this) == get\<I\>(x)` for all
              `I` values the source and destination have in common.
 */
template < typename CommutativeRing, std::size_t Rank >
template <
    typename U, std::size_t S,
    typename std::enable_if<
        !std::is_convertible<CommutativeRing, U>::value && (S < Rank),
        decltype( static_cast<U>(std::declval<CommutativeRing const &>()) )
    >::type*...
>
inline constexpr
cdh_complex_ar<CommutativeRing, Rank>::operator cdh_complex_ar<U, S>() const
 noexcept(
     noexcept( static_cast<U>(std::declval<value_type const &>()) )
     && std::is_nothrow_move_constructible<U>::value
 )
{ return b[0].operator cdh_complex_ar<U, S>(); }


}  // namespace math
}  // namespace boost


namespace std
{


//  Class-template specializations for std::tuple traits classes  ------------//

//! Meta-function for the number of component elements in `cdh_complex_ai`.
template < typename T, size_t R >
class tuple_size< ::boost::math::cdh_complex_ai<T,R> >
    : public integral_constant< size_t, ::boost::math::cdh_complex_ai<T,
      R>::dimensions >
{ };

//! Meta-function for the number of component elements in `cdh_complex_ar`.
template < typename T, size_t R >
class tuple_size< ::boost::math::cdh_complex_ar<T,R> >
    : public integral_constant< size_t, ::boost::math::cdh_complex_ar<T,
      R>::dimensions >
{ };

//! Meta-function for the component type within `cdh_complex_ai`.
template < size_t I, typename T, size_t R >
class tuple_element< I, ::boost::math::cdh_complex_ai<T,R> >
{
    static_assert( I < (::boost::math::cdh_complex_ai<T, R>::dimensions),
     "index out of bounds" );

public:
    typedef T  type;
};

//! Meta-function for the component type within `cdh_complex_ar`.
template < size_t I, typename T, size_t R >
class tuple_element< I, ::boost::math::cdh_complex_ar<T,R> >
{
    static_assert( I < (::boost::math::cdh_complex_ar<T, R>::dimensions),
     "index out of bounds" );

public:
    typedef T  type;
};


}  // namespace std


#endif  // BOOST_MATH_CD_HYPERCOMPLEX_CORE_HPP
