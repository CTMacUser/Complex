//  Boost.Complex unit test program file  ------------------------------------//

//  Copyright 2012 Daryle Walker.
//  Distributed under the Boost Software License, Version 1.0.  (See the
//  accompanying file LICENSE_1_0.txt or a copy at
//  <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

#define BOOST_TEST_MODULE Cayley-Dickson Hypercomplex Unit Test
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/output_test_stream.hpp>
#include "boost/math/cd_hypercomplex.hpp"

#include <boost/mpl/list.hpp>  // for boost::mpl::list

#include <cstddef>      // for std::size_t
#include <ctime>        // for std::time
#include <ios>          // for std::ios_base
#include <limits>       // for std::numeric_limits
#include <ostream>      // for std::basic_ostream
#include <random>       // for std::default_random_engine, etc.
#include <tuple>        // for std::tuple_size, tuple_element
#include <type_traits>  // for std::true_type, false_type, is_same, etc.

#ifndef CONTROL_USE_HW_RANDOM_ENGINE
#define CONTROL_USE_HW_RANDOM_ENGINE  0
#endif

#define SIMPLE_STATIC_ASSERT(...)  static_assert(__VA_ARGS__, #__VA_ARGS__)


// Common definitions  -------------------------------------------------------//

namespace {

// Elide need for full name qualifications
using boost::math::cdh_complex_ai;
using boost::math::cdh_complex_ar;

using boost::test_tools::output_test_stream;

// Use common names for each Cayley-Dickson level
template < typename T >  using real_ai_t       = cdh_complex_ai<T, 0>;
template < typename T >  using real_ar_t       = cdh_complex_ar<T, 0>;
template < typename T >  using complex_ai_t    = cdh_complex_ai<T, 1>;
template < typename T >  using complex_ar_t    = cdh_complex_ar<T, 1>;
template < typename T >  using quaternion_ai_t = cdh_complex_ai<T, 2>;
template < typename T >  using quaternion_ar_t = cdh_complex_ar<T, 2>;
template < typename T >  using octonion_ai_t   = cdh_complex_ai<T, 3>;
template < typename T >  using octonion_ar_t   = cdh_complex_ar<T, 3>;

// Lists of types to use in test runs
typedef boost::mpl::list<short, int, long>
 signed_types;
typedef boost::mpl::list<unsigned short, unsigned, unsigned long>
 unsigned_types;
typedef boost::mpl::list<float, double>
 floating_types;

typedef boost::mpl::list<short, unsigned short, int, unsigned, long,
 unsigned long>  integer_types;
typedef boost::mpl::list<short, unsigned short, int, unsigned, long,
 unsigned long, float, double>  numeric_types;
typedef boost::mpl::list<short, int, long, float, double>  signing_types;

// Random-number generator stuff
#if CONTROL_USE_HW_RANDOM_ENGINE
typedef std::random_device          re_t;
#else
typedef std::default_random_engine  re_t;
#endif

auto  make_re() -> re_t
{
#if CONTROL_USE_HW_RANDOM_ENGINE
    return {};
#else
    re_t  result;

    result.seed( static_cast<unsigned>(std::time( nullptr )) );
    return result;
#endif
}

auto  get_re() -> re_t &
{
    static auto  re = make_re();

    return re;
}

namespace detail
{
    template < typename T >
    auto  get_random_integer( std::true_type ) -> T
    { return std::uniform_int_distribution<T>{-100, +100}(get_re()); }

    template < typename T >
    auto  get_random_integer( std::false_type ) -> T
    { return std::uniform_int_distribution<T>{0u, 100u}(get_re()); }

    template < typename T >
    auto  get_random_number_impl( std::true_type ) -> T
    {
        return get_random_integer<T>( std::integral_constant<bool,
         std::numeric_limits<T>::is_signed>{} );
    }

    template < typename T >
    auto  get_random_number_impl( std::false_type ) -> T
    { return std::uniform_real_distribution<T>{-2.0f, +2.0f}(get_re()); }

    template <
        typename T,
        class Enable = typename std::enable_if<
            std::numeric_limits<T>::is_specialized
        >::type
    >
    auto  get_random_number() -> T
    {
        return get_random_number_impl<T>( std::integral_constant<bool,
         std::numeric_limits<T>::is_integer>{} );
    }
}

using detail::get_random_number;

template < typename T, std::size_t R >
auto  get_random_cdh_complex_ai() -> cdh_complex_ai<T, R>
{
    cdh_complex_ai<T, R>  result;

    for ( auto &x : result.c )
        x = get_random_number<T>();
    return result;
}

namespace detail
{
    template < typename T, std::size_t R >
    struct cdh_complex_ar_randomizer;

    template < typename T >
    struct cdh_complex_ar_randomizer<T, 0>
    {
        auto  operator ()() const -> cdh_complex_ar<T, 0>
        { return {{ get_random_number<T>() }}; }
    };

    template < typename T, std::size_t R >
    struct cdh_complex_ar_randomizer
    {
        auto  operator ()() const -> cdh_complex_ar<T, R>
        {
            cdh_complex_ar_randomizer<T, R - 1> const  rng;

            return { {rng(), rng()} };
        }
    };
}

template < typename T, std::size_t R >
auto  get_random_cdh_complex_ar() -> cdh_complex_ar<T, R>
{ return detail::cdh_complex_ar_randomizer<T, R>{}(); }

// Make sure an object never has a value of zero
// (Depends on two-way Boolean conversions and addition.)
template < typename T >
void  never_zero( T &t )  { t += !t; }

template < typename T, std::size_t R >
auto  get_random_nonzero_cdh_complex_ai() -> cdh_complex_ai<T, R>
{
    auto  result = get_random_cdh_complex_ai<T, R>();

    never_zero( result.c[std::uniform_int_distribution<size_t>{ 0,
     decltype(result)::dimensions - 1 }( get_re() )] );
    return result;
}

namespace detail
{
    template < typename T >
    void  never0_r( cdh_complex_ar<T, 0> &x )
    { never_zero(x.r[ 0 ]); }

    template < typename T, std::size_t R >
    void  never0_r( cdh_complex_ar<T, R> &x )
    { never0_r(x.b[ std::bernoulli_distribution{}(get_re()) ]); }
}

template < typename T, std::size_t R >
auto  get_random_nonzero_cdh_complex_ar() -> cdh_complex_ar<T, R>
{
   auto  result = get_random_cdh_complex_ar<T, R>();

   detail::never0_r( result );
   return result;
}

// Create a numeric type that doesn't implicitly convert to the built-in ones.
enum class two_bit_t
    : unsigned
{ zero, one, two, three };

template < typename Ch, class Tr >
std::basic_ostream<Ch, Tr> &
operator <<( std::basic_ostream<Ch, Tr> &o, two_bit_t x )
{ return o << static_cast<unsigned>(x); }

}  // anonymous namespace


// Unit tests  ---------------------------------------------------------------//

BOOST_AUTO_TEST_SUITE( core_cdhc_tests )

BOOST_AUTO_TEST_CASE_TEMPLATE( compile_time_attribute_tests, T, numeric_types )
{
    using std::is_same;

    // Confirm that the rank attribute works correctly
    SIMPLE_STATIC_ASSERT( real_ai_t<T>::rank == 0 );
    SIMPLE_STATIC_ASSERT( real_ar_t<T>::rank == 0 );
    SIMPLE_STATIC_ASSERT( complex_ai_t<T>::rank == 1 );
    SIMPLE_STATIC_ASSERT( complex_ar_t<T>::rank == 1 );
    SIMPLE_STATIC_ASSERT( quaternion_ai_t<T>::rank == 2 );
    SIMPLE_STATIC_ASSERT( quaternion_ar_t<T>::rank == 2 );
    SIMPLE_STATIC_ASSERT( octonion_ai_t<T>::rank == 3 );
    SIMPLE_STATIC_ASSERT( octonion_ar_t<T>::rank == 3 );

    // Confirm that the value_type attribute works correctly
    SIMPLE_STATIC_ASSERT( is_same<T,
     typename real_ai_t<T>::value_type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T,
     typename real_ar_t<T>::value_type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T,
     typename complex_ai_t<T>::value_type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T,
     typename complex_ar_t<T>::value_type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T,
     typename quaternion_ai_t<T>::value_type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T,
     typename quaternion_ar_t<T>::value_type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T,
     typename octonion_ai_t<T>::value_type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T,
     typename octonion_ar_t<T>::value_type>::value );

    // Confirm that the element_type attribute works correctly
    SIMPLE_STATIC_ASSERT( is_same<T,
     typename real_ai_t<T>::element_type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T,
     typename real_ar_t<T>::element_type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T,
     typename complex_ai_t<T>::element_type>::value );
    SIMPLE_STATIC_ASSERT( is_same<real_ar_t<T>,
     typename complex_ar_t<T>::element_type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T,
     typename quaternion_ai_t<T>::element_type>::value );
    SIMPLE_STATIC_ASSERT( is_same<complex_ar_t<T>,
     typename quaternion_ar_t<T>::element_type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T,
     typename octonion_ai_t<T>::element_type>::value );
    SIMPLE_STATIC_ASSERT( is_same<quaternion_ar_t<T>,
     typename octonion_ar_t<T>::element_type>::value );

    // Confirm that the element_count attribute works correctly
    SIMPLE_STATIC_ASSERT( real_ai_t<T>::element_count == 1 );
    SIMPLE_STATIC_ASSERT( real_ar_t<T>::element_count == 1 );
    SIMPLE_STATIC_ASSERT( complex_ai_t<T>::element_count == 2 );
    SIMPLE_STATIC_ASSERT( complex_ar_t<T>::element_count == 2 );
    SIMPLE_STATIC_ASSERT( quaternion_ai_t<T>::element_count == 4 );
    SIMPLE_STATIC_ASSERT( quaternion_ar_t<T>::element_count == 2 );
    SIMPLE_STATIC_ASSERT( octonion_ai_t<T>::element_count == 8 );
    SIMPLE_STATIC_ASSERT( octonion_ar_t<T>::element_count == 2 );

    // Confirm that the dimensions attribute works correctly
    SIMPLE_STATIC_ASSERT( real_ai_t<T>::dimensions == 1 );
    SIMPLE_STATIC_ASSERT( real_ar_t<T>::dimensions == 1 );
    SIMPLE_STATIC_ASSERT( complex_ai_t<T>::dimensions == 2 );
    SIMPLE_STATIC_ASSERT( complex_ar_t<T>::dimensions == 2 );
    SIMPLE_STATIC_ASSERT( quaternion_ai_t<T>::dimensions == 4 );
    SIMPLE_STATIC_ASSERT( quaternion_ar_t<T>::dimensions == 4 );
    SIMPLE_STATIC_ASSERT( octonion_ai_t<T>::dimensions == 8 );
    SIMPLE_STATIC_ASSERT( octonion_ar_t<T>::dimensions == 8 );

    // Size check
    SIMPLE_STATIC_ASSERT( sizeof(real_ai_t<T>) >= sizeof(T) );
    SIMPLE_STATIC_ASSERT( sizeof(real_ar_t<T>) >= sizeof(T) );
    SIMPLE_STATIC_ASSERT( sizeof(complex_ai_t<T>) >= 2u * sizeof(T) );
    SIMPLE_STATIC_ASSERT( sizeof(complex_ar_t<T>) >= 2u * sizeof(T) );
    SIMPLE_STATIC_ASSERT( sizeof(quaternion_ai_t<T>) >= 4u * sizeof(T) );
    SIMPLE_STATIC_ASSERT( sizeof(quaternion_ar_t<T>) >= 4u * sizeof(T) );
    SIMPLE_STATIC_ASSERT( sizeof(octonion_ai_t<T>) >= 8u * sizeof(T) );
    SIMPLE_STATIC_ASSERT( sizeof(octonion_ar_t<T>) >= 8u * sizeof(T) );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( bool_tests, T, numeric_types )
{
    real_ai_t<T> const        test1{};
    real_ar_t<T> const        test2{};
    complex_ai_t<T> const     test3{};
    complex_ar_t<T> const     test4{};
    quaternion_ai_t<T> const  test5{};
    quaternion_ar_t<T> const  test6{};
    octonion_ai_t<T> const    test7{};
    octonion_ar_t<T> const    test8{};

    // Direct boolean check
    BOOST_CHECK_EQUAL( static_cast<bool>(test1), false );
    BOOST_CHECK_EQUAL( static_cast<bool>(test2), false );
    BOOST_CHECK_EQUAL( static_cast<bool>(test3), false );
    BOOST_CHECK_EQUAL( static_cast<bool>(test4), false );
    BOOST_CHECK_EQUAL( static_cast<bool>(test5), false );
    BOOST_CHECK_EQUAL( static_cast<bool>(test6), false );
    BOOST_CHECK_EQUAL( static_cast<bool>(test7), false );
    BOOST_CHECK_EQUAL( static_cast<bool>(test8), false );

    // Indirect boolean check
    BOOST_CHECK( !test1 );
    BOOST_CHECK( !test2 );
    BOOST_CHECK( !test3 );
    BOOST_CHECK( !test4 );
    BOOST_CHECK( !test5 );
    BOOST_CHECK( !test6 );
    BOOST_CHECK( !test7 );
    BOOST_CHECK( !test8 );

    // Check when non-zero
    BOOST_CHECK( static_cast<bool>(get_random_nonzero_cdh_complex_ai<T,0>()) );
    BOOST_CHECK( static_cast<bool>(get_random_nonzero_cdh_complex_ar<T,0>()) );
    BOOST_CHECK( static_cast<bool>(get_random_nonzero_cdh_complex_ai<T,1>()) );
    BOOST_CHECK( static_cast<bool>(get_random_nonzero_cdh_complex_ar<T,1>()) );
    BOOST_CHECK( static_cast<bool>(get_random_nonzero_cdh_complex_ai<T,2>()) );
    BOOST_CHECK( static_cast<bool>(get_random_nonzero_cdh_complex_ar<T,2>()) );
    BOOST_CHECK( static_cast<bool>(get_random_nonzero_cdh_complex_ai<T,3>()) );
    BOOST_CHECK( static_cast<bool>(get_random_nonzero_cdh_complex_ar<T,3>()) );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( same_type_equality_tests, T, numeric_types )
{
    real_ai_t<T> const        test1a{};
    real_ai_t<T> const        test1b = get_random_nonzero_cdh_complex_ai<T,0>();
    real_ar_t<T> const        test2a{};
    real_ar_t<T> const        test2b = get_random_nonzero_cdh_complex_ar<T,0>();
    complex_ai_t<T> const     test3a{};
    complex_ai_t<T> const     test3b = get_random_nonzero_cdh_complex_ai<T,1>();
    complex_ar_t<T> const     test4a{};
    complex_ar_t<T> const     test4b = get_random_nonzero_cdh_complex_ar<T,1>();
    quaternion_ai_t<T> const  test5a{};
    quaternion_ai_t<T> const  test5b = get_random_nonzero_cdh_complex_ai<T,2>();
    quaternion_ar_t<T> const  test6a{};
    quaternion_ar_t<T> const  test6b = get_random_nonzero_cdh_complex_ar<T,2>();
    octonion_ai_t<T> const    test7a{};
    octonion_ai_t<T> const    test7b = get_random_nonzero_cdh_complex_ai<T,3>();
    octonion_ar_t<T> const    test8a{};
    octonion_ar_t<T> const    test8b = get_random_nonzero_cdh_complex_ar<T,3>();

    // Equality check
    BOOST_CHECK( test1a == test1a );
    BOOST_CHECK( test2a == test2a );
    BOOST_CHECK( test3a == test3a );
    BOOST_CHECK( test4a == test4a );
    BOOST_CHECK( test5a == test5a );
    BOOST_CHECK( test6a == test6a );
    BOOST_CHECK( test7a == test7a );
    BOOST_CHECK( test8a == test8a );

    BOOST_CHECK( test1b == test1b );
    BOOST_CHECK( test2b == test2b );
    BOOST_CHECK( test3b == test3b );
    BOOST_CHECK( test4b == test4b );
    BOOST_CHECK( test5b == test5b );
    BOOST_CHECK( test6b == test6b );
    BOOST_CHECK( test7b == test7b );
    BOOST_CHECK( test8b == test8b );

    BOOST_CHECK( !(test1a == test1b) );
    BOOST_CHECK( !(test2a == test2b) );
    BOOST_CHECK( !(test3a == test3b) );
    BOOST_CHECK( !(test4a == test4b) );
    BOOST_CHECK( !(test5a == test5b) );
    BOOST_CHECK( !(test6a == test6b) );
    BOOST_CHECK( !(test7a == test7b) );
    BOOST_CHECK( !(test8a == test8b) );

    // Inequality check
    BOOST_CHECK( test1a != test1b );
    BOOST_CHECK( test2a != test2b );
    BOOST_CHECK( test3a != test3b );
    BOOST_CHECK( test4a != test4b );
    BOOST_CHECK( test5a != test5b );
    BOOST_CHECK( test6a != test6b );
    BOOST_CHECK( test7a != test7b );
    BOOST_CHECK( test8a != test8b );

    // Boolean check comparison
    BOOST_CHECK_EQUAL( static_cast<bool>(test1b), test1b != test1a );
    BOOST_CHECK_EQUAL( static_cast<bool>(test2b), test2b != test2a );
    BOOST_CHECK_EQUAL( static_cast<bool>(test3b), test3b != test3a );
    BOOST_CHECK_EQUAL( static_cast<bool>(test4b), test4b != test4a );
    BOOST_CHECK_EQUAL( static_cast<bool>(test5b), test5b != test5a );
    BOOST_CHECK_EQUAL( static_cast<bool>(test6b), test6b != test6a );
    BOOST_CHECK_EQUAL( static_cast<bool>(test7b), test7b != test7a );
    BOOST_CHECK_EQUAL( static_cast<bool>(test8b), test8b != test8a );

    BOOST_CHECK_EQUAL( static_cast<bool>(test1a), test1a != test1a );
    BOOST_CHECK_EQUAL( static_cast<bool>(test2a), test2a != test2a );
    BOOST_CHECK_EQUAL( static_cast<bool>(test3a), test3a != test3a );
    BOOST_CHECK_EQUAL( static_cast<bool>(test4a), test4a != test4a );
    BOOST_CHECK_EQUAL( static_cast<bool>(test5a), test5a != test5a );
    BOOST_CHECK_EQUAL( static_cast<bool>(test6a), test6a != test6a );
    BOOST_CHECK_EQUAL( static_cast<bool>(test7a), test7a != test7a );
    BOOST_CHECK_EQUAL( static_cast<bool>(test8a), test8a != test8a );

    // Uneven-length checks
    complex_ai_t<T> const     test9{ {1, 2} };
    quaternion_ai_t<T> const  test10{ {1, 2, 3, 4} };
    quaternion_ai_t<T> const  test11{ {1, 2, 0, 0} };
    quaternion_ai_t<T> const  test12{ {6, 2, 0, 0} };

    BOOST_CHECK( test9 != test10 );
    BOOST_CHECK( !(test9 == test10) );
    BOOST_CHECK( test9 == test11 );
    BOOST_CHECK( !(test9 != test11) );
    BOOST_CHECK( test9 != test12 );
    BOOST_CHECK( !(test9 == test12) );

    BOOST_CHECK( test10 != test9 );
    BOOST_CHECK( !(test10 == test9) );
    BOOST_CHECK( test11 == test9 );
    BOOST_CHECK( !(test11 != test9) );
    BOOST_CHECK( test12 != test9 );
    BOOST_CHECK( !(test12 == test9) );

    real_ar_t<T> const        test13 = get_random_nonzero_cdh_complex_ar<T,0>();
    complex_ar_t<T> const     test14{ {test13, {{ 0 }}} };
    complex_ar_t<T> const     test15{ {test13, {{ 7 }}} };
    quaternion_ar_t<T> const  test16{ {test14, {}} };
    quaternion_ar_t<T> const  test17{ {test15, {}} };
    quaternion_ar_t<T> const  test18{ {test14, test15} };
    quaternion_ar_t<T> const  test19{{{{{{5}},{{9}}}},{{{{2}},{{8}}}}}};

    BOOST_CHECK( test13 == test14 );
    BOOST_CHECK( test14 == test13 );
    BOOST_CHECK( !(test13 != test14) );
    BOOST_CHECK( !(test14 != test13) );

    BOOST_CHECK( !(test13 == test15) );
    BOOST_CHECK( !(test15 == test13) );
    BOOST_CHECK( test13 != test15 );
    BOOST_CHECK( test15 != test13 );

    BOOST_CHECK( test13 == test16 );
    BOOST_CHECK( test16 == test13 );
    BOOST_CHECK( !(test13 != test16) );
    BOOST_CHECK( !(test16 != test13) );

    BOOST_CHECK( test14 == test16 );
    BOOST_CHECK( test16 == test14 );
    BOOST_CHECK( !(test14 != test16) );
    BOOST_CHECK( !(test16 != test14) );

    BOOST_CHECK( !(test13 == test17) );
    BOOST_CHECK( !(test17 == test13) );
    BOOST_CHECK( test13 != test17 );
    BOOST_CHECK( test17 != test13 );

    BOOST_CHECK( !(test14 == test17) );
    BOOST_CHECK( !(test17 == test14) );
    BOOST_CHECK( test14 != test17 );
    BOOST_CHECK( test17 != test14 );

    BOOST_CHECK( !(test13 == test18) );
    BOOST_CHECK( !(test18 == test13) );
    BOOST_CHECK( test13 != test18 );
    BOOST_CHECK( test18 != test13 );

    BOOST_CHECK( !(test14 == test18) );
    BOOST_CHECK( !(test18 == test14) );
    BOOST_CHECK( test14 != test18 );
    BOOST_CHECK( test18 != test14 );

    BOOST_CHECK( !(test13 == test19) );
    BOOST_CHECK( !(test19 == test13) );
    BOOST_CHECK( test13 != test19 );
    BOOST_CHECK( test19 != test13 );

    BOOST_CHECK( !(test14 == test19) );
    BOOST_CHECK( !(test19 == test14) );
    BOOST_CHECK( test14 != test19 );
    BOOST_CHECK( test19 != test14 );

    BOOST_CHECK( !(test15 == test19) );
    BOOST_CHECK( !(test19 == test15) );
    BOOST_CHECK( test15 != test19 );
    BOOST_CHECK( test19 != test15 );

    BOOST_CHECK( !(test16 == test19) );
    BOOST_CHECK( !(test19 == test16) );
    BOOST_CHECK( test16 != test19 );
    BOOST_CHECK( test19 != test16 );
}

BOOST_AUTO_TEST_CASE( cross_type_equality_tests )
{
    real_ai_t<int> const           test1a{ {3} };
    real_ai_t<unsigned> const      test2a{ {8} }, test3a{ {3} };
    complex_ai_t<int> const        test4a{ {8, 0} }, test5a{ {8, 2} };
    quaternion_ai_t<double> const  test6a{ {8, 0, 0, 0} }, test7a{{8, 0, 1, 0}};

    BOOST_CHECK( !(test1a == test2a) );
    BOOST_CHECK( !(test2a == test1a) );
    BOOST_CHECK( test1a != test2a );
    BOOST_CHECK( test2a != test1a );

    BOOST_CHECK( test1a == test3a );
    BOOST_CHECK( test3a == test1a );
    BOOST_CHECK( !(test1a != test3a) );
    BOOST_CHECK( !(test3a != test1a) );

    BOOST_CHECK( test2a == test4a );
    BOOST_CHECK( test4a == test2a );
    BOOST_CHECK( !(test2a != test4a) );
    BOOST_CHECK( !(test4a != test2a) );

    BOOST_CHECK( !(test2a == test5a) );
    BOOST_CHECK( !(test5a == test2a) );
    BOOST_CHECK( test2a != test5a );
    BOOST_CHECK( test5a != test2a );

    BOOST_CHECK( test2a == test6a );
    BOOST_CHECK( test6a == test2a );
    BOOST_CHECK( !(test2a != test6a) );
    BOOST_CHECK( !(test6a != test2a) );

    BOOST_CHECK( !(test2a == test7a) );
    BOOST_CHECK( !(test7a == test2a) );
    BOOST_CHECK( test2a != test7a );
    BOOST_CHECK( test7a != test2a );

    BOOST_CHECK( test4a == test6a );
    BOOST_CHECK( test6a == test4a );
    BOOST_CHECK( !(test4a != test6a) );
    BOOST_CHECK( !(test6a != test4a) );

    BOOST_CHECK( !(test4a == test7a) );
    BOOST_CHECK( !(test7a == test4a) );
    BOOST_CHECK( test4a != test7a );
    BOOST_CHECK( test7a != test4a );

    real_ar_t<int> const           test1b{ {3} };
    real_ar_t<unsigned> const      test2b{ {8} }, test3b{ {3} };
    complex_ar_t<int> const        test4b{ {{ {8} }, { {0} }} };
    complex_ar_t<int> const        test5b{ {{ {8} }, { {2} }} };
    quaternion_ar_t<double> const  test6b{{{{{{8}}, {{0}}}}, {{{{0}}, {{0}}}}}};
    quaternion_ar_t<double> const  test7b{{{{{{8}}, {{0}}}}, {{{{1}}, {{0}}}}}};

    BOOST_CHECK( !(test1b == test2b) );
    BOOST_CHECK( !(test2b == test1b) );
    BOOST_CHECK( test1b != test2b );
    BOOST_CHECK( test2b != test1b );

    BOOST_CHECK( test1b == test3b );
    BOOST_CHECK( test3b == test1b );
    BOOST_CHECK( !(test1b != test3b) );
    BOOST_CHECK( !(test3b != test1b) );

    BOOST_CHECK( test2b == test4b );
    BOOST_CHECK( test4b == test2b );
    BOOST_CHECK( !(test2b != test4b) );
    BOOST_CHECK( !(test4b != test2b) );

    BOOST_CHECK( !(test2b == test5b) );
    BOOST_CHECK( !(test5b == test2b) );
    BOOST_CHECK( test2b != test5b );
    BOOST_CHECK( test5b != test2b );

    BOOST_CHECK( test2b == test6b );
    BOOST_CHECK( test6b == test2b );
    BOOST_CHECK( !(test2b != test6b) );
    BOOST_CHECK( !(test6b != test2b) );

    BOOST_CHECK( !(test2b == test7b) );
    BOOST_CHECK( !(test7b == test2b) );
    BOOST_CHECK( test2b != test7b );
    BOOST_CHECK( test7b != test2b );

    BOOST_CHECK( test4b == test6b );
    BOOST_CHECK( test6b == test4b );
    BOOST_CHECK( !(test4b != test6b) );
    BOOST_CHECK( !(test6b != test4b) );

    BOOST_CHECK( !(test4b == test7b) );
    BOOST_CHECK( !(test7b == test4b) );
    BOOST_CHECK( test4b != test7b );
    BOOST_CHECK( test7b != test4b );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( dynamic_rank_tests, T, numeric_types )
{
    using boost::math::dynamic_rank;

    real_ai_t<T> const        test0a{}, test1a{ {1} };
    complex_ai_t<T> const     test2a{}, test3a{ {1} }, test4a{ {0, 1} },
     test5a{ {1, 1} };
    quaternion_ai_t<T> const  test6a{}, test7a{ {1} }, test8a{ {0, 1} },
     test9a{ {0, 0, 1} }, test10a{ {0, 0, 0, 1} };
    quaternion_ai_t<T> const  test11a{ {1, 1} }, test12a{ {1, 0, 1} },
     test13a{ {1, 0, 0, 1} };

    BOOST_CHECK_EQUAL( dynamic_rank(test0a), 0u );
    BOOST_CHECK_EQUAL( dynamic_rank(test1a), 0u );
    BOOST_CHECK_EQUAL( dynamic_rank(test2a), 0u );
    BOOST_CHECK_EQUAL( dynamic_rank(test3a), 0u );
    BOOST_CHECK_EQUAL( dynamic_rank(test4a), 1u );
    BOOST_CHECK_EQUAL( dynamic_rank(test5a), 1u );
    BOOST_CHECK_EQUAL( dynamic_rank(test6a), 0u );
    BOOST_CHECK_EQUAL( dynamic_rank(test7a), 0u );
    BOOST_CHECK_EQUAL( dynamic_rank(test8a), 1u );
    BOOST_CHECK_EQUAL( dynamic_rank(test9a), 2u );
    BOOST_CHECK_EQUAL( dynamic_rank(test10a), 2u );
    BOOST_CHECK_EQUAL( dynamic_rank(test11a), 1u );
    BOOST_CHECK_EQUAL( dynamic_rank(test12a), 2u );
    BOOST_CHECK_EQUAL( dynamic_rank(test13a), 2u );

    real_ar_t<T> const        test0b{}, test1b{ {1} };
    complex_ar_t<T> const     test2b{}, test3b{ {test1b} },
     test4b{ {test0b, test1b} }, test5b{ {test1b, test1b} };
    quaternion_ar_t<T> const  test6b{}, test7b{ {test3b} }, test8b{ {test4b} },
     test9b{ {test2b, test3b} }, test10b{ {test2b, test4b} };
    quaternion_ar_t<T> const  test11b{ {test5b} },
     test12b{ {test3b, test3b} }, test13b{ {test3b, test4b} };

    BOOST_CHECK_EQUAL( dynamic_rank(test0b), 0u );
    BOOST_CHECK_EQUAL( dynamic_rank(test1b), 0u );
    BOOST_CHECK_EQUAL( dynamic_rank(test2b), 0u );
    BOOST_CHECK_EQUAL( dynamic_rank(test3b), 0u );
    BOOST_CHECK_EQUAL( dynamic_rank(test4b), 1u );
    BOOST_CHECK_EQUAL( dynamic_rank(test5b), 1u );
    BOOST_CHECK_EQUAL( dynamic_rank(test6b), 0u );
    BOOST_CHECK_EQUAL( dynamic_rank(test7b), 0u );
    BOOST_CHECK_EQUAL( dynamic_rank(test8b), 1u );
    BOOST_CHECK_EQUAL( dynamic_rank(test9b), 2u );
    BOOST_CHECK_EQUAL( dynamic_rank(test10b), 2u );
    BOOST_CHECK_EQUAL( dynamic_rank(test11b), 1u );
    BOOST_CHECK_EQUAL( dynamic_rank(test12b), 2u );
    BOOST_CHECK_EQUAL( dynamic_rank(test13b), 2u );
}

BOOST_AUTO_TEST_SUITE( core_writing_cdhc_tests )

BOOST_AUTO_TEST_CASE_TEMPLATE( basic_output_tests, T, numeric_types )
{
    using std::ios_base;

    output_test_stream  o;

    // Real values don't get parenthesized.
    o << real_ai_t<T>{ {0} };
    BOOST_CHECK( o.is_equal("0") );
    o << real_ai_t<T>{ {7} };
    BOOST_CHECK( o.is_equal("7") );

    o << real_ar_t<T>{ {0} };
    BOOST_CHECK( o.is_equal("0") );
    o << real_ar_t<T>{ {11} };
    BOOST_CHECK( o.is_equal("11") );

    // Iterative- and recursive-designs differ in parenthesization.
    o << complex_ai_t<T>{ {9, 4} };
    BOOST_CHECK( o.is_equal("(9,4)") );
    o << complex_ar_t<T>{ {{ {4} }, { {9} }} };
    BOOST_CHECK( o.is_equal("(4,9)") );

    o << quaternion_ai_t<T>{ {16, 27, 25, 8} };
    BOOST_CHECK( o.is_equal("(16,27,25,8)") );
    o << quaternion_ar_t<T>{{{ {{ {8} }, { {27} }} },{ {{ {16} }, { {25} }} }}};
    BOOST_CHECK( o.is_equal("((8,27),(16,25))") );

    // Zero-valued (top-level) upper-barrages aren't printed.
    o << complex_ai_t<T>{ {5, 0} };
    BOOST_CHECK( o.is_equal("5") );
    o << complex_ar_t<T>{ {{ {6} }, { {0} }} };
    BOOST_CHECK( o.is_equal("6") );

    o << quaternion_ai_t<T>{ {12, 31, 0, 0} };
    BOOST_CHECK( o.is_equal("(12,31)") );
    o << quaternion_ar_t<T>{{{ {{ {33} }, { {28} }} },{ {{ {0} }, { {0} }} }}};
    BOOST_CHECK( o.is_equal("(33,28)") );

    o << quaternion_ai_t<T>{ {3, 0, 0, 0} };
    BOOST_CHECK( o.is_equal("3") );
    o << quaternion_ar_t<T>{{{ {{ {2} }, { {0} }} },{ {{ {0} }, { {0} }} }}};
    BOOST_CHECK( o.is_equal("2") );

    o << quaternion_ai_t<T>{ {13, 0, 15, 0} };
    BOOST_CHECK( o.is_equal("(13,0,15,0)") );
    o << quaternion_ar_t<T>{{{ {{ {0} }, { {0} }} },{ {{ {17} }, { {0} }} }}};
    BOOST_CHECK( o.is_equal("((0,0),(17,0))") );

    // Width attribute supported
    quaternion_ai_t<T> const  sample1{ {1, 10, 14, 18} };
    quaternion_ar_t<T> const  sample2{{{{{{1}},{{10}}}},{{{{14}},{{18}}}}}};

    o << sample1;
    BOOST_CHECK( o.is_equal("(1,10,14,18)") );
    o << sample2;
    BOOST_CHECK( o.is_equal("((1,10),(14,18))") );

    o.setf( ios_base::left, ios_base::adjustfield );
    o.width( 20 );
    o << sample1;
    BOOST_CHECK( o.is_equal("(1,10,14,18)        ") );
    o.width( 20 );
    o << sample2;
    BOOST_CHECK( o.is_equal("((1,10),(14,18))    ") );

    o.setf( ios_base::right, ios_base::adjustfield );
    o.width( 20 );
    o << sample1;
    BOOST_CHECK( o.is_equal("        (1,10,14,18)") );
    o.width( 20 );
    o << sample2;
    BOOST_CHECK( o.is_equal("    ((1,10),(14,18))") );

    o.setf( ios_base::internal, ios_base::adjustfield );
    o.width( 20 );
    o << sample1;
    BOOST_CHECK( o.is_equal("   (  1, 10, 14, 18)") );
     // Why?  Start with 20 spaces.  There are 4 components to be printed,
     // needing 4 + 1 = 5 punctuation marks, leaving 20 - 5 = 15 remaining
     // spaces.  Each component gets 15 / 4 = 3 spaces each.  Each of the
     // 4 components can be printed in 3 spaces, with the numbers right-
     // justified.  There are 15 % 4 = 3 spaces left over, being the extra
     // space when the tuple as a whole is right-justified.
    o.width( 20 );
    o << sample2;
    BOOST_CHECK( o.is_equal("   (( 1,10),(14,18))") );
     // Why?  Start with 20 spaces.  There are two barrages to be printed, with
     // 3 punctuation marks.  That's (20 - 3) / 2 = 8 spaces per barrage, with
     // (20 - 3) % 2 = 1 space left over.  For each barrage, there are two
     // numbers to be printed, along with 3 punctuation marks.  That's (8 - 3) /
     // 2 = 2 spaces per number, with (8 - 3) % 2 = 1 space left over.  Each of
     // the 4 components can be printed in 3 spaces, with the numbers right-
     // justified.  Note that we only use 2 * 2 + 3 = 7 spaces per barrage, and
     // 2 * 7 + 3 = 17 spaces for the entire tuple, so we really get 3 spaces
     // left over (instead of 1) as the spacing for right-justifying the tuple.
}

BOOST_AUTO_TEST_CASE_TEMPLATE( sign_and_width_output_tests, T, signing_types )
{
    using std::ios_base;

    output_test_stream        o;
    quaternion_ai_t<T> const  sample1{ {1, 10, 14, 18} };
    quaternion_ar_t<T> const  sample2{{{{{{1}},{{10}}}},{{{{14}},{{18}}}}}};

    o.setf( ios_base::showpos );
    o.setf( ios_base::left, ios_base::adjustfield );
    o.width( 20 );
    o << sample1;
    BOOST_CHECK( o.is_equal("(+1,+10,+14,+18)    ") );
    o.width( 20 );
    o << sample2;
    BOOST_CHECK( o.is_equal("((+1,+10),(+14,+18))") );

    o.setf( ios_base::right, ios_base::adjustfield );
    o.width( 20 );
    o << sample1;
    BOOST_CHECK( o.is_equal("    (+1,+10,+14,+18)") );
    o.width( 20 );
    o << sample2;
    BOOST_CHECK( o.is_equal("((+1,+10),(+14,+18))") );

    o.setf( ios_base::internal, ios_base::adjustfield );
    o.width( 20 );
    o << sample1;
    BOOST_CHECK( o.is_equal("   (+ 1,+10,+14,+18)") );
    o.width( 20 );
    o << sample2;
    BOOST_CHECK( o.is_equal("((+1,+10),(+14,+18))") );
}

BOOST_AUTO_TEST_SUITE_END()  // core_writing_cdhc_tests
BOOST_AUTO_TEST_SUITE_END()  // core_cdhc_tests

BOOST_AUTO_TEST_SUITE( tuple_cdhc_tests )

BOOST_AUTO_TEST_CASE_TEMPLATE( compile_time_tuple_tests, T, numeric_types )
{
    using std::tuple_size;
    using std::is_same;
    using std::tuple_element;

    // Confirm that the tuple-size attribute works correctly
    SIMPLE_STATIC_ASSERT( tuple_size<real_ai_t<T>>::value == 1 );
    SIMPLE_STATIC_ASSERT( tuple_size<real_ar_t<T>>::value == 1 );
    SIMPLE_STATIC_ASSERT( tuple_size<complex_ai_t<T>>::value == 2 );
    SIMPLE_STATIC_ASSERT( tuple_size<complex_ar_t<T>>::value == 2 );
    SIMPLE_STATIC_ASSERT( tuple_size<quaternion_ai_t<T>>::value == 4 );
    SIMPLE_STATIC_ASSERT( tuple_size<quaternion_ar_t<T>>::value == 4 );
    SIMPLE_STATIC_ASSERT( tuple_size<octonion_ai_t<T>>::value == 8 );
    SIMPLE_STATIC_ASSERT( tuple_size<octonion_ar_t<T>>::value == 8 );

    // Confirm that the tuple-element attribute works correctly
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<0,
     real_ai_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<0,
     real_ar_t<T>>::type>::value );

    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<0,
     complex_ai_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<1,
     complex_ai_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<0,
     complex_ar_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<1,
     complex_ar_t<T>>::type>::value );

    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<0,
     quaternion_ai_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<1,
     quaternion_ai_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<2,
     quaternion_ai_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<3,
     quaternion_ai_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<0,
     quaternion_ar_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<1,
     quaternion_ar_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<2,
     quaternion_ar_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<3,
     quaternion_ar_t<T>>::type>::value );

    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<0,
     octonion_ai_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<1,
     octonion_ai_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<2,
     octonion_ai_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<3,
     octonion_ai_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<4,
     octonion_ai_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<5,
     octonion_ai_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<6,
     octonion_ai_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<7,
     octonion_ai_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<0,
     octonion_ar_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<1,
     octonion_ar_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<2,
     octonion_ar_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<3,
     octonion_ar_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<4,
     octonion_ar_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<5,
     octonion_ar_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<6,
     octonion_ar_t<T>>::type>::value );
    SIMPLE_STATIC_ASSERT( is_same<T, typename tuple_element<7,
     octonion_ar_t<T>>::type>::value );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( tuple_access_tests, T, numeric_types )
{
    using boost::math::get;

    // Mutable l-value access
    real_ai_t<T>        test0a{ {1} };
    complex_ai_t<T>     test1a{ {2, 3} };
    quaternion_ai_t<T>  test2a{ {4, 5, 6, 7} };
    octonion_ai_t<T>    test3a{ {8, 9, 10, 11, 12, 13, 14, 15} };

    BOOST_CHECK_EQUAL( get<0>(test0a), (T)1 );
    BOOST_CHECK_EQUAL( get<0>(test1a), (T)2 );
    BOOST_CHECK_EQUAL( get<1>(test1a), (T)3 );
    BOOST_CHECK_EQUAL( get<0>(test2a), (T)4 );
    BOOST_CHECK_EQUAL( get<1>(test2a), (T)5 );
    BOOST_CHECK_EQUAL( get<2>(test2a), (T)6 );
    BOOST_CHECK_EQUAL( get<3>(test2a), (T)7 );
    BOOST_CHECK_EQUAL( get<0>(test3a), (T)8 );
    BOOST_CHECK_EQUAL( get<1>(test3a), (T)9 );
    BOOST_CHECK_EQUAL( get<2>(test3a), (T)10 );
    BOOST_CHECK_EQUAL( get<3>(test3a), (T)11 );
    BOOST_CHECK_EQUAL( get<4>(test3a), (T)12 );
    BOOST_CHECK_EQUAL( get<5>(test3a), (T)13 );
    BOOST_CHECK_EQUAL( get<6>(test3a), (T)14 );
    BOOST_CHECK_EQUAL( get<7>(test3a), (T)15 );

    // Immutable l-value access
    real_ai_t<T> const &        ctest0a = test0a;
    complex_ai_t<T> const &     ctest1a = test1a;
    quaternion_ai_t<T> const &  ctest2a = test2a;
    octonion_ai_t<T> const &    ctest3a = test3a;

    BOOST_CHECK_EQUAL( get<0>(ctest0a), (T)1 );
    BOOST_CHECK_EQUAL( get<0>(ctest1a), (T)2 );
    BOOST_CHECK_EQUAL( get<1>(ctest1a), (T)3 );
    BOOST_CHECK_EQUAL( get<0>(ctest2a), (T)4 );
    BOOST_CHECK_EQUAL( get<1>(ctest2a), (T)5 );
    BOOST_CHECK_EQUAL( get<2>(ctest2a), (T)6 );
    BOOST_CHECK_EQUAL( get<3>(ctest2a), (T)7 );
    BOOST_CHECK_EQUAL( get<0>(ctest3a), (T)8 );
    BOOST_CHECK_EQUAL( get<1>(ctest3a), (T)9 );
    BOOST_CHECK_EQUAL( get<2>(ctest3a), (T)10 );
    BOOST_CHECK_EQUAL( get<3>(ctest3a), (T)11 );
    BOOST_CHECK_EQUAL( get<4>(ctest3a), (T)12 );
    BOOST_CHECK_EQUAL( get<5>(ctest3a), (T)13 );
    BOOST_CHECK_EQUAL( get<6>(ctest3a), (T)14 );
    BOOST_CHECK_EQUAL( get<7>(ctest3a), (T)15 );

    // (Mutable) r-value access
    auto  rtest0a = [=] () -> decltype(test0a) { return test0a; };
    auto  rtest1a = [=] () -> decltype(test1a) { return test1a; };
    auto  rtest2a = [=] () -> decltype(test2a) { return test2a; };
    auto  rtest3a = [=] () -> decltype(test3a) { return test3a; };

    BOOST_CHECK_EQUAL( get<0>(rtest0a()), (T)1 );
    BOOST_CHECK_EQUAL( get<0>(rtest1a()), (T)2 );
    BOOST_CHECK_EQUAL( get<1>(rtest1a()), (T)3 );
    BOOST_CHECK_EQUAL( get<0>(rtest2a()), (T)4 );
    BOOST_CHECK_EQUAL( get<1>(rtest2a()), (T)5 );
    BOOST_CHECK_EQUAL( get<2>(rtest2a()), (T)6 );
    BOOST_CHECK_EQUAL( get<3>(rtest2a()), (T)7 );
    BOOST_CHECK_EQUAL( get<0>(rtest3a()), (T)8 );
    BOOST_CHECK_EQUAL( get<1>(rtest3a()), (T)9 );
    BOOST_CHECK_EQUAL( get<2>(rtest3a()), (T)10 );
    BOOST_CHECK_EQUAL( get<3>(rtest3a()), (T)11 );
    BOOST_CHECK_EQUAL( get<4>(rtest3a()), (T)12 );
    BOOST_CHECK_EQUAL( get<5>(rtest3a()), (T)13 );
    BOOST_CHECK_EQUAL( get<6>(rtest3a()), (T)14 );
    BOOST_CHECK_EQUAL( get<7>(rtest3a()), (T)15 );

    // Repeat
    real_ar_t<T>        test0b{ {16} };
    complex_ar_t<T>     test1b{ {{ {17} }, { {18} }} };
    quaternion_ar_t<T>  test2b{ {{{{{19}}, {{20}}}}, {{{{21}}, {{22}}}}} };
    octonion_ar_t<T>    test3b{ {{{{{{{23}}, {{24}}}}, {{{{25}}, {{26}}}}}},
     {{{{{{27}}, {{28}}}}, {{{{29}}, {{30}}}}}}} };

    real_ar_t<T> const &        ctest0b = test0b;
    complex_ar_t<T> const &     ctest1b = test1b;
    quaternion_ar_t<T> const &  ctest2b = test2b;
    octonion_ar_t<T> const &    ctest3b = test3b;

    auto  rtest0b = [=] () -> decltype(test0b) { return test0b; };
    auto  rtest1b = [=] () -> decltype(test1b) { return test1b; };
    auto  rtest2b = [=] () -> decltype(test2b) { return test2b; };
    auto  rtest3b = [=] () -> decltype(test3b) { return test3b; };

    BOOST_CHECK_EQUAL( get<0>(test0b), (T)16 );
    BOOST_CHECK_EQUAL( get<0>(test1b), (T)17 );
    BOOST_CHECK_EQUAL( get<1>(test1b), (T)18 );
    BOOST_CHECK_EQUAL( get<0>(test2b), (T)19 );
    BOOST_CHECK_EQUAL( get<1>(test2b), (T)20 );
    BOOST_CHECK_EQUAL( get<2>(test2b), (T)21 );
    BOOST_CHECK_EQUAL( get<3>(test2b), (T)22 );
    BOOST_CHECK_EQUAL( get<0>(test3b), (T)23 );
    BOOST_CHECK_EQUAL( get<1>(test3b), (T)24 );
    BOOST_CHECK_EQUAL( get<2>(test3b), (T)25 );
    BOOST_CHECK_EQUAL( get<3>(test3b), (T)26 );
    BOOST_CHECK_EQUAL( get<4>(test3b), (T)27 );
    BOOST_CHECK_EQUAL( get<5>(test3b), (T)28 );
    BOOST_CHECK_EQUAL( get<6>(test3b), (T)29 );
    BOOST_CHECK_EQUAL( get<7>(test3b), (T)30 );

    BOOST_CHECK_EQUAL( get<0>(ctest0b), (T)16 );
    BOOST_CHECK_EQUAL( get<0>(ctest1b), (T)17 );
    BOOST_CHECK_EQUAL( get<1>(ctest1b), (T)18 );
    BOOST_CHECK_EQUAL( get<0>(ctest2b), (T)19 );
    BOOST_CHECK_EQUAL( get<1>(ctest2b), (T)20 );
    BOOST_CHECK_EQUAL( get<2>(ctest2b), (T)21 );
    BOOST_CHECK_EQUAL( get<3>(ctest2b), (T)22 );
    BOOST_CHECK_EQUAL( get<0>(ctest3b), (T)23 );
    BOOST_CHECK_EQUAL( get<1>(ctest3b), (T)24 );
    BOOST_CHECK_EQUAL( get<2>(ctest3b), (T)25 );
    BOOST_CHECK_EQUAL( get<3>(ctest3b), (T)26 );
    BOOST_CHECK_EQUAL( get<4>(ctest3b), (T)27 );
    BOOST_CHECK_EQUAL( get<5>(ctest3b), (T)28 );
    BOOST_CHECK_EQUAL( get<6>(ctest3b), (T)29 );
    BOOST_CHECK_EQUAL( get<7>(ctest3b), (T)30 );

    BOOST_CHECK_EQUAL( get<0>(rtest0b()), (T)16 );
    BOOST_CHECK_EQUAL( get<0>(rtest1b()), (T)17 );
    BOOST_CHECK_EQUAL( get<1>(rtest1b()), (T)18 );
    BOOST_CHECK_EQUAL( get<0>(rtest2b()), (T)19 );
    BOOST_CHECK_EQUAL( get<1>(rtest2b()), (T)20 );
    BOOST_CHECK_EQUAL( get<2>(rtest2b()), (T)21 );
    BOOST_CHECK_EQUAL( get<3>(rtest2b()), (T)22 );
    BOOST_CHECK_EQUAL( get<0>(rtest3b()), (T)23 );
    BOOST_CHECK_EQUAL( get<1>(rtest3b()), (T)24 );
    BOOST_CHECK_EQUAL( get<2>(rtest3b()), (T)25 );
    BOOST_CHECK_EQUAL( get<3>(rtest3b()), (T)26 );
    BOOST_CHECK_EQUAL( get<4>(rtest3b()), (T)27 );
    BOOST_CHECK_EQUAL( get<5>(rtest3b()), (T)28 );
    BOOST_CHECK_EQUAL( get<6>(rtest3b()), (T)29 );
    BOOST_CHECK_EQUAL( get<7>(rtest3b()), (T)30 );
}

BOOST_AUTO_TEST_SUITE_END()  // tuple_cdhc_tests

BOOST_AUTO_TEST_CASE_TEMPLATE( iteration_tests, T, numeric_types )
{
    using boost::math::get;

    octonion_ar_t<T>          test1{};
    octonion_ar_t<T> const &  ctest1 = test1;
    unsigned                  counter = 0u;

    ctest1.iterate( [&counter](T const &x){counter += !x;} );
    BOOST_CHECK_EQUAL( counter, ctest1.size() );

    counter = 1u;
    test1.iterate( [&counter](T &x){x = counter * counter; ++counter;} );
    BOOST_CHECK_EQUAL( get<0>(ctest1), (T)1 );
    BOOST_CHECK_EQUAL( get<1>(ctest1), (T)4 );
    BOOST_CHECK_EQUAL( get<2>(ctest1), (T)9 );
    BOOST_CHECK_EQUAL( get<3>(ctest1), (T)16 );
    BOOST_CHECK_EQUAL( get<4>(ctest1), (T)25 );
    BOOST_CHECK_EQUAL( get<5>(ctest1), (T)36 );
    BOOST_CHECK_EQUAL( get<6>(ctest1), (T)49 );
    BOOST_CHECK_EQUAL( get<7>(ctest1), (T)64 );
}

BOOST_AUTO_TEST_SUITE( conversion_tests )

BOOST_AUTO_TEST_CASE_TEMPLATE( same_element_type_upsize_conversion_test, T,
 numeric_types )
{
    // Check if a higher-rung conversion zero-fills the extra spots
    T const                   v = get_random_number<T>();
    real_ai_t<T> const        test1a{ {v} };
    complex_ai_t<T> const     test2a = test1a;
    quaternion_ai_t<T> const  test3a_1 = test1a, test3a_2 = test2a;
    octonion_ai_t<T> const    test4a_1 = test1a, test4a_2 = test2a,
     test4a_3 = test3a_1;

    BOOST_CHECK_EQUAL( test2a, complex_ai_t<T>{{ v }} );
    BOOST_CHECK_EQUAL( test3a_1, quaternion_ai_t<T>{{ v }} );
    BOOST_CHECK_EQUAL( test3a_2, quaternion_ai_t<T>{{ v }} );
    BOOST_CHECK_EQUAL( test4a_1, octonion_ai_t<T>{{ v }} );
    BOOST_CHECK_EQUAL( test4a_2, octonion_ai_t<T>{{ v }} );
    BOOST_CHECK_EQUAL( test4a_3, octonion_ai_t<T>{{ v }} );

    // Repeat
    real_ar_t<T> const        test1b{ {v} };
    complex_ar_t<T> const     test2b = test1b;
    quaternion_ar_t<T> const  test3b_1 = test1b, test3b_2 = test2b;
    octonion_ar_t<T> const    test4b_1 = test1b, test4b_2 = test2b,
     test4b_3 = test3b_2;

    BOOST_CHECK_EQUAL( test2b, complex_ar_t<T>{{ {{ v }} }} );
    BOOST_CHECK_EQUAL( test3b_1, quaternion_ar_t<T>{{ {{ {{ v }} }} }} );
    BOOST_CHECK_EQUAL( test3b_2, quaternion_ar_t<T>{{ {{ {{ v }} }} }} );
    BOOST_CHECK_EQUAL( test4b_1, octonion_ar_t<T>{{ {{ {{ {{ v }} }} }} }} );
    BOOST_CHECK_EQUAL( test4b_2, octonion_ar_t<T>{{ {{ {{ {{ v }} }} }} }} );
    BOOST_CHECK_EQUAL( test4b_3, octonion_ar_t<T>{{ {{ {{ {{ v }} }} }} }} );
}

BOOST_AUTO_TEST_CASE( diff_elements_same_size_conversion_test )
{
    using boost::math::get;

    // Check if elements convert
    real_ai_t<unsigned short> const  test1a{ {6} };
    real_ai_t<unsigned> const        test2a = test1a;
    real_ai_t<int> const             test3a = test1a;
    real_ai_t<long> const            test4a = test1a;
    real_ai_t<float> const           test5a = test1a;
    real_ai_t<double> const          test6a = test1a;

    BOOST_CHECK_EQUAL( get<0>(test2a), 6u );
    BOOST_CHECK_EQUAL( get<0>(test3a), 6 );
    BOOST_CHECK_EQUAL( get<0>(test4a), 6L );
    BOOST_CHECK_EQUAL( get<0>(test5a), 6.0f );
    BOOST_CHECK_EQUAL( get<0>(test6a), 6.0 );

    complex_ai_t<unsigned short> const  test7a{ {2, 5} };
    complex_ai_t<unsigned> const        test8a = test7a;
    complex_ai_t<int> const             test9a = test7a;
    complex_ai_t<long> const            test10a = test7a;
    complex_ai_t<float> const           test11a = test7a;
    complex_ai_t<double> const          test12a = test7a;

    BOOST_CHECK_EQUAL( get<0>(test8a), 2u );
    BOOST_CHECK_EQUAL( get<0>(test9a), 2 );
    BOOST_CHECK_EQUAL( get<0>(test10a), 2L );
    BOOST_CHECK_EQUAL( get<0>(test11a), 2.0f );
    BOOST_CHECK_EQUAL( get<0>(test12a), 2.0 );
    BOOST_CHECK_EQUAL( get<1>(test8a), 5u );
    BOOST_CHECK_EQUAL( get<1>(test9a), 5 );
    BOOST_CHECK_EQUAL( get<1>(test10a), 5L );
    BOOST_CHECK_EQUAL( get<1>(test11a), 5.0f );
    BOOST_CHECK_EQUAL( get<1>(test12a), 5.0 );

    // Repeat
    real_ar_t<unsigned short> const  test1b{ {7} };
    real_ar_t<unsigned> const        test2b = test1b;
    real_ar_t<int> const             test3b = test1b;
    real_ar_t<long> const            test4b = test1b;
    real_ar_t<float> const           test5b = test1b;
    real_ar_t<double> const          test6b = test1b;

    BOOST_CHECK_EQUAL( get<0>(test2b), 7u );
    BOOST_CHECK_EQUAL( get<0>(test3b), 7 );
    BOOST_CHECK_EQUAL( get<0>(test4b), 7L );
    BOOST_CHECK_EQUAL( get<0>(test5b), 7.0f );
    BOOST_CHECK_EQUAL( get<0>(test6b), 7.0 );

    complex_ar_t<unsigned short> const  test7b{ {8, 15} };
    complex_ar_t<unsigned> const        test8b = test7b;
    complex_ar_t<int> const             test9b = test7b;
    complex_ar_t<long> const            test10b = test7b;
    complex_ar_t<float> const           test11b = test7b;
    complex_ar_t<double> const          test12b = test7b;

    BOOST_CHECK_EQUAL( get<0>(test8b), 8u );
    BOOST_CHECK_EQUAL( get<0>(test9b), 8 );
    BOOST_CHECK_EQUAL( get<0>(test10b), 8L );
    BOOST_CHECK_EQUAL( get<0>(test11b), 8.0f );
    BOOST_CHECK_EQUAL( get<0>(test12b), 8.0 );
    BOOST_CHECK_EQUAL( get<1>(test8b), 15u );
    BOOST_CHECK_EQUAL( get<1>(test9b), 15 );
    BOOST_CHECK_EQUAL( get<1>(test10b), 15L );
    BOOST_CHECK_EQUAL( get<1>(test11b), 15.0f );
    BOOST_CHECK_EQUAL( get<1>(test12b), 15.0 );

    // Check again, with larger tuples
    quaternion_ai_t<int> const    test13a{ {-4, 3, 10, -21} };
    quaternion_ai_t<long> const   test14a = test13a;
    quaternion_ai_t<float> const  test15a = test13a;

    BOOST_CHECK_EQUAL( get<0>(test14a), -4L );
    BOOST_CHECK_EQUAL( get<1>(test14a), 3L );
    BOOST_CHECK_EQUAL( get<2>(test14a), 10L );
    BOOST_CHECK_EQUAL( get<3>(test14a), -21L );
    BOOST_CHECK_EQUAL( get<0>(test15a), -4.0f );
    BOOST_CHECK_EQUAL( get<1>(test15a), 3.0f );
    BOOST_CHECK_EQUAL( get<2>(test15a), 10.0f );
    BOOST_CHECK_EQUAL( get<3>(test15a), -21.0f );

    octonion_ar_t<float> const   test16b{{ {{ {{ {{ 30.5f }}, {{ -31.25f }} }},
     {{ {{ 32.0f }}, {{ -33.75f }} }} }}, {{ {{ {{ 34.625f }}, {{ -35.5f }} }},
     {{ {{ 36.125f }}, {{ -37.25f }} }} }} }};
    octonion_ar_t<double> const  test17b = test16b;

    BOOST_CHECK_EQUAL( get<0>(test17b), +30.5 );
    BOOST_CHECK_EQUAL( get<1>(test17b), -31.25 );
    BOOST_CHECK_EQUAL( get<2>(test17b), +32.0 );
    BOOST_CHECK_EQUAL( get<3>(test17b), -33.75 );
    BOOST_CHECK_EQUAL( get<4>(test17b), +34.625 );
    BOOST_CHECK_EQUAL( get<5>(test17b), -35.5 );
    BOOST_CHECK_EQUAL( get<6>(test17b), +36.125 );
    BOOST_CHECK_EQUAL( get<7>(test17b), -37.25 );

    // Check with possibly-narrowing conversions, which are still implicit
    // (The ones with "v" prevent constexpr elision of narrowing.)
    double const             v = get_random_number<double>();
    real_ai_t<double> const  test18a{ {-6.25} }, test19a{ {v} };
    real_ai_t<float> const   test20a = test18a, test21a = test19a;
    real_ar_t<double> const  test18b{ {-6.25} }, test19b{ {v} };
    real_ar_t<float> const   test20b = test18b, test21b = test19b;

    BOOST_CHECK_EQUAL( get<0>(test20a), -6.25f );
    BOOST_CHECK_EQUAL( get<0>(test21a), static_cast<float>(v) );
    BOOST_CHECK_EQUAL( get<0>(test20b), -6.25f );
    BOOST_CHECK_EQUAL( get<0>(test21b), static_cast<float>(v) );
}

BOOST_AUTO_TEST_CASE( diff_elements_and_size_conversion_test )
{
    using boost::math::get;

    // Check if conversion and zero-fill work, with non-adjacent rungs
    real_ai_t<long> const        test1a{ {-5L} };
    octonion_ai_t<double> const  test2a = test1a;

    BOOST_CHECK_EQUAL( get<0>(test2a), -5.0 );
    BOOST_CHECK_EQUAL( get<1>(test2a),  0.0 );
    BOOST_CHECK_EQUAL( get<2>(test2a),  0.0 );
    BOOST_CHECK_EQUAL( get<3>(test2a),  0.0 );
    BOOST_CHECK_EQUAL( get<4>(test2a),  0.0 );
    BOOST_CHECK_EQUAL( get<5>(test2a),  0.0 );
    BOOST_CHECK_EQUAL( get<6>(test2a),  0.0 );
    BOOST_CHECK_EQUAL( get<7>(test2a),  0.0 );

    // Check with two adjacent non-real rungs
    complex_ai_t<int> const       test3a{ {8, -9} };
    quaternion_ai_t<float> const  test4a = test3a;

    BOOST_CHECK_EQUAL( get<0>(test4a), +8.0f );
    BOOST_CHECK_EQUAL( get<1>(test4a), -9.0f );
    BOOST_CHECK_EQUAL( get<2>(test4a),  0.0f );
    BOOST_CHECK_EQUAL( get<3>(test4a),  0.0f );

    // Repeat
    real_ar_t<float> const       test1b{ {-5.5f} };
    octonion_ar_t<double> const  test2b = test1b;

    BOOST_CHECK_EQUAL( get<0>(test2b), -5.5 );
    BOOST_CHECK_EQUAL( get<1>(test2b),  0.0 );
    BOOST_CHECK_EQUAL( get<2>(test2b),  0.0 );
    BOOST_CHECK_EQUAL( get<3>(test2b),  0.0 );
    BOOST_CHECK_EQUAL( get<4>(test2b),  0.0 );
    BOOST_CHECK_EQUAL( get<5>(test2b),  0.0 );
    BOOST_CHECK_EQUAL( get<6>(test2b),  0.0 );
    BOOST_CHECK_EQUAL( get<7>(test2b),  0.0 );

    complex_ar_t<unsigned> const  test3b{ {{ {7u} }, { {12u} }} };
    quaternion_ar_t<long> const   test4b = test3b;

    BOOST_CHECK_EQUAL( get<0>(test4b), 7L );
    BOOST_CHECK_EQUAL( get<1>(test4b), 12L );
    BOOST_CHECK_EQUAL( get<2>(test4b), 0L );
    BOOST_CHECK_EQUAL( get<3>(test4b), 0L );

    // Check with possibly-narrowing conversions, which are still implicit
    // (The ones with "v" prevent constexpr elision of narrowing.)
    double const             v1 = get_random_number<double>(),
     v2 = get_random_number<double>();
    real_ai_t<double> const        test5a{ {-7.0} }, test6a{ {v1} };
    complex_ai_t<float> const      test7a = test5a, test8a = test6a;
    complex_ar_t<double> const     test5b{ {{ {+3.5} }, { {-6.25} }} },
     test6b{ {{ {v1} }, { {v2} }} };
    quaternion_ar_t<float> const   test7b = test5b, test8b = test6b;

    BOOST_CHECK_EQUAL( get<0>(test7a), -7.0f );
    BOOST_CHECK_EQUAL( get<1>(test7a),  0.0f );
    BOOST_CHECK_EQUAL( get<0>(test8a), static_cast<float>(v1) );
    BOOST_CHECK_EQUAL( get<1>(test8a), 0.0f );
    BOOST_CHECK_EQUAL( get<0>(test7b), +3.50f );
    BOOST_CHECK_EQUAL( get<1>(test7b), -6.25f );
    BOOST_CHECK_EQUAL( get<2>(test7b),  0.00f );
    BOOST_CHECK_EQUAL( get<3>(test7b),  0.00f );
    BOOST_CHECK_EQUAL( get<0>(test8b), static_cast<float>(v1) );
    BOOST_CHECK_EQUAL( get<1>(test8b), static_cast<float>(v2) );
    BOOST_CHECK_EQUAL( get<2>(test8b), 0.0f );
    BOOST_CHECK_EQUAL( get<3>(test8b), 0.0f );
}

BOOST_AUTO_TEST_CASE( explicit_conversion_test )
{
    using boost::math::get;

    // Use a conversion that needs to be explicit
    two_bit_t const  v = two_bit_t::two;
    auto const       test1a = real_ai_t<two_bit_t>{ {v} };
    auto const       test2a = static_cast<real_ai_t<unsigned>>( test1a );
     // not "{ test1a }" nor "= test1a"

    BOOST_REQUIRE_EQUAL( static_cast<unsigned>(v), 2u );
    BOOST_CHECK_EQUAL( get<0>(test2a), 2u );

    // Convert to a smaller tuple
    quaternion_ai_t<int> const  test3a{ {1, 2, 3, 4} };
    auto const                  test4a = static_cast<complex_ai_t<int>>(test3a);

    BOOST_CHECK_EQUAL( get<0>(test4a), 1 );
    BOOST_CHECK_EQUAL( get<1>(test4a), 2 );

    // Convert to smaller tuple with different type
    complex_ai_t<int> const  test5a{ {-7, -2} };
    auto const               test6a = static_cast<real_ai_t<double>>( test5a );

    BOOST_CHECK_EQUAL( get<0>(test6a), -7.0 );

    // Same, but with explicit element conversion
    octonion_ai_t<two_bit_t> const  test7a{ {two_bit_t::three, two_bit_t::two,
     two_bit_t::one, two_bit_t::zero, two_bit_t::zero, two_bit_t::one,
     two_bit_t::two, two_bit_t::three} };
    auto const         test8a = static_cast<complex_ai_t<unsigned>>( test7a );

    BOOST_CHECK_EQUAL( get<0>(test8a), 3u );
    BOOST_CHECK_EQUAL( get<1>(test8a), 2u );

    // Use explicit element conversion with longer tuple
    complex_ai_t<two_bit_t> const  test9a{ {two_bit_t::zero, two_bit_t::one} };
    auto const      test10a = static_cast<quaternion_ai_t<unsigned>>( test9a );

    BOOST_CHECK_EQUAL( get<0>(test10a), 0u );
    BOOST_CHECK_EQUAL( get<1>(test10a), 1u );
    BOOST_CHECK_EQUAL( get<2>(test10a), 0u );
    BOOST_CHECK_EQUAL( get<3>(test10a), 0u );

    // Repeat w/ recursively-defined types; base case
    auto const  test1b = real_ar_t<two_bit_t>{ {v} };
    auto const  test2b = static_cast<real_ar_t<unsigned>>( test1b );

    BOOST_CHECK_EQUAL( get<0>(test2b), 2u );

    // Explicit conversion, base case to non-base case
    auto const  test3b = static_cast<quaternion_ar_t<unsigned>>( test1b );

    BOOST_CHECK_EQUAL( get<0>(test3b), 2u );
    BOOST_CHECK_EQUAL( get<1>(test3b), 0u );
    BOOST_CHECK_EQUAL( get<2>(test3b), 0u );
    BOOST_CHECK_EQUAL( get<3>(test3b), 0u );

    // Explicit conversion, two same-size non-base cases
    octonion_ar_t<two_bit_t> const  test4b{ {{ {{ {{ {two_bit_t::three} },
     { {two_bit_t::two} }} }, { {{ {two_bit_t::one} },
     { {two_bit_t::zero} }} }} }, { {{ {{ {two_bit_t::zero} },
     { {two_bit_t::one} }} }, { {{ {two_bit_t::two} },
     { {two_bit_t::three} }} }} }} };
    auto const    test5b = static_cast<octonion_ar_t<unsigned>>( test4b );

    BOOST_CHECK_EQUAL( get<0>(test5b), 3u );
    BOOST_CHECK_EQUAL( get<1>(test5b), 2u );
    BOOST_CHECK_EQUAL( get<2>(test5b), 1u );
    BOOST_CHECK_EQUAL( get<3>(test5b), 0u );
    BOOST_CHECK_EQUAL( get<4>(test5b), 0u );
    BOOST_CHECK_EQUAL( get<5>(test5b), 1u );
    BOOST_CHECK_EQUAL( get<6>(test5b), 2u );
    BOOST_CHECK_EQUAL( get<7>(test5b), 3u );

    // Explicit conversion, non-base case to larger one
    auto const  test6b = static_cast<octonion_ar_t<unsigned>>( test4b.b[1] );

    BOOST_CHECK_EQUAL( get<0>(test6b), 0u );
    BOOST_CHECK_EQUAL( get<1>(test6b), 1u );
    BOOST_CHECK_EQUAL( get<2>(test6b), 2u );
    BOOST_CHECK_EQUAL( get<3>(test6b), 3u );
    BOOST_CHECK_EQUAL( get<4>(test6b), 0u );
    BOOST_CHECK_EQUAL( get<5>(test6b), 0u );
    BOOST_CHECK_EQUAL( get<6>(test6b), 0u );
    BOOST_CHECK_EQUAL( get<7>(test6b), 0u );

    // Conversion, same element type but smaller sizes
    auto const  test7b = static_cast<quaternion_ar_t<unsigned>>( test6b );
    auto const  test8b = static_cast<complex_ar_t<unsigned>>( test6b );
    auto const  test9b = static_cast<real_ar_t<unsigned>>( test6b );

    BOOST_CHECK_EQUAL( get<0>(test7b), 0u );
    BOOST_CHECK_EQUAL( get<1>(test7b), 1u );
    BOOST_CHECK_EQUAL( get<2>(test7b), 2u );
    BOOST_CHECK_EQUAL( get<3>(test7b), 3u );

    BOOST_CHECK_EQUAL( get<0>(test8b), 0u );
    BOOST_CHECK_EQUAL( get<1>(test8b), 1u );
    BOOST_CHECK_EQUAL( get<0>(test9b), 0u );

    // Conversion, implicitly-convertible element type but smaller sizes
    auto const  test10b = static_cast<quaternion_ar_t<int>>( test6b );
    auto const  test11b = static_cast<complex_ar_t<int>>( test6b );
    auto const  test12b = static_cast<real_ar_t<int>>( test6b );

    BOOST_CHECK_EQUAL( get<0>(test10b), 0 );
    BOOST_CHECK_EQUAL( get<1>(test10b), 1 );
    BOOST_CHECK_EQUAL( get<2>(test10b), 2 );
    BOOST_CHECK_EQUAL( get<3>(test10b), 3 );
    BOOST_CHECK_EQUAL( get<0>(test11b), 0 );
    BOOST_CHECK_EQUAL( get<1>(test11b), 1 );
    BOOST_CHECK_EQUAL( get<0>(test12b), 0 );

    // Explicit conversion, to smaller sizes
    auto const  test13b = static_cast<quaternion_ar_t<unsigned>>( test4b );
    auto const  test14b = static_cast<complex_ar_t<unsigned>>( test4b );
    auto const  test15b = static_cast<real_ar_t<unsigned>>( test4b );

    BOOST_CHECK_EQUAL( get<0>(test13b), 3u );
    BOOST_CHECK_EQUAL( get<1>(test13b), 2u );
    BOOST_CHECK_EQUAL( get<2>(test13b), 1u );
    BOOST_CHECK_EQUAL( get<3>(test13b), 0u );
    BOOST_CHECK_EQUAL( get<0>(test14b), 3u );
    BOOST_CHECK_EQUAL( get<1>(test14b), 2u );
    BOOST_CHECK_EQUAL( get<0>(test15b), 3u );
}

BOOST_AUTO_TEST_SUITE_END()  // conversion_tests
