//  Boost Complex Numbers, recursive-mode unit test program file  ------------//

//  Copyright 2013 Daryle Walker.
//  Distributed under the Boost Software License, Version 1.0.  (See the
//  accompanying file LICENSE_1_0.txt or a copy at
//  <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/output_test_stream.hpp>
#include <boost/mpl/list.hpp>

#include "boost/math/complex_rt.hpp"
#include "boost/math/complex_it.hpp"

#include <cstddef>
#include <ios>
#include <tuple>
#include <type_traits>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>


// Common definitions  -------------------------------------------------------//

namespace {

    // Save time writing out long types
    namespace mp = boost::multiprecision;
    using boost::mpl::list;
    using boost::test_tools::output_test_stream;
    using boost::math::complex_rt;

    // Sample testing types for components
    typedef mp::number<mp::cpp_dec_float<50>, mp::et_off>          my_float;
    typedef list<int, unsigned, double, mp::int512_t, my_float>  test_types;
    typedef list<int, unsigned, mp::int512_t>            test_integer_types;
    typedef list<double, my_float>                      test_floating_types;
    typedef list<int, unsigned, double>                  test_builtin_types;

}

// Flag un-printable types here.


BOOST_AUTO_TEST_SUITE( complex_rt_tests )

BOOST_AUTO_TEST_SUITE( core_tests )

// Check the various compile-time attributes.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_complex_compile_time, T, test_types )
{
    using std::is_same;
    using std::size_t;

    typedef complex_rt<T, 0>  real_type;
    typedef complex_rt<T, 1>  complex_type;
    typedef complex_rt<T, 2>  quaternion_type;
    typedef complex_rt<T, 3>  octonion_type;

    // Template parameters
    BOOST_REQUIRE( (is_same<typename real_type::value_type, T>::value) );
    BOOST_REQUIRE( (is_same<typename complex_type::value_type, T>::value) );
    BOOST_REQUIRE( (is_same<typename quaternion_type::value_type, T>::value) );
    BOOST_REQUIRE( (is_same<typename octonion_type::value_type, T>::value) );

    BOOST_REQUIRE_EQUAL( real_type::rank, 0u );
    BOOST_REQUIRE_EQUAL( complex_type::rank, 1u );
    BOOST_REQUIRE_EQUAL( quaternion_type::rank, 2u );
    BOOST_REQUIRE_EQUAL( octonion_type::rank, 3u );

    // Support types and values
    BOOST_REQUIRE( (is_same<typename real_type::size_type, size_t>::value) );
    BOOST_REQUIRE( (is_same<typename complex_type::size_type, size_t>::value) );
    BOOST_REQUIRE( (is_same<typename quaternion_type::size_type,
     size_t>::value) );
    BOOST_REQUIRE( (is_same<typename octonion_type::size_type,
     size_t>::value) );

    BOOST_REQUIRE_EQUAL( real_type::static_size, 1u );
    BOOST_REQUIRE_EQUAL( complex_type::static_size, 2u );
    BOOST_REQUIRE_EQUAL( quaternion_type::static_size, 4u );
    BOOST_REQUIRE_EQUAL( octonion_type::static_size, 8u );

    BOOST_REQUIRE( (is_same<typename complex_type::barrage_type,
     real_type>::value) );
    BOOST_REQUIRE( (is_same<typename quaternion_type::barrage_type,
     complex_type>::value) );
    BOOST_REQUIRE( (is_same<typename octonion_type::barrage_type,
     quaternion_type>::value) );
}

// Check the most basic operation, component-level access.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_complex_component_access1, T,
 test_integer_types )
{
    typedef boost::rational<T>            rational_type;
    typedef complex_rt<T, 0>                  real_type;
    typedef complex_rt<T, 1>               complex_type;
    typedef complex_rt<T, 2>            quaternion_type;
    typedef complex_rt<rational_type, 3>  octonion_type;

    real_type               a;
    real_type const &       aa = a;
    complex_type            b;
    complex_type const &    bb = b;
    quaternion_type         c;
    quaternion_type const & cc = c;
    octonion_type           d;
    octonion_type const &   dd = d;

    a[ 0 ] = 6;
    BOOST_CHECK_EQUAL( a[0], T(6) );
    BOOST_CHECK_EQUAL( T(6), aa[0] );

    b[ 0 ] = 5;
    b[ 1 ] = 7;
    BOOST_CHECK_EQUAL( b[0], T(5) );
    BOOST_CHECK_EQUAL( T(5), bb[0] );
    BOOST_CHECK_EQUAL( b[1], T(7) );
    BOOST_CHECK_EQUAL( T(7), bb[1] );

    c[ 0 ] = 10;
    c[ 1 ] = 11;
    c[ 2 ] = 12;
    c[ 3 ] = 13;
    BOOST_CHECK_EQUAL( c[0], T(10) );
    BOOST_CHECK_EQUAL( T(10), cc[0] );
    BOOST_CHECK_EQUAL( c[1], T(11) );
    BOOST_CHECK_EQUAL( T(11), cc[1] );
    BOOST_CHECK_EQUAL( c[2], T(12) );
    BOOST_CHECK_EQUAL( T(12), cc[2] );
    BOOST_CHECK_EQUAL( c[3], T(13) );
    BOOST_CHECK_EQUAL( T(13), cc[3] );

    d[ 0 ] = rational_type( 1 );
    d[ 1 ] = rational_type{};
    d[ 2 ] = rational_type( 2, 3 );
    d[ 3 ] = rational_type( 3, 5 );
    d[ 4 ] = rational_type( 8, 13 );
    d[ 5 ] = rational_type( 5, 8 );
    d[ 6 ] = rational_type( 21, 34 );
    d[ 7 ] = rational_type( 89, 55);
    BOOST_CHECK_EQUAL( d[0].numerator(), T(1) );
    BOOST_CHECK_EQUAL( T(1), dd[0].denominator() );
    BOOST_CHECK_EQUAL( d[1].numerator(), T(0) );
    BOOST_CHECK_EQUAL( T(1), dd[1].denominator() );
    BOOST_CHECK_EQUAL( d[2].numerator(), T(2) );
    BOOST_CHECK_EQUAL( T(3), dd[2].denominator() );
    BOOST_CHECK_EQUAL( d[3].numerator(), T(3) );
    BOOST_CHECK_EQUAL( T(5), dd[3].denominator() );
    BOOST_CHECK_EQUAL( d[4].numerator(), T(8) );
    BOOST_CHECK_EQUAL( T(13), dd[4].denominator() );
    BOOST_CHECK_EQUAL( d[5].numerator(), T(5) );
    BOOST_CHECK_EQUAL( T(8), dd[5].denominator() );
    BOOST_CHECK_EQUAL( d[6].numerator(), T(21) );
    BOOST_CHECK_EQUAL( T(34), dd[6].denominator() );
    BOOST_CHECK_EQUAL( d[7].numerator(), T(89) );
    BOOST_CHECK_EQUAL( T(55), dd[7].denominator() );
}

// Check the most basic operation, component-level access, floating-point.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_complex_component_access2, T,
 test_floating_types )
{
    typedef complex_rt<T, 0>     real_type;
    typedef complex_rt<T, 1>  complex_type;

    real_type               a;
    real_type const &       aa = a;
    complex_type            b;
    complex_type const &    bb = b;

    a[ 0 ] = 6.0;
    BOOST_CHECK_CLOSE( a[0], T(6.0), 0.1 );
    BOOST_CHECK_CLOSE( T(6.0), aa[0], 0.1 );

    b[ 0 ] = 5.5;
    b[ 1 ] = -7.0;
    BOOST_CHECK_CLOSE( b[0], T(5.5), 0.1 );
    BOOST_CHECK_CLOSE( T(5.5), bb[0], 0.1 );
    BOOST_CHECK_CLOSE( b[1], T(-7.0), 0.1 );
    BOOST_CHECK_CLOSE( T(-7.0), bb[1], 0.1 );
}

// Check Boolean conversion.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_complex_to_boolean, T, test_types )
{
    // Real support
    complex_rt<T, 0>  r;

    r[ 0 ] = T{};
    BOOST_CHECK( !r );
    r[ 0 ] = (T)2;
    BOOST_CHECK( (bool)r );

    // Check with multi-scalar
    complex_rt<T, 2>  q;

    q[ 0 ] = q[ 1 ] = q[ 2 ] = q[ 3 ] = T{};
    BOOST_CHECK( !q );
    q[ 2 ] = (T)3;
    BOOST_CHECK( (bool)q );
    q[ 3 ] = (T)5;
    BOOST_CHECK( (bool)q );
    q[ 2 ] = q[ 3 ] = (T)0;
    BOOST_CHECK( !q );
}

// Check barrage-level access.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_complex_barrages, T, test_types )
{
    // Degenerate case
    complex_rt<T, 0>  r;
    auto const &      rr = r;

    r[ 0 ] = (T)6;
    BOOST_CHECK_EQUAL( rr[0], rr.lower_barrage()[0] );
    BOOST_CHECK_EQUAL( rr[0], rr.upper_barrage()[0] );

    // (Regular) complex
    complex_rt<T, 1>  c;
    auto const &      cc = c;

    c[ 0 ] = (T)7;
    c[ 1 ] = (T)18;

    auto const  cc1 = cc.lower_barrage();
    auto const  cc2 = cc.upper_barrage();

    BOOST_CHECK_EQUAL( cc1[0], cc[0] );
    BOOST_CHECK_EQUAL( cc2[0], cc[1] );

    // Extreme
    complex_rt<T, 3>  o;
    auto const &      oo = o;

    o[ 0 ] = (T)0;
    o[ 1 ] = (T)1;
    o[ 2 ] = (T)4;
    o[ 3 ] = (T)9;
    o[ 4 ] = (T)16;
    o[ 5 ] = (T)25;
    o[ 6 ] = (T)36;
    o[ 7 ] = (T)49;

    auto const  oo1 = oo.lower_barrage();
    auto const  oo2 = oo.upper_barrage();

    BOOST_CHECK_EQUAL( oo1[0], oo[0] );
    BOOST_CHECK_EQUAL( oo1[1], oo[1] );
    BOOST_CHECK_EQUAL( oo1[2], oo[2] );
    BOOST_CHECK_EQUAL( oo1[3], oo[3] );
    BOOST_CHECK_EQUAL( oo2[0], oo[4] );
    BOOST_CHECK_EQUAL( oo2[1], oo[5] );
    BOOST_CHECK_EQUAL( oo2[2], oo[6] );
    BOOST_CHECK_EQUAL( oo2[3], oo[7] );

    // Mutability
    o.lower_barrage() = oo2;
    o.upper_barrage() = oo1;
    BOOST_CHECK_EQUAL( oo1[0], oo[4] );
    BOOST_CHECK_EQUAL( oo1[1], oo[5] );
    BOOST_CHECK_EQUAL( oo1[2], oo[6] );
    BOOST_CHECK_EQUAL( oo1[3], oo[7] );
    BOOST_CHECK_EQUAL( oo2[0], oo[0] );
    BOOST_CHECK_EQUAL( oo2[1], oo[1] );
    BOOST_CHECK_EQUAL( oo2[2], oo[2] );
    BOOST_CHECK_EQUAL( oo2[3], oo[3] );

    // Degenerate mutability
    complex_rt<T, 0>  s;

    s[ 0 ] = (T)63;
    r.lower_barrage() = s;
    BOOST_CHECK_EQUAL( r[0], r.lower_barrage()[0] );
    s[ 0 ] = (T)65;
    r.upper_barrage() = s;
    BOOST_CHECK_EQUAL( r[0], r.upper_barrage()[0] );
}

// Check comparisons between hypercomplex and scalars.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_complex_real_equality, T, test_types )
{
    T                 sample = 6;
    complex_rt<T, 0>  r;
    complex_rt<T, 1>  c;
    complex_rt<T, 2>  q;
    complex_rt<T, 3>  o;

    // Zero check
    r[ 0 ] = T{};
    BOOST_CHECK( r != sample );
    BOOST_CHECK( sample != r );
    BOOST_CHECK( !(r == sample) );
    BOOST_CHECK( !(sample == r) );

    c[ 0 ] = c[ 1 ] = T{};
    BOOST_CHECK( c != sample );
    BOOST_CHECK( sample != c );
    BOOST_CHECK( !(c == sample) );
    BOOST_CHECK( !(sample == c) );

    q[ 0 ] = q[ 1 ] = q[ 2 ] = q[ 3 ] = T{};
    BOOST_CHECK( q != sample );
    BOOST_CHECK( sample != q );
    BOOST_CHECK( !(q == sample) );
    BOOST_CHECK( !(sample == q) );

    o[ 0 ] = o[ 1 ] = o[ 2 ] = o[ 3 ] = T{};
    o[ 4 ] = o[ 5 ] = o[ 6 ] = o[ 7 ] = T{};
    BOOST_CHECK( o != sample );
    BOOST_CHECK( sample != o );
    BOOST_CHECK( !(o == sample) );
    BOOST_CHECK( !(sample == o) );

    // Equality check
    r[ 0 ] = c[ 0 ] = q[ 0 ] = o[ 0 ] = sample;

    BOOST_CHECK( r == sample );
    BOOST_CHECK( sample == r );
    BOOST_CHECK( !(r != sample) );
    BOOST_CHECK( !(sample != r) );

    BOOST_CHECK( c == sample );
    BOOST_CHECK( sample == c );
    BOOST_CHECK( !(c != sample) );
    BOOST_CHECK( !(sample != c) );

    BOOST_CHECK( q == sample );
    BOOST_CHECK( sample == q );
    BOOST_CHECK( !(q != sample) );
    BOOST_CHECK( !(sample != q) );

    BOOST_CHECK( o == sample );
    BOOST_CHECK( sample == o );
    BOOST_CHECK( !(o != sample) );
    BOOST_CHECK( !(sample != o) );

    // Change amoung non-real components
    c[ 1 ] = q[ 1 ] = o[ 1 ] = (T)1;

    BOOST_CHECK( c != sample );
    BOOST_CHECK( sample != c );
    BOOST_CHECK( !(c == sample) );
    BOOST_CHECK( !(sample == c) );

    BOOST_CHECK( q != sample );
    BOOST_CHECK( sample != q );
    BOOST_CHECK( !(q == sample) );
    BOOST_CHECK( !(sample == q) );

    BOOST_CHECK( o != sample );
    BOOST_CHECK( sample != o );
    BOOST_CHECK( !(o == sample) );
    BOOST_CHECK( !(sample == o) );

    q[ 1 ] = o[ 1 ] = T{};
    q[ 3 ] = o[ 3 ] = (T)2;

    BOOST_CHECK( q != sample );
    BOOST_CHECK( sample != q );
    BOOST_CHECK( !(q == sample) );
    BOOST_CHECK( !(sample == q) );

    BOOST_CHECK( o != sample );
    BOOST_CHECK( sample != o );
    BOOST_CHECK( !(o == sample) );
    BOOST_CHECK( !(sample == o) );
}

// Check comparisons between two hypercomplex numbers.
BOOST_AUTO_TEST_CASE( test_complex_equality )
{
    // Two reals, same component type
    complex_rt<int, 0>  a, b;

    a[ 0 ] = b[ 0 ] = 0;
    BOOST_CHECK( a == b );
    BOOST_CHECK( !(a != b) );
    a[ 0 ] = -2;
    b[ 0 ] = +3;
    BOOST_CHECK( a != b );
    BOOST_CHECK( !(a == b) );

    // Two reals, differing component types
    complex_rt<long, 0>  c;

    c[ 0 ] = -2L;
    BOOST_CHECK( a == c );
    BOOST_CHECK( !(a != c) );
    BOOST_CHECK( b != c );
    BOOST_CHECK( !(b == c) );

    // Two (regular) complex, same component type
    complex_rt<int, 1>  d, e;

    d[ 0 ] = d[ 1 ] = e[ 0 ] = e[ 1 ] = 0;
    BOOST_CHECK( d == e );
    BOOST_CHECK( !(d != e) );
    d[ 0 ] = -3;
    BOOST_CHECK( d != e );
    BOOST_CHECK( !(d == e) );
    e[ 0 ] = -3;
    e[ 1 ] = +2;
    BOOST_CHECK( d != e );
    BOOST_CHECK( !(d == e) );
    d[ 1 ] = +2;
    BOOST_CHECK( d == e );
    BOOST_CHECK( !(d != e) );

    // Two (regular) complex, differing component types
    complex_rt<long, 1>  f;

    f[ 0 ] = f[ 1 ] = 0L;
    BOOST_CHECK( d != f );
    BOOST_CHECK( !(d == f) );
    f[ 0 ] = d[ 0 ];
    BOOST_CHECK( d != f );
    BOOST_CHECK( !(d == f) );
    f[ 0 ] = f[ 1 ] = d[ 1 ];
    BOOST_CHECK( d != f );
    BOOST_CHECK( !(d == f) );
    f[ 0 ] = d[ 0 ];
    BOOST_CHECK( d == f );
    BOOST_CHECK( !(d != f) );

    // Mixed levels
    complex_rt<int, 2>  g;

    g[ 0 ] = g[ 1 ] = g[ 2 ] = g[ 3 ] = 0;
    BOOST_CHECK( e != g );
    BOOST_CHECK( !(e == g) );
    BOOST_CHECK( g != e );
    BOOST_CHECK( !(g == e) );
    g[ 0 ] = e[ 0 ];
    g[ 1 ] = e[ 1 ];
    BOOST_CHECK( e == g );
    BOOST_CHECK( !(e != g) );
    BOOST_CHECK( g == e );
    BOOST_CHECK( !(g != e) );
    ++g[ 1 ];
    BOOST_CHECK( e != g );
    BOOST_CHECK( !(e == g) );
    BOOST_CHECK( g != e );
    BOOST_CHECK( !(g == e) );
    --g[ 1 ];
    BOOST_CHECK( e == g );
    BOOST_CHECK( !(e != g) );
    BOOST_CHECK( g == e );
    BOOST_CHECK( !(g != e) );
    g[ 3 ] = 5;
    BOOST_CHECK( e != g );
    BOOST_CHECK( !(e == g) );
    BOOST_CHECK( g != e );
    BOOST_CHECK( !(g == e) );

    // Mixed levels, differing component types
    complex_rt<long, 3>  h;

    h[ 0 ] = h[ 1 ] = h[ 2 ] = h[ 3 ] = 0L;
    h[ 4 ] = h[ 5 ] = h[ 6 ] = h[ 7 ] = 0L;
    BOOST_CHECK( e != h );
    BOOST_CHECK( !(e == h) );
    BOOST_CHECK( h != e );
    BOOST_CHECK( !(h == e) );
    h[ 0 ] = g[ 0 ];
    h[ 1 ] = g[ 1 ];
    BOOST_CHECK( e == h );
    BOOST_CHECK( !(e != h) );
    BOOST_CHECK( h == e );
    BOOST_CHECK( !(h != e) );
    BOOST_CHECK( g != h );
    BOOST_CHECK( !(g == h) );
    BOOST_CHECK( h != g );
    BOOST_CHECK( !(h == g) );
    h[ 3 ] = g[ 3 ];
    BOOST_CHECK( e != h );
    BOOST_CHECK( !(e == h) );
    BOOST_CHECK( h != e );
    BOOST_CHECK( !(h == e) );
    BOOST_CHECK( g == h );
    BOOST_CHECK( !(g != h) );
    BOOST_CHECK( h == g );
    BOOST_CHECK( !(h != g) );
    h[ 6 ] = -7L;
    BOOST_CHECK( e != h );
    BOOST_CHECK( !(e == h) );
    BOOST_CHECK( h != e );
    BOOST_CHECK( !(h == e) );
    BOOST_CHECK( g != h );
    BOOST_CHECK( !(g == h) );
    BOOST_CHECK( h != g );
    BOOST_CHECK( !(h == g) );
}

// Check output routines.
BOOST_AUTO_TEST_CASE( test_complex_output )
{
    output_test_stream  ots;
    complex_rt<int, 0>  r;
    complex_rt<int, 1>  c;
    complex_rt<int, 2>  q;
    complex_rt<int, 3>  o;

    r[ 0 ] = 1;
    ots << r;
    BOOST_CHECK( ots.is_equal("1") );

    c[ 0 ] = 2;
    c[ 1 ] = 3;
    ots << c;
    BOOST_CHECK( ots.is_equal("(2,3)") );

    q[ 0 ] = -4;
    q[ 1 ] = +5;
    q[ 2 ] = -6;
    q[ 3 ] = +7;
    ots << std::showpos << q;
    BOOST_CHECK( ots.is_equal("((-4,+5),(-6,+7))") );

    o[ 0 ] = -10;
    o[ 1 ] = 11;
    o[ 2 ] = 12;
    o[ 3 ] = -13;
    o[ 4 ] = 14;
    o[ 5 ] = 15;
    o[ 6 ] = -16;
    o[ 7 ] = 101;
    ots << std::noshowpos << o;
    BOOST_CHECK( ots.is_equal("(((-10,11),(12,-13)),((14,15),(-16,101)))") );
}

BOOST_AUTO_TEST_SUITE_END()  // core_tests

BOOST_AUTO_TEST_SUITE( constructor_tests )

// Check the results of default- and scalar-construction.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_default_real_construction, T, test_types )
{
    // Real
    complex_rt<T, 0>  a = {}, b = { (T)2 };

    BOOST_CHECK_EQUAL( a[0], T{} );
    BOOST_CHECK_EQUAL( b[0], T(2) );

    // (Regular) complex
    complex_rt<T, 1>  c = {}, d = { (T)7 };

    BOOST_CHECK_EQUAL( c[0], T{} );
    BOOST_CHECK_EQUAL( c[1], T{} );
    BOOST_CHECK_EQUAL( d[0], T(7) );
    BOOST_CHECK_EQUAL( d[1], T{} );

    // Quaternions
    complex_rt<T, 2>  e = {}, f = { (T)19 };

    BOOST_CHECK_EQUAL( e[0], T{} );
    BOOST_CHECK_EQUAL( e[1], T{} );
    BOOST_CHECK_EQUAL( e[2], T{} );
    BOOST_CHECK_EQUAL( e[3], T{} );
    BOOST_CHECK_EQUAL( f[0], T(19) );
    BOOST_CHECK_EQUAL( f[1], T{} );
    BOOST_CHECK_EQUAL( f[2], T{} );
    BOOST_CHECK_EQUAL( f[3], T{} );

    // Octonions
    complex_rt<T, 3>  g = {}, h = { (T)101 };

    BOOST_CHECK_EQUAL( g[0], T{} );
    BOOST_CHECK_EQUAL( g[1], T{} );
    BOOST_CHECK_EQUAL( g[2], T{} );
    BOOST_CHECK_EQUAL( g[3], T{} );
    BOOST_CHECK_EQUAL( g[4], T{} );
    BOOST_CHECK_EQUAL( g[5], T{} );
    BOOST_CHECK_EQUAL( g[6], T{} );
    BOOST_CHECK_EQUAL( g[7], T{} );
    BOOST_CHECK_EQUAL( h[0], T(101) );
    BOOST_CHECK_EQUAL( h[1], T{} );
    BOOST_CHECK_EQUAL( h[2], T{} );
    BOOST_CHECK_EQUAL( h[3], T{} );
    BOOST_CHECK_EQUAL( h[4], T{} );
    BOOST_CHECK_EQUAL( h[5], T{} );
    BOOST_CHECK_EQUAL( h[6], T{} );
    BOOST_CHECK_EQUAL( h[7], T{} );
}

// Check conversion from complex_it objects.
BOOST_AUTO_TEST_CASE( test_cross_philosophy_construction )
{
    using boost::math::complex_it;

    // Real-to-real
    complex_rt<int, 0> const  a{ complex_it<int, 0>{5} };

    BOOST_CHECK_EQUAL( a[0], 5 );

    // Real-to-real, different component types
    complex_rt<long, 0> const  b{ complex_it<int, 0>{2} };

    BOOST_CHECK_EQUAL( b[0], 2L );

    // Downgrade to real
    complex_rt<int, 0> const  c{ complex_it<int, 1>{3, -7} };

    BOOST_CHECK_EQUAL( c[0], 3 );

    // Downgrade to real, different component types
    complex_rt<int, 0> const  d{ complex_it<long, 1>{-11, 13} };

    BOOST_CHECK_EQUAL( d[0], -11 );

    // Upgrade from real
    complex_rt<int, 1> const  e{ complex_it<int, 0>{17} };

    BOOST_CHECK_EQUAL( e[0], 17 );
    BOOST_CHECK_EQUAL( e[1], 0 );

    // Upgrade from real, different component types
    complex_rt<long, 2> const  f{ complex_it<int, 0>{-19} };

    BOOST_CHECK_EQUAL( f[0], -19L );
    BOOST_CHECK_EQUAL( f[1], 0L );
    BOOST_CHECK_EQUAL( f[2], 0L );
    BOOST_CHECK_EQUAL( f[3], 0L );

    // Same post-real rank
    complex_rt<int, 2> const  g{ complex_it<int, 2>{23, 29, 31, 37} };

    BOOST_CHECK_EQUAL( g[0], 23 );
    BOOST_CHECK_EQUAL( g[1], 29 );
    BOOST_CHECK_EQUAL( g[2], 31 );
    BOOST_CHECK_EQUAL( g[3], 37 );

    // Same post-real rank, different component types
    complex_rt<unsigned long, 2> const  h{ complex_it<int, 2>{41, 43, 47, 53} };

    BOOST_CHECK_EQUAL( h[0], 41UL );
    BOOST_CHECK_EQUAL( h[1], 43UL );
    BOOST_CHECK_EQUAL( h[2], 47UL );
    BOOST_CHECK_EQUAL( h[3], 53UL );

    // Downgrade post-real ranks
    complex_rt<int, 1> const  k{ complex_it<int, 2>{57, 59, 61, 67} };

    BOOST_CHECK_EQUAL( k[0], 57 );
    BOOST_CHECK_EQUAL( k[1], 59 );

    // Downgrade post-real ranks, different component types
    complex_rt<long, 1> const  m{ complex_it<unsigned, 2>{71u, 73u, 79u, 83u} };

    BOOST_CHECK_EQUAL( m[0], 71L );
    BOOST_CHECK_EQUAL( m[1], 73L );

    // Upgrade post-real ranks
    complex_rt<int, 2> const  n{ complex_it<int, 1>{87, 89} };

    BOOST_CHECK_EQUAL( n[0], 87 );
    BOOST_CHECK_EQUAL( n[1], 89 );
    BOOST_CHECK_EQUAL( n[2], 0 );
    BOOST_CHECK_EQUAL( n[3], 0 );

    // Upgrade post-real ranks, different component types
    complex_rt<unsigned, 3> const  p{ complex_it<long, 1>{93L, 97L} };

    BOOST_CHECK_EQUAL( p[0], 93u );
    BOOST_CHECK_EQUAL( p[1], 97u );
    BOOST_CHECK_EQUAL( p[2], 0u );
    BOOST_CHECK_EQUAL( p[3], 0u );
    BOOST_CHECK_EQUAL( p[4], 0u );
    BOOST_CHECK_EQUAL( p[5], 0u );
    BOOST_CHECK_EQUAL( p[6], 0u );
    BOOST_CHECK_EQUAL( p[7], 0u );
}

// Check conversion from a list of scalars.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_multireal_construction, T, test_types )
{
    // (Regular) complex
    complex_rt<T, 1> const  a = { (T)2, (T)3 };

    BOOST_CHECK_EQUAL( a[0], T(2) );
    BOOST_CHECK_EQUAL( a[1], T(3) );

    // Quaternions
    complex_rt<T, 2> const  b = { (T)5, (T)7 }, c = { (T)11, (T)13, (T)17,
     (T)19 };

    BOOST_CHECK_EQUAL( b[0], T(5) );
    BOOST_CHECK_EQUAL( b[1], T(7) );
    BOOST_CHECK_EQUAL( b[2], T{} );
    BOOST_CHECK_EQUAL( b[3], T{} );
    BOOST_CHECK_EQUAL( c[0], T(11) );
    BOOST_CHECK_EQUAL( c[1], T(13) );
    BOOST_CHECK_EQUAL( c[2], T(17) );
    BOOST_CHECK_EQUAL( c[3], T(19) );
}

// Check conversions keeping the same length, but different component types.
BOOST_AUTO_TEST_CASE( test_same_size_diff_type_conversion )
{
    // Note that same-type/same-size is covered by the automatically-defined
    // copy constructor, which (usually) takes priority over constructor
    // templates (such as the ones used here).

    // Reals
    complex_rt<unsigned, 0> const  a = { complex_rt<unsigned char, 0>{'\0'} };

    BOOST_CHECK_EQUAL( a[0], 0u );

    // (Regular) Complexes
    complex_rt<long, 1> const  b = { complex_rt<int, 1>{-2, +3} };

    BOOST_CHECK_EQUAL( b[0], -2L );
    BOOST_CHECK_EQUAL( b[1], +3L );

    // Quaternions
    complex_rt<double,2> const  c = {complex_rt<float,2>{+5.5f, -7.0f, +11.0f}};

    BOOST_CHECK_CLOSE( c[0], +5.5, 0.1 );
    BOOST_CHECK_CLOSE( c[1], -7.0, 0.1 );
    BOOST_CHECK_CLOSE( c[2], 11.0, 0.1 );
    BOOST_CHECK_CLOSE( c[3],  0.0, 0.1 );

    // Since brace-initialization is used in the implementation, narrowing
    // conversions can get flagged as warnings.
}

// Check conversions with immediately-lower rank, any component type.
BOOST_AUTO_TEST_CASE( test_barrage_conversion )
{
    // Same type between barrages and composite
    complex_rt<int, 0> const     a1 = { 2 }, a2 = { -3 };
    complex_rt<int, 1> const     a = { a1, a2 };
    complex_rt<double, 1> const  b1 = { -5.5 }, b2 = { +7.1, -11.3 };
    complex_rt<double, 2> const  b = { b1, b2 };

    BOOST_CHECK_EQUAL( a[0], a1[0] );
    BOOST_CHECK_EQUAL( a[1], a2[0] );

    BOOST_CHECK_CLOSE( b[0], b1[0], 0.1 );
    BOOST_CHECK_CLOSE( b[1], b1[1], 0.1 );
    BOOST_CHECK_CLOSE( b[2], b2[0], 0.1 );
    BOOST_CHECK_CLOSE( b[3], b2[1], 0.1 );

    // Only one barrage
    complex_rt<int, 1> const     aa = { a2 };
    complex_rt<double, 2> const  bb = { b1 };

    BOOST_CHECK_EQUAL( aa[0], a2[0] );
    BOOST_CHECK_EQUAL( aa[1], 0 );

    BOOST_CHECK_CLOSE( bb[0], b1[0], 0.1 );
    BOOST_CHECK_CLOSE( bb[1], b1[1], 0.1 );
    BOOST_CHECK_CLOSE( bb[2], 0.0, 0.1 );
    BOOST_CHECK_CLOSE( bb[3], 0.0, 0.1 );

    // Mixed types
    complex_rt<long, 1> const         c = { complex_rt<int, 0>{-13},
     complex_rt<long, 0>{17L} };
    complex_rt<long double, 2> const  d = { complex_rt<float, 1>{-19.4f},
     complex_rt<double, 1>{23.0, -29.8} };

    BOOST_CHECK_EQUAL( c[0], -13L );
    BOOST_CHECK_EQUAL( c[1], 17L );

    BOOST_CHECK_CLOSE( d[0], -19.4L, 0.1 );
    BOOST_CHECK_CLOSE( d[1], 0.0L, 0.1 );
    BOOST_CHECK_CLOSE( d[2], 23.0L, 0.1 );
    BOOST_CHECK_CLOSE( d[3], -29.8L, 0.1 );

    // One mixed barrage
    complex_rt<short, 1> const  e = { complex_rt<char, 0>{125} };

    BOOST_CHECK_EQUAL( e[0], 125 );
    BOOST_CHECK_EQUAL( e[1], 0 );

    // Since brace-initialization is used in the implementation, narrowing
    // conversions can get flagged as warnings.
}

// Check conversions with severely-lower rank, any component type.
BOOST_AUTO_TEST_CASE( test_subbarrage_conversion )
{
    // Same type between pieces and whole
    complex_rt<int, 0> const  qc[] = { {2}, {-3}, {5}, {-7} };
    complex_rt<int, 2> const  q = { qc[0], qc[1], qc[2], qc[3] };
    complex_rt<int, 0> const  oc[] = { {11}, {-13}, {17}, {-19} };
    complex_rt<int, 3> const  o = { qc[0], qc[1], qc[2], qc[3], oc[0], oc[1],
     oc[2], oc[3] };

    BOOST_CHECK_EQUAL( q, (decltype(q){2, -3, 5, -7}) );
    BOOST_CHECK_EQUAL( o, (decltype(o){2, -3, 5, -7, 11, -13, 17, -19}) );

    // Differing types
    complex_rt<long, 1> const       xc1[] = { {23L, -29L}, {31L, -37L} };
    complex_rt<int, 1> const        xc2[] = { {-41, 43}, {47, -53} };
    complex_rt<long long, 3> const  x = { xc1[0], xc2[1], xc2[0], xc1[1] };
    complex_rt<long, 3> const       y = { xc2[1] };

    BOOST_CHECK_EQUAL( x, (decltype(x){23LL, -29LL, 47LL, -53LL, -41LL, 43LL,
     31LL, -37LL}) );
    BOOST_CHECK_EQUAL( y, (decltype(y){47L, -53L, 0L}) );
}

// Check conversions with a greater rank, any component type.
BOOST_AUTO_TEST_CASE( test_supersize_conversion )
{
    // Integer
    complex_rt<int, 3> const        o = { -2, 3, -5, 7, -11, 13, -17, 19 };
    complex_rt<int, 2> const        q1{ o };
    complex_rt<long, 2> const       q2{ o };
    complex_rt<long, 1> const       c1{ o };
    complex_rt<int, 1> const        c2{ q1 };
    complex_rt<long long, 0> const  r1{ o };
    complex_rt<int, 0> const        r2{ o };

    BOOST_CHECK_EQUAL( q1[0], o[0] );
    BOOST_CHECK_EQUAL( q1[1], o[1] );
    BOOST_CHECK_EQUAL( q1[2], o[2] );
    BOOST_CHECK_EQUAL( q1[3], o[3] );
    BOOST_CHECK_EQUAL( q2, (decltype(q2){-2L, 3L, -5L, 7L}) );
    BOOST_CHECK_EQUAL( c1, (decltype(c1){-2L, 3L}) );
    BOOST_CHECK_EQUAL( c2, (decltype(c2){-2, 3}) );
    BOOST_CHECK_EQUAL( r1[0], -2LL );
    BOOST_CHECK_EQUAL( r2[0], -2 );

    // Floating
    complex_rt<float, 2> const        q = { -23.3f, +29.9f, -31.1f };
    complex_rt<double, 1> const       c3{ q };
    complex_rt<float, 1> const        c4{ q };
    complex_rt<float, 0> const        r3{ q };
    complex_rt<long double, 0> const  r4{ q };
    complex_rt<double, 0> const       r5{ c3 };

    BOOST_CHECK_CLOSE( c3[0], -23.3, 0.1 );
    BOOST_CHECK_CLOSE( c3[1], +29.9, 0.1 );
    BOOST_CHECK_CLOSE( c4[0], -23.3f, 0.1 );
    BOOST_CHECK_CLOSE( c4[1], +29.9f, 0.1 );
    BOOST_CHECK_CLOSE( r3[0], -23.3f, 0.1 );
    BOOST_CHECK_CLOSE( r4[0], -23.3L, 0.1 );
    BOOST_CHECK_CLOSE( r5[0], -23.3, 0.1 );
}

// Check explicit-conversions with component type needing explicit conversion.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_explicit_conversions, T, test_types )
{
    // Sample type needing explicit conversions
    struct proxy
    {
        T  value;

        proxy() = default;
        proxy( T const &t )  : value{ t }  {}

        explicit  operator T() const  { return value; }
    };

    // All reals, but one longer
    complex_rt<T, 0> const      a = { (T)2 };
    complex_rt<proxy, 0> const  b = { a };
    complex_rt<T, 0> const      c{ b };
    complex_rt<T, 3> const      o{ b };

    BOOST_CHECK_EQUAL( a, c );
    BOOST_CHECK_EQUAL( a, o );

    // (Regular) complex, different-length destinations
    complex_rt<proxy, 1> const  d = { proxy((T)3), proxy((T)5) };
    complex_rt<T, 0> const      e{ d };
    complex_rt<T, 1> const      f{ d };
    complex_rt<T, 2> const      g{ d };

    BOOST_CHECK_EQUAL( e[0], (T)3 );
    BOOST_CHECK_EQUAL( f[0], (T)3 );
    BOOST_CHECK_EQUAL( f[1], (T)5 );
    BOOST_CHECK_EQUAL( g[0], (T)3 );
    BOOST_CHECK_EQUAL( g[1], (T)5 );
    BOOST_CHECK_EQUAL( g[2], T{} );
    BOOST_CHECK_EQUAL( g[3], T{} );
}

BOOST_AUTO_TEST_SUITE_END()  // constructor_tests

BOOST_AUTO_TEST_SUITE( operation_tests )

// Check conversion to complex_it objects.
BOOST_AUTO_TEST_CASE( test_cross_philosophy_conversion )
{
    using boost::math::complex_it;

    // Real-to-real
    complex_it<int, 0> const  a{ complex_rt<int, 0>{5} };

    BOOST_CHECK_EQUAL( a[0], 5 );

    // Real-to-real, different component types
    complex_it<long, 0> const  b{ complex_rt<int, 0>{2} };

    BOOST_CHECK_EQUAL( b[0], 2L );

    // Downgrade to real
    complex_it<int, 0> const  c{ complex_rt<int, 1>{3, -7} };

    BOOST_CHECK_EQUAL( c[0], 3 );

    // Downgrade to real, different component types
    complex_it<int, 0> const  d{ complex_rt<long, 1>{-11, 13} };

    BOOST_CHECK_EQUAL( d[0], -11 );

    // Upgrade from real
    complex_it<int, 1> const  e{ complex_rt<int, 0>{17} };

    BOOST_CHECK_EQUAL( e[0], 17 );
    BOOST_CHECK_EQUAL( e[1], 0 );

    // Upgrade from real, different component types
    complex_it<long, 2> const  f{ complex_rt<int, 0>{-19} };

    BOOST_CHECK_EQUAL( f[0], -19L );
    BOOST_CHECK_EQUAL( f[1], 0L );
    BOOST_CHECK_EQUAL( f[2], 0L );
    BOOST_CHECK_EQUAL( f[3], 0L );

    // Same post-real rank
    complex_it<int, 2> const  g{ complex_rt<int, 2>{23, 29, 31, 37} };

    BOOST_CHECK_EQUAL( g[0], 23 );
    BOOST_CHECK_EQUAL( g[1], 29 );
    BOOST_CHECK_EQUAL( g[2], 31 );
    BOOST_CHECK_EQUAL( g[3], 37 );

    // Same post-real rank, different component types
    complex_it<unsigned long, 2> const  h{ complex_rt<int, 2>{41, 43, 47, 53} };

    BOOST_CHECK_EQUAL( h[0], 41UL );
    BOOST_CHECK_EQUAL( h[1], 43UL );
    BOOST_CHECK_EQUAL( h[2], 47UL );
    BOOST_CHECK_EQUAL( h[3], 53UL );

    // Downgrade post-real ranks
    complex_it<int, 1> const  k{ complex_rt<int, 2>{57, 59, 61, 67} };

    BOOST_CHECK_EQUAL( k[0], 57 );
    BOOST_CHECK_EQUAL( k[1], 59 );

    // Downgrade post-real ranks, different component types
    complex_it<long, 1> const  m{ complex_rt<unsigned, 2>{71u, 73u, 79u, 83u} };

    BOOST_CHECK_EQUAL( m[0], 71L );
    BOOST_CHECK_EQUAL( m[1], 73L );

    // Upgrade post-real ranks
    complex_it<int, 2> const  n{ complex_rt<int, 1>{87, 89} };

    BOOST_CHECK_EQUAL( n[0], 87 );
    BOOST_CHECK_EQUAL( n[1], 89 );
    BOOST_CHECK_EQUAL( n[2], 0 );
    BOOST_CHECK_EQUAL( n[3], 0 );

    // Upgrade post-real ranks, different component types
    complex_it<unsigned, 3> const  p{ complex_rt<long, 1>{93L, 97L} };

    BOOST_CHECK_EQUAL( p[0], 93u );
    BOOST_CHECK_EQUAL( p[1], 97u );
    BOOST_CHECK_EQUAL( p[2], 0u );
    BOOST_CHECK_EQUAL( p[3], 0u );
    BOOST_CHECK_EQUAL( p[4], 0u );
    BOOST_CHECK_EQUAL( p[5], 0u );
    BOOST_CHECK_EQUAL( p[6], 0u );
    BOOST_CHECK_EQUAL( p[7], 0u );
}

// Check the swapping of states.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_swap, T, test_types )
{
    // Reals
    complex_rt<T, 0>  a = { (T)2 }, b = { (T)3 };

    BOOST_CHECK_EQUAL( a[0], T(2) );
    BOOST_CHECK_EQUAL( b[0], T(3) );
    swap( a, b );
    BOOST_CHECK_EQUAL( a[0], T(3) );
    BOOST_CHECK_EQUAL( b[0], T(2) );

    // Quaternions
    complex_rt<T, 2>  c = { (T)5, (T)7, (T)11, (T)13 };
    complex_rt<T, 2>  d = { (T)17, (T)19, (T)23 };

    swap( c, d );
    BOOST_CHECK_EQUAL( c[0], T(17) );
    BOOST_CHECK_EQUAL( c[1], T(19) );
    BOOST_CHECK_EQUAL( c[2], T(23) );
    BOOST_CHECK_EQUAL( c[3], T{} );
    BOOST_CHECK_EQUAL( d[0], T(5) );
    BOOST_CHECK_EQUAL( d[1], T(7) );
    BOOST_CHECK_EQUAL( d[2], T(11) );
    BOOST_CHECK_EQUAL( d[3], T(13) );
}

// Check conjugation.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_conj, T, test_builtin_types )
{
    // Reals
    complex_rt<T, 0> const  a = {},         b = { (T)2 };
    auto const             aa = conj( a ), bb = conj( b );

    BOOST_CHECK_EQUAL( aa[0], +T{} );
    BOOST_CHECK_EQUAL( bb[0], +T(2) );

    // (Regular) complex
    complex_rt<T, 1> const  c = { (T)3, (T)5 }, d = { (T)7 };
    auto const             cc = conj( c ),     dd = conj( d );

    BOOST_CHECK_EQUAL( cc[0], +T(3) );
    BOOST_CHECK_EQUAL( cc[1], -T(5) );
    BOOST_CHECK_EQUAL( dd[0], +T(7) );
    BOOST_CHECK_EQUAL( dd[1], -T{} );

    // Quaternions
    complex_rt<T, 2> const  e = { (T)11, (T)13, -(T)17, (T)19 };
    auto const             ee = conj( e );

    BOOST_CHECK_EQUAL( ee[0], +T(11) );
    BOOST_CHECK_EQUAL( ee[1], -T(13) );
    BOOST_CHECK_EQUAL( ee[2], +T(17) );  // double negative!
    BOOST_CHECK_EQUAL( ee[3], -T(19) );
}

// Check the real- and imaginary-component member functions.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_member_real_imag, T, test_types )
{
    // Reals
    complex_rt<T, 0>  a = {}, b = { (T)2 };

    BOOST_CHECK_EQUAL( a.real(), T{} );
    BOOST_CHECK_EQUAL( a.imag(), T{} );
    BOOST_CHECK_EQUAL( b.real(), T(2) );
    BOOST_CHECK_EQUAL( b.imag(), T{} );
    a.real( (T)3 );
    b.real( (T)5 );
    BOOST_CHECK_EQUAL( a.real(), T(3) );
    BOOST_CHECK_EQUAL( a.imag(), T{} );
    BOOST_CHECK_EQUAL( b.real(), T(5) );
    BOOST_CHECK_EQUAL( b.imag(), T{} );

    // (Regular) complexes
    complex_rt<T, 1>  c = { (T)7, (T)11 }, d = { (T)13 };

    BOOST_CHECK_EQUAL( c.real(), T(7) );
    BOOST_CHECK_EQUAL( c.imag(), T(11) );
    BOOST_CHECK_EQUAL( d.real(), T(13) );
    BOOST_CHECK_EQUAL( d.imag(), T{} );
    c.real( (T)17 );
    c.imag( (T)19 );
    d.real( (T)23 );
    d.imag( (T)29 );
    BOOST_CHECK_EQUAL( c.real(), T(17) );
    BOOST_CHECK_EQUAL( c.imag(), T(19) );
    BOOST_CHECK_EQUAL( d.real(), T(23) );
    BOOST_CHECK_EQUAL( d.imag(), T(29) );

    // Quaternions
    complex_rt<T, 2>  e = { (T)31, (T)37, (T)41, (T)43 }, f = { (T)47 };

    BOOST_CHECK_EQUAL( e.real(), T(31) );
    BOOST_CHECK_EQUAL( e.imag(), T(37) );
    BOOST_CHECK_EQUAL( f.real(), T(47) );
    BOOST_CHECK_EQUAL( f.imag(), T{} );
    e.real( (T)53 );
    e.imag( (T)59 );
    f.real( (T)61 );
    f.imag( (T)67 );
    BOOST_CHECK_EQUAL( e.real(), T(53) );
    BOOST_CHECK_EQUAL( e.imag(), T(59) );
    BOOST_CHECK_EQUAL( f.real(), T(61) );
    BOOST_CHECK_EQUAL( f.imag(), T(67) );

    BOOST_CHECK_EQUAL( e[2], T(41) );
    BOOST_CHECK_EQUAL( e[3], T(43) );
    BOOST_CHECK_EQUAL( f[2], T{} );
    BOOST_CHECK_EQUAL( f[3], T{} );
}

// Check the unreal member functions.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_member_unreal, T, test_types )
{
    // Reals
    complex_rt<T, 0>  a = {}, b = { (T)2 };

    BOOST_CHECK_EQUAL( a.unreal(), decltype(a){} );
    BOOST_CHECK_EQUAL( b.unreal(), decltype(b){} );
    a.unreal( decltype(a){(T)3} );
    b.unreal( decltype(b){(T)5} );
    BOOST_CHECK_EQUAL( a.unreal(), decltype(a){} );
    BOOST_CHECK_EQUAL( b.unreal(), decltype(b){} );
    BOOST_CHECK_EQUAL( a.real(), T{} );
    BOOST_CHECK_EQUAL( b.real(), T(2) );

    // (Regular) complexes
    complex_rt<T, 1>  c = { (T)7, (T)11 }, d = { (T)13, (T)17 };

    BOOST_CHECK_EQUAL( c.unreal(), (decltype( c ){ T{}, T(11) }) );
    BOOST_CHECK_EQUAL( d.unreal(), (decltype( d ){ T{}, T(17) }) );
    c.unreal( decltype(c){(T)19, (T)23} );
    d.unreal( decltype(d){(T)29, (T)31} );
    BOOST_CHECK_EQUAL( c.unreal(), (decltype( c ){ T{}, T(23) }) );
    BOOST_CHECK_EQUAL( d.unreal(), (decltype( d ){ T{}, T(31) }) );
    BOOST_CHECK_EQUAL( c.real(), T(7) );
    BOOST_CHECK_EQUAL( d.real(), T(13) );

    // Quaternions
    complex_rt<T, 2>  e = { (T)37, (T)41, T(43), T(47) }, f = { (T)53, (T)59 };

    BOOST_CHECK_EQUAL( e.unreal(),(decltype( e ){ T{}, T(41), T(43), T(47) }) );
    BOOST_CHECK_EQUAL( f.unreal(),(decltype( f ){ T{}, T(59) }) );
    e.unreal( decltype(e){(T)61, (T)67, (T)71, (T)73} );
    f.unreal( decltype(f){(T)79, (T)83, (T)87, (T)89} );
    BOOST_CHECK_EQUAL( e.unreal(),(decltype( e ){ T{}, T(67), T(71), T(73) }) );
    BOOST_CHECK_EQUAL( f.unreal(),(decltype( f ){ T{}, T(83), T(87), T(89) }) );
    BOOST_CHECK_EQUAL( e.real(), T(37) );
    BOOST_CHECK_EQUAL( f.real(), T(53) );
}

// Check (Cayley) norm, with integer types.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_norm1, T, test_integer_types )
{
    // Reals
    complex_rt<T, 0> const  a = {}, b = { (T)2 };

    BOOST_CHECK_EQUAL( norm(a), T{} );
    BOOST_CHECK_EQUAL( norm(b), T(4) );

    // (Regular) complex
    complex_rt<T, 1> const  c = {}, d = { (T)3 }, e = { (T)5, (T)7 };

    BOOST_CHECK_EQUAL( norm(c), T{} );
    BOOST_CHECK_EQUAL( norm(d), T(9) );
    BOOST_CHECK_EQUAL( norm(e), T(74) );

    // Quaternions
    complex_rt<T, 2> const  f = {}, g = { (T)11, (T)13, (T)17, (T)19 };

    BOOST_CHECK_EQUAL( norm(f), T{} );
    BOOST_CHECK_EQUAL( norm(g), T(940) );
}

// Check (Cayley) norm, with floating-point types.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_norm2, T, test_floating_types )
{
    // Reals
    complex_rt<T, 0> const  a = {}, b = { (T)2 };

    BOOST_CHECK_CLOSE( norm(a), T{}, 0.1 );
    BOOST_CHECK_CLOSE( norm(b), T(4), 0.1 );

    // (Regular) complex
    complex_rt<T, 1> const  c = {}, d = { -(T)3 }, e = { (T)5, (T)-7 };

    BOOST_CHECK_CLOSE( norm(c), T{}, 0.1 );
    BOOST_CHECK_CLOSE( norm(d), T(9), 0.1 );
    BOOST_CHECK_CLOSE( norm(e), T(74), 0.1 );

    // Quaternions
    complex_rt<T, 2> const  f = {}, g = { (T)11, (T)13, -(T)17, (T)19 };

    BOOST_CHECK_CLOSE( norm(f), T{}, 0.1 );
    BOOST_CHECK_CLOSE( norm(g), T(940), 0.1 );
}

BOOST_AUTO_TEST_SUITE_END()  // operation_tests

BOOST_AUTO_TEST_SUITE( tuple_tests )

// Check the various compile-time tuple descriptors.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_tuple_support, T, test_types )
{
    using std::tuple_size;
    using std::is_same;
    using std::tuple_element;

    typedef complex_rt<T, 0>  real_type;
    typedef complex_rt<T, 1>  complex_type;
    typedef complex_rt<T, 2>  quaternion_type;
    typedef complex_rt<T, 3>  octonion_type;

    // Tuple size
    BOOST_REQUIRE_EQUAL( tuple_size<real_type>::value, 1u );
    BOOST_REQUIRE_EQUAL( tuple_size<complex_type>::value, 2u );
    BOOST_REQUIRE_EQUAL( tuple_size<quaternion_type>::value, 4u );
    BOOST_REQUIRE_EQUAL( tuple_size<octonion_type>::value, 8u );

    // Tuple element types
    BOOST_REQUIRE( (is_same<typename tuple_element<0, real_type>::type,
     T>::value) );

    BOOST_REQUIRE( (is_same<typename tuple_element<0, complex_type>::type,
     T>::value) );
    BOOST_REQUIRE( (is_same<typename tuple_element<1, complex_type>::type,
     T>::value) );

    BOOST_REQUIRE( (is_same<typename tuple_element<0, quaternion_type>::type,
     T>::value) );
    BOOST_REQUIRE( (is_same<typename tuple_element<1, quaternion_type>::type,
     T>::value) );
    BOOST_REQUIRE( (is_same<typename tuple_element<2, quaternion_type>::type,
     T>::value) );
    BOOST_REQUIRE( (is_same<typename tuple_element<3, quaternion_type>::type,
     T>::value) );

    BOOST_REQUIRE( (is_same<typename tuple_element<0, octonion_type>::type,
     T>::value) );
    BOOST_REQUIRE( (is_same<typename tuple_element<1, octonion_type>::type,
     T>::value) );
    BOOST_REQUIRE( (is_same<typename tuple_element<2, octonion_type>::type,
     T>::value) );
    BOOST_REQUIRE( (is_same<typename tuple_element<3, octonion_type>::type,
     T>::value) );
    BOOST_REQUIRE( (is_same<typename tuple_element<4, octonion_type>::type,
     T>::value) );
    BOOST_REQUIRE( (is_same<typename tuple_element<5, octonion_type>::type,
     T>::value) );
    BOOST_REQUIRE( (is_same<typename tuple_element<6, octonion_type>::type,
     T>::value) );
    BOOST_REQUIRE( (is_same<typename tuple_element<7, octonion_type>::type,
     T>::value) );

#if 0
    // These would give errors at compile-time, since the indices are too large!
    BOOST_REQUIRE( (is_same<typename tuple_element<10, real_type>::type,
     T>::value) );
    BOOST_REQUIRE( (is_same<typename tuple_element<10, complex_type>::type,
     T>::value) );
    BOOST_REQUIRE( (is_same<typename tuple_element<10, quaternion_type>::type,
     T>::value) );
    BOOST_REQUIRE( (is_same<typename tuple_element<10, octonion_type>::type,
     T>::value) );
#endif
}

// Check the tuple access functions.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_tuple_get, T, test_types )
{
    // Need to bring this into scope.
    using boost::math::get;

    // Reals
    complex_rt<T, 0>           a = { (T)2 };
    complex_rt<T, 0> const &  aa = a;

    BOOST_CHECK_EQUAL( get<0>(a), T(2) );
    get<0>( a ) = get<0>( complex_rt<T, 0>{(T)3} );
    BOOST_CHECK_EQUAL( T(3), get<0>(aa) );

    // Quaternions
    complex_rt<T, 2>           b = { (T)5, (T)7, (T)11, (T)13 };
    complex_rt<T, 2> const &  bb = b;

    BOOST_CHECK_EQUAL( get<0>(b), T(5) );
    BOOST_CHECK_EQUAL( get<1>(b), T(7) );
    BOOST_CHECK_EQUAL( get<2>(b), T(11) );
    BOOST_CHECK_EQUAL( get<3>(b), T(13) );
    get<0>( b ) = get<0>( complex_rt<T, 1>{(T)17, (T)8} );
    get<1>( b ) = get<1>( complex_rt<T, 1>{(T)9, (T)19} );
    get<2>( b ) = get<2>( complex_rt<T, 2>{(T)1, (T)6, (T)23} );
    get<3>( b ) = get<3>( complex_rt<T, 3>{a} );
    BOOST_CHECK_EQUAL( T(17), get<0>(bb) );
    BOOST_CHECK_EQUAL( T(19), get<1>(bb) );
    BOOST_CHECK_EQUAL( T(23), get<2>(bb) );
    BOOST_CHECK_EQUAL( T{}, get<3>(bb) );
}

BOOST_AUTO_TEST_SUITE_END()  // tuple_tests

BOOST_AUTO_TEST_SUITE( operator_tests )

// Check the identity operator.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_identity, T, test_builtin_types )
{
    // Reals
    complex_rt<T, 0> const   a = {},  b = { (T)2 };
    auto const              aa = +a, bb = +b;

    BOOST_CHECK_EQUAL( aa[0], +T{} );
    BOOST_CHECK_EQUAL( bb[0], +T(2) );

    // Quaternions
    complex_rt<T, 2> const  c = { (T)3, (T)5, (T)7 };
    auto const             cc = +c;

    BOOST_CHECK_EQUAL( cc[0], +T(3) );
    BOOST_CHECK_EQUAL( cc[1], +T(5) );
    BOOST_CHECK_EQUAL( cc[2], +T(7) );
    BOOST_CHECK_EQUAL( cc[3], +T{} );
}

// Check the addition operators.
BOOST_AUTO_TEST_CASE( test_addition )
{
    // Type-aliases
    typedef complex_rt<int, 0>        real_type;
    typedef complex_rt<int, 1>     complex_type;
    typedef complex_rt<int, 2>  quaternion_type;

    // Reals
    BOOST_CHECK_EQUAL( real_type(2) + real_type(3), real_type(5) );
    BOOST_CHECK_EQUAL( real_type(7) + 11, real_type(18) );
    BOOST_CHECK_EQUAL( 11 + real_type(13), real_type(24) );

    // (Regular) complexes
    BOOST_CHECK_EQUAL( complex_type(17, 19) + complex_type(23, 29),
     complex_type(40, 48) );
    BOOST_CHECK_EQUAL( complex_type(31, 37) + 41, complex_type(72, 37) );
    BOOST_CHECK_EQUAL( 43 + complex_type(47, 53), complex_type(90, 53) );
    BOOST_CHECK_EQUAL( real_type(59) + complex_type(61, 67), complex_type(120,
     67) );
    BOOST_CHECK_EQUAL( complex_type(71, 73) + real_type(79), complex_type(150,
     73) );

    // Quaternions
    BOOST_CHECK_EQUAL( quaternion_type(83, 89, 97, 101) + quaternion_type(103,
     107, 109, 113), quaternion_type(186, 196, 206, 214) );
    BOOST_CHECK_EQUAL( quaternion_type(127, 131, 137, 139) + 149,
     quaternion_type(276, 131, 137, 139) );
    BOOST_CHECK_EQUAL( 151 + quaternion_type(157, 163, 167, 173),
     quaternion_type(308, 163, 167, 173) );
    BOOST_CHECK_EQUAL( real_type(179) + quaternion_type(181, 191, 193, 197),
     quaternion_type(360, 191, 193, 197) );
    BOOST_CHECK_EQUAL( quaternion_type(199, 211, 223, 227) + real_type(229),
     quaternion_type(428, 211, 223, 227) );
    BOOST_CHECK_EQUAL( complex_type(233, 239) + quaternion_type(241, 251, 257,
     263), quaternion_type(474, 490, 257, 263) );
    BOOST_CHECK_EQUAL( quaternion_type(269, 271, 277, 281) + complex_type(283,
     293), quaternion_type(552, 564, 277, 281) );

    // Add-assignment
    real_type        a = { 1 };
    complex_type     b = { 2, 3 };
    quaternion_type  c = { 4, 5, 6, 7 };

    a += real_type{ 8 };
    BOOST_CHECK_EQUAL( a, real_type(9) );
    a += 10;
    BOOST_CHECK_EQUAL( a, real_type(19) );
    b += complex_type{ 11, 12 };
    BOOST_CHECK_EQUAL( b, complex_type(13, 15) );
    b += a;
    BOOST_CHECK_EQUAL( b, complex_type(32, 15) );
    b += 14;
    BOOST_CHECK_EQUAL( b, complex_type(46, 15) );
    c += quaternion_type{ 16, 17, 18, 20 };
    BOOST_CHECK_EQUAL( c, quaternion_type(20, 22, 24, 27) );
    c += b;
    BOOST_CHECK_EQUAL( c, quaternion_type(66, 37, 24, 27) );
    c += a;
    BOOST_CHECK_EQUAL( c, quaternion_type(85, 37, 24, 27) );
    c += 21;
    BOOST_CHECK_EQUAL( c, quaternion_type(106, 37, 24, 27) );

    // Successor
    BOOST_CHECK_EQUAL( ++a, real_type(20) );
    BOOST_CHECK_EQUAL( a++, real_type(20) );
    BOOST_CHECK_EQUAL( a, real_type(21) );

    BOOST_CHECK_EQUAL( ++b, complex_type(47, 15) );
    BOOST_CHECK_EQUAL( b++, complex_type(47, 15) );
    BOOST_CHECK_EQUAL( b, complex_type(48, 15) );

    BOOST_CHECK_EQUAL( ++c, quaternion_type(107, 37, 24, 27) );
    BOOST_CHECK_EQUAL( c++, quaternion_type(107, 37, 24, 27) );
    BOOST_CHECK_EQUAL( c, quaternion_type(108, 37, 24, 27) );
}

// Check the negation operator.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_negation, T, test_builtin_types )
{
    // Reals
    complex_rt<T, 0> const   a = {},  b = { (T)2 };
    auto const              aa = -a, bb = -b;

    BOOST_CHECK_EQUAL( aa[0], -T{} );
    BOOST_CHECK_EQUAL( bb[0], -T(2) );

    // Quaternions
    complex_rt<T, 2> const  c = { (T)3, (T)5, (T)7 };
    auto const             cc = -c;

    BOOST_CHECK_EQUAL( cc[0], -T(3) );
    BOOST_CHECK_EQUAL( cc[1], -T(5) );
    BOOST_CHECK_EQUAL( cc[2], -T(7) );
    BOOST_CHECK_EQUAL( cc[3], -T{} );
}

// Check the subtraction operators.
BOOST_AUTO_TEST_CASE( test_subtraction )
{
    // Type-aliases
    typedef complex_rt<int, 0>        real_type;
    typedef complex_rt<int, 1>     complex_type;
    typedef complex_rt<int, 2>  quaternion_type;

    // Reals
    BOOST_CHECK_EQUAL( real_type(3) - real_type(2), real_type(1) );
    BOOST_CHECK_EQUAL( real_type(5) - 7, real_type(-2) );
    BOOST_CHECK_EQUAL( -11 - real_type(-13), real_type(2) );

    // (Regular) complexes
    BOOST_CHECK_EQUAL( complex_type(23, 19) - complex_type(17, 29),
     complex_type(6, -10) );
    BOOST_CHECK_EQUAL( complex_type(37, -31) - -41, complex_type(78, -31) );
    BOOST_CHECK_EQUAL( 43 - complex_type(47, 53), complex_type(-4, -53) );
    BOOST_CHECK_EQUAL( real_type(61) - complex_type(59, 67), complex_type(2,
     -67) );
    BOOST_CHECK_EQUAL( complex_type(71, 73) - real_type(79), complex_type(-8,
     73) );

    // Quaternions
    BOOST_CHECK_EQUAL( quaternion_type(103, 107, 109, 113) - quaternion_type(83,
     89, 97, 101), quaternion_type(20, 18, 12, 12) );
    BOOST_CHECK_EQUAL( quaternion_type(127, 137, 139, 149) - 131,
     quaternion_type(-4, 137, 139, 149) );
    BOOST_CHECK_EQUAL( 157 - quaternion_type(151, 163, 167, 173),
     quaternion_type(6, -163, -167, -173) );
    BOOST_CHECK_EQUAL( real_type(-179) - quaternion_type(-181, 191, -193, 197),
     quaternion_type(2, -191, 193, -197) );
    BOOST_CHECK_EQUAL( quaternion_type(199, 211, 227, 229) - real_type(223),
     quaternion_type(-24, 211, 227, 229) );
    BOOST_CHECK_EQUAL( complex_type(241, 239) - quaternion_type(233, 251, -257,
     263), quaternion_type(8, -12, 257, -263) );
    BOOST_CHECK_EQUAL( quaternion_type(271, 281, 283, 293) - complex_type(269,
     277), quaternion_type(2, 4, 283, 293) );

    // Subtract-assignment
    real_type        a = { 1 };
    complex_type     b = { 2, 3 };
    quaternion_type  c = { 4, 5, 6, 7 };

    a -= real_type{ 8 };
    BOOST_CHECK_EQUAL( a, real_type(-7) );
    a -= -10;
    BOOST_CHECK_EQUAL( a, real_type(3) );
    b -= complex_type{ -1, 4 };
    BOOST_CHECK_EQUAL( b, complex_type(3, -1) );
    b -= a;
    BOOST_CHECK_EQUAL( b, complex_type(0, -1) );
    b -= 14;
    BOOST_CHECK_EQUAL( b, complex_type(-14, -1) );
    c -= quaternion_type{ 16, 17, 18, 20 };
    BOOST_CHECK_EQUAL( c, quaternion_type(-12, -12, -12, -13) );
    c -= b;
    BOOST_CHECK_EQUAL( c, quaternion_type(2, -11, -12, -13) );
    c -= a;
    BOOST_CHECK_EQUAL( c, quaternion_type(-1, -11, -12, -13) );
    c -= -21;
    BOOST_CHECK_EQUAL( c, quaternion_type(20, -11, -12, -13) );

    // Predecessor
    BOOST_CHECK_EQUAL( --a, real_type(2) );
    BOOST_CHECK_EQUAL( a--, real_type(2) );
    BOOST_CHECK_EQUAL( a, real_type(1) );

    BOOST_CHECK_EQUAL( --b, complex_type(-15, -1) );
    BOOST_CHECK_EQUAL( b--, complex_type(-15, -1) );
    BOOST_CHECK_EQUAL( b, complex_type(-16, -1) );

    BOOST_CHECK_EQUAL( --c, quaternion_type(19, -11, -12, -13) );
    BOOST_CHECK_EQUAL( c--, quaternion_type(19, -11, -12, -13) );
    BOOST_CHECK_EQUAL( c, quaternion_type(18, -11, -12, -13) );
}

// Check the (new) conjugation operator.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_conjugation, T, test_builtin_types )
{
    // Reals
    complex_rt<T, 0> const   a = {},  b = { (T)2 };
    auto const              aa = ~a, bb = ~b;

    BOOST_CHECK_EQUAL( aa[0], +T{} );
    BOOST_CHECK_EQUAL( bb[0], +T(2) );

    // Quaternions
    complex_rt<T, 2> const  c = { (T)3, (T)5, (T)7 };
    auto const             cc = ~c;

    BOOST_CHECK_EQUAL( cc[0], +T(3) );
    BOOST_CHECK_EQUAL( cc[1], -T(5) );
    BOOST_CHECK_EQUAL( cc[2], -T(7) );
    BOOST_CHECK_EQUAL( cc[3], -T{} );
}

// Check the multiplication-with-scalar operators, integer edition.
BOOST_AUTO_TEST_CASE_TEMPLATE(test_scalar_multiplication1,T,test_integer_types)
{
    // Type-aliases
    typedef complex_rt<T, 0>        real_type;
    typedef complex_rt<T, 1>     complex_type;
    typedef complex_rt<T, 2>  quaternion_type;

    // Reals
    auto const  z = -T( 2 ) * real_type( T(1) );
    auto const  y = real_type( T(7) ) * T( 5 );

    BOOST_CHECK_EQUAL( z[0], -T(2) );
    BOOST_CHECK_EQUAL( y[0], T(35) );

    // (Regular) complexes
    auto const  a = T( 3 ) * complex_type{ -T(2), T(4) };
    auto const  b = complex_type{ T{}, -T(5) } * -T( 4 );

    BOOST_CHECK_EQUAL( a[0], -T(6) );
    BOOST_CHECK_EQUAL( a[1], T(12) );

    BOOST_CHECK_EQUAL( b[0], T{} );
    BOOST_CHECK_EQUAL( b[1], T(20) );

    // Quaternions
    auto const  c = T( 3 ) * quaternion_type{ T(4), T{}, -T(10), T(6) };
    auto const  d = quaternion_type{ -T(3), T(11), -T(5) } * -T( 1 );

    BOOST_CHECK_EQUAL( c[0], T(12) );
    BOOST_CHECK_EQUAL( c[1], T{} );
    BOOST_CHECK_EQUAL( c[2], -T(30) );
    BOOST_CHECK_EQUAL( c[3], T(18) );

    BOOST_CHECK_EQUAL( d[0], T(3) );
    BOOST_CHECK_EQUAL( d[1], -T(11) );
    BOOST_CHECK_EQUAL( d[2], T(5) );
    BOOST_CHECK_EQUAL( d[3], T{} );

    // Multiply-assignment
    real_type        e = { T(1) };
    complex_type     f = { T(2), T(3) };
    quaternion_type  g = { T(4), T(5), T(6), T(7) };

    e *= T( 10 );
    BOOST_CHECK_EQUAL( e[0], T(10) );
    f *= -T( 3 );
    BOOST_CHECK_EQUAL( f[0], -T(6) );
    BOOST_CHECK_EQUAL( f[1], -T(9) );
    g *= -T( 5 );
    BOOST_CHECK_EQUAL( g[0], -T(20) );
    BOOST_CHECK_EQUAL( g[1], -T(25) );
    BOOST_CHECK_EQUAL( g[2], -T(30) );
    BOOST_CHECK_EQUAL( g[3], -T(35) );
}

// Check the multiplication-with-scalar operators, floating-point edition.
BOOST_AUTO_TEST_CASE_TEMPLATE(test_scalar_multiplication2,T,test_floating_types)
{
    // Type-aliases
    typedef complex_rt<T, 0>        real_type;
    typedef complex_rt<T, 1>     complex_type;
    typedef complex_rt<T, 2>  quaternion_type;

    // Reals
    auto const  z = T( -2.5 ) * real_type( T(1.25) );
    auto const  y = real_type( T(7) ) * T( 5 );

    BOOST_CHECK_CLOSE( z[0], T(-3.125), 0.0001);
    BOOST_CHECK_CLOSE( y[0], T(35.0), 0.0001 );

    // (Regular) complexes
    auto const  a = T( 3 ) * complex_type{ T(-2), T(4) };
    auto const  b = complex_type{ T{}, T(-5) } * T( -4 );

    BOOST_CHECK_CLOSE( a[0], T(-6), 0.0001 );
    BOOST_CHECK_CLOSE( a[1], T(12), 0.0001 );

    BOOST_CHECK_CLOSE( b[0], T{}, 0.0001 );
    BOOST_CHECK_CLOSE( b[1], T(20), 0.0001 );

    // Quaternions
    auto const  c = T( 3.12 ) * quaternion_type{ T(4.4), T{}, T(-10), T(6) };
    auto const  d = quaternion_type{ T(-3), T(11), T(-5.1) } * T( -0.39 );

    BOOST_CHECK_CLOSE( c[0], T(13.728), 0.0001 );
    BOOST_CHECK_CLOSE( c[1], T{}, 0.0001 );
    BOOST_CHECK_CLOSE( c[2], T(-31.2), 0.0001 );
    BOOST_CHECK_CLOSE( c[3], T(18.72), 0.0001 );

    BOOST_CHECK_CLOSE( d[0], T(1.17), 0.0001 );
    BOOST_CHECK_CLOSE( d[1], T(-4.29), 0.0001 );
    BOOST_CHECK_CLOSE( d[2], T(1.989), 0.0001 );
    BOOST_CHECK_CLOSE( d[3], T{}, 0.0001 );

    // Multiply-assignment
    real_type        e = { T(1) };
    complex_type     f = { T(2), T(3) };
    quaternion_type  g = { T(4), T(5), T(6), T(7) };

    e *= T( 10.1 );
    BOOST_CHECK_CLOSE( e[0], T(10.1), 0.0001 );
    f *= T( -3 );
    BOOST_CHECK_CLOSE( f[0], T(-6), 0.0001 );
    BOOST_CHECK_CLOSE( f[1], T(-9), 0.0001 );
    g *= T( -2.5 );
    BOOST_CHECK_CLOSE( g[0], T(-10), 0.0001 );
    BOOST_CHECK_CLOSE( g[1], T(-12.5), 0.0001 );
    BOOST_CHECK_CLOSE( g[2], T(-15), 0.0001 );
    BOOST_CHECK_CLOSE( g[3], T(-17.5), 0.0001 );
}

// Check the division-with-scalar (including modulus) operators.
BOOST_AUTO_TEST_CASE( test_scalar_division_and_modulus )
{
    // Integer
    complex_rt<int, 0>  a = { 9 };
    complex_rt<int, 1>  b = { 8, 3 };
    complex_rt<int, 2>  c = { 0, 4, 8, 12 };

    BOOST_CHECK_EQUAL( real(a / 5), 1 );
    BOOST_CHECK_EQUAL( real(a % 5), 4 );

    BOOST_CHECK_EQUAL( (b / 4), decltype(b)(2, 0) );
    BOOST_CHECK_EQUAL( (b % 4), decltype(b)(0, 3) );

    BOOST_CHECK_EQUAL( (c / 3), decltype(c)(0, 1, 2, 4) );
    BOOST_CHECK_EQUAL( (c % 3), decltype(c)(0, 1, 2, 0) );
    BOOST_CHECK_EQUAL( (c / 3) * 3 + (c % 3), c );

    a /= 2;
    BOOST_CHECK_EQUAL( a, decltype(a)(4) );
    a = 9;
    a %= 2;
    BOOST_CHECK_EQUAL( a, decltype(a)(1) );
    b /= 1;
    BOOST_CHECK_EQUAL( b, decltype(b)(8, 3) );
    b %= 1;
    BOOST_CHECK_EQUAL( b, decltype(b){} );
    c *= 7;
    BOOST_REQUIRE_EQUAL( c, decltype(c)(0, 28, 56, 84) );
    c /= 5;
    BOOST_CHECK_EQUAL( c, decltype(c)(0, 5, 11, 16) );
    c = { 0, 28, 56, 84 };
    c %= 5;
    BOOST_CHECK_EQUAL( c, decltype(c)(0, 3, 1, 4) );

    // Floating
    complex_rt<double, 0>  d = { -16.5 };
    complex_rt<double, 1>  e = { 7.0, 0.0 };
    complex_rt<double, 2>  f = { -3.0, 2.1, 1.21, -100.7 };

    BOOST_CHECK_CLOSE( (d / -4.1)[0], 4.024, 0.1 );

    BOOST_CHECK_CLOSE( (e / 0.5)[0], 14.0, 0.1 );
    BOOST_CHECK_CLOSE( (e / 0.5)[1], 0.0, 0.1 );

    BOOST_CHECK_CLOSE( (f / -0.02)[0],  150.0, 0.1 );
    BOOST_CHECK_CLOSE( (f / -0.02)[1], -105.0, 0.1 );
    BOOST_CHECK_CLOSE( (f / -0.02)[2],  -60.5, 0.1 );
    BOOST_CHECK_CLOSE( (f / -0.02)[3], 5035.0, 0.1 );

    d /= 0.25;
    BOOST_CHECK_CLOSE( real(d), -66.0, 0.1 );
    e /= -3.0;
    BOOST_CHECK_CLOSE( real(e), -2.3333, 0.1 );
    BOOST_CHECK_CLOSE( imag(e),  0.0,    0.1 );
    f /= 12.1;
    BOOST_CHECK_CLOSE( f[0], -0.2479, 0.1 );
    BOOST_CHECK_CLOSE( f[1], +0.1735, 0.1 );
    BOOST_CHECK_CLOSE( f[2],  0.1,    0.1 );
    BOOST_CHECK_CLOSE( f[3], -8.3223, 0.1 );
}

BOOST_AUTO_TEST_SUITE_END()  // operator_tests

BOOST_AUTO_TEST_SUITE( function_tests )

// Check free-function real, imag, and unreal.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_real_imag_unreal, T, test_types )
{
    // Reals
    complex_rt<T, 0> const  a = {}, b = { (T)2 };

    BOOST_CHECK_EQUAL( real(a), T{} );
    BOOST_CHECK_EQUAL( imag(a), T{} );
    BOOST_CHECK_EQUAL( unreal(a), decltype(a){} );
    BOOST_CHECK_EQUAL( real(b), T(2) );
    BOOST_CHECK_EQUAL( imag(b), T{} );
    BOOST_CHECK_EQUAL( unreal(b), decltype(b){} );

    // (Regular) complexes
    complex_rt<T, 1> const  c = {}, d = { (T)3, (T)5 };

    BOOST_CHECK_EQUAL( real(c), T{} );
    BOOST_CHECK_EQUAL( imag(c), T{} );
    BOOST_CHECK_EQUAL( unreal(c), decltype(c){} );
    BOOST_CHECK_EQUAL( real(d), T(3) );
    BOOST_CHECK_EQUAL( imag(d), T(5) );
    BOOST_CHECK_EQUAL( unreal(d), (decltype(d){ T{}, T(5) }) );

    // Quaternions
    complex_rt<T, 2> const  e = {}, f = { (T)7, (T)11, (T)13 };

    BOOST_CHECK_EQUAL( real(e), T{} );
    BOOST_CHECK_EQUAL( imag(e), T{} );
    BOOST_CHECK_EQUAL( unreal(e), decltype(e){} );
    BOOST_CHECK_EQUAL( real(f), T(7) );
    BOOST_CHECK_EQUAL( imag(f), T(11) );
    BOOST_CHECK_EQUAL( unreal(f), (decltype(f){ T{}, T(11), T(13), T{} }) );
}

// Check the (non-Cayley) norm functions.
BOOST_AUTO_TEST_CASE( test_norm )
{
    // Integer
    complex_rt<int, 0>  a = { +4 }, b = { -3 };
    complex_rt<int, 1>  c = { -4, +3 }, d = { a, b }, e = { -a, +b },
                        f = { +a, -b };
    complex_rt<int, 2>  g = { 5, 0, -2, 1 };

    BOOST_CHECK_EQUAL( taxi(a), 4 );
    BOOST_CHECK_CLOSE( abs(a), 4.0, 0.1 );
    BOOST_CHECK_EQUAL( sup(a), 4 );
    BOOST_CHECK_EQUAL( taxi(b), 3 );
    BOOST_CHECK_CLOSE( abs(b), 3.0, 0.1 );
    BOOST_CHECK_EQUAL( sup(b), 3 );

    BOOST_CHECK_EQUAL( taxi(c), 7 );
    BOOST_CHECK_EQUAL( taxi(d), 7 );
    BOOST_CHECK_EQUAL( taxi(e), 7 );
    BOOST_CHECK_EQUAL( taxi(f), 7 );
    BOOST_CHECK_CLOSE( abs(c), 5.0, 0.1 );
    BOOST_CHECK_CLOSE( abs(d), 5.0, 0.1 );
    BOOST_CHECK_CLOSE( abs(e), 5.0, 0.1 );
    BOOST_CHECK_CLOSE( abs(f), 5.0, 0.1 );
    BOOST_CHECK_EQUAL( sup(c), 4 );
    BOOST_CHECK_EQUAL( sup(d), 4 );
    BOOST_CHECK_EQUAL( sup(e), 4 );
    BOOST_CHECK_EQUAL( sup(f), 4 );

    BOOST_CHECK_EQUAL( taxi(g), 8 );
    BOOST_CHECK_CLOSE( abs(g), 5.4772, 0.1 );
    BOOST_CHECK_EQUAL( sup(g), 5 );

    // Floating
    complex_rt<double, 0>  aa = { +4.0 }, bb = { -3.0 };
    complex_rt<double, 1>  cc = { -4.0, +3.0 }, dd = { aa, bb },
                           ee = { -aa, +bb }, ff = { +aa, -bb };
    complex_rt<double, 2>  gg = { 5.0, 0.0, -2.0, 1.0 },
                            h = { +3.0, -4.0, +12.0, -84.0 };
    complex_rt<double, 3>   k = { 6.7, -0.9, -11.2, 0.01, 4.33, -8.25, 255.5 };

    BOOST_CHECK_CLOSE( taxi(aa), 4.0, 0.1 );
    BOOST_CHECK_CLOSE( abs(aa), 4.0, 0.1 );
    BOOST_CHECK_CLOSE( sup(aa), 4.0, 0.1 );
    BOOST_CHECK_CLOSE( taxi(bb), 3.0, 0.1 );
    BOOST_CHECK_CLOSE( abs(bb), 3.0, 0.1 );
    BOOST_CHECK_CLOSE( sup(bb), 3.0, 0.1 );

    BOOST_CHECK_CLOSE( taxi(cc), 7.0, 0.1 );
    BOOST_CHECK_CLOSE( taxi(dd), 7.0, 0.1 );
    BOOST_CHECK_CLOSE( taxi(ee), 7.0, 0.1 );
    BOOST_CHECK_CLOSE( taxi(ff), 7.0, 0.1 );
    BOOST_CHECK_CLOSE( abs(cc), 5.0, 0.1 );
    BOOST_CHECK_CLOSE( abs(dd), 5.0, 0.1 );
    BOOST_CHECK_CLOSE( abs(ee), 5.0, 0.1 );
    BOOST_CHECK_CLOSE( abs(ff), 5.0, 0.1 );
    BOOST_CHECK_CLOSE( sup(cc), 4.0, 0.1 );
    BOOST_CHECK_CLOSE( sup(dd), 4.0, 0.1 );
    BOOST_CHECK_CLOSE( sup(ee), 4.0, 0.1 );
    BOOST_CHECK_CLOSE( sup(ff), 4.0, 0.1 );

    BOOST_CHECK_CLOSE( taxi(gg), 8.0, 0.1 );
    BOOST_CHECK_CLOSE( abs(gg), 5.4772, 0.1 );
    BOOST_CHECK_CLOSE( sup(gg), 5.0, 0.1 );

    BOOST_CHECK_CLOSE( taxi(h), 103.0, 0.1 );
    BOOST_CHECK_CLOSE( abs(h), 85.0, 0.1 );
    BOOST_CHECK_CLOSE( sup(h), 84.0, 0.1 );

    BOOST_CHECK_CLOSE( taxi(k), 286.89, 0.1 );
    BOOST_CHECK_CLOSE( abs(k), 256.00, 0.1 );
    BOOST_CHECK_CLOSE( sup(k), 255.5, 0.1 );
}

// Check the sign function.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_sgn, T, test_floating_types )
{
    // Reals
    complex_rt<T, 0> const  a = { T(+4) }, b = {}, c = { T(-3) },
                            d = { T(0.1) }, e = { T(-0.5) };

    BOOST_CHECK_CLOSE( sgn(a)[0], T(+1.0), 0.1 );
    BOOST_CHECK_CLOSE( sgn(b)[0], T( 0.0), 0.1 );
    BOOST_CHECK_CLOSE( sgn(c)[0], T(-1.0), 0.1 );
    BOOST_CHECK_CLOSE( sgn(d)[0], T(+1.0), 0.1 );
    BOOST_CHECK_CLOSE( sgn(e)[0], T(-1.0), 0.1 );

    // (Regular) complexes
    complex_rt<T, 1> const  f = {}, g = { T(3.0), T(-4.0) };
    auto const              f_sgn = sgn( f ), g_sgn = sgn( g );

    BOOST_CHECK_CLOSE( f_sgn[0], T{}, 0.1 );
    BOOST_CHECK_CLOSE( f_sgn[1], T{}, 0.1 );
    BOOST_CHECK_CLOSE( g_sgn[0], T(+0.6), 0.1 );
    BOOST_CHECK_CLOSE( g_sgn[1], T(-0.8), 0.1 );

    // Quaternions
    complex_rt<T, 2> const  h = { T{}, T{}, T(0.1), T{} }, k = {},
                            m = { T(+0.6), T(-0.8), T(+2.4), T(-16.8) };
    auto const          h_sgn = sgn( h ), k_sgn = sgn( k ), m_sgn = sgn( m );

    BOOST_CHECK_CLOSE( h_sgn[0], T( 0.0), 0.1 );
    BOOST_CHECK_CLOSE( h_sgn[1], T( 0.0), 0.1 );
    BOOST_CHECK_CLOSE( h_sgn[2], T(+1.0), 0.1 );
    BOOST_CHECK_CLOSE( h_sgn[3], T( 0.0), 0.1 );
    BOOST_CHECK_CLOSE( k_sgn[0], T{}, 0.1 );
    BOOST_CHECK_CLOSE( k_sgn[1], T{}, 0.1 );
    BOOST_CHECK_CLOSE( k_sgn[2], T{}, 0.1 );
    BOOST_CHECK_CLOSE( k_sgn[3], T{}, 0.1 );
    BOOST_CHECK_CLOSE( m_sgn[0], T(+0.03529), 0.1 );
    BOOST_CHECK_CLOSE( m_sgn[1], T(-0.04706), 0.1 );
    BOOST_CHECK_CLOSE( m_sgn[2], T(+0.14118), 0.1 );
    BOOST_CHECK_CLOSE( m_sgn[3], T(-0.98824), 0.1 );
}

BOOST_AUTO_TEST_SUITE_END()  // function_tests

BOOST_AUTO_TEST_SUITE_END()  // complex_rt_tests
