//  Boost Complex Numbers, iterative-mode unit test program file  ------------//

//  Copyright 2013 Daryle Walker.
//  Distributed under the Boost Software License, Version 1.0.  (See the
//  accompanying file LICENSE_1_0.txt or a copy at
//  <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/output_test_stream.hpp>
#include <boost/mpl/list.hpp>

#include "boost/math/complex_it.hpp"

#include <cstddef>
#include <ios>
#include <tuple>
#include <type_traits>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/rational.hpp>


// Common definitions  -------------------------------------------------------//

namespace {

    // Save time writing out long types
    namespace mp = boost::multiprecision;
    using boost::mpl::list;
    using boost::test_tools::output_test_stream;
    using boost::math::complex_it;

    // Sample testing types for components
    typedef mp::number<mp::cpp_dec_float<50>, mp::et_off>          my_float;
    typedef list<int, unsigned, double, mp::int512_t, my_float>  test_types;
    typedef list<int, unsigned, mp::int512_t>            test_integer_types;
    typedef list<double, my_float>                      test_floating_types;

}

// Flag un-printable types here.


BOOST_AUTO_TEST_SUITE( complex_it_tests )

BOOST_AUTO_TEST_SUITE( core_tests )

// Check the various compile-time attributes.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_complex_compile_time, T, test_types )
{
    using std::is_same;
    using std::size_t;

    typedef complex_it<T, 0>  real_type;
    typedef complex_it<T, 1>  complex_type;
    typedef complex_it<T, 2>  quaternion_type;
    typedef complex_it<T, 3>  octonion_type;

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
    typedef complex_it<T, 0>                  real_type;
    typedef complex_it<T, 1>               complex_type;
    typedef complex_it<T, 2>            quaternion_type;
    typedef complex_it<rational_type, 3>  octonion_type;

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
    typedef complex_it<T, 0>     real_type;
    typedef complex_it<T, 1>  complex_type;

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

// Check iterator access, for support of range-for.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_complex_iteration, T, test_types )
{
    // Real support
    complex_it<T, 0>     a;
    decltype(a) const &  aa = a;
    auto                 ab = begin( a );
    auto const           ae = end( a );
    auto                 aab = begin( aa );
    auto const           aae = end( aa );

    BOOST_CHECK_EQUAL( unsigned(ae - ab), decltype(a)::static_size );
    *ab++ = (T)6;
    BOOST_CHECK_EQUAL( aa[0], T(6) );
    BOOST_CHECK_EQUAL( ab, ae );
    a[ 0 ] = (T)7;
    BOOST_CHECK_EQUAL( *aab++, T(7) );
    BOOST_CHECK_EQUAL( aab, aae );

    // Check with multi-scalar
    complex_it<T, 2>     q;
    decltype(q) const &  qq = q;
    auto                 qb = begin( q );
    auto const           qe = end( q );
    auto                 qqb = begin( qq );
    auto const           qqe = end( qq );

    BOOST_CHECK_EQUAL( unsigned(qqe - qqb), decltype(q)::static_size );
    *qb++ = (T)10;
    BOOST_CHECK_EQUAL( qq[0], T(10) );
    *qb++ = (T)11;
    BOOST_CHECK_EQUAL( qq[1], T(11) );
    *qb++ = (T)12;
    BOOST_CHECK_EQUAL( qq[2], T(12) );
    *qb++ = (T)13;
    BOOST_CHECK_EQUAL( qq[3], T(13) );
    BOOST_CHECK_EQUAL( qb, qe );
    q[ 0 ] = (T)100;
    q[ 1 ] = (T)101;
    q[ 2 ] = (T)102;
    q[ 3 ] = (T)103;
    BOOST_CHECK_EQUAL( *qqb++, T(100) );
    BOOST_CHECK_EQUAL( *qqb++, T(101) );
    BOOST_CHECK_EQUAL( *qqb++, T(102) );
    BOOST_CHECK_EQUAL( *qqb++, T(103) );
    BOOST_CHECK_EQUAL( qqb, qqe );
}

// Check Boolean conversion.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_complex_to_boolean, T, test_types )
{
    // Real support
    complex_it<T, 0>  r;

    r[ 0 ] = T{};
    BOOST_CHECK( !r );
    r[ 0 ] = (T)2;
    BOOST_CHECK( (bool)r );

    // Check with multi-scalar
    complex_it<T, 2>  q;

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
    complex_it<T, 0>  r;
    auto const &      rr = r;

    r[ 0 ] = (T)6;
    BOOST_CHECK_EQUAL( rr[0], rr.lower_barrage()[0] );
    BOOST_CHECK_EQUAL( rr[0], rr.upper_barrage()[0] );

    // (Regular) complex
    complex_it<T, 1>  c;
    auto const &      cc = c;

    c[ 0 ] = (T)7;
    c[ 1 ] = (T)18;

    auto const  cc1 = cc.lower_barrage();
    auto const  cc2 = cc.upper_barrage();

    BOOST_CHECK_EQUAL( cc1[0], cc[0] );
    BOOST_CHECK_EQUAL( cc2[0], cc[1] );

    // Extreme
    complex_it<T, 3>  o;
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
}

// Check comparisons between hypercomplex and scalars.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_complex_real_equality, T, test_types )
{
    T                 sample = 6;
    complex_it<T, 0>  r;
    complex_it<T, 1>  c;
    complex_it<T, 2>  q;
    complex_it<T, 3>  o;

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
    complex_it<int, 0>  a, b;

    a[ 0 ] = b[ 0 ] = 0;
    BOOST_CHECK( a == b );
    BOOST_CHECK( !(a != b) );
    a[ 0 ] = -2;
    b[ 0 ] = +3;
    BOOST_CHECK( a != b );
    BOOST_CHECK( !(a == b) );

    // Two reals, differing component types
    complex_it<long, 0>  c;

    c[ 0 ] = -2L;
    BOOST_CHECK( a == c );
    BOOST_CHECK( !(a != c) );
    BOOST_CHECK( b != c );
    BOOST_CHECK( !(b == c) );

    // Two (regular) complex, same component type
    complex_it<int, 1>  d, e;

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
    complex_it<long, 1>  f;

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
    complex_it<int, 2>  g;

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
    complex_it<long, 3>  h;

    h[ 0 ] = h[ 1 ] = h[ 2 ] = h[ 3 ] = 0;
    h[ 4 ] = h[ 5 ] = h[ 6 ] = h[ 7 ] = 0;
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
    complex_it<int, 0>  r;
    complex_it<int, 1>  c;
    complex_it<int, 2>  q;
    complex_it<int, 3>  o;

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
    BOOST_CHECK( ots.is_equal("(-4,+5,-6,+7)") );

    o[ 0 ] = -10;
    o[ 1 ] = 11;
    o[ 2 ] = 12;
    o[ 3 ] = -13;
    o[ 4 ] = 14;
    o[ 5 ] = 15;
    o[ 6 ] = -16;
    o[ 7 ] = 101;
    ots << std::noshowpos << o;
    BOOST_CHECK( ots.is_equal("(-10,11,12,-13,14,15,-16,101)") );
}

BOOST_AUTO_TEST_SUITE_END()  // core_tests

BOOST_AUTO_TEST_SUITE( constructor_tests )

// Check the results of default- and scalar-construction.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_default_real_construction, T, test_types )
{
    // Real
    complex_it<T, 0>  a{}, b{ (T)2 };

    BOOST_CHECK_EQUAL( a[0], T{} );
    BOOST_CHECK_EQUAL( b[0], T(2) );

    // (Regular) complex
    complex_it<T, 1>  c{}, d{ (T)7 };

    BOOST_CHECK_EQUAL( c[0], T{} );
    BOOST_CHECK_EQUAL( c[1], T{} );
    BOOST_CHECK_EQUAL( d[0], T(7) );
    BOOST_CHECK_EQUAL( d[1], T{} );

    // Quaternions
    complex_it<T, 2>  e{}, f{ (T)19 };

    BOOST_CHECK_EQUAL( e[0], T{} );
    BOOST_CHECK_EQUAL( e[1], T{} );
    BOOST_CHECK_EQUAL( e[2], T{} );
    BOOST_CHECK_EQUAL( e[3], T{} );
    BOOST_CHECK_EQUAL( f[0], T(19) );
    BOOST_CHECK_EQUAL( f[1], T{} );
    BOOST_CHECK_EQUAL( f[2], T{} );
    BOOST_CHECK_EQUAL( f[3], T{} );

    // Octonions
    complex_it<T, 3>  g{}, h{ (T)101 };

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

BOOST_AUTO_TEST_CASE_TEMPLATE( test_multireal_construction, T, test_types )
{
    // (Regular) complex
    complex_it<T, 1> const  a{ (T)2, (T)3 };

    BOOST_CHECK_EQUAL( a[0], T(2) );
    BOOST_CHECK_EQUAL( a[1], T(3) );

    // Quaternions
    complex_it<T, 2> const  b{ (T)5, (T)7 }, c{ (T)11, (T)13, (T)17, (T)19 };

    BOOST_CHECK_EQUAL( b[0], T(5) );
    BOOST_CHECK_EQUAL( b[1], T(7) );
    BOOST_CHECK_EQUAL( b[2], T{} );
    BOOST_CHECK_EQUAL( b[3], T{} );
    BOOST_CHECK_EQUAL( c[0], T(11) );
    BOOST_CHECK_EQUAL( c[1], T(13) );
    BOOST_CHECK_EQUAL( c[2], T(17) );
    BOOST_CHECK_EQUAL( c[3], T(19) );
}

BOOST_AUTO_TEST_CASE( test_same_size_diff_type_conversion )
{
    // Note that same-type/same-size is covered by the automatically-defined
    // copy constructor, which (usually) takes priority over constructor
    // templates (such as the ones used here).

    // Reals
    complex_it<unsigned, 0> const  a{ complex_it<unsigned char, 0>{'\0'} };

    BOOST_CHECK_EQUAL( a[0], 0u );

    // (Regular) Complexes
    complex_it<long, 1> const  b{ complex_it<int, 1>{-2, +3} };

    BOOST_CHECK_EQUAL( b[0], -2L );
    BOOST_CHECK_EQUAL( b[1], +3L );

    // Quaternions
    complex_it<double, 2> const  c{ complex_it<float,2>{+5.5f, -7.0f, +11.0f} };

    BOOST_CHECK_CLOSE( c[0], +5.5, 0.1 );
    BOOST_CHECK_CLOSE( c[1], -7.0, 0.1 );
    BOOST_CHECK_CLOSE( c[2], 11.0, 0.1 );
    BOOST_CHECK_CLOSE( c[3],  0.0, 0.1 );
}

BOOST_AUTO_TEST_CASE( test_barrage_conversion )
{
    // Same type between barrages and composite
    complex_it<int, 0> const     a1{ 2 }, a2{ -3 };
    complex_it<int, 1> const     a{ a1, a2 };
    complex_it<double, 1> const  b1{ -5.5 }, b2{ +7.1, -11.3 };
    complex_it<double, 2> const  b{ b1, b2 };

    BOOST_CHECK_EQUAL( a[0], a1[0] );
    BOOST_CHECK_EQUAL( a[1], a2[0] );

    BOOST_CHECK_CLOSE( b[0], b1[0], 0.1 );
    BOOST_CHECK_CLOSE( b[1], b1[1], 0.1 );
    BOOST_CHECK_CLOSE( b[2], b2[0], 0.1 );
    BOOST_CHECK_CLOSE( b[3], b2[1], 0.1 );

    // Only one barrage
    complex_it<int, 1> const     aa{ a2 };
    complex_it<double, 2> const  bb{ b1 };

    BOOST_CHECK_EQUAL( aa[0], a2[0] );
    BOOST_CHECK_EQUAL( aa[1], 0 );

    BOOST_CHECK_CLOSE( bb[0], b1[0], 0.1 );
    BOOST_CHECK_CLOSE( bb[1], b1[1], 0.1 );
    BOOST_CHECK_CLOSE( bb[2], 0.0, 0.1 );
    BOOST_CHECK_CLOSE( bb[3], 0.0, 0.1 );

    // Mixed types
    complex_it<long, 1> const         c{ complex_it<int, 0>{-13},
     complex_it<long, 0>{17L} };
    complex_it<long double, 2> const  d{ complex_it<float, 1>{-19.4f},
     complex_it<double, 1>{23.0, -29.8} };

    BOOST_CHECK_EQUAL( c[0], -13L );
    BOOST_CHECK_EQUAL( c[1], 17L );

    BOOST_CHECK_CLOSE( d[0], -19.4L, 0.1 );
    BOOST_CHECK_CLOSE( d[1], 0.0L, 0.1 );
    BOOST_CHECK_CLOSE( d[2], 23.0L, 0.1 );
    BOOST_CHECK_CLOSE( d[3], -29.8L, 0.1 );

    // One mixed barrage
    complex_it<short, 1> const  e{ complex_it<char, 0>{125} };

    BOOST_CHECK_EQUAL( e[0], 125 );
    BOOST_CHECK_EQUAL( e[1], 0 );
}

BOOST_AUTO_TEST_SUITE_END()  // constructor_tests

BOOST_AUTO_TEST_SUITE( tuple_tests )

// Check the various compile-time tuple descriptors.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_tuple_support, T, test_types )
{
    using std::tuple_size;
    using std::is_same;
    using std::tuple_element;

    typedef complex_it<T, 0>  real_type;
    typedef complex_it<T, 1>  complex_type;
    typedef complex_it<T, 2>  quaternion_type;
    typedef complex_it<T, 3>  octonion_type;

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

BOOST_AUTO_TEST_SUITE_END()  // tuple_tests

BOOST_AUTO_TEST_SUITE_END()  // complex_it_tests
