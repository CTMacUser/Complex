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
#include <cstdint>
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
    typedef list<int, unsigned, double>                  test_builtin_types;
    typedef list<int, std::intmax_t, mp::int512_t>        test_signed_types;

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

    // Mutability
    o.lower_barrage( oo2 );
    o.upper_barrage( oo1 );
    BOOST_CHECK_EQUAL( oo1[0], oo[4] );
    BOOST_CHECK_EQUAL( oo1[1], oo[5] );
    BOOST_CHECK_EQUAL( oo1[2], oo[6] );
    BOOST_CHECK_EQUAL( oo1[3], oo[7] );
    BOOST_CHECK_EQUAL( oo2[0], oo[0] );
    BOOST_CHECK_EQUAL( oo2[1], oo[1] );
    BOOST_CHECK_EQUAL( oo2[2], oo[2] );
    BOOST_CHECK_EQUAL( oo2[3], oo[3] );

    // Degenerate mutability
    complex_it<T, 0>  s;

    s[ 0 ] = (T)63;
    r.lower_barrage( s );
    BOOST_CHECK_EQUAL( r[0], r.lower_barrage()[0] );
    s[ 0 ] = (T)65;
    r.upper_barrage( s );
    BOOST_CHECK_EQUAL( r[0], r.upper_barrage()[0] );
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
    complex_it<T, 0>  a = {}, b = { (T)2 };

    BOOST_CHECK_EQUAL( a[0], T{} );
    BOOST_CHECK_EQUAL( b[0], T(2) );

    // (Regular) complex
    complex_it<T, 1>  c = {}, d = { (T)7 };

    BOOST_CHECK_EQUAL( c[0], T{} );
    BOOST_CHECK_EQUAL( c[1], T{} );
    BOOST_CHECK_EQUAL( d[0], T(7) );
    BOOST_CHECK_EQUAL( d[1], T{} );

    // Quaternions
    complex_it<T, 2>  e = {}, f = { (T)19 };

    BOOST_CHECK_EQUAL( e[0], T{} );
    BOOST_CHECK_EQUAL( e[1], T{} );
    BOOST_CHECK_EQUAL( e[2], T{} );
    BOOST_CHECK_EQUAL( e[3], T{} );
    BOOST_CHECK_EQUAL( f[0], T(19) );
    BOOST_CHECK_EQUAL( f[1], T{} );
    BOOST_CHECK_EQUAL( f[2], T{} );
    BOOST_CHECK_EQUAL( f[3], T{} );

    // Octonions
    complex_it<T, 3>  g = {}, h = { (T)101 };

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

// Check conversion from a list of scalars.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_multireal_construction, T, test_types )
{
    // (Regular) complex
    complex_it<T, 1> const  a = { (T)2, (T)3 };

    BOOST_CHECK_EQUAL( a[0], T(2) );
    BOOST_CHECK_EQUAL( a[1], T(3) );

    // Quaternions
    complex_it<T, 2> const  b = { (T)5, (T)7 }, c = { (T)11, (T)13, (T)17,
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
    complex_it<unsigned, 0> const  a = { complex_it<unsigned char, 0>{'\0'} };

    BOOST_CHECK_EQUAL( a[0], 0u );

    // (Regular) Complexes
    complex_it<long, 1> const  b = { complex_it<int, 1>{-2, +3} };

    BOOST_CHECK_EQUAL( b[0], -2L );
    BOOST_CHECK_EQUAL( b[1], +3L );

    // Quaternions
    complex_it<double,2> const  c = {complex_it<float,2>{+5.5f, -7.0f, +11.0f}};

    BOOST_CHECK_CLOSE( c[0], +5.5, 0.1 );
    BOOST_CHECK_CLOSE( c[1], -7.0, 0.1 );
    BOOST_CHECK_CLOSE( c[2], 11.0, 0.1 );
    BOOST_CHECK_CLOSE( c[3],  0.0, 0.1 );
}

// Check conversions with immediately-lower rank, any component type.
BOOST_AUTO_TEST_CASE( test_barrage_conversion )
{
    // Same type between barrages and composite
    complex_it<int, 0> const     a1 = { 2 }, a2 = { -3 };
    complex_it<int, 1> const     a = { a1, a2 };
    complex_it<double, 1> const  b1 = { -5.5 }, b2 = { +7.1, -11.3 };
    complex_it<double, 2> const  b = { b1, b2 };

    BOOST_CHECK_EQUAL( a[0], a1[0] );
    BOOST_CHECK_EQUAL( a[1], a2[0] );

    BOOST_CHECK_CLOSE( b[0], b1[0], 0.1 );
    BOOST_CHECK_CLOSE( b[1], b1[1], 0.1 );
    BOOST_CHECK_CLOSE( b[2], b2[0], 0.1 );
    BOOST_CHECK_CLOSE( b[3], b2[1], 0.1 );

    // Only one barrage
    complex_it<int, 1> const     aa = { a2 };
    complex_it<double, 2> const  bb = { b1 };

    BOOST_CHECK_EQUAL( aa[0], a2[0] );
    BOOST_CHECK_EQUAL( aa[1], 0 );

    BOOST_CHECK_CLOSE( bb[0], b1[0], 0.1 );
    BOOST_CHECK_CLOSE( bb[1], b1[1], 0.1 );
    BOOST_CHECK_CLOSE( bb[2], 0.0, 0.1 );
    BOOST_CHECK_CLOSE( bb[3], 0.0, 0.1 );

    // Mixed types
    complex_it<long, 1> const         c = { complex_it<int, 0>{-13},
     complex_it<long, 0>{17L} };
    complex_it<long double, 2> const  d = { complex_it<float, 1>{-19.4f},
     complex_it<double, 1>{23.0, -29.8} };

    BOOST_CHECK_EQUAL( c[0], -13L );
    BOOST_CHECK_EQUAL( c[1], 17L );

    BOOST_CHECK_CLOSE( d[0], -19.4L, 0.1 );
    BOOST_CHECK_CLOSE( d[1], 0.0L, 0.1 );
    BOOST_CHECK_CLOSE( d[2], 23.0L, 0.1 );
    BOOST_CHECK_CLOSE( d[3], -29.8L, 0.1 );

    // One mixed barrage
    complex_it<short, 1> const  e = { complex_it<char, 0>{125} };

    BOOST_CHECK_EQUAL( e[0], 125 );
    BOOST_CHECK_EQUAL( e[1], 0 );
}

// Check conversions with severely-lower rank, any component type.
BOOST_AUTO_TEST_CASE( test_subbarrage_conversion )
{
    // Same type between pieces and whole
    complex_it<int, 0> const  qc[] = { {2}, {-3}, {5}, {-7} };
    complex_it<int, 2> const  q = { qc[0], qc[1], qc[2], qc[3] };
    complex_it<int, 0> const  oc[] = { {11}, {-13}, {17}, {-19} };
    complex_it<int, 3> const  o = { qc[0], qc[1], qc[2], qc[3], oc[0], oc[1],
     oc[2], oc[3] };

    BOOST_CHECK_EQUAL( q, (decltype(q){2, -3, 5, -7}) );
    BOOST_CHECK_EQUAL( o, (decltype(o){2, -3, 5, -7, 11, -13, 17, -19}) );

    // Differing types
    complex_it<long, 1> const       xc1[] = { {23L, -29L}, {31L, -37L} };
    complex_it<int, 1> const        xc2[] = { {-41, 43}, {47, -53} };
    complex_it<long long, 3> const  x = { xc1[0], xc2[1], xc2[0], xc1[1] };
    complex_it<long, 3> const       y = { xc2[1] };

    BOOST_CHECK_EQUAL( x, (decltype(x){23LL, -29LL, 47LL, -53LL, -41LL, 43LL,
     31LL, -37LL}) );
    BOOST_CHECK_EQUAL( y, (decltype(y){47L, -53L, 0L}) );
}

// Check conversions with a greater rank, any component type.
BOOST_AUTO_TEST_CASE( test_supersize_conversion )
{
    // Integer
    complex_it<int, 3> const        o = { -2, 3, -5, 7, -11, 13, -17, 19 };
    complex_it<int, 2> const        q1{ o };
    complex_it<long, 2> const       q2{ o };
    complex_it<long, 1> const       c1{ o };
    complex_it<int, 1> const        c2{ q1 };
    complex_it<long long, 0> const  r1{ o };
    complex_it<int, 0> const        r2{ o };

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
    complex_it<float, 2> const        q = { -23.3f, +29.9f, -31.1f };
    complex_it<double, 1> const       c3{ q };
    complex_it<float, 1> const        c4{ q };
    complex_it<float, 0> const        r3{ q };
    complex_it<long double, 0> const  r4{ q };
    complex_it<double, 0> const       r5{ c3 };

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
    complex_it<T, 0> const      a = { (T)2 };
    complex_it<proxy, 0> const  b = { a };
    complex_it<T, 0> const      c{ b };
    complex_it<T, 3> const      o{ b };

    BOOST_CHECK_EQUAL( a, c );
    BOOST_CHECK_EQUAL( a, o );

    // (Regular) complex, different-length destinations
    complex_it<proxy, 1> const  d = { proxy((T)3), proxy((T)5) };
    complex_it<T, 0> const      e{ d };
    complex_it<T, 1> const      f{ d };
    complex_it<T, 2> const      g{ d };

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

// Check the swapping of states.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_swap, T, test_types )
{
    // Reals
    complex_it<T, 0>  a = { (T)2 }, b = { (T)3 };

    BOOST_CHECK_EQUAL( a[0], T(2) );
    BOOST_CHECK_EQUAL( b[0], T(3) );
    swap( a, b );
    BOOST_CHECK_EQUAL( a[0], T(3) );
    BOOST_CHECK_EQUAL( b[0], T(2) );

    // Quaternions
    complex_it<T, 2>  c = { (T)5, (T)7, (T)11, (T)13 };
    complex_it<T, 2>  d = { (T)17, (T)19, (T)23 };

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
BOOST_AUTO_TEST_CASE_TEMPLATE( test_conj, T, test_types )
{
    // Reals
    complex_it<T, 0> const  a = {},         b = { (T)2 };
    auto const             aa = conj( a ), bb = conj( b );

    BOOST_CHECK_EQUAL( aa[0], +T{} );
    BOOST_CHECK_EQUAL( bb[0], +T(2) );

    // (Regular) complex
    complex_it<T, 1> const  c = { (T)3, (T)5 }, d = { (T)7 };
    auto const             cc = conj( c ),     dd = conj( d );

    BOOST_CHECK_EQUAL( cc[0], +T(3) );
    BOOST_CHECK_EQUAL( cc[1], -T(5) );
    BOOST_CHECK_EQUAL( dd[0], +T(7) );
    BOOST_CHECK_EQUAL( dd[1], -T{} );

    // Quaternions
    complex_it<T, 2> const  e = { (T)11, (T)13, -(T)17, (T)19 };
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
    complex_it<T, 0>  a = {}, b = { (T)2 };

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
    complex_it<T, 1>  c = { (T)7, (T)11 }, d = { (T)13 };

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
    complex_it<T, 2>  e = { (T)31, (T)37, (T)41, (T)43 }, f = { (T)47 };

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
    complex_it<T, 0>  a = {}, b = { (T)2 };

    BOOST_CHECK_EQUAL( a.unreal(), decltype(a){} );
    BOOST_CHECK_EQUAL( b.unreal(), decltype(b){} );
    a.unreal( decltype(a){(T)3} );
    b.unreal( decltype(b){(T)5} );
    BOOST_CHECK_EQUAL( a.unreal(), decltype(a){} );
    BOOST_CHECK_EQUAL( b.unreal(), decltype(b){} );
    BOOST_CHECK_EQUAL( a.real(), T{} );
    BOOST_CHECK_EQUAL( b.real(), T(2) );

    // (Regular) complexes
    complex_it<T, 1>  c = { (T)7, (T)11 }, d = { (T)13, (T)17 };

    BOOST_CHECK_EQUAL( c.unreal(), (decltype( c ){ T{}, T(11) }) );
    BOOST_CHECK_EQUAL( d.unreal(), (decltype( d ){ T{}, T(17) }) );
    c.unreal( decltype(c){(T)19, (T)23} );
    d.unreal( decltype(d){(T)29, (T)31} );
    BOOST_CHECK_EQUAL( c.unreal(), (decltype( c ){ T{}, T(23) }) );
    BOOST_CHECK_EQUAL( d.unreal(), (decltype( d ){ T{}, T(31) }) );
    BOOST_CHECK_EQUAL( c.real(), T(7) );
    BOOST_CHECK_EQUAL( d.real(), T(13) );

    // Quaternions
    complex_it<T, 2>  e = { (T)37, (T)41, T(43), T(47) }, f = { (T)53, (T)59 };

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
    complex_it<T, 0> const  a = {}, b = { (T)2 };

    BOOST_CHECK_EQUAL( norm(a), T{} );
    BOOST_CHECK_EQUAL( norm(b), T(4) );

    // (Regular) complex
    complex_it<T, 1> const  c = {}, d = { (T)3 }, e = { (T)5, (T)7 };

    BOOST_CHECK_EQUAL( norm(c), T{} );
    BOOST_CHECK_EQUAL( norm(d), T(9) );
    BOOST_CHECK_EQUAL( norm(e), T(74) );

    // Quaternions
    complex_it<T, 2> const  f = {}, g = { (T)11, (T)13, (T)17, (T)19 };

    BOOST_CHECK_EQUAL( norm(f), T{} );
    BOOST_CHECK_EQUAL( norm(g), T(940) );
}

// Check (Cayley) norm, with floating-point types.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_norm2, T, test_floating_types )
{
    // Reals
    complex_it<T, 0> const  a = {}, b = { (T)2 };

    BOOST_CHECK_CLOSE( norm(a), T{}, 0.1 );
    BOOST_CHECK_CLOSE( norm(b), T(4), 0.1 );

    // (Regular) complex
    complex_it<T, 1> const  c = {}, d = { -(T)3 }, e = { (T)5, (T)-7 };

    BOOST_CHECK_CLOSE( norm(c), T{}, 0.1 );
    BOOST_CHECK_CLOSE( norm(d), T(9), 0.1 );
    BOOST_CHECK_CLOSE( norm(e), T(74), 0.1 );

    // Quaternions
    complex_it<T, 2> const  f = {}, g = { (T)11, (T)13, -(T)17, (T)19 };

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

// Check the tuple access functions.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_tuple_get, T, test_types )
{
    // Need to bring this into scope.
    using boost::math::get;

    // Reals
    complex_it<T, 0>           a = { (T)2 };
    complex_it<T, 0> const &  aa = a;

    BOOST_CHECK_EQUAL( get<0>(a), T(2) );
    get<0>( a ) = get<0>( complex_it<T, 0>{(T)3} );
    BOOST_CHECK_EQUAL( T(3), get<0>(aa) );

    // Quaternions
    complex_it<T, 2>           b = { (T)5, (T)7, (T)11, (T)13 };
    complex_it<T, 2> const &  bb = b;

    BOOST_CHECK_EQUAL( get<0>(b), T(5) );
    BOOST_CHECK_EQUAL( get<1>(b), T(7) );
    BOOST_CHECK_EQUAL( get<2>(b), T(11) );
    BOOST_CHECK_EQUAL( get<3>(b), T(13) );
    get<0>( b ) = get<0>( complex_it<T, 1>{(T)17, (T)8} );
    get<1>( b ) = get<1>( complex_it<T, 1>{(T)9, (T)19} );
    get<2>( b ) = get<2>( complex_it<T, 2>{(T)1, (T)6, (T)23} );
    get<3>( b ) = get<3>( complex_it<T, 3>{a} );
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
    complex_it<T, 0> const   a = {},  b = { (T)2 };
    auto const              aa = +a, bb = +b;

    BOOST_CHECK_EQUAL( aa[0], +T{} );
    BOOST_CHECK_EQUAL( bb[0], +T(2) );

    // Quaternions
    complex_it<T, 2> const  c = { (T)3, (T)5, (T)7 };
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
    typedef complex_it<int, 0>        real_type;
    typedef complex_it<int, 1>     complex_type;
    typedef complex_it<int, 2>  quaternion_type;

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
BOOST_AUTO_TEST_CASE_TEMPLATE( test_negation, T, test_types )
{
    // Reals
    complex_it<T, 0> const   a = {},  b = { (T)2 };
    auto const              aa = -a, bb = -b;

    BOOST_CHECK_EQUAL( aa[0], -T{} );
    BOOST_CHECK_EQUAL( bb[0], -T(2) );

    // Quaternions
    complex_it<T, 2> const  c = { (T)3, (T)5, (T)7 };
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
    typedef complex_it<int, 0>        real_type;
    typedef complex_it<int, 1>     complex_type;
    typedef complex_it<int, 2>  quaternion_type;

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
BOOST_AUTO_TEST_CASE_TEMPLATE( test_conjugation, T, test_types )
{
    // Reals
    complex_it<T, 0> const   a = {},  b = { (T)2 };
    auto const              aa = ~a, bb = ~b;

    BOOST_CHECK_EQUAL( aa[0], +T{} );
    BOOST_CHECK_EQUAL( bb[0], +T(2) );

    // Quaternions
    complex_it<T, 2> const  c = { (T)3, (T)5, (T)7 };
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
    typedef complex_it<T, 0>        real_type;
    typedef complex_it<T, 1>     complex_type;
    typedef complex_it<T, 2>  quaternion_type;

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
    typedef complex_it<T, 0>        real_type;
    typedef complex_it<T, 1>     complex_type;
    typedef complex_it<T, 2>  quaternion_type;

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

// Check the multiplication operators.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_cayley_multiplication,T,test_signed_types )
{
    // Type-aliases
    typedef complex_it<T, 0>        real_type;
    typedef complex_it<T, 1>     complex_type;
    typedef complex_it<T, 2>  quaternion_type;
    typedef complex_it<T, 3>    octonion_type;

    // Reals
    real_type const  a = { T(8) }, b = { T(-9) }, c = { T(7) };

    BOOST_CHECK_EQUAL( a * b, real_type(T( -72 )) );
    BOOST_CHECK_EQUAL( b * a, real_type(T( -72 )) );

    BOOST_CHECK_EQUAL( (a * b) * c, real_type(T( -504 )) );
    BOOST_CHECK_EQUAL( a * (b * c), real_type(T( -504 )) );

    // (Regular) complexes
    complex_type const  d = { T(2), T(5) }, e = { T(4), T(-6) }, f = { T(-9),
     T(6) }, g = { T(10) };

    BOOST_CHECK_EQUAL( d * e, complex_type(T( 38 ), T( 8 )) );
    BOOST_CHECK_EQUAL( e * d, complex_type(T( 38 ), T( 8 )) );

    BOOST_CHECK_EQUAL( (d * e) * f, complex_type(T( -390 ), T( 156 )) );
    BOOST_CHECK_EQUAL( d * (e * f), complex_type(T( -390 ), T( 156 )) );

    BOOST_CHECK_EQUAL( e * f, complex_type(T( 0 ), T( 78 )) );
    BOOST_CHECK_EQUAL( f * g, complex_type(T( -90 ), T( 60 )) );
    BOOST_CHECK_EQUAL( g * f, complex_type(T( -90 ), T( 60 )) );

    BOOST_CHECK_EQUAL( a * f, complex_type(T( -72 ), T( 48 )) );
    BOOST_CHECK_EQUAL( f * a, complex_type(T( -72 ), T( 48 )) );

    // Quaternions
    quaternion_type const  h = { T(2), T(13), T(-5), T(17) }, k = { T(11), T(3),
     T(-7), T(19) },       m = { T(-1), T(4), T{}, T(-9) },   n = { T(-6) };

    BOOST_CHECK_EQUAL( h * k, quaternion_type(T( -375 ), T( 173 ), T( -265 ),
     T( 149 )) );
    BOOST_CHECK_EQUAL( k * h, quaternion_type(T( -375 ), T( 125 ), T( +127 ),
     T( 301 )) );

    BOOST_CHECK_EQUAL( (h * k) * m, quaternion_type(T( 1024 ), T( 712 ),
     T( 2418 ), T( 4286 )) );
    BOOST_CHECK_EQUAL( h * (k * m), quaternion_type(T( 1024 ), T( 712 ),
     T( 2418 ), T( 4286 )) );

    BOOST_CHECK_EQUAL( m * n, quaternion_type(T( 6 ), T( -24 ), T{}, T( 54 )) );
    BOOST_CHECK_EQUAL( n * m, quaternion_type(T( 6 ), T( -24 ), T{}, T( 54 )) );

    BOOST_CHECK_EQUAL( m * c, quaternion_type(T( -7 ), T( 28 ),T{}, T( -63 )) );
    BOOST_CHECK_EQUAL( c * m, quaternion_type(T( -7 ), T( 28 ),T{}, T( -63 )) );
    BOOST_CHECK_EQUAL( d * m, quaternion_type(T(-22), T(3), T(+45), T(-18)) );
    BOOST_CHECK_EQUAL( m * d, quaternion_type(T(-22), T(3), T(-45), T(-18)) );

    // Octonions
    octonion_type const p = {T(7), T(-2), T(0), T(7), T(-8), T(6), T(1), T(-6)},
     q = { T(3), T(3), T(-13), T(-8), T(11), T(12), T(-4), T(-11) },
     r = { T(-5), T(9), T(10), T(6), T(-10), T(-11), T(9), T(9) };

    BOOST_CHECK_EQUAL( p * q, octonion_type(T( 37 ), T( -21 ), T( -59 ),
     T( 181 ), T( 207 ), T( 162 ), T( -221 ), T( -9 )) );
    BOOST_CHECK_EQUAL( q * p, octonion_type(T( 37 ), T( 51 ), T( -123 ),
     T( -251 ), T( -101 ), T( 42 ), T( 171 ), T( -181 )) );

    BOOST_CHECK_EQUAL( (p * q) * r, octonion_type(T( 5430 ), T(  -475 ),
     T(  3432 ), T( 2384 ), T( -3540 ), T(  526 ), T( 2813 ), T( -5445 )) );
    BOOST_CHECK_EQUAL( p * (q * r), octonion_type(T( 5430 ), T( -1581 ),
     T( -2530 ), T( 3460 ), T( -5078 ), T( 1362 ), T( 4369 ), T(  -675 )) );

    BOOST_CHECK_EQUAL( (p * p) * p, octonion_type(T( -3647 ), T( 86 ), T( 0 ),
     T( -301 ), T( 344 ), T( -258 ), T( -43 ), T( 258 )) );
    BOOST_CHECK_EQUAL( p * (p * p), octonion_type(T( -3647 ), T( 86 ), T( 0 ),
     T( -301 ), T( 344 ), T( -258 ), T( -43 ), T( 258 )) );

    BOOST_CHECK_EQUAL( b * q, octonion_type(T( -27 ), T( -27 ), T( 117 ),
     T( 72 ), T( -99 ), T( -108 ), T( 36 ), T( 99 )) );
    BOOST_CHECK_EQUAL( q * b, octonion_type(T( -27 ), T( -27 ), T( 117 ),
     T( 72 ), T( -99 ), T( -108 ), T( 36 ), T( 99 )) );
    BOOST_CHECK_EQUAL( e * r, octonion_type(T( 34 ), T( 66 ), T( 76 ), T( -36 ),
     T( -106 ), T(   16 ), T( -18 ), T(  90 )) );
    BOOST_CHECK_EQUAL( r * e, octonion_type(T( 34 ), T( 66 ), T(  4 ), T(  84 ),
     T(   26 ), T( -104 ), T(  90 ), T( -18 )) );
    BOOST_CHECK_EQUAL( h * r, octonion_type(T( -179 ), T( -247 ), T( 120 ),
     T(  102 ), T(  15 ), T(  46 ), T(  372 ), T( -214 )) );
    BOOST_CHECK_EQUAL( r * h, octonion_type(T( -179 ), T(  153 ), T( -30 ),
     T( -248 ), T( -55 ), T( -90 ), T( -336 ), T(  250 )) );

    // Multiply-assignment
    real_type        aa = a;
    complex_type     dd = d;
    quaternion_type  hh = h;
    octonion_type    pp = p;

    aa *= b;
    BOOST_CHECK_EQUAL( aa, a * b );
    dd *= e;
    BOOST_CHECK_EQUAL( dd, d * e );
    dd = f;
    dd *= a;
    BOOST_CHECK_EQUAL( dd, f * a );
    hh *= k;
    BOOST_CHECK_EQUAL( hh, h * k );
    hh = m;
    hh *= c;
    BOOST_CHECK_EQUAL( hh, m * c );
    hh = m;
    hh *= d;
    BOOST_CHECK_EQUAL( hh, m * d );
    pp *= q;
    BOOST_CHECK_EQUAL( pp, p * q );
    pp = q;
    pp *= b;
    BOOST_CHECK_EQUAL( pp, q * b );
    pp = r;
    pp *= e;
    BOOST_CHECK_EQUAL( pp, r * e );
    pp = r;
    pp *= h;
    BOOST_CHECK_EQUAL( pp, r * h );
}

// Check the division-with-scalar (including modulus) operators.
BOOST_AUTO_TEST_CASE( test_scalar_division_and_modulus )
{
    // Integer
    complex_it<int, 0>  a = { 9 };
    complex_it<int, 1>  b = { 8, 3 };
    complex_it<int, 2>  c = { 0, 4, 8, 12 };

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
    complex_it<double, 0>  d = { -16.5 };
    complex_it<double, 1>  e = { 7.0, 0.0 };
    complex_it<double, 2>  f = { -3.0, 2.1, 1.21, -100.7 };

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

// Check the division (including modulus) operators.
BOOST_AUTO_TEST_CASE( test_division_and_modulus )
{
    // Rational
    {
    typedef boost::rational<long long>      rational_type;
    typedef complex_it<rational_type, 0>        real_type;
    typedef complex_it<rational_type, 1>     complex_type;
    typedef complex_it<rational_type, 2>  quaternion_type;
    typedef complex_it<rational_type, 3>    octonion_type;

    real_type const           two = { 2 },five = { 5 }, six = { 6 };
    complex_type const  eight_ten = { 8, 10 },      one_one = { 1, 1 };
    quaternion_type const  primed = { 2, 3, 5, 7 }, squared = { 4, 9, 25, 49 };
    octonion_type const     super = { rational_type{57,180}, rational_type{-222,
     180}, rational_type{287, 180}, rational_type{177, 180}, rational_type{-105,
     180}, rational_type{-112, 180}, rational_type{-168, 180},
     rational_type{-78, 180} },                       omega = {
     rational_type{-571, 143}, rational_type{-34, 143}, rational_type{703, 143},
     rational_type{676, 143}, rational_type{572, 143}, rational_type{-285, 143},
     rational_type{263, 143}, rational_type{828, 143} };
    auto  temp0 = six;
    auto  temp1 = eight_ten;
    auto  temp2 = primed;
    auto  temp3 = omega;

    BOOST_CHECK_EQUAL( six / two, real_type(rational_type{ 3 }) );
    BOOST_CHECK_EQUAL( five / six, real_type(rational_type( 5, 6 )) );
    temp0 /= two;
    BOOST_CHECK_EQUAL( temp0, real_type(rational_type{ 3 }) );

    BOOST_CHECK_EQUAL( eight_ten / two, complex_type(4, 5) );
    BOOST_CHECK_EQUAL( eight_ten / five, complex_type(rational_type(8,5), 2) );
    BOOST_CHECK_EQUAL( two / eight_ten, complex_type(rational_type( 4, 41 ),
     rational_type( -5, 41 )) );
    BOOST_CHECK_EQUAL( six / one_one, complex_type(3, -3) );
    temp1 /= two;
    BOOST_CHECK_EQUAL( temp1, complex_type(4, 5) );

    BOOST_CHECK_EQUAL( eight_ten / one_one, complex_type(9, 1) );
    BOOST_CHECK_EQUAL( (eight_ten + rational_type{ 1 }) / one_one,
     complex_type(rational_type( 19, 2 ), rational_type( 1, 2 )) );
    temp1 = eight_ten;
    temp1 /= one_one;
    BOOST_CHECK_EQUAL( temp1, complex_type(9, 1) );

    BOOST_CHECK_EQUAL( primed / two, quaternion_type(rational_type{ 1 },
     rational_type{ 3, 2 }, rational_type{ 5, 2 }, rational_type{ 7, 2 }) );
    BOOST_CHECK_EQUAL( primed / -two, quaternion_type(rational_type{ -1 },
     rational_type{ -3, 2 }, rational_type{ -5, 2 }, rational_type{ -7, 2 }) );
    BOOST_CHECK_EQUAL( two / primed, quaternion_type(rational_type{ 4, 87 },
     rational_type{ -6,87 },rational_type{ -10,87 },rational_type{ -14,87 }) );
    BOOST_CHECK_EQUAL( (two / primed) * (primed / two), rational_type{1} );
    BOOST_CHECK_EQUAL( (rational_type{25} * two) / primed, rational_type(25, 87)
     * quaternion_type(4, -6, -10, -14) );
    BOOST_CHECK_EQUAL( primed / eight_ten, quaternion_type(46, 4, -30, 106) /
     rational_type{164} );
    BOOST_CHECK_EQUAL( eight_ten / primed, quaternion_type(46, -4, 30, -106) /
     rational_type{87} );
    BOOST_CHECK_EQUAL( squared / primed, quaternion_type(503, 76, -54, 100) /
     rational_type{87} );
    BOOST_CHECK_EQUAL( primed / squared, quaternion_type(503, -76, 54, -100) /
     rational_type{3123} );
    temp2 /= two;
    BOOST_CHECK_EQUAL( temp2, quaternion_type(rational_type{ 1 },
     rational_type{ 3, 2 }, rational_type{ 5, 2 }, rational_type{ 7, 2 }) );
    temp2 = primed;
    temp2 /= eight_ten;
    BOOST_CHECK_EQUAL( temp2, quaternion_type(46, 4, -30, 106) /
     rational_type{164} );
    temp2 = primed;
    temp2 /= squared;
    BOOST_CHECK_EQUAL( temp2, quaternion_type(503, -76, 54, -100) /
     rational_type{3123} );

    // (The above calculations can fit within "int," but the next two require
    // "long," and the third one needs "long long" for the intermediate
    // calculations.)
    BOOST_CHECK_EQUAL( super / omega, octonion_type(22809358L, -21944780L,
     -43116931L, -4047329L, 68594526L, 49063586L, 7821099L, -33736703L) /
     rational_type{439477920L} );
    BOOST_CHECK_EQUAL( omega / super, octonion_type(7177770L, 6905700L,
     13568265L, 1273635L, -21585690L, -15439590L, -2461185L, 10616445L) /
     rational_type{8011861L} );
    BOOST_CHECK_EQUAL( (super / omega) * (omega / super), rational_type{1LL} );
    temp3 /= super;
    BOOST_CHECK_EQUAL( temp3, octonion_type(7177770L, 6905700L,
     13568265L, 1273635L, -21585690L, -15439590L, -2461185L, 10616445L) /
     rational_type{8011861L} );
    }

    // Floating
    {
    complex_it<double, 0> const  two = { 2.0 }, five = { 5.0 }, six = { 6.0 };
    complex_it<double, 1> const  eight_ten = { 8.0,10.0 },one_one = { 1.0,1.0 };
    complex_it<double, 2> const  primed = { 2.0, 3.0, 5.0, 7.0 },
     squared = { 4.0, 9.0, 25.0, 49.0 };
    complex_it<double, 3> const  super = { 0.316667, -1.233333, 1.594444,
     0.983333, -0.583333, -0.622222, -0.933333, -0.433333 },
     omega = { -3.993007, -0.237762, 4.916084, 4.727273, 4.0, -1.993007,
     1.839161, 5.79021 };
    auto  temp0 = six;
    auto  temp1 = eight_ten;
    auto  temp2 = primed;
    auto  temp3 = omega;

    BOOST_CHECK_CLOSE( (six / two)[0], 3.0, 0.1 );
    BOOST_CHECK_CLOSE( (five / six)[0], 0.8333, 0.1 );

    temp0 /= two;
    BOOST_CHECK_CLOSE( temp0[0], 3.0, 0.1 );

    auto const  ett = eight_ten / two, etf = eight_ten / five;
    auto const  tet = two / eight_ten, soo = six / one_one;

    BOOST_CHECK_CLOSE( ett[0], 4.0, 0.1 );
    BOOST_CHECK_CLOSE( ett[1], 5.0, 0.1 );

    BOOST_CHECK_CLOSE( etf[0], 1.6, 0.1 );
    BOOST_CHECK_CLOSE( etf[1], 2.0, 0.1 );

    BOOST_CHECK_CLOSE( tet[0],  4.0 / 41.0, 0.1 );  //  0.09756
    BOOST_CHECK_CLOSE( tet[1], -5.0 / 41.0, 0.1 );  // -0.12195

    BOOST_CHECK_CLOSE( soo[0],  3.0, 0.1 );
    BOOST_CHECK_CLOSE( soo[1], -3.0, 0.1 );

    temp1 /= two;
    BOOST_CHECK_CLOSE( temp1[0], 4.0, 0.1 );
    BOOST_CHECK_CLOSE( temp1[1], 5.0, 0.1 );

    auto const  etoo = eight_ten / one_one, etooo = (eight_ten + 1.0) / one_one;

    BOOST_CHECK_CLOSE( etoo[0], 9.0, 0.1 );
    BOOST_CHECK_CLOSE( etoo[1], 1.0, 0.1 );

    BOOST_CHECK_CLOSE( etooo[0], 9.5, 0.1 );
    BOOST_CHECK_CLOSE( etooo[1], 0.5, 0.1 );

    temp1 = eight_ten;
    temp1 /= one_one;
    BOOST_CHECK_CLOSE( temp1[0], 9.0, 0.1 );
    BOOST_CHECK_CLOSE( temp1[1], 1.0, 0.1 );

    auto const  pt = primed / two, pnt = primed / -two;
    auto const  tp = two / primed, tppt = tp * pt;
    auto const  tftp = (25.0 * two) / primed;

    BOOST_CHECK_CLOSE( pt[0], 1.0, 0.1 );
    BOOST_CHECK_CLOSE( pt[1], 1.5, 0.1 );
    BOOST_CHECK_CLOSE( pt[2], 2.5, 0.1 );
    BOOST_CHECK_CLOSE( pt[3], 3.5, 0.1 );

    BOOST_CHECK_CLOSE( pnt[0], -1.0, 0.1 );
    BOOST_CHECK_CLOSE( pnt[1], -1.5, 0.1 );
    BOOST_CHECK_CLOSE( pnt[2], -2.5, 0.1 );
    BOOST_CHECK_CLOSE( pnt[3], -3.5, 0.1 );

    BOOST_CHECK_CLOSE( tp[0],   4.0 / 87.0, 0.1 );  //  0.045977
    BOOST_CHECK_CLOSE( tp[1],  -6.0 / 87.0, 0.1 );  // -0.068966
    BOOST_CHECK_CLOSE( tp[2], -10.0 / 87.0, 0.1 );  // -0.114943
    BOOST_CHECK_CLOSE( tp[3], -14.0 / 87.0, 0.1 );  // -0.16092

    BOOST_CHECK_CLOSE( tppt[0], 1.0, 0.1 );
    BOOST_CHECK_SMALL( tppt[1], 0.001 );
    BOOST_CHECK_SMALL( tppt[2], 0.001 );
    BOOST_CHECK_SMALL( tppt[3], 0.001 );

    BOOST_CHECK_CLOSE( tftp[0], 25.0 / 87.0 *   4.0, 0.1 );  //  1.149425
    BOOST_CHECK_CLOSE( tftp[1], 25.0 / 87.0 *  -6.0, 0.1 );  // -1.724138
    BOOST_CHECK_CLOSE( tftp[2], 25.0 / 87.0 * -10.0, 0.1 );  // -2.873563
    BOOST_CHECK_CLOSE( tftp[3], 25.0 / 87.0 * -14.0, 0.1 );  // -4.022989

    temp2 /= two;
    BOOST_CHECK_CLOSE( temp2[0], 1.0, 0.1 );
    BOOST_CHECK_CLOSE( temp2[1], 1.5, 0.1 );
    BOOST_CHECK_CLOSE( temp2[2], 2.5, 0.1 );
    BOOST_CHECK_CLOSE( temp2[3], 3.5, 0.1 );

    auto const  pet = primed / eight_ten, etp = eight_ten / primed;
    auto const  sp = squared / primed, ps = primed / squared;

    BOOST_CHECK_CLOSE( pet[0],  46.0 / 164.0, 0.1 );  //  0.280488
    BOOST_CHECK_CLOSE( pet[1],   4.0 / 164.0, 0.1 );  //  0.024390
    BOOST_CHECK_CLOSE( pet[2], -30.0 / 164.0, 0.1 );  // -0.182927
    BOOST_CHECK_CLOSE( pet[3], 106.0 / 164.0, 0.1 );  //  0.646341

    BOOST_CHECK_CLOSE( etp[0],   46.0 / 87.0, 0.1 );  //  0.528736
    BOOST_CHECK_CLOSE( etp[1],   -4.0 / 87.0, 0.1 );  // -0.045977
    BOOST_CHECK_CLOSE( etp[2],   30.0 / 87.0, 0.1 );  //  0.344828
    BOOST_CHECK_CLOSE( etp[3], -106.0 / 87.0, 0.1 );  // -1.218391

    BOOST_CHECK_CLOSE( sp[0], 503.0 / 87.0, 0.1 );  //  5.781609
    BOOST_CHECK_CLOSE( sp[1],  76.0 / 87.0, 0.1 );  //  0.873563
    BOOST_CHECK_CLOSE( sp[2], -54.0 / 87.0, 0.1 );  // -0.62069
    BOOST_CHECK_CLOSE( sp[3], 100.0 / 87.0, 0.1 );  //  1.149425

    BOOST_CHECK_CLOSE( ps[0],  503.0 / 3123.0, 0.1 );  //  0.161063
    BOOST_CHECK_CLOSE( ps[1],  -76.0 / 3123.0, 0.1 );  // -0.024336
    BOOST_CHECK_CLOSE( ps[2],   54.0 / 3123.0, 0.1 );  //  0.017291
    BOOST_CHECK_CLOSE( ps[3], -100.0 / 3123.0, 0.1 );  // -0.032020

    temp2 = primed;
    temp2 /= eight_ten;
    BOOST_CHECK_CLOSE( temp2[0],  46.0 / 164.0, 0.1 );
    BOOST_CHECK_CLOSE( temp2[1],   4.0 / 164.0, 0.1 );
    BOOST_CHECK_CLOSE( temp2[2], -30.0 / 164.0, 0.1 );
    BOOST_CHECK_CLOSE( temp2[3], 106.0 / 164.0, 0.1 );

    temp2 = primed;
    temp2 /= squared;
    BOOST_CHECK_CLOSE( temp2[0],  503.0 / 3123.0, 0.1 );
    BOOST_CHECK_CLOSE( temp2[1],  -76.0 / 3123.0, 0.1 );
    BOOST_CHECK_CLOSE( temp2[2],   54.0 / 3123.0, 0.1 );
    BOOST_CHECK_CLOSE( temp2[3], -100.0 / 3123.0, 0.1 );

    auto const  so = super / omega, os = omega / super, soos = so * os;

    BOOST_CHECK_CLOSE( so[0],  22809358.0 / 439477920.0, 0.1 );  //  0.051901
    BOOST_CHECK_CLOSE( so[1], -21944780.0 / 439477920.0, 0.1 );  // -0.049934
    BOOST_CHECK_CLOSE( so[2], -43116931.0 / 439477920.0, 0.1 );  // -0.098109
    BOOST_CHECK_CLOSE( so[3],  -4047329.0 / 439477920.0, 0.1 );  // -0.009209
    BOOST_CHECK_CLOSE( so[4],  68594526.0 / 439477920.0, 0.1 );  //  0.156082
    BOOST_CHECK_CLOSE( so[5],  49063586.0 / 439477920.0, 0.1 );  //  0.111641
    BOOST_CHECK_CLOSE( so[6],   7821099.0 / 439477920.0, 0.1 );  //  0.017796
    BOOST_CHECK_CLOSE( so[7], -33736703.0 / 439477920.0, 0.1 );  // -0.076765

    BOOST_CHECK_CLOSE( os[0],   7177770.0 / 8011861.0, 0.1 );  //  0.895893
    BOOST_CHECK_CLOSE( os[1],   6905700.0 / 8011861.0, 0.1 );  //  0.861935
    BOOST_CHECK_CLOSE( os[2],  13568265.0 / 8011861.0, 0.1 );  //  1.693522
    BOOST_CHECK_CLOSE( os[3],   1273635.0 / 8011861.0, 0.1 );  //  0.158969
    BOOST_CHECK_CLOSE( os[4], -21585690.0 / 8011861.0, 0.1 );  // -2.694217
    BOOST_CHECK_CLOSE( os[5], -15439590.0 / 8011861.0, 0.1 );  // -1.927092
    BOOST_CHECK_CLOSE( os[6],  -2461185.0 / 8011861.0, 0.1 );  // -0.307193
    BOOST_CHECK_CLOSE( os[7],  10616445.0 / 8011861.0, 0.1 );  //  1.325091

    BOOST_CHECK_CLOSE( soos[0], 1.0, 0.1 );
    BOOST_CHECK_SMALL( soos[1], 0.001 );
    BOOST_CHECK_SMALL( soos[2], 0.001 );
    BOOST_CHECK_SMALL( soos[3], 0.001 );
    BOOST_CHECK_SMALL( soos[4], 0.001 );
    BOOST_CHECK_SMALL( soos[5], 0.001 );
    BOOST_CHECK_SMALL( soos[6], 0.001 );
    BOOST_CHECK_SMALL( soos[7], 0.001 );

    temp3 /= super;
    BOOST_CHECK_CLOSE( temp3[0],   7177770.0 / 8011861.0, 0.1 );
    BOOST_CHECK_CLOSE( temp3[1],   6905700.0 / 8011861.0, 0.1 );
    BOOST_CHECK_CLOSE( temp3[2],  13568265.0 / 8011861.0, 0.1 );
    BOOST_CHECK_CLOSE( temp3[3],   1273635.0 / 8011861.0, 0.1 );
    BOOST_CHECK_CLOSE( temp3[4], -21585690.0 / 8011861.0, 0.1 );
    BOOST_CHECK_CLOSE( temp3[5], -15439590.0 / 8011861.0, 0.1 );
    BOOST_CHECK_CLOSE( temp3[6],  -2461185.0 / 8011861.0, 0.1 );
    BOOST_CHECK_CLOSE( temp3[7],  10616445.0 / 8011861.0, 0.1 );
    }

    // Integer
    {
    typedef complex_it<int, 0>        real_type;
    typedef complex_it<int, 1>     complex_type;
    typedef complex_it<int, 2>  quaternion_type;
    typedef complex_it<int, 3>    octonion_type;

    complex_it<int, 0> const  two = { 2 }, five = { 5 }, six = { 6 };
    complex_it<int, 1> const  eight_ten = { 8, 10 }, one_one = { 1, 1 };
    complex_it<int, 2> const  primed = { 2, 3, 5, 7 }, squared = {4, 9, 25, 49};
    complex_it<long,3> const  super = { 39, -81, 88, 60, -75, -94, -73, 68 },
     omega = { -2030, -8258, 5636, 4613, -24813, -11991, 4862, -53 };
    auto  temp0 = six;
    auto  temp1 = eight_ten;
    auto  temp2 = primed;
    auto  temp3 = omega;

    BOOST_CHECK_EQUAL( six / two, real_type(3) );
    BOOST_CHECK_EQUAL( six % two, real_type{} );
    BOOST_CHECK_EQUAL( five / six, real_type{} );
    BOOST_CHECK_EQUAL( five % six, real_type(5) );
    BOOST_CHECK_EQUAL( six / five, real_type(1) );
    BOOST_CHECK_EQUAL( six % five, real_type(1) );

    BOOST_CHECK_EQUAL( temp0 /= two, real_type(3) );
    temp0 = six;
    BOOST_CHECK_EQUAL( temp0 %= two, real_type{} );
    temp0 = five;
    BOOST_CHECK_EQUAL( temp0 /= six, real_type{} );
    temp0 = five;
    BOOST_CHECK_EQUAL( temp0 %= six, real_type(5) );
    temp0 = six;
    BOOST_CHECK_EQUAL( temp0 /= five, real_type(1) );
    temp0 = six;
    BOOST_CHECK_EQUAL( temp0 %= five, real_type(1) );

    BOOST_CHECK_EQUAL( eight_ten / two, complex_type(4, 5) );
    BOOST_CHECK_EQUAL( eight_ten % two, complex_type{} );
    BOOST_CHECK_EQUAL( eight_ten / five, complex_type(1, 2) );
     // = trunc( {8/5, 2} )
    BOOST_CHECK_EQUAL( eight_ten % five, complex_type(3, 0) );
    BOOST_CHECK_EQUAL( two / eight_ten, complex_type{} );
     // = trunc( {4/41, -5/41} )
    BOOST_CHECK_EQUAL( two % eight_ten, two );
    BOOST_CHECK_EQUAL( six / one_one, complex_type(3, -3) );
    BOOST_CHECK_EQUAL( six % one_one, complex_type{} );

    BOOST_CHECK_EQUAL( temp1 /= two, complex_type(4, 5) );
    temp1 = eight_ten;
    BOOST_CHECK_EQUAL( temp1 %= two, complex_type{} );

    BOOST_CHECK_EQUAL( eight_ten / one_one, complex_type(9, 1) );
    BOOST_CHECK_EQUAL( eight_ten % one_one, complex_type{} );
    BOOST_CHECK_EQUAL( (eight_ten + 1) / one_one, complex_type(9, 0) );
     // = trunc( {9.5, 1/2} )
    BOOST_CHECK_EQUAL( (eight_ten + 1) % one_one, complex_type(0, 1) );

    temp1 = eight_ten;
    BOOST_CHECK_EQUAL( temp1 /= one_one, complex_type(9, 1) );
    temp1 = eight_ten;
    BOOST_CHECK_EQUAL( temp1 %= one_one, complex_type{} );

    BOOST_CHECK_EQUAL( primed / two, quaternion_type(1, 1, 2, 3) );
     // = trunc( {1, 3/2, 5/2, 7/2} )
    BOOST_CHECK_EQUAL( primed % two, quaternion_type(0, 1, 1, 1) );
    BOOST_CHECK_EQUAL( primed / -two, quaternion_type(-1, -1, -2, -3) );
     // = trunc( {-1, -3/2, -5/2, -7/2} )
    BOOST_CHECK_EQUAL( primed % -two, quaternion_type( 0, +1, +1, +1) );
    BOOST_CHECK_EQUAL( two / primed, quaternion_type{} );
     // = trunc( {4/87, -6/87, -10/87, -14/87} )
    BOOST_CHECK_EQUAL( two % primed, two );
    BOOST_CHECK_EQUAL( (25 * two) / primed, quaternion_type(1, -1, -2, -4) );
     // = trunc( {100/87, -150/87, -250/87, -350/87} )
    BOOST_CHECK_EQUAL( (25 * two) % primed, quaternion_type(7, -7, 4) );

    BOOST_CHECK_EQUAL( temp2 /= two, quaternion_type(1, 1, 2, 3) );
    temp2 = primed;
    BOOST_CHECK_EQUAL( temp2 %= two, quaternion_type(0, 1, 1, 1) );

    BOOST_CHECK_EQUAL( primed / eight_ten, quaternion_type{} );
     // = trunc( {23/82, 1/41, -15/82, 53/82} )
    BOOST_CHECK_EQUAL( primed % eight_ten, primed );
    BOOST_CHECK_EQUAL( eight_ten / primed, quaternion_type(0, 0, 0, -1) );
     // = trunc( {46/87, -4/87, 30/87, -106/87} )
    BOOST_CHECK_EQUAL( eight_ten % primed, quaternion_type(1, 5, 3, 2) );
    BOOST_CHECK_EQUAL( squared / primed, quaternion_type(5, 0, 0, 1) );
     // = trunc( {503/87, 76/87, -54/87, 100/87} )
    BOOST_CHECK_EQUAL( squared % primed, quaternion_type(1, -1, -3, 12) );
    BOOST_CHECK_EQUAL( primed / squared, quaternion_type{} );
     // = trunc( {503/3123, -76/3123, 6/347, -100/3123} )
    BOOST_CHECK_EQUAL( primed % squared, primed );

    temp2 = primed;
    BOOST_CHECK_EQUAL( temp2 /= eight_ten, quaternion_type{} );
    temp2 = primed;
    BOOST_CHECK_EQUAL( temp2 %= eight_ten, primed );
    temp2 = primed;
    BOOST_CHECK_EQUAL( temp2 /= squared, quaternion_type{} );
    temp2 = primed;
    BOOST_CHECK_EQUAL( temp2 %= squared, primed );

    BOOST_CHECK_EQUAL( super / omega, octonion_type{} );
     // = trunc( {3992075/908470632, 762529/454235316, 269753/227117658,
     // -3595525/908470632, 459349/227117658, -582805/302823544,
     // 289247/302823544, 66950/113558829} )
    BOOST_CHECK_EQUAL( super % omega, super );
    BOOST_CHECK_EQUAL( omega / super, octonion_type(91, -34, -24, 82, -41, 39,
     -19, -12) );
     // = trunc( {798415/8768, -762529/21920, -1969/80, 719105/8768,
     // -459349/10960, 349683/8768, -867741/43840, -6695/548} )
    BOOST_CHECK_EQUAL( omega % super, octonion_type(-37, 148, 74, -54, -314, 83,
     74, 36) );
 
    BOOST_CHECK_EQUAL( temp3 /= super, octonion_type(91, -34, -24, 82, -41, 39,
     -19, -12) );
    temp3 = omega;
    BOOST_CHECK_EQUAL( temp3 %= super, octonion_type(-37, 148, 74, -54, -314,
     83, 74, 36) );
    }
}

BOOST_AUTO_TEST_SUITE_END()  // operator_tests

BOOST_AUTO_TEST_SUITE( function_tests )

// Check free-function real, imag, and unreal.
BOOST_AUTO_TEST_CASE_TEMPLATE( test_real_imag_unreal, T, test_types )
{
    // Reals
    complex_it<T, 0> const  a = {}, b = { (T)2 };

    BOOST_CHECK_EQUAL( real(a), T{} );
    BOOST_CHECK_EQUAL( imag(a), T{} );
    BOOST_CHECK_EQUAL( unreal(a), decltype(a){} );
    BOOST_CHECK_EQUAL( real(b), T(2) );
    BOOST_CHECK_EQUAL( imag(b), T{} );
    BOOST_CHECK_EQUAL( unreal(b), decltype(b){} );

    // (Regular) complexes
    complex_it<T, 1> const  c = {}, d = { (T)3, (T)5 };

    BOOST_CHECK_EQUAL( real(c), T{} );
    BOOST_CHECK_EQUAL( imag(c), T{} );
    BOOST_CHECK_EQUAL( unreal(c), decltype(c){} );
    BOOST_CHECK_EQUAL( real(d), T(3) );
    BOOST_CHECK_EQUAL( imag(d), T(5) );
    BOOST_CHECK_EQUAL( unreal(d), (decltype(d){ T{}, T(5) }) );

    // Quaternions
    complex_it<T, 2> const  e = {}, f = { (T)7, (T)11, (T)13 };

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
    complex_it<int, 0>  a = { +4 }, b = { -3 };
    complex_it<int, 1>  c = { -4, +3 }, d = { a, b }, e = { -a, +b },
                        f = { +a, -b };
    complex_it<int, 2>  g = { 5, 0, -2, 1 };

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
    complex_it<double, 0>  aa = a, bb = b;
    complex_it<double, 1>  cc = c, dd = d, ee = e, ff = f;
    complex_it<double, 2>  gg = g, h = { +3.0, -4.0, +12.0, -84.0 };
    complex_it<double, 3>  k = { 6.7, -0.9, -11.2, 0.01, 4.33, -8.25, 255.5 };

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
    complex_it<T, 0> const  a = { T(+4) }, b = {}, c = { T(-3) },
                            d = { T(0.1) }, e = { T(-0.5) };

    BOOST_CHECK_CLOSE( sgn(a)[0], T(+1.0), 0.1 );
    BOOST_CHECK_CLOSE( sgn(b)[0], T( 0.0), 0.1 );
    BOOST_CHECK_CLOSE( sgn(c)[0], T(-1.0), 0.1 );
    BOOST_CHECK_CLOSE( sgn(d)[0], T(+1.0), 0.1 );
    BOOST_CHECK_CLOSE( sgn(e)[0], T(-1.0), 0.1 );

    // (Regular) complexes
    complex_it<T, 1> const  f = {}, g = { T(3.0), T(-4.0) };
    auto const              f_sgn = sgn( f ), g_sgn = sgn( g );

    BOOST_CHECK_CLOSE( f_sgn[0], T{}, 0.1 );
    BOOST_CHECK_CLOSE( f_sgn[1], T{}, 0.1 );
    BOOST_CHECK_CLOSE( g_sgn[0], T(+0.6), 0.1 );
    BOOST_CHECK_CLOSE( g_sgn[1], T(-0.8), 0.1 );

    // Quaternions
    complex_it<T, 2> const  h = { T{}, T{}, T(0.1), T{} }, k = {},
                            m = { T(+6.0), T(-8.0), T(+24.0), T(-168.0) };
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

BOOST_AUTO_TEST_SUITE_END()  // complex_it_tests
