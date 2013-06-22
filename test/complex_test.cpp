//  Boost Complex Numbers, general ops unit test program file  ---------------//

//  Copyright 2013 Daryle Walker.
//  Distributed under the Boost Software License, Version 1.0.  (See the
//  accompanying file LICENSE_1_0.txt or a copy at
//  <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

#define BOOST_TEST_MAIN  "Complex Number Unit Tests"

#include <boost/test/unit_test.hpp>
//#include <boost/test/floating_point_comparison.hpp>
#include <boost/mpl/list.hpp>

#include "boost/math/complex.hpp"

// Put standard includes here.

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>


// Common definitions  -------------------------------------------------------//

namespace {

    // Save time writing out long types
    using boost::multiprecision::cpp_int;
    using boost::multiprecision::cpp_dec_float_50;
    using boost::math::complex_rt;
    using boost::math::complex_it;

    // Sample testing types for components
    typedef boost::mpl::list<int, unsigned, double, cpp_int, cpp_dec_float_50>
      test_types;

}

// Flag un-printable types here.


BOOST_AUTO_TEST_SUITE( general_complex_tests )

// Complex-number structures should have minimal padding.
BOOST_AUTO_TEST_CASE_TEMPLATE( complex_size_demo, T, test_types )
{
    BOOST_REQUIRE( sizeof(complex_rt<T, 0>) >= sizeof(T) );
    BOOST_REQUIRE( sizeof(complex_rt<T, 1>) >= sizeof(T) * 2 );
    BOOST_REQUIRE( sizeof(complex_rt<T, 2>) >= sizeof(T) * 4 );
    BOOST_REQUIRE( sizeof(complex_rt<T, 3>) >= sizeof(T) * 8 );

    BOOST_REQUIRE( sizeof(complex_it<T, 0>) >= sizeof(T) );
    BOOST_REQUIRE( sizeof(complex_it<T, 1>) >= sizeof(T) * 2 );
    BOOST_REQUIRE( sizeof(complex_it<T, 2>) >= sizeof(T) * 4 );
    BOOST_REQUIRE( sizeof(complex_it<T, 3>) >= sizeof(T) * 8 );

    BOOST_WARN( not (complex_rt<T, 0>::has_padding) );
    BOOST_WARN( not (complex_rt<T, 1>::has_padding) );
    BOOST_WARN( not (complex_rt<T, 2>::has_padding) );
    BOOST_WARN( not (complex_rt<T, 3>::has_padding) );

    BOOST_WARN( not (complex_it<T, 0>::has_padding) );
    BOOST_WARN( not (complex_it<T, 1>::has_padding) );
    BOOST_WARN( not (complex_it<T, 2>::has_padding) );
    BOOST_WARN( not (complex_it<T, 3>::has_padding) );

#if BOOST_MATH_COMPLEX_IS_PACKED
    BOOST_REQUIRE_EQUAL( sizeof(complex_rt<T, 0>), sizeof(T) );
    BOOST_REQUIRE_EQUAL( sizeof(complex_rt<T, 1>), sizeof(T) * 2 );
    BOOST_REQUIRE_EQUAL( sizeof(complex_rt<T, 2>), sizeof(T) * 4 );
    BOOST_REQUIRE_EQUAL( sizeof(complex_rt<T, 3>), sizeof(T) * 8 );

    BOOST_REQUIRE_EQUAL( sizeof(complex_it<T, 0>), sizeof(T) );
    BOOST_REQUIRE_EQUAL( sizeof(complex_it<T, 1>), sizeof(T) * 2 );
    BOOST_REQUIRE_EQUAL( sizeof(complex_it<T, 2>), sizeof(T) * 4 );
    BOOST_REQUIRE_EQUAL( sizeof(complex_it<T, 3>), sizeof(T) * 8 );
#endif
}

// The complex-number class templates can substitute for std::complex.
BOOST_AUTO_TEST_CASE_TEMPLATE( complex_substitution_demo, T, test_types )
{
    BOOST_REQUIRE_EQUAL( complex_rt<T>::rank, 1u );
    BOOST_REQUIRE_EQUAL( complex_it<T>::rank, 1u );
}

BOOST_AUTO_TEST_SUITE_END()  // general_complex_tests
