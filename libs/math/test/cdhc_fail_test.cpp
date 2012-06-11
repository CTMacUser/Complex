//  Boost.Complex Library, cd_hypercomplex compile-time-fail test file  ------//

//  Copyright 2012 Daryle Walker.  Use, modification, and distribution are
//  subject to the Boost Software License, Version 1.0.  (See accompanying file
//  LICENSE_1_0.txt or a copy at <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

#include "boost/math/cd_hypercomplex.hpp"

#include <cassert>
#include <tuple>


// Pre-processor controls for each error type.  Change the control values within
// this file or use pre-#definitions from the command-line to override.  This
// enables doing subsets of all the tests.  (This helps when the compiler stops
// searching after finding some of the errors, so later ones get missed.)
#ifndef CONTROL_TEST_FAIL_AI_INDIRECT_BOOL_CONVERT
#define CONTROL_TEST_FAIL_AI_INDIRECT_BOOL_CONVERT  1
#endif
#ifndef CONTROL_TEST_FAIL_AI_DIRECT_BOOL_CONVERT
#define CONTROL_TEST_FAIL_AI_DIRECT_BOOL_CONVERT  1
#endif
#ifndef CONTROL_TEST_FAIL_AR_INDIRECT_BOOL_CONVERT
#define CONTROL_TEST_FAIL_AR_INDIRECT_BOOL_CONVERT  1
#endif
#ifndef CONTROL_TEST_FAIL_AR_DIRECT_BOOL_CONVERT
#define CONTROL_TEST_FAIL_AR_DIRECT_BOOL_CONVERT  1
#endif
#ifndef CONTROL_TEST_FAIL_AR_MULTILEVEL_BOOL_CONVERT
#define CONTROL_TEST_FAIL_AR_MULTILEVEL_BOOL_CONVERT  1
#endif
#ifndef CONTROL_TEST_FAIL_AI_OP_EQUALS
#define CONTROL_TEST_FAIL_AI_OP_EQUALS  1
#endif
#ifndef CONTROL_TEST_FAIL_AI_OP_NONEQUAL
#define CONTROL_TEST_FAIL_AI_OP_NONEQUAL  1
#endif
#ifndef CONTROL_TEST_FAIL_AR_OP_EQUALS_SAME_LENGTH
#define CONTROL_TEST_FAIL_AR_OP_EQUALS_SAME_LENGTH  1
#endif
#ifndef CONTROL_TEST_FAIL_AR_OP_NONEQUAL_SAME_LENGTH
#define CONTROL_TEST_FAIL_AR_OP_NONEQUAL_SAME_LENGTH  1
#endif
#ifndef CONTROL_TEST_FAIL_AR_OP_NONEQUAL_S2L
#define CONTROL_TEST_FAIL_AR_OP_NONEQUAL_S2L  1
#endif
#ifndef CONTROL_TEST_FAIL_AR_OP_EQUALS_S2L
#define CONTROL_TEST_FAIL_AR_OP_EQUALS_S2L  1
#endif
#ifndef CONTROL_TEST_FAIL_AR_OP_NONEQUAL_L2S
#define CONTROL_TEST_FAIL_AR_OP_NONEQUAL_L2S  1
#endif
#ifndef CONTROL_TEST_FAIL_AR_OP_EQUALS_L2S
#define CONTROL_TEST_FAIL_AR_OP_EQUALS_L2S  1
#endif
#ifndef CONTROL_TEST_FAIL_AI_TUPLE_ELEMENT
#define CONTROL_TEST_FAIL_AI_TUPLE_ELEMENT  1
#endif
#ifndef CONTROL_TEST_FAIL_AR_TUPLE_ELEMENT
#define CONTROL_TEST_FAIL_AR_TUPLE_ELEMENT  1
#endif
#ifndef CONTROL_TEST_FAIL_AI_GET
#define CONTROL_TEST_FAIL_AI_GET  1
#endif
#ifndef CONTROL_TEST_FAIL_AR_BASE_GET
#define CONTROL_TEST_FAIL_AR_BASE_GET  1
#endif
#ifndef CONTROL_TEST_FAIL_AR_MULTILEVEL_GET
#define CONTROL_TEST_FAIL_AR_MULTILEVEL_GET  1
#endif
#ifndef CONTROL_TEST_FAIL_AI_CONVERSIONS
#define CONTROL_TEST_FAIL_AI_CONVERSIONS  1
#endif
#ifndef CONTROL_TEST_FAIL_AR_CONVERSIONS
#define CONTROL_TEST_FAIL_AR_CONVERSIONS  1
#endif


// Compact type expressions
using boost::math::cdh_complex_ai;
using boost::math::cdh_complex_ar;


// Classes without some sort of feature
template < typename T >
struct no_boolean_convert_t
{
    T  v;
};

template < typename T >
struct no_implicit_boolean_convert_t
{
    T  v;
    explicit operator bool() const  { return static_cast<bool>(v); }
};


// Use "assert," but make sure implicit Boolean conversion is done
bool  delayed_assert( bool x )  { assert(x); return x; }


// Main program
int  main()
{
    // Sample values
    no_implicit_boolean_convert_t<int>  v1{ 0 }, v2{ 5 };
    no_boolean_convert_t<int>           w1{ 0 }, w2{ -4 };

    // Check if these compile
    cdh_complex_ai<no_implicit_boolean_convert_t<int>, 0>  test1{ {v1} };
    cdh_complex_ai<no_implicit_boolean_convert_t<int>, 0>  test2{ {v2} };
    cdh_complex_ai<no_boolean_convert_t<int>, 0>           test3{ {w1} };
    cdh_complex_ai<no_boolean_convert_t<int>, 0>           test4{ {w2} };
    cdh_complex_ar<no_implicit_boolean_convert_t<int>, 0>  test5{ {v1} };
    cdh_complex_ar<no_implicit_boolean_convert_t<int>, 0>  test6{ {v2} };
    cdh_complex_ar<no_boolean_convert_t<int>, 0>           test7{ {w1} };
    cdh_complex_ar<no_boolean_convert_t<int>, 0>           test8{ {w2} };
    cdh_complex_ar<no_boolean_convert_t<int>, 1> test9{ {test7, test8} };

    assert( !v1 );          // OK: op! can use explicit Boolean conversion
#if 0
    //assert( v2 );         // workability depends on "assert"'s implmentation
    delayed_assert( v2 );   // no implicit Boolean conversion
    assert( !w1 );          // no Boolean conversion
    assert( w2 );           // no Boolean conversion
#endif
    assert( !test1 );
    assert( static_cast<bool>(test2) );
#if CONTROL_TEST_FAIL_AI_INDIRECT_BOOL_CONVERT
    assert( !test3 );
#endif
#if CONTROL_TEST_FAIL_AI_DIRECT_BOOL_CONVERT
    assert( static_cast<bool>(test4) );
#endif
    assert( !test5 );
    assert( static_cast<bool>(test6) );
#if CONTROL_TEST_FAIL_AR_INDIRECT_BOOL_CONVERT
    assert( !test7 );
#endif
#if CONTROL_TEST_FAIL_AR_DIRECT_BOOL_CONVERT
    assert( static_cast<bool>(test8) );
#endif
#if CONTROL_TEST_FAIL_AR_MULTILEVEL_BOOL_CONVERT
    assert( static_cast<bool>(test9) );
#endif

    // The previously used types don't have equality operators either....
#if CONTROL_TEST_FAIL_AI_OP_EQUALS
    assert( !(test1 == test2) );
#endif
#if CONTROL_TEST_FAIL_AI_OP_NONEQUAL
    assert( test1 != test2 );
#endif

#if CONTROL_TEST_FAIL_AR_OP_EQUALS_SAME_LENGTH
    assert( test7 == test7 );
    assert( !(test7 == test8) );
#endif

#if CONTROL_TEST_FAIL_AR_OP_NONEQUAL_SAME_LENGTH
    assert( !(test7 != test7) );
    assert( test7 != test8 );
#endif

#if CONTROL_TEST_FAIL_AR_OP_NONEQUAL_S2L
    assert( test7 != test9 );
#endif
#if CONTROL_TEST_FAIL_AR_OP_EQUALS_S2L
    assert( !(test7 == test9) );
#endif
#if CONTROL_TEST_FAIL_AR_OP_NONEQUAL_L2S
    assert( test9 != test7 );
#endif
#if CONTROL_TEST_FAIL_AR_OP_EQUALS_L2S
    assert( !(test9 == test7) );
#endif

    // Errors happen if the tuple_element index meets or exceeds the value
    // from tuple_size.
    using std::tuple_element;

#if CONTROL_TEST_FAIL_AI_TUPLE_ELEMENT
    typename tuple_element<1, cdh_complex_ai<no_boolean_convert_t<int>,
     0>>::type  test10;
    typename tuple_element<2, cdh_complex_ai<no_boolean_convert_t<int>,
     1>>::type  test11;
    typename tuple_element<4, cdh_complex_ai<no_boolean_convert_t<int>,
     2>>::type  test12;
    typename tuple_element<8, cdh_complex_ai<no_boolean_convert_t<int>,
     3>>::type  test13;

    assert( sizeof(test10) >= sizeof(int) );
    assert( sizeof(test11) >= sizeof(int) );
    assert( sizeof(test12) >= sizeof(int) );
    assert( sizeof(test13) >= sizeof(int) );
#endif

#if CONTROL_TEST_FAIL_AR_TUPLE_ELEMENT
    typename tuple_element<3, cdh_complex_ar<no_boolean_convert_t<int>,
     0>>::type  test14;
    typename tuple_element<5, cdh_complex_ar<no_boolean_convert_t<int>,
     1>>::type  test15;
    typename tuple_element<7, cdh_complex_ar<no_boolean_convert_t<int>,
     2>>::type  test16;
    typename tuple_element<9, cdh_complex_ar<no_boolean_convert_t<int>,
     3>>::type  test17;

    assert( sizeof(test14) >= sizeof(int) );
    assert( sizeof(test15) >= sizeof(int) );
    assert( sizeof(test16) >= sizeof(int) );
    assert( sizeof(test17) >= sizeof(int) );
#endif

    // The index-out-of-bounds errors can happen with "get" too.
    using boost::math::get;

#if CONTROL_TEST_FAIL_AI_GET
    cdh_complex_ai<no_boolean_convert_t<int>, 0>          test18;
    cdh_complex_ai<no_boolean_convert_t<int>, 1>          test19;
    cdh_complex_ai<no_boolean_convert_t<int>, 2>          test20;
    cdh_complex_ai<no_boolean_convert_t<int>, 3>          test21;
    cdh_complex_ai<no_boolean_convert_t<int>, 0> const &  ctest18 = test18;
    cdh_complex_ai<no_boolean_convert_t<int>, 1> const &  ctest19 = test19;
    cdh_complex_ai<no_boolean_convert_t<int>, 2> const &  ctest20 = test20;
    cdh_complex_ai<no_boolean_convert_t<int>, 3> const &  ctest21 = test21;

    assert( sizeof(get<3>( test18 )) >= sizeof(int) );
    assert( sizeof(get<5>( test19 )) >= sizeof(int) );
    assert( sizeof(get<7>( test20 )) >= sizeof(int) );
    assert( sizeof(get<9>( test21 )) >= sizeof(int) );

    assert( sizeof(get<4>( ctest18 )) >= sizeof(int) );
    assert( sizeof(get<6>( ctest19 )) >= sizeof(int) );
    assert( sizeof(get<9>( ctest20 )) >= sizeof(int) );
    assert( sizeof(get<11>( ctest21 )) >= sizeof(int) );

    assert( sizeof(get<5>( decltype(test18){} )) >= sizeof(int) );
    assert( sizeof(get<7>( decltype(test19){} )) >= sizeof(int) );
    assert( sizeof(get<11>( decltype(test20){} )) >= sizeof(int) );
    assert( sizeof(get<13>( decltype(test21){} )) >= sizeof(int) );
#endif

#if CONTROL_TEST_FAIL_AR_BASE_GET
    cdh_complex_ar<no_boolean_convert_t<int>, 0>          test22;
    cdh_complex_ar<no_boolean_convert_t<int>, 0> const &  ctest22 = test22;

    assert( sizeof(get<1>( test22 )) >= sizeof(int) );
    assert( sizeof(get<2>( ctest22 )) >= sizeof(int) );
    assert( sizeof(get<3>( decltype(test22){} )) >= sizeof(int) );
#endif

#if CONTROL_TEST_FAIL_AR_MULTILEVEL_GET
    cdh_complex_ar<no_boolean_convert_t<int>, 1>          test23;
    cdh_complex_ar<no_boolean_convert_t<int>, 2>          test24;
    cdh_complex_ar<no_boolean_convert_t<int>, 3>          test25;
    cdh_complex_ar<no_boolean_convert_t<int>, 1> const &  ctest23 = test23;
    cdh_complex_ar<no_boolean_convert_t<int>, 2> const &  ctest24 = test24;
    cdh_complex_ar<no_boolean_convert_t<int>, 3> const &  ctest25 = test25;

    assert( sizeof(get<2>( test23 )) >= sizeof(int) );
    assert( sizeof(get<3>( ctest23 )) >= sizeof(int) );
    assert( sizeof(get<4>( decltype(test23){} )) >= sizeof(int) );

    assert( sizeof(get<4>( test24 )) >= sizeof(int) );
    assert( sizeof(get<5>( ctest24 )) >= sizeof(int) );
    assert( sizeof(get<6>( decltype(test24){} )) >= sizeof(int) );

    assert( sizeof(get<8>( test25 )) >= sizeof(int) );
    assert( sizeof(get<9>( ctest25 )) >= sizeof(int) );
    assert( sizeof(get<10>( decltype(test25){} )) >= sizeof(int) );
#endif

    // Complex conversions can only be done if the element types are
    // convertible.
#if CONTROL_TEST_FAIL_AI_CONVERSIONS
    cdh_complex_ai<no_boolean_convert_t<int>, 1>  test26{ {w1, w2} };
    auto const  test27 =
     static_cast<cdh_complex_ai<no_implicit_boolean_convert_t<int>, 0>>( test26
     );
    auto const  test28 =
     static_cast<cdh_complex_ai<no_implicit_boolean_convert_t<int>, 1>>( test26
     );
    auto const  test29 =
     static_cast<cdh_complex_ai<no_implicit_boolean_convert_t<int>, 2>>( test26
     );

    assert( test27.c[0].v == 0 );
    assert( test28.c[0].v == 0 );
    assert( test28.c[1].v == -4 );
    assert( test29.c[0].v == 0 );
    assert( test29.c[1].v == -4 );
    assert( test29.c[2].v == 0 );
    assert( test29.c[3].v == 0 );
#endif

#if CONTROL_TEST_FAIL_AR_CONVERSIONS
    auto const  test30 =
     static_cast<cdh_complex_ar<no_implicit_boolean_convert_t<int>, 0>>(test9);
    auto const  test31 =
     static_cast<cdh_complex_ar<no_implicit_boolean_convert_t<int>, 1>>(test9);
    auto const  test32 =
     static_cast<cdh_complex_ar<no_implicit_boolean_convert_t<int>, 2>>(test9);

    assert( test30.r[0].v == test7.r[0].v );
    assert( test31.b[0].r[0].v == test7.r[0].v );
    assert( test31.b[1].r[0].v == test8.r[0].v );
    assert( test32.b[0].b[0].r[0].v == test7.r[0].v );
    assert( test32.b[0].b[1].r[0].v == test8.r[0].v );
    assert( !test32.b[1] );
#endif

    return 0;
}
