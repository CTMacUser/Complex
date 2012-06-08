//  Boost.Complex Library, cd_hypercomplex compile-time-fail test file  ------//

//  Copyright 2012 Daryle Walker.  Use, modification, and distribution are
//  subject to the Boost Software License, Version 1.0.  (See accompanying file
//  LICENSE_1_0.txt or a copy at <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

#include "boost/math/cd_hypercomplex.hpp"

#include <cassert>


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

    return 0;
}
