//  Boost.Complex Library, Cayley-Dickson hypercomplex example program file  -//

//  Copyright 2012 Daryle Walker.  Use, modification, and distribution are
//  subject to the Boost Software License, Version 1.0.  (See accompanying file
//  LICENSE_1_0.txt or a copy at <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

#include "boost/math/cd_hypercomplex.hpp"  // for boost::math::cdh_complex...

#include <cassert>   // for assert
#include <ctime>     // for std::time
#include <iostream>  // for std::cout
#include <ostream>   // for std::endl
#include <random>    // for std::default_random_engine, uniform_int_distribution


// The types we'll we working with...
using boost::math::cdh_complex_ai;
using boost::math::cdh_complex_ar;

template < typename T >  using real_ai_t = cdh_complex_ai<T, 0>;
template < typename U >  using real_ar_t = cdh_complex_ar<U, 0>;
template < typename T >  using complex_ai_t = cdh_complex_ai<T, 1>;
template < typename U >  using complex_ar_t = cdh_complex_ar<U, 1>;
template < typename T >  using quaternion_ai_t = cdh_complex_ai<T, 2>;
template < typename U >  using quaternion_ar_t = cdh_complex_ar<U, 2>;
template < typename T >  using octonion_ai_t = cdh_complex_ai<T, 3>;
template < typename U >  using octonion_ar_t = cdh_complex_ar<U, 3>;

typedef real_ai_t<int>        real_aii_t;
typedef real_ar_t<int>        real_ari_t;
typedef complex_ai_t<int>     complex_aii_t;
typedef complex_ar_t<int>     complex_ari_t;
typedef quaternion_ai_t<int>  quaternion_aii_t;
typedef quaternion_ar_t<int>  quaternion_ari_t;
typedef octonion_ai_t<int>    octonion_aii_t;
typedef octonion_ar_t<int>    octonion_ari_t;


// Main program
int  main()
{
    using std::cout;
    using std::endl;
    using boost::math::dynamic_rank;

    // Initialization demonstration
    real_aii_t        t1 = { {1} }, t1a = { 1 };
    real_ari_t        t2 = { {2} }, t2a = { 2 };
    complex_aii_t     t3 = { {3, 4} }, t3a = { 3, 4 };
    complex_ari_t     t4 = { {{ {5} }, { {6} }} }, t4a = { 5, 6 };
    quaternion_aii_t  t5 = { {7, 8, 9, 10} }, t5a = { 7, 8, 9, 10 };
    quaternion_ari_t  t6 = { {{ {{ {11} }, { {12} }} }, { {{ {13} }, { {14} }}
                           }} },
                     t6a = { 11, 12, 13, 14 };
    octonion_aii_t    t7 = { {15, 16, 17, 18, 19, 20, 21, 22} },
                     t7a = { 15, 16, 17, 18, 19, 20, 21, 22 };
    octonion_ari_t    t8 = { {{ {{ {{ {23} }, { {24} }} }, { {{ {25} }, { {26}
                           }} }} }, { {{ {{ {27} }, { {28} }} }, { {{ {29} }, {
                           {30} }} }} }} },
                     t8a = { 23, 24, 25, 26, 27, 28, 29, 30 };

    real_aii_t const &  ct1 = t1;

    // Show "size"
    for ( unsigned  i = 0u ; i < ct1.size() ; ++i )
        cout << ct1.c[ i ];
    cout << endl;
    assert( t2.size() == (1u << 0) );
    assert( t4a.size() == (1u << 1) );

    // Show "begin" and "end" via range-for
    std::default_random_engine          re;
    std::uniform_int_distribution<int>  di{ -100, 100 };
    octonion_aii_t                      demo1;
    octonion_aii_t const &              c_demo1 = demo1;
    real_ari_t                          demo2;
    real_ari_t const &                  c_demo2 = demo2;

    re.seed( (unsigned) std::time(nullptr) );

    for ( auto &x : demo1 )
        x = di( re );
    for ( auto const &x : c_demo1 )
        cout << x << ' ';
    cout << endl;

    for ( auto &x : demo2 )
        x = di( re );
    for ( auto const &x : c_demo2 )
        cout << x << ' ';
    cout << endl;

    // Show "dynamic_rank"
    assert( dynamic_rank(t1) == 0 );
    assert( dynamic_rank(t1a) == 0 );
    assert( dynamic_rank(t2) == 0 );
    assert( dynamic_rank(t2a) == 0 );
    assert( dynamic_rank(t3) == 1 );
    assert( dynamic_rank(t3a) == 1 );
    assert( dynamic_rank(t4) == 1 );
    assert( dynamic_rank(t4a) == 1 );
    assert( dynamic_rank(t5) == 2 );
    assert( dynamic_rank(t5a) == 2 );
    assert( dynamic_rank(t6) == 2 );
    assert( dynamic_rank(t6a) == 2 );
    assert( dynamic_rank(t7) == 3 );
    assert( dynamic_rank(t7a) == 3 );
    assert( dynamic_rank(t8) == 3 );
    assert( dynamic_rank(t8a) == 3 );

    return 0;
}
