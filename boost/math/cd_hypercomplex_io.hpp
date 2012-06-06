//  Boost.Complex Library, cd_hypercomplex_io.hpp header file  ---------------//

//  Copyright 2012 Daryle Walker.
//  Distributed under the Boost Software License, Version 1.0.  (See the
//  accompanying file LICENSE_1_0.txt or a copy at
//  <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

/** \file
    \brief  I/O operator templates for Cayley-Dickson hypercomplex numbers.

    \author  Daryle Walker

    \version  0.5

    \copyright  Boost Software License, version 1.0

    Contains the definitions of operator templates for performing standard
    input and output with objects of types from the `cdh_complex...` family.
 */

#ifndef BOOST_MATH_CD_HYPERCOMPLEX_IO_HPP
#define BOOST_MATH_CD_HYPERCOMPLEX_IO_HPP

#include "boost/math/cd_hypercomplex_core.hpp"

#include <algorithm>
#include <cstddef>
#include <ios>
#include <ostream>
#include <sstream>


namespace boost
{
namespace math
{


//  Implementation details  --------------------------------------------------//

//! \cond
namespace detail
{
    // Core output routine for `cdh_complex_ai`.
    template < typename Ch, class Tr, class Al, typename T, std::size_t R >
    void
    write_cdh_complex( std::basic_ostringstream<Ch, Tr, Al> &o,
     cdh_complex_ai<T, R> const &x, std::streamsize w )
    {
        using std::ios_base;

        auto const  L = 1UL << dynamic_rank( x );
        auto const  ww = ( w < (L + 1u) ) ? 0 : ( w - L - 1u );
        Ch          prefix = o.widen( '(' );

        std::for_each( x.begin(), x.begin() + L, [=,&prefix,&o] ( T const &v ) {
            o.put( prefix );
            o.width( ww / L );
            o << v;
            prefix = o.widen( ',' );
        } );
        o.put( o.widen(')') );
    }

    // Core output routine for `cdh_complex_ar`, base case.
    template < typename Ch, class Tr, class Al, typename T >
    void
    write_cdh_complex( std::basic_ostringstream<Ch, Tr, Al> &o,
     cdh_complex_ar<T, 0> const &x, std::streamsize w )
    {
        o.width( w );
        o << x.r[0];
    }

    // Core output routine for `cdh_complex_ar`, recursive case.
    template < typename Ch, class Tr, class Al, typename T, std::size_t R >
    void
    write_cdh_complex( std::basic_ostringstream<Ch, Tr, Al> &o,
     cdh_complex_ar<T, R> const &x, std::streamsize w )
    {
        auto const  ww = ( w < 3 ) ? 0 : ( w - 3 ) ;

        o.put( o.widen('(') );
        write_cdh_complex( o, x.b[0], ww / 2 );
        o.put( o.widen(',') );
        write_cdh_complex( o, x.b[1], ww / 2 );
        o.put( o.widen(')') );
    }

}  // namespace detail
//! \endcond


//  C.D.-hypercomplex output operator template definitions  ------------------//

/** \brief  Standard stream-insertion operator for `cdh_complex_ai`.

    Writes a hypercomplex object to the given output-stream.  It follows the
    suggested algorithm to print `std::complex` values from the C++11 standard,
    section 26.4.6, paragraphs 16 through 17.

    If the upper barrage of the value is zero, then only the lower barrage
    (applied recursively) is written.  The real part is always written, even if
    the entire value is zero.  If more than one part is written, then the parts
    will be surrounded by parentheses and comma-separated.

    \param[in,out] o  The stream to insert the text representation to.
    \param[in]     c  The number to have its text representation inserted.

    \throws  Anything stream-insertion of `T` objects might do, and/or general
             I/O actions.  Most of the time, error flags in `o` are set instead.

    \returns  `o`.
 */
template < typename Ch, class Tr, typename T, std::size_t R >
std::basic_ostream<Ch, Tr> &
operator <<( std::basic_ostream<Ch, Tr> &o, cdh_complex_ai<T, R> const &c )
{
    if( dynamic_rank(c) )
    {
        using std::ios_base;

        std::basic_ostringstream<Ch, Tr>  s;
        bool const          use_width = ( o.flags() & ios_base::adjustfield )
         == ios_base::internal;

        s.flags( o.flags() );
        s.imbue( o.getloc() );
        s.precision( o.precision() );
        s.fill( o.fill() );
        detail::write_cdh_complex( s, c, use_width ? o.width() : 0 );
        return o << s.str();

        // Extra space provided by the width setting is split among the
        // components.  Any remainder is applied to the whole text
        // representation.  The extra space will generally be partially
        // internal, seeming to violate a setting of "left" or "right" for the
        // "adjustfield."  Therefore, we turn off the extra space at the
        // component level except when "adjustfield" is "internal."
    }
    else
        return o << c.c[ 0 ];
}

/** \brief  Standard stream-insertion operator for `cdh_complex_ar`, zero-rank
            specialization.

    Since `cdh_complex_ar` is (partially) specialized when `rank` is zero, so
    is the text-writing function.  It applies the real part directly to the
    stream.

    \param[in,out] o  The stream to insert the text representation to.
    \param[in]     c  The number to have its text representation inserted.

    \throws  Anything stream-insertion of `T` objects might do, and/or general
             I/O actions.  Most of the time, error flags in `o` are set instead.

    \returns  `o`.
 */
template < typename Ch, class Tr, typename T >
std::basic_ostream<Ch, Tr> &
operator <<( std::basic_ostream<Ch, Tr> &o, cdh_complex_ar<T, 0> const &c )
{ return o << c.r[0]; }

/** \brief  Standard stream-insertion operator for `cdh_complex_ar`.

    Writes a hypercomplex object to the given output-stream.  It follows the
    suggested algorithm to print `std::complex` values from the C++11 standard,
    section 26.4.6, paragraphs 16 through 17.

    If the upper barrage of the value is zero, then only the lower barrage
    (applied recursively) is written.  The real part is always written, even if
    the entire value is zero.  If more than one part is written, then the two
    barrages are written (via recursive call) in order, surrounded by
    parentheses and comma-separated.

    \param[in,out] o  The stream to insert the text representation to.
    \param[in]     c  The number to have its text representation inserted.

    \throws  Anything stream-insertion of `T` objects might do, and/or general
             I/O actions.  Most of the time, error flags in `o` are set instead.

    \returns  `o`.
 */
template < typename Ch, class Tr, typename T, std::size_t R >
std::basic_ostream<Ch, Tr> &
operator <<( std::basic_ostream<Ch, Tr> &o, cdh_complex_ar<T, R> const &c )
{
    if( c.b[1] )
    {
        using std::ios_base;

        std::basic_ostringstream<Ch, Tr>  s;
        bool const          use_width = ( o.flags() & ios_base::adjustfield )
         == ios_base::internal;

        s.flags( o.flags() );
        s.imbue( o.getloc() );
        s.precision( o.precision() );
        s.fill( o.fill() );
        detail::write_cdh_complex( s, c, use_width ? o.width() : 0 );
        return o << s.str();

        // Extra space provided by the width setting is split among the
        // components.  Any remainder is applied to the whole text
        // representation.  The extra space will generally be partially
        // internal, seeming to violate a setting of "left" or "right" for the
        // "adjustfield."  Therefore, we turn off the extra space at the
        // component level except when "adjustfield" is "internal."
    }
    else
        return o << c.b[ 0 ];
}


}  // namespace math
}  // namespace boost


#endif  // BOOST_MATH_CD_HYPERCOMPLEX_IO_HPP
