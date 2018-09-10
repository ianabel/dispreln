
#include "RootFinder.h"

/*
 * These routines take a bounding box or simplex that is guaranteed to contain one root
 * ( possibly of high multiplicity roots are handled explicitly ) and then solve down to
 * the required tolerance.
 *
 * This should backend off of an algorithm for finding and bracketing the roots efficiently.
 * 
 * We do not use the available boost algorithms as they do not have derivative-free complex
 * root solvers.
 */

namespace RootFinder {
Complex DirectSolve( RootBoundingBox box, Func const & f, Real tol )
{
	// If box contains multiple roots, you are doing it all wrong
	if ( box.Index == 0 )
		throw std::invalid_argument( "No roots in box" );
	if ( box.Index < 0 )
		throw std::invalid_argument( "Box contains poles. Aborting." );
	// We will assume there is *one* root, possibly with a high multiplicity

	// We just do a secant solve from the centre, 

	Complex z = box.centre(),u;
	Complex fprimez = ( f( box.upper ) - f( box.lower ) )/( box.upper - box.lower );
	Complex fz = f( z );

	// Shouldn't need more than 10 steps.
	unsigned int MAX_ITER = 40;
	unsigned int i;
	for ( i=0; i<MAX_ITER; i++ )
	{
		// Compute next point
		u = z - fz/fprimez;
		// Check for convergence
		// if( |u-z| < tol * |z| ) we're done
		if ( std::abs( u - z ) < tol*std::abs( z ) )
			break;
		// Compute an approximation to f'(u)
		Complex fu = f( u );
		fprimez = ( fu - fz )/( u - z );

		// move from z -> u
		z = u;
		fz = fu;
	}

	if ( i == MAX_ITER )
	{
		// Something probably went wrong, this shouldn't really happen.
		throw std::logic_error( "Maximum secant iterations exceeded" );
	}

	return u;
}


Complex DirectSolve( Simplex T, unsigned int Index, Func const & F, Real tol )
{
	Complex r[ 3 ];
	Complex f[ 3 ],fc;
	for ( unsigned int i=0; i<3; ++i )
	{
		r[ i ] = T[ i ] - T.centre();
		f[ i ] = F( T[ i ] );
	}
	fc = F( T.centre() );
	Complex fprime_c;
	Complex lambda[ 3 ];

	lambda[ 0 ] = 1.0;
	lambda[ 1 ] = ( r[ 0 ]*r[ 0 ] )*( r[ 0 ] - r[ 2 ] )/( r[ 1 ]*r[ 1 ] * ( r[ 2 ] - r[ 1 ] ) );
	lambda[ 2 ] = ( r[ 0 ]*r[ 0 ] )*( r[ 1 ] - r[ 0 ] )/( r[ 2 ]*r[ 2 ] * ( r[ 2 ] - r[ 1 ] ) );

	Complex l_sum = lambda[ 0 ] + lambda[ 1 ] + lambda[ 2 ];
	Complex lr_sum = r[ 0 ]*lambda[ 0 ] + r[ 1 ]*lambda[ 1 ] + r[ 2 ]*lambda[ 2 ];
	Complex lf_sum = f[ 0 ]*lambda[ 0 ] + f[ 1 ]*lambda[ 1 ] + f[ 2 ]*lambda[ 2 ];

	fprime_c = ( lf_sum - l_sum*fc )/lr_sum;
	// fprime2_c = f''(c)/2 
	Complex fprime2_c = ( f[ 0 ] + f[ 1 ] + f[ 2 ] - 3.*fc )/( r[ 0 ]*r[ 0 ] + r[ 1 ]*r[ 1 ] + r[ 2 ]*r[ 2 ] );

	Complex fn2 = fc,xn2 = T.centre();
	Complex xn1 = T.centre() - ( fc/fprime_c )*( 1.0 + ( fc/fprime_c )*fprime2_c );
	Complex fn1 = F( xn1 );

	// f(x_i) = f(u) + (x_i - u)*f'(u) + (x_i-u)^2 f''(u)/2 + ...
	// sum_i f_i = 3 f(u) + 3( c - u) f'(u) +   [ sum_i ( x_i - u )^2 ] f''(u)/2
	// f(c) = f(u) + (c - u)*f'(u) + (c-u)^2 f''(u)/2
	// f''(u)/2 = (sum_i f_i - 3*f(c)) / ( [ sum_i ( x_i - u )^2 ] - (c - u)^2 )
	
	for ( unsigned int i=0; i<3; ++i )
		r[ i ] = T[ i ] - xn1;
	Complex r2sum = r[ 0 ]*r[ 0 ] + r[ 1 ]*r[ 1 ] + r[ 2 ]*r[ 2 ];
	Complex fprime2_u = ( f[ 0 ] + f[ 1 ] + f[ 2 ]  - 3.*fc )/( r2sum  - ( fc/fprime_c )*( fc/fprime_c ) );
	Complex fprime_u = ( fc - fn1 - r2sum*fprime2_u )/( fc/ fprime_c );

	Complex xn = xn1 - ( fn1/ fprime_u )*( 1.0 + ( fn1/fprime_u )*fprime2_u );

	Complex fn = F( xn );
	

	unsigned int MAX_ITER = 40;
	Complex u = 0.0;
	for ( unsigned int i=0; i<MAX_ITER; i++ )
	{
		// fn2 =  f(x_(n-2)) ; fn1 = f(x_(n-1)) ; fn = f(x_n)
		// xn2 = x_(n-2) etc
		// Do inverse quadratic interpolation
		u = fn1*fn*xn2/( ( fn2 - fn1 )*( fn2 - fn ) ) + fn2*fn*xn1/( ( fn2 - fn1 )*( fn - fn1 ) ) + fn2*fn1*xn/( ( fn - fn2 )*( fn - fn1 ) );
		if( std::abs( u - xn2 ) < tol*std::abs( u ) )
			return u;

		fn2 = fn1;
		xn2 = xn1;
		fn1 = fn;
		xn1 = xn;
		xn = u;
		fn = F( u );
	}
	return u;
}




}
