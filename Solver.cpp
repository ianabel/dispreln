
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

struct Simplex {
	Complex vertices[ 3 ];
	Simplex( Complex a, Complex b, Complex c )
	{
		vertices[ 0 ] = a;
		vertices[ 1 ] = b;
		vertices[ 2 ] = c;
	}
	Simplex( Simplex const& other )
	{
		for ( unsigned int i=0; i<3; i++ )
			vertices[ i ] = other.vertices[ i ];
	}
	Complex centre() const
	{
		return ( vertices[ 0 ] + vertices[ 1 ] + vertices[ 2 ] )/3.0;
	}
};

Complex DirectSolve( Simplex region, unsigned int Index, Func const & f, Real tol )
{
	throw std::logic_error( "Unimplemented." );
	Complex u( 0.0, 0.0 );
	return u;
}
}
