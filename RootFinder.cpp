
#include "RootFinder.h"

namespace RootFinder {
/*
 * Computes the winding number of a curve C around the origin.
 * As usual, this is counted counter-clockwise.
 */
int WindingNumber( Path C, unsigned int N )
{
	// Check C is in fact a closed curve on [0,1]
	Complex z0,z1;
	Real eps = 1e-6;
	z0 = C( 0.0 );
	z1 = C( 1.0 );

	if ( std::abs( z0 - z1 ) > eps )
	{
		throw std::invalid_argument( "Path is not closed." );
	}

	// Divide [0,1] into N pieces as a first guess
	Real increment = 1.0/N;



	// Now compute winding number as # of crossings of negative real axis
	
	int Windings = 0;

	for ( Real s = 0.0;  s < 1.0; )
	{
		Real r = s + increment;

		Complex Cs = C( s ),Cr = C( r );

		Real delta = ::abs( Cr - Cs ) / ::abs( Cr );
		if ( delta > 0.25  )
		{
			increment *= .5;
			continue;
		}
		else if ( delta < .005 )
		{
			increment *= 1.5;
			r = s + increment;
			Cr = C( r );
		}

		if ( std::abs( Cs ) < 1e-100 )
		{
			throw std::invalid_argument( "Curve gets too close to zero" );
		}
		if ( Cs.real() < 0.0 && Cr.real() < 0.0 )
		{
			if ( Cs.imag() <= 0.0 && Cr.imag() > 0.0 )
				--Windings;
			else if ( Cs.imag() > 0.0 && Cr.imag() <= 0.0 ) 
				++Windings;
		}
		s = r;
	}
	return Windings;
}

Path Rectangle( Complex a, Complex b )
{
	return [ = ]( Real s ) {
		if ( 0.0 < s && s <= 0.25 )
			return std::complex<Real>( a.real() + 4.0*s*( b.real() - a.real() ), a.imag() );
		else if ( 0.25 < s && s <= 0.5 )
			return std::complex<Real>( b.real(), a.imag() + ( 4.0*s - 1.0 )*( b.imag()-a.imag() ) );
		else if ( 0.5 < s && s <= 0.75 )
			return std::complex<Real>( b.real() - ( 4.0*s - 2.0 )*( b.real() - a.real() ), b.imag() );
		else if ( 0.75 < s && s <= 1.0 )
			return std::complex<Real>( a.real(), b.imag() - ( 4.0*s - 3.0 )*( b.imag()-a.imag() ) );
		else
			return a;
	};
}

Path RectangleImage( Complex a, Complex b, Func const & f )
{
	// generate a path by taking f( Rectangle( s ) )
	return [ = ]( Real s ) {
		if ( 0.0 < s && s <= 0.25 )
			return f( std::complex<Real>( a.real() + 4.0*s*( b.real() - a.real() ), a.imag() ) );
		else if ( 0.25 < s && s <= 0.5 )
			return f( std::complex<Real>( b.real(), a.imag() + ( 4.0*s - 1.0 )*( b.imag()-a.imag() ) ) );
		else if ( 0.5 < s && s <= 0.75 )
			return f( std::complex<Real>( b.real() - ( 4.0*s - 2.0 )*( b.real() - a.real() ), b.imag() ) );
		else if ( 0.75 < s && s <= 1.0 )
			return f( std::complex<Real>( a.real(), b.imag() - ( 4.0*s - 3.0 )*( b.imag()-a.imag() ) ) );
		else
			return f( a );
	};
}

Path RectangleImage( RootBoundingBox const & b, Func const& f )
{
	return RectangleImage( b.lower, b.upper, f );
}

std::list< RootBoundingBox > Refine( RootBoundingBox outerBox, Func const & f, unsigned int N )
{
	std::list<RootBoundingBox> results,fine_grid;
	Complex outer_a = outerBox.lower;
	Complex outer_b = outerBox.upper;
	unsigned int toFind = outerBox.Index;
	Complex d1( ( outer_b.real() - outer_a.real() )/N, 0.0 ),d2( 0.0, ( outer_b.imag() - outer_a.imag() )/N );
	for ( unsigned int i=0; i < N; i++ )
	{
		for ( unsigned int j=0; j < N; j++ )
		{
			// i + 0.0 promotes i to a double to multiply the complex number
			fine_grid.emplace_back( outer_a + ( i + 0.0 )*d1 + ( j + 0.0 )*d2, outer_a + ( i+1.0 )*d1 + ( j + 1.0 )*d2, -1 );
		}
	}

	for ( auto& box : fine_grid )
	{
		box.Index = WindingNumber( RectangleImage( box.lower, box.upper, f ) );
		if ( box.Index > 0 )
		{
			results.emplace_back( box.lower, box.upper, box.Index );
			toFind -= box.Index;
			if ( toFind == 0 )
				break;
		}
	}
	return results;
}

bool Resolved( RootBoundingBox b, Real tol )
{
	return ( ( std::abs( b.upper - b.lower )/std::abs( b.upper + b.lower ) ) < tol );
}


std::list< RootBoundingBox > RefineAll( std::list<RootBoundingBox> &Boxes, Func const & f, Real tol )
{
	for ( auto it = Boxes.begin() ; it != Boxes.end() ; )
	{
		if ( Resolved( *it, tol ) )
			++it;
		else
		{
			auto subdivisions = Refine( *it, f );
			auto refined = RefineAll( subdivisions, f, tol );
			Boxes.splice( it, refined );
			it = Boxes.erase( it );
		}
	}
	return Boxes;
}

std::list<Complex> FindWithin( RootBoundingBox box, Func const& f, Real tol )
{
	std::list<RootBoundingBox> grid;
	box.Index = WindingNumber( RectangleImage( box, f ) );

	grid.push_back( box );

	// Bisect to get most of the work done
	RefineAll( grid, f, ::sqrt( tol ) );

	std::list<Complex> roots;

	for ( auto& patch : grid )
	{
		roots.push_back( DirectSolve( patch, f, tol ) );
	}

	return roots;
}

}
