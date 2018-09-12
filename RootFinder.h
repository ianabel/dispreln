#ifndef ROOTFINDER_H
#define ROOTFINDER_H

#include <complex>
#include <functional>
#include <vector>

#include <deque>
#include <list>
#include <iostream>
#include <boost/math/constants/constants.hpp>
#include <algorithm>


namespace RootFinder {
	using Real = double;
	using Complex = std::complex<Real>;
	static const Complex I( 0.0, 1.0 );
	struct Simplex;

	using Path = std::function< Complex( Real )>;
	using Func = std::function< Complex( Complex )>;

	int WindingNumber( Path, unsigned int N=400 );
	Path Rectangle( Complex lower, Complex upper );
	Path RectangleImage( Complex a, Complex b, Func const & f );
	Path Image( Simplex const &T, Func const & f );

	struct RootBoundingBox {
		Complex lower;
		Complex upper; 
		int Index;

		RootBoundingBox() : lower( 0.0 ), upper( 0.0 ), Index( 0 ) {};
		RootBoundingBox( RootBoundingBox const& rb ) : lower( rb.lower ), upper( rb.upper ), Index( rb.Index ) {};
		RootBoundingBox( Complex a, Complex b, unsigned int i ) : Index( i ) {
			if ( a.real() < b.real() )
			{
				lower.real(a.real());
				upper.real(b.real());

			}
			else
			{
				lower.real(b.real());
				upper.real(a.real());
			}

			if ( a.imag() < b.imag() )
			{
				lower.imag(a.imag());
				upper.imag(b.imag());

			}
			else
			{
				lower.imag(b.imag());
				upper.imag(a.imag());
			}
		};
		friend std::ostream& operator<<( std::ostream& os, const RootBoundingBox& Box ) {
			os << "[" << Box.lower << "," << Box.upper << "]";
			return os;
		};

		Complex centre() const { return ( lower + upper )/2.0; };
		Complex diag() const { return ( lower - upper )/2.0; };

		bool contains( Complex x ) {
			return ( x.real() > lower.real() && x.real() < upper.real() ) && ( x.imag() > lower.imag() && x.imag() < upper.imag() );
		}

		void scale( Real fac )
		{
			// Keep Centre fixed, and multiply dimensions by fac
			Complex new_diag = diag() * fac;
			Complex _centre = centre();
			upper = ( _centre + new_diag );
			lower = ( _centre - new_diag );
		}

		void recentre( Complex x )
		{
			// Move centre to x
			Complex _diag = diag();
			upper = ( x + _diag );
			lower = ( x - _diag );
		}
	};


	void scale( RootBoundingBox &box, Real fac );
	void recentre( RootBoundingBox &box, Complex x );

	Path RectangleImage( RootBoundingBox const& b, Func const & f );

	std::list< RootBoundingBox > Refine( RootBoundingBox outerBox, Func const & f, unsigned int N=2 );
	std::list< RootBoundingBox > RefineAll( std::list<RootBoundingBox> & Boxes, Func const & f, Real tol = 0.01 );

	Complex DirectSolve( RootBoundingBox box, Func const & f, Real tol = 1e-10 );
	std::list<Complex> FindWithin( RootBoundingBox box, Func const& f, Real tol = 1e-6 );

	template<typename T> std::deque< std::pair<Complex,Real> > TrackRoot( RootBoundingBox Initial, T const& G, std::list<Real> params, Real tol = 1e-3 )
	{
		std::deque< std::pair<Complex,Real> > RootList;
		RootBoundingBox bounds = Initial;
		static unsigned int MAX_REFINE = 20;
		static Real ExpandFac = 1.25;
		Real size_lim = 2.0*std::abs( Initial.diag() );
		RootList.clear();

		for ( auto alpha : params )
		{
			T F( G );
			F.set_param( alpha );

			int N_ROOT = -1;
			// Move box to be centred on the last root.
			if ( RootList.size() > 0 )
			{
				bounds.recentre( RootList.rbegin()->first );
				if ( abs( bounds.diag() ) > size_lim )
					bounds.scale( 0.5 );
			}

			Complex delta_guess,root_guess;
			if ( RootList.size() > 2 )
			{
				auto rn1 = RootList[ RootList.size() - 1 ];
				auto rn2 = RootList[ RootList.size() - 2 ];

				delta_guess = ( alpha - rn1.second )*( rn1.first - rn2.first )/( rn1.second - rn2.second );
				// recentre?
				root_guess =bounds.centre() + delta_guess;
				bounds.recentre( root_guess );
			}
				
			// std::cerr << "New alpha is " << alpha << " and we are now looking in " << bounds; 

			unsigned int n_points = 400;
			for ( unsigned i=0; N_ROOT != 1 && i < MAX_REFINE; i++ )
			{
				N_ROOT = WindingNumber( RectangleImage( bounds, F ), n_points );
				if ( N_ROOT == 0 )
					bounds.scale( ExpandFac );
				else if ( N_ROOT < 0 )
				{
					n_points *= 2;
					if ( n_points > 16000 )
						break;
				}
				else if ( N_ROOT >= 1 )
					break;
			}
			
			if ( N_ROOT <= 0 )
			{
				// throw std::logic_error( "Could not expand the box enough to find a root. Has it disappeared?" );
				std::cerr << "Scan terminated prematurely at alpha = " << alpha << std::endl;
				break;
			}


			auto roots = FindWithin( bounds, F, tol );

			roots.sort( [=]( Complex const& a, Complex const& b ){ return ( std::abs( a - root_guess )  < std::abs( b - root_guess ) ); } );

			Complex root = roots.front();

			if ( std::abs( N_ROOT ) != roots.size() )
			{
				std::cerr << "Possible root crossing at alpha = " << alpha << " and root location " << root << std::endl;
			}

			RootList.emplace_back( root, alpha );
		}

		return RootList;
	}

	bool inline MostUnstable( Complex const& a, Complex const& b )
	{
		if ( a.imag() >= b.imag() )
			return true;
		else
			return false;
	}

	template<typename T> std::deque< std::pair< Complex, Real > > MostUnstableModes( RootBoundingBox Initial, T const& G, std::list<Real> parameters, Real tol = 1e-3 )
	{
		RootBoundingBox box = Initial;
		Real centre_in_omega = box.centre().real();
		std::deque< std::pair< Complex, Real > > RootList;

		for ( auto alpha : parameters ) 
		{
			// Copy the dispersion object
			T DispersionObject( G );
			// Adjust it
			DispersionObject.set_param( alpha );
			auto roots = FindWithin( box, DispersionObject, 1e-3 );
			if ( roots.size() == 0 )
				throw std::logic_error( "Bad Box" );
			roots.sort( MostUnstable );
			Complex root = roots.front();
			RootList.emplace_back( root, alpha );
			box.recentre( centre_in_omega + root.imag() * I );
		}
		return RootList;
	}

	struct Simplex {
		Simplex( Complex a, Complex b, Complex c ) {
			vertices.reserve( 3 );
			vertices[ 0 ] = a;
			vertices[ 1 ] = b;
			vertices[ 2 ] = c;
			Complex bc = ( a + b + c )/3.0;
			std::sort( vertices.begin(), vertices.end(), [=]( Complex x, Complex y ){ return ( std::arg( x - bc ) < std::arg( y - bc ) );} );
		}
		Simplex( Simplex const& other )
		{
			for ( unsigned int i=0; i<3; i++ )
				vertices[ i ] = other.vertices[ i ];
		}

		Simplex Extend( Simplex, Complex );
		Complex centre() { return ( vertices[ 0 ] + vertices[ 1 ] + vertices[ 2 ] )/ 3.0;}
		bool inside( Complex x ) {
			bool Side1,Side2,Side3;
			Side1 = ( ( ( x - vertices[ 0 ] )*( vertices[ 1 ] - vertices[ 0 ] ) ).real() >=0 );
			Side2 = ( ( ( x - vertices[ 0 ] )*( vertices[ 1 ] - vertices[ 0 ] ) ).real() >=0 );
			Side3 = ( ( ( x - vertices[ 0 ] )*( vertices[ 1 ] - vertices[ 0 ] ) ).real() >=0 );
			return Side1 && Side2 && Side3;
		}

		Complex operator[]( unsigned int i ) const
		{
			return vertices[ i ];
		}

		Complex& operator[]( unsigned int i )
		{
			return vertices[ i ];
		}
		private:
			std::vector<Complex> vertices;
	};


}
#endif // ROOTFINDER_H
