#ifndef ROOTFINDER_H
#define ROOTFINDER_H

#include <complex>
#include <functional>

#include <list>
#include <iostream>
#include <boost/math/constants/constants.hpp>


namespace RootFinder {
	using Real = double;
	using Complex = std::complex<Real>;
	static const Complex I( 0.0, 1.0 );

	using Path = std::function< Complex( Real )>;
	using Func = std::function< Complex( Complex )>;

	int WindingNumber( Path, unsigned int N=400 );
	Path Rectangle( Complex lower, Complex upper );
	Path RectangleImage( Complex a, Complex b, Func const & f );

	struct RootBoundingBox {
		Complex lower;
		Complex upper; 
		int Index;

		RootBoundingBox() : lower( 0.0 ), upper( 0.0 ), Index( 0 ) {};
		RootBoundingBox( RootBoundingBox const& rb ) : lower( rb.lower ), upper( rb.upper ), Index( rb.Index ) {};
		RootBoundingBox( Complex a, Complex b, unsigned int i ) : lower( a ), upper( b ), Index( i ) {};
		friend std::ostream& operator<<( std::ostream& os, const RootBoundingBox& Box ) {
			os << "[" << Box.lower << "," << Box.upper << "]";
			return os;
		};

		Complex centre() const { return ( lower + upper )/2.0; };
		Complex diag() const { return ( lower - upper )/2.0; };


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

	template<typename T> std::list< std::pair<Complex,Real> > TrackRoot( RootBoundingBox Initial, T const& G, std::list<Real> params, Real tol = 1e-3 )
	{
		std::list< std::pair<Complex,Real> > RootList;
		RootBoundingBox bounds = Initial;
		static unsigned int MAX_REFINE = 20;
		static Real RefineFac = 0.75;
		static Real ExpandFac = 1.25;
		RootList.clear();

		for ( auto alpha : params )
		{
			T F( G );
			F.set_param( alpha );

			int N_ROOT = -1;
			// Move box to be centred on a guess at the next root
			if ( RootList.size() >= 2 ) {
				/*
				// Linear extrapolation
				auto old_soln = RootList.rbegin();
				auto older_soln = --RootList.rbegin();
				Complex old_root = old_soln->first;
				Complex older_root = older_soln->first;
				Real old_alpha = old_soln->second;
				Real older_alpha = older_soln->second;
				Complex rootGuess = old_root + ( alpha - old_alpha )*( old_root - older_root )/( old_alpha - older_alpha );
				bounds.recentre( rootGuess );
				*/
				bounds.recentre( RootList.rbegin()->first );
				if ( abs( bounds.diag() ) > abs( Initial.diag() )*2.0 )
					bounds.scale( 0.5 );
			} else if ( RootList.size() == 1 ) {
				// Only have one data point, so just centre on that
				bounds.recentre( RootList.rbegin()->first );
			} else {
				bounds = Initial;
			}
				
			// std::cerr << "New alpha is " << alpha << " and we are now looking in " << bounds; 

			for ( unsigned i=0; N_ROOT != 1 && i < MAX_REFINE; i++ )
			{
				N_ROOT = WindingNumber( RectangleImage( bounds, F ) );
				if ( N_ROOT == 0 )
					bounds.scale( ExpandFac );
				else if ( N_ROOT > 1 )
					bounds.scale( RefineFac );
				else if ( N_ROOT < 0 )
					throw std::logic_error( "Cannot find roots of non-holomorphic functions" );
				else if ( N_ROOT == 1 )
					break;
			}

			if ( N_ROOT == 1 )
			{
				auto roots = FindWithin( bounds, F, tol );
				Complex root = roots.front();
				RootList.emplace_back( root, alpha );
			}
			else
			{
				std::cerr << "Old root was " << RootList.rbegin()->first << " at alpha = " << RootList.rbegin()->second << ";" << std::endl;
				std::cerr << "New alpha is " << alpha << " and we are now looking in " << bounds; 
				throw std::logic_error( "Cannot adjust box to only contain new root -- parameter step is likely too large" );
			}

			// std::cerr << " where we found the root " << RootList.back().first << std::endl;
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

	template<typename T> std::list< std::pair< Complex, Real > > MostUnstableModes( RootBoundingBox Initial, T const& G, std::list<Real> parameters, Real tol = 1e-3 )
	{
		RootBoundingBox box = Initial;
		Real centre_in_omega = box.centre().real();
		std::list< std::pair< Complex, Real > > RootList;

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

}
#endif // ROOTFINDER_H
