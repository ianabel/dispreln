
#define BOOST_TEST_MODULE gk_dispersion_tests
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/math/constants/constants.hpp>

#include "RootFinder.h"



BOOST_AUTO_TEST_SUITE( root_finding_test_suite, * boost::unit_test::tolerance( 1e-14 ) )

using namespace RootFinder;

BOOST_AUTO_TEST_CASE( winding_tests )
{
	Real M_2PI = boost::math::double_constants::two_pi;
	Path Circle = [=]( Real x ) { return std::complex<Real>( std::sin( M_2PI * x ), std::cos( M_2PI * x ) );};
	Path DoubleCircle = [=]( Real x ) { return std::complex<Real>( std::sin( M_2PI * x * 2 ), std::cos( M_2PI * x * 2) );};
	Path Squircle = [=]( Real x ) { return ( ( x - .5 )*( x-.5 ) + 0.5 )*std::complex<Real>( std::cos( M_2PI * x * 2 ), std::sin( M_2PI * x * 2) );};

	BOOST_TEST( WindingNumber( Circle ) == -1.0 );
	BOOST_TEST( WindingNumber( DoubleCircle ) == -2.0 );
	BOOST_TEST( WindingNumber( Squircle ) == 2.0 );

	Complex a( -1,-1 ),b( 1,1 );
	BOOST_TEST( WindingNumber( Rectangle( a, b ) ) == 1.0 );
	BOOST_TEST( WindingNumber( Rectangle( b, a ) ) == 1.0 );

	Complex c( 2,2 );
	BOOST_TEST( WindingNumber( Rectangle( b, c ) ) == 0.0 );
}

BOOST_AUTO_TEST_CASE( polynomial_root_tests, * boost::unit_test::tolerance( 1e-14 ) )
{
	// z^2 - 4i == 0, with roots +/- sqrt(2)*(1+i)
	Func P = []( Complex z ){ return z*z - std::complex<Real>( 0.0, 4.0 );};

	Complex a( 0.0, 0.0 ), b( 1.0, 1.0 ), c( 2.0, 2.0 ),d( -2.0, -2.0 );

	BOOST_TEST( WindingNumber( RectangleImage( a, b, P ) ) == 0 );
	BOOST_TEST( WindingNumber( RectangleImage( b, c, P ) ) == 1 );
	BOOST_TEST( WindingNumber( RectangleImage( d, c, P ) ) == 2 );

	std::list<Complex> roots = FindWithin( RootBoundingBox( b, c, 1 ), P, 1e-13 );
	Complex exact = ::sqrt( 2.0 )*std::complex<double>( 1.0, 1.0 );

	BOOST_TEST( roots.size() == 1 );
	BOOST_TEST( std::abs( roots.front() - exact ) == 0.0 );

	// And with a double root 
	Func Q = [=]( Complex z ){ return ( z-exact )*( z*z - std::complex<Real>( 0.0, 4.0 ) );};
	roots = FindWithin( RootBoundingBox( b, c, 1 ), Q, 1e-15 );

	BOOST_TEST( roots.size() == 1 );
	BOOST_TEST( std::abs( roots.front() - exact ) == 0.0 );
}

BOOST_AUTO_TEST_SUITE_END()

