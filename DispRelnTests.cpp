
#define BOOST_TEST_MODULE gk_dispersion_tests
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <algorithm>
#include <complex>
#include <vector>
#include <utility>

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/math/constants/constants.hpp>

#include "DispReln.h"
#include "RootFinder.h"
#include "Config.h"



BOOST_AUTO_TEST_SUITE( special_function_test_suite, * boost::unit_test::tolerance( 1e-14 ) )

static const std::complex<double> I( 0.0, 1.0 );

BOOST_AUTO_TEST_CASE( z_function_test )
{
	// Set of test points computed with Mathematica
	std::list< std::pair< std::complex<double>, std::complex<double> > > Z_Function_Values{ 
		{0.0, I*boost::math::double_constants::root_pi},
		{1.0, -1.076159013825537 + 0.6520493321732922 * I}
	};
	for ( auto x : Z_Function_Values ) 
		BOOST_TEST( std::abs( DispReln::Z( x.first ) - x.second  ) == 0.0 );
}

BOOST_AUTO_TEST_CASE( gamma_function_test )
{
	// Arbitary set of test points
	std::list< double > Test_Values{ 0.0, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6};
	// Check we're returning things we expect
	for ( auto x : Test_Values ) 
	{
		double g1 = std::exp( -x ) * std::cyl_bessel_i( 0.0, x );
		double g2 = std::exp( -x ) * ( std::cyl_bessel_i( 0.0, x ) - std::cyl_bessel_i( 1.0, x ) );
		double g3 = std::exp( -x ) * ( std::cyl_bessel_i( 0.0, x ) + x * ( std::cyl_bessel_i( 1.0, x ) - std::cyl_bessel_i( 0.0, x ) ) );
		double g1_a = DispReln::Gamma<1,0,0>( x );
		double g2_a = DispReln::Gamma<2,1,0>( x );
		double g3_a = DispReln::Gamma<3,0,0>( x );
		BOOST_TEST( g1 == g1_a );
		BOOST_TEST( g2 == g2_a );
		BOOST_TEST( g3 == g3_a );
	}

	// k rho ~ = 45 => alpha ~= 1000. so check we're good up to alpha of 5000 (and ky rho of 100)
	std::vector< double > LargeValues{ 250.0, 500.0, 750.0, 1000.0, 5000.0 };
	// ExpBessel<0> of these values, as determined by mathematica 
	std::vector< double > EB0Vals{0.025243969387054753633,0.017845706500153167237,0.014569742116743979078,0.012617240455891256586,0.0056420368987445886570};
	std::vector< double > EB1Vals{0.025193430757117305262,0.017827851852898056461,0.014560025713286366714,0.012610930256928629470,0.0056414726668388859036};
	for ( unsigned int i=0; i<LargeValues.size(); ++i )
	{
		double x = LargeValues[ i ];
		double eb0 = DispReln::ExpBessel<0>( x );
		double eb1 = DispReln::ExpBessel<1>( x );
		
		BOOST_TEST( eb0 == EB0Vals[ i ] );
		BOOST_TEST( eb1 == EB1Vals[ i ] );
	}
}

BOOST_AUTO_TEST_CASE( zeta_function_test )
{
	// Arbitary set of test points
	std::list< std::complex<double> > Test_Values{ 0.0, 0.1 + 0.2*I, 0.2 + 0.3*I, 0.4 + 0.4*I, 0.8 + 0.5*I, 1.6 + 0.6*I, 3.2 - I, 6.4 - 2.0*I, 12.8 - 3.0*I, 25.6 - 5.0*I};
	// Check we're returning things we expect
	for ( auto x : Test_Values ) 
	{
		std::complex<double> Z = DispReln::Z( x );
		std::complex<double> Z0 = Z - DispReln::Zeta<0>( x );
		std::complex<double> Z1 = ( 1.0 + x*Z ) - DispReln::Zeta<1>( x );
		std::complex<double> Z2 = x*( 1.0 + x*Z ) - DispReln::Zeta<2>( x );
		std::complex<double> Z3 = 0.5 + x*x*( 1.0 + x*Z ) - DispReln::Zeta<3>( x );
		BOOST_TEST( std::abs( Z0 ) == 0.0 );
		BOOST_TEST( std::abs( Z1 ) == 0.0 );
		BOOST_TEST( std::abs( Z2 ) == 0.0 );
		BOOST_TEST( std::abs( Z3 ) == 0.0 );
	}
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( disp_reln_solver_test, * boost::unit_test::tolerance( 1e-11 ) )

using namespace RootFinder;
static const std::complex<double> I( 0.0, 1.0 );

BOOST_AUTO_TEST_CASE( acoustic_wave_test )
{
	RootBoundingBox FirstTwo( 0.0 - 2.0*I, 10.0, 0.0 );
	Func DispFn = std::bind( DispReln::AcousticDisp, std::placeholders::_1, 1.0 );
	std::list<Complex> roots = FindWithin< std::list< Complex > >( FirstTwo, DispFn, 1e-11 );
	roots.sort( []( const Complex &x, const Complex &y ){ return ( x.imag() > y.imag() );} );

	BOOST_TEST( roots.size() == 2 );
	Complex r1 = 1.446673204192393 - 0.6019815403718067*I;
	Complex r2 = 2.358299486175531 - 1.790624079933201*I;

	BOOST_TEST( std::abs( roots.front() - r1 ) == 0.0 );
	BOOST_TEST( std::abs( roots.back() - r2 ) == 0.0 );


	// Z/tau = 0.01, so cold ions, which gives rise to the undamped ion acoustic mode
	DispFn = std::bind( DispReln::AcousticDisp, std::placeholders::_1, 0.01 );
	roots = FindWithin< std::list< Complex > >( FirstTwo, DispFn, 1e-11 );
	roots.sort( []( const Complex &x, const Complex &y ){ return ( x.imag() > y.imag() );} );

	BOOST_TEST( roots.size() == 2 );
	r1  = 7.178530545559893 - 6.189514943424657e-14*I; // Basically undamped 
	r2  = 2.5729233681035124- 1.2398216368657875*I;
	BOOST_TEST( std::abs( roots.front() - r1 ) == 0.0 );
	BOOST_TEST( std::abs( roots.back() - r2 ) == 0.0 );
}

BOOST_AUTO_TEST_CASE( acoustic_tracking_test )
{
	RootBoundingBox box( 0.0, 0.0, 0 );

	std::list<Real> tau_scan{ 1.0, 0.8, 0.6, 0.4, 0.2, 0.1 };
	
	box.lower = std::complex<double>( 1.0, -1.0 );
	box.upper = std::complex<double>( 2.0, 0.0 );

	struct AcousticDR_t {
		double Z_over_tau;
		Complex operator()( Complex xi ) { return DispReln::AcousticDisp( xi, Z_over_tau );};
		AcousticDR_t( AcousticDR_t const& other ) : Z_over_tau( other.Z_over_tau ) {};
		AcousticDR_t() : Z_over_tau( 0.0 ) {};
		void set_param( Real x ) { Z_over_tau = x; };
	} AcousticDR;

	auto scan = TrackRoot( box, AcousticDR, tau_scan, 1e-8 );
	
	Complex answers[] = { 1.4466732041923935 - I*0.6019815403718068, 1.515528614457207 - 0.5370210294765353*I, 1.6140317173734795 - 0.4539510799509855*I, 1.7749127224445853 - 0.3399231863870913*I, 2.1291713454169376 - 0.164215945462315*I, 2.6366843740563994 - 0.0412507862131464*I};
	int j=0;
	for ( auto i = scan.begin(); i != scan.end(); ++i )
	{

		Complex answer = answers[ j ];
		BOOST_TEST( std::abs( i->first - answer ) == 0.0 );
		j++;
	}

	box.lower = std::complex<double>( 2.0, -2.0 );
	box.upper = std::complex<double>( 3.0, -1.0 );

	scan = TrackRoot(  box, AcousticDR, tau_scan, 1e-8 );

	Complex answers2[] = {2.3582994861755306 - 1.7906240799332012*I,2.3862789863171483 - 1.7639099893198387*I,2.4224152344165897 - 1.7292115194660806*I,2.4731713608236823 - 1.679218893327795*I,2.557186342255047 - 1.587013942433162*I,2.627188266700257 - 1.4768051228022558*I};
	j=0;
		
	for ( auto i = scan.begin(); i != scan.end(); ++i )
	{
		Complex answer = answers2[ j ];
		BOOST_TEST( std::abs( i->first - answer ) == 0.0 );
		j++;
	}

}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( config_tests )

BOOST_AUTO_TEST_CASE( example1_test )
{
	static std::string config_file( "test-data/DriftWaveConfig.xml" );
	DispReln::ElectrostaticSlab DriftSlab =  DispReln::Config::ReadConfig<DispReln::ElectrostaticSlab>( config_file );

	BOOST_TEST( DriftSlab.SpeciesList.size() == 2 );
	DispReln::Species ions = DriftSlab.SpeciesList.front().s;
	DispReln::Species electrons = DriftSlab.SpeciesList.back().s;

	BOOST_TEST( ions.Temperature == 1.0 );
	BOOST_TEST( ions.Density == 1.0 );
	BOOST_TEST( ions.Z == 1.0 );
	BOOST_TEST( ions.mass == 1.0 );
	BOOST_TEST( ions.fprim == 10.0 );
	BOOST_TEST( ions.tprim == 0.0 );

	BOOST_TEST( electrons.Temperature == 1.0 );
	BOOST_TEST( electrons.Density == 1.0 );
	BOOST_TEST( electrons.Z == -1.0 );
	BOOST_TEST( electrons.mass == 0.0002724 );
	BOOST_TEST( electrons.fprim == 10.0 );
	BOOST_TEST( electrons.tprim == 0.0 );

	auto scans = DispReln::Config::GenerateScans( config_file );
	BOOST_TEST( scans.size() == 1 );
	auto MainScan = scans.front();
	BOOST_TEST( MainScan.fixed.size() == 2 );
	BOOST_TEST( ( int )MainScan.fixed.front().first == ( int )DispReln::Config::ScanTypes::kpar );
	BOOST_TEST( ( int )MainScan.fixed.back().first == ( int )DispReln::Config::ScanTypes::kx );

	BOOST_TEST( MainScan.fixed.front().second == 1.0 );
	BOOST_TEST( MainScan.fixed.back().second == 0.0 );

	BOOST_TEST( ( int )MainScan.parameter == ( int )DispReln::Config::ScanTypes::ky );
		
	
	std::list<double> RefValues{0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0};

	BOOST_TEST( MainScan.values == RefValues, boost::test_tools::per_element() );

	std::complex<double> l( -5.0, -2.0 ),u( 5.0, 2.0 );

	BOOST_TEST( MainScan.box.lower == l );
	BOOST_TEST( MainScan.box.upper == u );

}

BOOST_AUTO_TEST_SUITE_END()
