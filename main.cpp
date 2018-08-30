
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <complex>
#include <utility>

#include <boost/math/constants/constants.hpp>
#include "RootFinder.h"
#include "DispReln.h"
#include "ParameterScan.h"
#include "Config.h"

using namespace RootFinder;


std::ostream& operator<<( std::ostream& os, std::pair< Complex, Real > p )
{
	os << p.first << " " << p.second;
	return os;
}

std::ostream& operator<<( std::ostream& os, std::list<std::pair< Complex, Real >> l )
{
	for ( auto x : l )
		os << x.second << "\t" << x.first.real() << "\t" << x.first.imag() << std::endl;
	return os;
}

std::ostream& operator<<( std::ostream& os, std::list<Complex> l )
{
	for ( auto x : l )
		os << x << std::endl;
	return os;
}

int main( int argc, char** argv )
{
	RootBoundingBox box( 0.0, 0.0, 0 );
	std::list<Complex> roots;

	/*
	// Masses in amu
	double electron_mass = 5.48579909e-4;
	double deuteron_mass = 2.01355321;
	double proton_mass   = 1.00727647;
	DispReln::Species Ions( 1.0, 1.0, 1.0, 1.0, 7.0, 0.0 ),Electrons( 1.0, 1.0, -1.0, electron_mass/deuteron_mass, 7.0, 0.0 );
	DispReln::ElectrostaticSlab DriftSlab( {Ions,Electrons} );
	*/

	DispReln::ElectrostaticSlab DriftSlab =  DispReln::Config::ReadConfig<DispReln::ElectrostaticSlab>( "tmp-test.xml" );

	// omega_star / k_|| v_th  = .5 * (k_y rho_i)/(k_|| L_n) = .5*(k_y * rho)*(fprim / kpar)
	// for v_th_i , we want om_star / k_|| v_th >> 1 @ k_perp rho_i ~ 0.5 - 1.5 ish.
	// Let fprim = 8, so xi_star_i ~= 2-6 ; xi_star_e ~= .033 - .099
	//
	// Drift-wave root ~ omega_star_e , which means xi ~ 1.4 - 4.2 
	
	auto scans = DispReln::Config::GenerateScans( "tmp-test.xml" );
	auto MainScan = scans.front();
	for ( auto &f : MainScan.fixed )
	{
		switch ( f.first )
		{
			case DispReln::Config::ScanTypes::kpar:
				DriftSlab.set_kpar( f.second );
				break;
			case DispReln::Config::ScanTypes::ky:
				DriftSlab.set_ky( f.second );
				break;
			case DispReln::Config::ScanTypes::kx:
				DriftSlab.set_kx( f.second );
				break;
			default:
				break;
		}
	}

	std::list< std::pair< std::complex<double>, double > scan;

	switch ( MainScan.parameter )
	{
		case DispReln::Config::ScanTypes::kpar:
			DispReln::kpar_scanner kparScan( DriftSlab );
			scan = MostUnstableModes( MainScan.box, kparScan, MainScan.values );
			break;
		case DispReln::Config::ScanTypes::ky:
			DispReln::ky_scanner kyScan( DriftSlab );
			scan = MostUnstableModes( MainScan.box, kyScan, MainScan.values );
			break;
		case DispReln::Config::ScanTypes::kx:
			DispReln::kx_scanner kxScan( DriftSlab );
			scan = MostUnstableModes( MainScan.box, kxScan, MainScan.values );
			break;
		default:
			break;
	}
	
	//{0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0}
	auto scan = MostUnstableModes( MainScan.box, kyScan, MainScan.values );

	std::cout << scan << std::endl;
	
	return 0;
}

