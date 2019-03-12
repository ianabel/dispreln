#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <complex>
#include <utility>
#include <deque>
#include <list>

#include <boost/program_options.hpp>
#include <boost/math/constants/constants.hpp>
#include "RootFinder.h"
#include "DispReln.h"
#include "ParameterScan.h"
#include "Config.h"

using Complex = std::complex<double>;
using Real = double;
using RootPair = std::pair< Complex, Real >;
using RootList_t = std::deque<RootPair>;

std::ostream& operator<<( std::ostream& os, RootPair p )
{
	os << p.first << " " << p.second;
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
	
	DispReln::Species ions( 1.0, 1.0, 1.0, 1.0, 0.0, 20.0 );
	DispReln::Species electrons( 1.0, 1.0, -1.0, 0.0, 0.0, 0.0 );
	
	DispReln::ElectrostaticSlab ESSlab{ ions, electrons };

	std::cout << std::setprecision( 10 );

	ESSlab.set_kpar( 1.0 );
	ESSlab.set_ky( 1.05 );
	ESSlab.set_kx( 0.0  );

	std::complex<double> z(-1.974275369,0.5043842457 );

	std::cout << " A( " << z << " ) == " << ESSlab( z ) << std::endl;

	std::cout << " Zeta<0>( " << z << " ) == " << DispReln::Zeta<0>( z ) << std::endl;
	std::cout << " Gamma<1,0,0>( " << 1.05 << " ) == " << DispReln::Gamma<1,0,0>( 1.05*1.05/2 ) << std::endl;

	return 0;
	
}


