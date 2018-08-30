
#include "DispReln.h"

#include "Faddeeva.hh"


namespace DispReln {

	std::complex<double> Z( std::complex<double> xi )
	{
		return std::complex<double>( 0.0, 1.0 )*( boost::math::double_constants::root_pi )*Faddeeva::w( xi );
	}

	std::complex<double> AcousticDisp( std::complex<double> xi, double ZoverTau )
	{
		return 1.0 + ZoverTau + xi*Z( xi );
	}

	std::complex<double> FLRAcousticDisp( std::complex<double> xi, double ZoverTau, double kperp )
	{
		return 1.0 + ZoverTau + xi*Z( xi )*Gamma<1,0,0>( kperp*kperp/2.0 );
	}	

}

