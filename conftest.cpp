
#include "Config.h"
#include <iostream>

int main( int argc, char **argv )
{

	DispReln::GKSlab foo = DispReln::Config::ReadConfig<DispReln::GKSlab>( "DriftWaveConfig.xml" );
	std::cout << "Reference Beta is " << foo.beta_ref << std::endl;
	for ( auto &x : foo.SpeciesList )
	{
		std::cout << "Species with Temperature = " << x.s.Temperature << ", Density = " << x.s.Density << ", Charge = " << x.s.Z << "e, Mass = " << x.s.mass << " amu" << std::endl;
	}
	for ( auto &y : DispReln::Config::GenerateScans( "DriftWaveConfig.xml" ) )
	{
		switch ( y.parameter )
		{
			case DispReln::Config::ScanTypes::ky:
				std::cout << "Scanning in ky" << std::endl;
			default:
				break;
		}
		std::cout << "Holding ";
		for ( auto &f : y.fixed )
		{
			switch ( f.first )
			{
				case DispReln::Config::ScanTypes::kpar:
					std::cout << "k_|| = " << f.second << " ";
					break;
				case DispReln::Config::ScanTypes::ky:
					std::cout << "k_y = " << f.second << " ";
					break;
				case DispReln::Config::ScanTypes::kx:
					std::cout << "k_x = " << f.second << " ";
					break;
				default:
					break;
			}
		}
		std::cout << "fixed." << std::endl;
		std::cout << "Values to be scanned : ";
		for ( auto d : y.values )
			std::cout << d << " ";
		std::cout << std::endl;
	}
	return 0;
}
