
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <complex>
#include <utility>

#include <boost/program_options.hpp>
#include <boost/math/constants/constants.hpp>
#include "RootFinder.h"
#include "DispReln.h"
#include "ParameterScan.h"
#include "Config.h"

using Complex = std::complex<double>;
using Real = double;

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

std::string Header( DispReln::Config::Scan Scan )
{
	std::ostringstream out( "" );
	std::string col_header;
	switch ( Scan.parameter )
	{
		case DispReln::Config::ScanTypes::kpar:
			out << "# Scanning in k_|| ";
			col_header = "# k_||\tomega\tgamma";
			break;
		case DispReln::Config::ScanTypes::ky:
			out << "# Scanning in k_y ";
			col_header = "# k_y\tomega\tgamma";
			break;
		case DispReln::Config::ScanTypes::kx:
			out << "# Scanning in k_x ";
			col_header = "# k_x\tomega\tgamma";
			break;
		default:
			break;
	}
	switch ( Scan.mode )
	{
		case DispReln::Config::ScanMode::TrackRoot:
			out << "and tracking a single root of the dispersion relation" << std::endl;
			break;
		case DispReln::Config::ScanMode::MostUnstableMode:
			out << "and tracking the most unstable mode" << std::endl;
			break;
	}

	out << col_header << std::endl;
	return out.str();
}

template<typename T> std::list<std::pair<Complex,Real> > DoScan( DispReln::Config::Scan MainScan, T scanner )
{
	switch ( MainScan.mode )
	{
		case DispReln::Config::ScanMode::TrackRoot:
			return TrackRoot( MainScan.box, scanner, MainScan.values );
			break;
		case DispReln::Config::ScanMode::MostUnstableMode:
			return MostUnstableModes( MainScan.box, scanner, MainScan.values );
			break;
	}
	return std::list<std::pair<Complex,Real> >();
}

template<typename T> std::list<std::pair<Complex,Real> > PerformScan( DispReln::Config::Scan MainScan, T physics )
{
	for ( auto &f : MainScan.fixed )
	{
		switch ( f.first )
		{
			case DispReln::Config::ScanTypes::kpar:
				physics.set_kpar( f.second );
				break;
			case DispReln::Config::ScanTypes::ky:
				physics.set_ky( f.second );
				break;
			case DispReln::Config::ScanTypes::kx:
				physics.set_kx( f.second );
				break;
			default:
				break;
		}
	}

	std::list< std::pair< Complex, Real > > scan;

	DispReln::ky_scanner kyScan( physics );
	DispReln::kpar_scanner kparScan( physics );
	DispReln::kx_scanner kxScan( physics );
	DispReln::fprim_scanner fprimScan( physics, MainScan.sIndex );
	DispReln::tprim_scanner tprimScan( physics, MainScan.sIndex );

	switch ( MainScan.parameter )
	{
		case DispReln::Config::ScanTypes::kpar:
			scan = DoScan( MainScan, kparScan );
			break;
		case DispReln::Config::ScanTypes::ky:
			scan = DoScan( MainScan, kyScan );
			break;
		case DispReln::Config::ScanTypes::kx:
			scan = DoScan( MainScan, kxScan );
			break;
		case DispReln::Config::ScanTypes::fprim:
			scan = DoScan( MainScan, fprimScan );
			break;
		case DispReln::Config::ScanTypes::tprim:
			scan = DoScan( MainScan, tprimScan );
			break;

		default:
			scan.clear();
			break;
	}
	return scan;
}

template<typename T> int RunAllScans( std::list<DispReln::Config::Scan> ScanList, T physics )
{
	for ( auto &S : ScanList )
	{
		std::list< std::pair<Complex,Real> > scan_results;
		scan_results = PerformScan( S, physics );
		std::cout << Header( S ) << std::endl;
		std::cout << scan_results << std::endl;
	}
	return 0;
}




namespace po = boost::program_options;

int main( int argc, char** argv )
{
	po::options_description desc( "Allowed Options" );
	desc.add_options()( "help", "Produce this help message" )( "config-file,c", po::value<std::string>(), "Configuration file" );

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm); 
	
	std::string config_file( "DispSolver.xml" );
	if ( vm.count( "config-file" ) )
	{
		config_file = vm[ "config-file" ].as<std::string>();
	}


	DispReln::Config::PhysicsTypes physics = DispReln::Config::GetPhysicsModel( config_file );
	auto scans = DispReln::Config::GenerateScans( config_file );
	DispReln::GKSlab FullSlab;
	DispReln::ElectrostaticSlab ESSlab;

	switch ( physics )
	{
		case DispReln::Config::PhysicsTypes::GKSlab:
			FullSlab =  DispReln::Config::ReadConfig<DispReln::GKSlab>( config_file );
			std::cout << "# Solutions to the full gyrokinetic dispersion relation in a slab" << std::endl;
			return RunAllScans( scans, FullSlab );
			break;
		case DispReln::Config::PhysicsTypes::ElectrostaticSlab:
			ESSlab =  DispReln::Config::ReadConfig<DispReln::ElectrostaticSlab>( config_file );
			std::cout << "# Solutions to the electrostatic gyrokinetic dispersion relation in a slab" << std::endl;
			return RunAllScans( scans, ESSlab );
			break;
		default:
			std::cerr << "Config file specifies an invalid physics model" << std::endl;
			return -1;
			break;
	}
	
	return -2;
}


