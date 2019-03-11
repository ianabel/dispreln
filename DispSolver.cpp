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

std::string Header( DispReln::Config::Scan const& Scan, bool want_col_header = true )
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
		case DispReln::Config::ScanTypes::fprim:
			out << "# Scanning in fprim of species " << Scan.sIndex << " ";
			col_header = "# L_ref/L_n\tomega\tgamma";
			break;
		case DispReln::Config::ScanTypes::fprim_adj:
			out << "# Scanning in fprim of species " << Scan.sIndex << " and adjusting fprim of species " << Scan.adjIndex << " to maintain quasineutrality";
			out << std::endl << "# ";
			col_header = "# L_ref/L_n\tomega\tgamma";
			break;
		default:
			break;
	}
	switch ( Scan.mode )
	{
		case DispReln::Config::ScanMode::AllRoots:
			out << "and finding all roots of the dispersion relation inside the region " << Scan.box << std::endl;
			break;
		case DispReln::Config::ScanMode::TrackRoot:
			out << "and tracking a single root of the dispersion relation" << std::endl;
			break;
		case DispReln::Config::ScanMode::MostUnstableMode:
			out << "and tracking the most unstable mode" << std::endl;
			break;
		default:
			break;
	}

	out << "# The scan will be performed holding the following values fixed:" << std::endl;
	std::map< DispReln::Config::ScanTypes, std::string > outfixed{
		{ DispReln::Config::ScanTypes::kpar, "k_||" },
		{ DispReln::Config::ScanTypes::kx, "k_x" },
		{ DispReln::Config::ScanTypes::ky, "k_y" },
		{ DispReln::Config::ScanTypes::beta, "beta" },
		{ DispReln::Config::ScanTypes::fprim, "L_ref/L_n" },
		{ DispReln::Config::ScanTypes::tprim, "L_ref/L_T" },
	};

	for ( auto &x : Scan.fixed )
		out << "# " << outfixed[ x.first ] << " = " << x.second << std::endl;



	if ( want_col_header ) 
		out << col_header << std::endl;
	return out.str();
}

std::string Footer( DispReln::Config::Scan const& Scan )
{
	return "\n\n";
}

template<typename T> std::deque<RootPair> DoScan( DispReln::Config::Scan MainScan, T scanner )
{
	switch ( MainScan.mode )
	{
		case DispReln::Config::ScanMode::TrackRoot:
			return TrackRoot( MainScan.box, scanner, MainScan.values, MainScan.tolerance );
			break;
		case DispReln::Config::ScanMode::MostUnstableMode:
			return MostUnstableModes( MainScan.box, scanner, MainScan.values, MainScan.tolerance );
			break;
		case DispReln::Config::ScanMode::AllRoots:
			return AllRoots( MainScan.box, scanner, MainScan.values, MainScan.tolerance );
			break;
	}
	return std::deque<RootPair>();
}

template<typename T> std::deque<RootPair> PerformScan( DispReln::Config::Scan MainScan, T physics )
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
			case DispReln::Config::ScanTypes::beta:
				physics.set_beta( f.second );
				break;
			default:
				break;
		}
	}

	std::deque<RootPair> scan;

	DispReln::ky_scanner<T> kyScan( physics );
	DispReln::kpar_scanner<T> kparScan( physics );
	DispReln::kx_scanner<T> kxScan( physics );
	DispReln::fprim_scanner<T> fprimScan( physics, MainScan.sIndex );
	DispReln::tprim_scanner<T> tprimScan( physics, MainScan.sIndex );
	DispReln::fprim_adj_scanner<T> fprimAdjScan( physics, MainScan.sIndex, MainScan.adjIndex );


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
		case DispReln::Config::ScanTypes::fprim_adj:
			scan = DoScan( MainScan, fprimAdjScan );
			break;

		default:
			scan.clear();
			break;
	}
	return scan;
}

template<typename T> void OutputScan( DispReln::Config::Scan const& scan, T physics, RootList_t scan_results );
template<> void OutputScan<DispReln::ElectrostaticSlab>( DispReln::Config::Scan const& scan, DispReln::ElectrostaticSlab physics, RootList_t scan_results )
{
	std::cout << Header( scan ) << std::endl;

	double norm;

	switch ( scan.normalization )
	{
		case DispReln::Config::Normalization::Default:
			norm = 1.0;
			break;
		case DispReln::Config::Normalization::kRef:
			norm = physics.get_kpar();
			break;
		default:
			std::cerr << ( unsigned int )scan.normalization << std::endl;
			throw std::invalid_argument( "Attempting to use an Alfvenic normalization for an electrostatic run!" );
			break;
	}
	for ( auto &x : scan_results )
	{
		std::cout << x.second << "\t" << x.first.real()/norm << "\t" << x.first.imag()/norm << std::endl;
	}

	std::cout << Footer( scan ) << std::endl;

}

template<> void OutputScan<DispReln::GKSlab>( DispReln::Config::Scan const& scan, DispReln::GKSlab physics, RootList_t scan_results )
{
	std::cout << Header( scan ) << std::endl;

	double norm,beta;
	norm = 1.0;

	switch ( scan.normalization )
	{
		case DispReln::Config::Normalization::Default:
			break;
		case DispReln::Config::Normalization::kRef:
			norm = physics.get_kpar();
			break;
		case DispReln::Config::Normalization::kAlfven:
			norm = physics.get_kpar();
		case DispReln::Config::Normalization::Alfven:
			beta = physics.beta_ref;
			if ( scan.beta != 0.0 )
				beta = scan.beta;
			norm /= ::sqrt( beta );
			break;
		default:
			throw std::invalid_argument( "Unsupported Normalization" );
			break;
	}
	for ( auto &x : scan_results )
	{
		std::cout << x.second << "\t" << x.first.real()/norm << "\t" << x.first.imag()/norm << std::endl;
	}

	std::cout << Footer( scan ) << std::endl;

}

template<typename T> int RunAllScans( std::list<DispReln::Config::Scan> ScanList, T physics )
{
	for ( auto &S : ScanList )
	{
		RootList_t scan_results;
		scan_results = PerformScan( S, physics );
		OutputScan( S, physics, scan_results );
	}
	return 0;
}


template<typename T> void ListScan( DispReln::Config::Scan const& scan, T const& physics )
{
	std::cout << Header( scan ) << std::endl;


	switch ( scan.normalization )
	{
		case DispReln::Config::Normalization::Default:
			std::cout << "Frequencies and growth rates will be normalized to v_tref / L_ref" << std::endl;
			break;
		case DispReln::Config::Normalization::kRef:
			std::cout << "Frequencies and growth rates will be normalized to k_|| v_tref" << std::endl;
			break;
		case DispReln::Config::Normalization::kAlfven:
			std::cout << "Frequencies and growth rates will be normalized to k_|| v_A" << std::endl;
		case DispReln::Config::Normalization::Alfven:
			std::cout << "Frequencies and growth rates will be normalized to v_A / L_ref" << std::endl;
			break;
		default:
			throw std::invalid_argument( "Unsupported Normalization" );
			break;
	}

	std::cout << "The scan will be performed holding the following values fixed:" << std::endl;
	std::map< DispReln::Config::ScanTypes, std::string > outfixed{
		{ DispReln::Config::ScanTypes::kpar, "k_||" },
		{ DispReln::Config::ScanTypes::kx, "k_x" },
		{ DispReln::Config::ScanTypes::ky, "k_y" },
		{ DispReln::Config::ScanTypes::beta, "beta" },
		{ DispReln::Config::ScanTypes::fprim, "L_ref/L_n" },
		{ DispReln::Config::ScanTypes::tprim, "L_ref/L_T" },
	};

	for ( auto &x : scan.fixed )
		std::cout << outfixed[ x.first ] << " = " << x.second << std::endl;

	std::cout << "The scan will range over the following values of " << outfixed[ scan.parameter ] << std::endl;
	std::cout << "\t{ ";
	for ( auto &x : scan.values )
		std::cout << x << ", ";
	std::cout << "}" << std::endl;

	std::cout << Footer( scan ) << std::endl;
}


template<typename T> int ListScans( std::list<DispReln::Config::Scan> ScanList, T physics )
{
	for ( auto &S : ScanList )
	{
		ListScan( S, physics );
	}
	return 0;
}


namespace po = boost::program_options;

int main( int argc, char** argv )
{
	po::options_description desc( "Allowed Options" );
	desc.add_options()
		( "help", "Produce this help message" )
		( "config-file,c", po::value<std::string>(), "Configuration file" )
		( "simulate,s", "Parse the configuration, and produce an output detailing what parameter scans will be done, without actually performing the scans." )
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm); 
	
	std::string config_file( "DispSolver.xml" );
	if ( vm.count( "config-file" ) )
	{
		config_file = vm[ "config-file" ].as<std::string>();
	}

	bool Simulate = false;
	if ( vm.count( "simulate" ) )
	{
		Simulate = true;
	}

	DispReln::Config::PhysicsTypes physics = DispReln::Config::GetPhysicsModel( config_file );
	auto scans = DispReln::Config::GenerateScans( config_file );
	DispReln::GKSlab FullSlab;
	DispReln::ElectrostaticSlab ESSlab;

	std::cout << std::setprecision( 10 );

	if ( Simulate )
	{
		switch ( physics )
		{
			case DispReln::Config::PhysicsTypes::GKSlab:
				FullSlab =  DispReln::Config::ReadConfig<DispReln::GKSlab>( config_file );
				std::cout << "# Solutions to the full gyrokinetic dispersion relation in a slab" << std::endl;
				return ListScans( scans, FullSlab );
				break;
			case DispReln::Config::PhysicsTypes::ElectrostaticSlab:
				ESSlab =  DispReln::Config::ReadConfig<DispReln::ElectrostaticSlab>( config_file );
				std::cout << "# Solutions to the electrostatic gyrokinetic dispersion relation in a slab" << std::endl;
				return ListScans( scans, ESSlab );
				break;
			default:
				std::cerr << "Config file specifies an invalid physics model" << std::endl;
				return -1;
				break;
		}
	}

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


