

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Config.h"

#include <map>
#include <iostream>

// Masses in amu
static const double electron_mass = 5.48579909e-4;
static const double deuteron_mass = 2.01355321;
static const double proton_mass   = 1.00727647;

/*
static const std::map< std::string, DispReln::Species > PredefinedSpecies{
	{ "Deuterium", DispReln::Species( 1.0, 1.0, 1.0, deuterium_mass, 0.0, 0.0 )},
	{ "Hydrogen", DispReln::Species( 1.0, 1.0, 1.0, proton_mass, 0.0, 0.0 )},
	{ "Electron", DispReln::Species( 1.0, 1.0, -1.0, electron_mass, 0.0, 0.0 )}
};
*/


namespace pt = boost::property_tree;

namespace DispReln {
	namespace Config {
		using namespace std::string_literals;
		static const std::map< std::string, ScanTypes > ScanMap{
			{ "kpar", ScanTypes::kpar },
			{ "kx", ScanTypes::kx },
			{ "ky", ScanTypes::ky },
			{ "fprim", ScanTypes::fprim },
			{ "tprim", ScanTypes::tprim },
			{ "beta_ref", ScanTypes::beta }
		};

		static const std::map< std::string, Normalization > NormMap{
			{ "default", Normalization::Default },
			{ "ref", Normalization::Default },
			{ "Alfven", Normalization::Alfven },
			{ "alfven", Normalization::Alfven },
			{ "kpar", Normalization::kRef },
			{ "kpar_va", Normalization::kAlfven }
		};



		static const std::map< std::string, PhysicsTypes> PhysicsMap{
			{ "gkslab", PhysicsTypes::GKSlab },
			{ "electrostatic_slab", PhysicsTypes::ElectrostaticSlab },
			{ "edge_slab", PhysicsTypes::EdgeSlab }
		};

		std::list<DispReln::Species> GetSpeciesList( pt::ptree const &tree )
		{
			std::list<DispReln::Species> spec_list;
			static const std::string spec_tag( "species" );
			for ( auto &subtag : tree )
			{
				double T,n,m,q,fp,tp;
				if ( subtag.first != spec_tag )
					continue;
				T = subtag.second.get<double>( "<xmlattr>.Temperature" );
				n = subtag.second.get<double>( "<xmlattr>.Density" );
				m = subtag.second.get<double>( "<xmlattr>.Mass" );
				q = subtag.second.get<double>( "<xmlattr>.Charge" );
				fp = subtag.second.get<double>( "<xmlattr>.fprim" );
				tp = subtag.second.get<double>( "<xmlattr>.tprim" );
				spec_list.emplace_back( T, n, q, m, fp, tp );
			}
			return spec_list;
		}

		PhysicsTypes GetPhysicsModel( std::string const& filename )
		{
			pt::ptree root;
			pt::read_xml( filename, root );
			PhysicsTypes type = PhysicsTypes::Invalid;
			std::string scan_tag( "scan" );

			for ( auto &tags : root )
			{
				// There can be multiple scan tags
				if ( tags.first == scan_tag )
					continue;
				else
				{
					if ( type != PhysicsTypes::Invalid )
						throw std::invalid_argument( "There can be only one top-level tag that is not a scan tag" );
					else
					{
						auto it = PhysicsMap.find( tags.first );
						if ( it == PhysicsMap.end() )
							throw std::invalid_argument( "Top level tags that aren't a scan tag should be a physics tag!" );
						else
							type = it->second;
					}
				}
			}
			if ( type == PhysicsTypes::Invalid )
				throw std::invalid_argument( "You need a physics tag!" );
			return type;
		}


		std::list<Scan> GenerateScans( std::string const& filename )
		{
			pt::ptree root;
			pt::read_xml( filename, root );
			std::string scan_tag( "scan" );
			std::list<Scan> Scans;

			for ( auto &tags : root )
			{
				// There can be multiple scan tags
				if ( tags.first == scan_tag )
				{
					ScanTypes param;
					ScanMode mode;
					unsigned int speciesIndex = 0;

					std::string parameter = tags.second.get<std::string>( "<xmlattr>.parameter" );
					auto it = ScanMap.find( parameter );

					if ( it != ScanMap.end() )
						param = it->second;
					else
						throw std::invalid_argument( "Unknown Paramter to scan in: " + parameter  );

					if ( param == ScanTypes::fprim || param == ScanTypes::tprim )
					{
						// Needs a species index
						try {
							speciesIndex = tags.second.get<unsigned int>( "<xmlattr>.species" );
						} catch ( ... )
						{
							throw std::invalid_argument( "When you scan in density or temperature gradient you need to specify which species to scan in." );
						}
					}
					
					std::string type = tags.second.get<std::string>( "<xmlattr>.type" );
					if ( type == "MostUnstableMode" )
						mode = ScanMode::MostUnstableMode;
					else if ( type == "TrackRoot" )
						mode = ScanMode::TrackRoot;
					else
						throw std::invalid_argument( "Unknown Scanning Mode: " + type );

					Scan badger_Scan( param, mode );
					badger_Scan.beta = 0.0;
					badger_Scan.sIndex = speciesIndex;
					for ( auto &subtags : tags.second )
					{
						if ( subtags.first == "val" )
							badger_Scan.values.emplace_back( std::stod( subtags.second.data() ) );
						else if ( subtags.first == "fixed" )
						{
							for ( auto &fix : subtags.second.get_child( "<xmlattr>" ) )
							{
								auto it = ScanMap.find( fix.first );
								if ( it == ScanMap.end() )
									throw std::invalid_argument( "Can't fix the following variable -- " + fix.first );
								else 
								{
									badger_Scan.fixed.emplace_back( it->second, std::stod( fix.second.data() ) );
									if ( it->second == ScanTypes::beta )
										badger_Scan.beta = std::stod( fix.second.data() );
								}
							}
						}
						else if ( subtags.first == "range" )
						{
							double llim = subtags.second.get<double>( "llim" );
							double ulim = subtags.second.get<double>( "ulim" );
							double increment = subtags.second.get<double>( "increment" );
							for ( double x = llim; x <= ulim; x += increment )
								badger_Scan.values.emplace_back( x );
						}
						else if ( subtags.first == "box" )
						{
							std::istringstream lower_s( subtags.second.get<std::string>( "lower" ) );
							std::istringstream upper_s( subtags.second.get<std::string>( "upper" ) );
							lower_s >> badger_Scan.box.lower;
							upper_s >> badger_Scan.box.upper;
						}
						else if ( subtags.first == "output" )
						{
							std::string normalize_str = subtags.second.get( "<xmlattr>.normalization", "default" );
							auto it = NormMap.find( normalize_str );
							if ( it == NormMap.end() )
								throw std::invalid_argument( "Unknown output normalization " + normalize_str );
							else
								badger_Scan.normalization = it->second;
						}
						else if ( subtags.first == "<xmlattr>" )
						{
							continue;
						}
						else
						{
							throw std::invalid_argument( "Unknown tag " + subtags.first + " in configuration file." );
						}
					}

					badger_Scan.values.sort();
					double tol=1e-14;
					badger_Scan.values.unique( [ = ]( double x, double y ){ return ( std::abs( x - y ) < tol );});

					Scans.push_back( badger_Scan );

				}
			}

			return Scans;
		}

		template <> DispReln::GKSlab ReadConfig( std::string const& filename )
		{
			pt::ptree root;
			pt::read_xml( filename, root );
			pt::ptree gkslab;
			try {
				gkslab = root.get_child( "gkslab" );
			} catch (std::exception &bar)
			{
				throw std::invalid_argument( "XML File should contain a top-level <gkslab>" );
			}
			double beta_ref;
			try {
				beta_ref = gkslab.get<double>( "<xmlattr>.beta_ref" );
			} catch ( std::exception &foo )
			{
				throw std::invalid_argument( "Reference beta must be set for a gk slab" );
			}
			auto speclist = GetSpeciesList( gkslab );
			return DispReln::GKSlab( speclist, beta_ref );
		}

		template <> DispReln::ElectrostaticSlab ReadConfig( std::string const& filename )
		{
			pt::ptree root;
			pt::read_xml( filename, root );
			pt::ptree gkslab;
			try {
				gkslab = root.get_child( "electrostatic_slab" );
			} catch (std::exception &bar)
			{
				throw std::invalid_argument( "XML File should contain a top-level <electrostatic_slab>" );
			}
			auto speclist = GetSpeciesList( gkslab );
			return DispReln::ElectrostaticSlab( speclist );
		}
	}
}
