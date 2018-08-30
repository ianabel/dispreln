#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include "DispReln.h"
#include "RootFinder.h"



namespace DispReln {
	namespace Config {

		enum struct ScanTypes { kpar, kx, ky, fprim, tprim, Invalid = 255 };
		enum struct PhysicsTypes { GKSlab, ElectrostaticSlab, EdgeSlab, Invalid = 255 };
		enum struct ScanMode { TrackRoot, MostUnstableMode };

		template<typename T> T ReadConfig( std::string const& );
		template<> DispReln::GKSlab ReadConfig( std::string const& );
		template<> DispReln::ElectrostaticSlab ReadConfig( std::string const& );
		PhysicsTypes GetPhysicsModel( std::string const& filename );
	
		using Fixture = std::pair< ScanTypes, double >;
		struct Scan {
			ScanTypes parameter;
			ScanMode  mode;
			std::list<double> values;
			std::list<Fixture> fixed;
			RootFinder::RootBoundingBox box;
			Scan( ScanTypes x, ScanMode y ) { parameter = x; mode = y; values.clear(); fixed.clear();};
			Scan( Scan const& o ) : parameter( o.parameter ), mode( o.mode ), values( o.values ), fixed( o.fixed ), box( o.box ) {};
		};

		std::list<Scan> GenerateScans( std::string const& filename );
	}
}

#endif // CONFIG_H
