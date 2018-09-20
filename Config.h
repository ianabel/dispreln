#ifndef CONFIG_H
#define CONFIG_H

#include "DispReln.h"
#include "RootFinder.h"

#include <string>


namespace DispReln {
	namespace Config {

		enum struct ScanTypes { kpar, kx, ky, fprim, tprim, beta, Invalid = 255 };
		enum struct PhysicsTypes { GKSlab, ElectrostaticSlab, EdgeSlab, Invalid = 255 };
		enum struct ScanMode { TrackRoot, MostUnstableMode, AllRoots };
		enum struct Normalization { Default = 1, Alfven, kRef, kAlfven, Invalid = 255 };

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
			unsigned int sIndex;
			Normalization normalization;
			double beta,tolerance;
			Scan( ScanTypes x, ScanMode y ) { parameter = x; mode = y; values.clear(); fixed.clear(); beta = 0.0; tolerance = 1e-3;};
			Scan( Scan const& o ) : parameter( o.parameter ), mode( o.mode ), values( o.values ), fixed( o.fixed ), 
											box( o.box ), sIndex( o.sIndex ), normalization( o.normalization ), beta( o.beta ), tolerance( o.tolerance ) {};
		};

		std::list<Scan> GenerateScans( std::string const& filename );
	}
}

#endif // CONFIG_H
