#ifndef DISPRELN_H
#define DISPRELN_H

#include <boost/math/constants/constants.hpp>
#include <complex>
#include <list>


namespace DispReln {
	using Real = double;
	using Complex = std::complex<double>;

	std::complex<double> Z( std::complex<double> xi );

	std::complex<double> AcousticDisp( std::complex<double> xi, double ZoverTau );
	
	template<unsigned int> std::complex<double> Zeta( std::complex<double> xi );
	template<unsigned int,unsigned int,unsigned int> double Gamma( double );

	template<> double inline Gamma<1,0,0>( double alpha )
	{
		return std::exp( -alpha ) * std::cyl_bessel_i( 0.0, alpha );
	}

	template<> double inline Gamma<2,1,0>( double alpha )
	{
		return std::exp( -alpha ) * ( std::cyl_bessel_i( 0.0, alpha ) - std::cyl_bessel_i( 1.0, alpha ) );
	}

	template<> double inline Gamma<3,0,0>( double alpha )
	{
		return std::exp( -alpha ) * ( std::cyl_bessel_i( 0.0, alpha ) + alpha * ( std::cyl_bessel_i( 1.0, alpha ) - std::cyl_bessel_i( 0.0, alpha ) ) );
	}

	// General Gamma<l,m,n> can be done if you have a generalized hypergeometric pFq
	// evaluator.
	
	template<unsigned int l,unsigned int m,unsigned int n> double Gamma( double alpha ) { throw std::logic_error( "Unimplemented!" );}

	template<> std::complex<double> inline Zeta<0>( std::complex<double> xi ) 
	{
		return Z( xi );
	}

	template<> std::complex<double> inline Zeta<1>( std::complex<double> xi ) 
	{
		return 1.0 + xi*Z( xi );
	}

	template<> std::complex<double> inline Zeta<2>( std::complex<double> xi ) 
	{
		return xi*( 1.0 + xi*Z( xi ) );
	}

	template<unsigned int l> inline std::complex<double> Zeta( std::complex<double> xi )
	{
		if ( l % 2 == 0 )
			return xi*Zeta<l-1>( xi );
		else
			return xi*Zeta<l-1>( xi ) + std::tgamma( ( l )/2.0 )/boost::math::double_constants::root_pi;
	}
	
	struct Species {
		double Temperature; // Relative to a fiducial reference temperature.
		double Density;     // Relative to a fiducial reference density.
		double Z; // Charge in units of the elemental charge (electrons have Z = -1)
		double mass; // Relative to a fiducial mass.
		double fprim; // a/L_n , positive for centrally peaked profiles
		double tprim; // a/L_T , positive for centrally peaked profiles
		double rho; // rho_s in terms of rho_ref
		Species( double T_in, double Dens_in, double Z_in, double mass_in, double fprim_in, double tprim_in ) :
			Temperature( T_in ), Density( Dens_in ), Z( Z_in ), mass( mass_in ), fprim( fprim_in ), tprim( tprim_in )
		{
			rho = ::sqrt( Temperature * mass )/Z;
		}
		// Copy constructor
		Species( Species const& other ) :
			Temperature( other.Temperature ), Density( other.Density ), Z( other.Z ), mass( other.mass ), fprim( other.fprim ), 
			tprim( other.tprim ), rho( other.rho ) 
		{
			rho = ::sqrt( Temperature * mass )/::abs( Z );
		};

	};


	class ElectrostaticSlab {
		public:
			ElectrostaticSlab( std::list<Species> const & spec_list )
			{
				SpeciesList.clear();
				for ( auto i : spec_list )
					SpeciesList.emplace_back( i );
				_kpar = 1.0; _kx=0; _ky=0;
				recalculate();
			}

			ElectrostaticSlab( std::initializer_list<Species> const & spec_list )
			{
				SpeciesList.clear();
				for ( auto i : spec_list )
					SpeciesList.emplace_back( i );
				_kpar = 1.0; _kx=0; _ky=0;
				recalculate();
			}


			Complex operator()( Complex xi )
			{
				// _ky & _kx come normalized to rho_ref
				// k_|| normalized to a;
				Complex D( 0.0, 0.0 );
				for ( auto x : SpeciesList )
				{
					if ( x.s.mass == 0.0 )
					{
						D += x.boltz;
						continue;
					}

					Complex xi_s = x.vt*xi;
					
					D += ( x.boltz )*( 1.0 + ( xi_s - x.om_star )*Zeta<0>( xi_s )*Gamma<1,0,0>( x.alpha ) - x.om_star_t*( Zeta<0>( xi_s )*Gamma<3,0,0>( x.alpha ) + Zeta<2>( xi_s )*Gamma<1,0,0>( x.alpha ) ) );

				}
				return D;
			}
			void set_kpar( double kp ){_kpar = kp;recalculate();};
			void set_kx( double kx ){_kx= kx;recalculate();};
			void set_ky( double ky ){_ky= ky;recalculate();};

			struct _species {
				DispReln::Species s;
				Real om_star,om_star_t,alpha;
				Real vt,boltz;
				_species( Species const& s_in ) : s( s_in ),om_star( 0.0 ),om_star_t( 0.0 ),alpha( 0.0 )
				{
					vt = ::sqrt( s_in.mass/s_in.Temperature );
					boltz = s_in.Z*s_in.Z*s_in.Density / s_in.Temperature;
				};
				_species( _species const& _s_in ) : s( _s_in.s ), om_star( _s_in.om_star ), om_star_t( _s_in.om_star_t ), alpha( _s_in.alpha )
				{
					vt = ::sqrt( _s_in.s.mass/_s_in.s.Temperature );
					boltz = _s_in.s.Z*_s_in.s.Z*_s_in.s.Density / _s_in.s.Temperature;
				};
			};

			void recalculate()
			{
				for ( auto &x : SpeciesList ) 
				{
					x.s.rho = ::sqrt( x.s.Temperature * x.s.mass )/::abs( x.s.Z );
					x.om_star = ( x.s.Z / ::abs ( x.s.Z ) )*0.5 * _ky * x.s.rho * x.s.fprim /  _kpar;
					x.om_star_t = ( x.s.Z / ::abs ( x.s.Z ) )*0.5 * _ky * x.s.rho * x.s.tprim /  _kpar;
					x.alpha = ( _ky*_ky + _kx*_kx )*( x.s.rho * x.s.rho )/2.0;
					x.vt = ::sqrt( x.s.mass / x.s.Temperature );
					x.boltz = x.s.Z * x.s.Z * x.s.Density / x.s.Temperature;
				}
			}


			std::list<_species> SpeciesList;
			ElectrostaticSlab( ElectrostaticSlab const &other )
				: SpeciesList( other.SpeciesList )
			{
				_kpar = other._kpar; _kx = other._kx; _ky = other._ky;
				recalculate();
			};
			ElectrostaticSlab()
			{
				SpeciesList.clear();
				_kpar = 0.0;
				_ky = 0.0;
				_kx = 0.0;
			};
		protected:
			double _kpar,_kx,_ky;

	};

	class GKSlab {
		public:
			GKSlab( std::list<Species> const & spec_list, double beta )
			{
				SpeciesList.clear();
				for ( auto i : spec_list )
					SpeciesList.emplace_back( i );
				_kpar = 1.0; _kx=0; _ky=0;
				beta_ref = beta;
				recalculate();
			}

			GKSlab( std::initializer_list<Species> const & spec_list, double beta )
			{
				SpeciesList.clear();
				for ( auto i : spec_list )
					SpeciesList.emplace_back( i );
				_kpar = 1.0; _kx=0; _ky=0;
				beta_ref = beta;
				recalculate();
			}


			Complex operator()( Complex xi )
			{
				// _ky & _kx come normalized to rho_ref
				// k_|| normalized to a;
				Complex A( 0.0, 0.0 ),B( 0.0, 0.0 ),C( 0.0, 0.0 ),D( 0.0,0.0 ),E( 0.0, 0.0 );
				for ( auto x : SpeciesList )
				{
					Complex xi_s = x.vt*xi;

					A += ( x.boltz )*( 1.0 + ( xi_s - x.om_star )*Zeta<0>( xi_s )*Gamma<1,0,0>( x.alpha ) - x.om_star_t*( Zeta<0>( xi_s )*Gamma<3,0,0>( x.alpha ) + Zeta<2>( xi_s )*Gamma<1,0,0>( x.alpha ) ) );
					B += ( x.boltz/x.vt ) * ( xi_s  - ( xi_s - x.om_star )*Gamma<1,0,0>( x.alpha ) + x.om_star_t*( Gamma<3,0,0>( x.alpha ) + 0.5*Gamma<1,0,0>( x.alpha ) ) );
					C += ( x.s.Z * x.s.Density ) * ( ( xi_s - x.om_star )*Zeta<0>( xi_s )*Gamma<2,1,0>( x.alpha ) - x.om_star_t*( Zeta<0>( xi_s )*Gamma<4,1,0>( x.alpha ) + Zeta<2>( xi_s )*Gamma<2,1,0>( x.alpha ) ) );
					D += ( x.s.Temperature * x.s.Density ) * ( -( xi_s - x.om_star )*Zeta<0>( xi_s )*Gamma<3,1,1>( x.alpha ) + x.om_star_t*( Zeta<0>( xi_s )*Gamma<5,1,1>( x.alpha ) + Zeta<2>( xi_s )*Gamma<3,1,1>( x.alpha ) ) );
					E += ( x.boltz/x.vt ) * ( ( xi_s - x.om_star )*Gamma<2,1,0>( x.alpha ) - x.om_star_t*( Gamma<4,1,0>( x.alpha ) + 0.5*Gamma<2,1,0>( x.alpha ) ) );

				}

				return  ( A*alpha_ref/beta_ref - A*B + B*B )*( 2.0*A/beta_ref - A*D + C*C) - (A*E + B*C)*(A*E + B*C);
			}
			void set_kpar( double kp ){_kpar = kp;recalculate();};
			void set_kx( double kx ){_kx= kx;recalculate();};
			void set_ky( double ky ){_ky= ky;recalculate();};

			struct _species {
				DispReln::Species s;
				Real om_star,om_star_t,alpha;
				Real vt,boltz;
				_species( Species const& s_in ) : s( s_in ),om_star( 0.0 ),om_star_t( 0.0 ),alpha( 0.0 )
				{
					vt = ::sqrt( s_in.mass/s_in.Temperature );
					boltz = s_in.Z*s_in.Z*s_in.Density / s_in.Temperature;
				};
				_species( _species const& _s_in ) : s( _s_in.s ), om_star( _s_in.om_star ), om_star_t( _s_in.om_star_t ), alpha( _s_in.alpha )
				{
					vt = ::sqrt( _s_in.s.mass/_s_in.s.Temperature );
					boltz = _s_in.s.Z*_s_in.s.Z*_s_in.s.Density / _s_in.s.Temperature;
				};
			};

			void recalculate()
			{
				alpha_ref = ( _ky*_ky + _kx*_kx )/2.0;
				for ( auto &x : SpeciesList ) 
				{
					x.s.rho = ::sqrt( x.s.Temperature * x.s.mass )/::abs( x.s.Z );
					x.om_star = ( x.s.Z / ::abs ( x.s.Z ) )*0.5 * _ky * x.s.rho * x.s.fprim /  _kpar;
					x.om_star_t = ( x.s.Z / ::abs ( x.s.Z ) )*0.5 * _ky * x.s.rho * x.s.tprim /  _kpar;
					x.alpha = ( alpha_ref )*( x.s.rho * x.s.rho );
					x.vt = ::sqrt( x.s.mass / x.s.Temperature );
					x.boltz = x.s.Z * x.s.Z * x.s.Density / x.s.Temperature;
				}
			}


			std::list<_species> SpeciesList;
			GKSlab( GKSlab const &other )
				: SpeciesList( other.SpeciesList )
			{
				_kpar = other._kpar; _kx = other._kx; _ky = other._ky;
				beta_ref = other.beta_ref;
				recalculate();
			};
			GKSlab()
			{
				SpeciesList.clear();
				_kpar = 0.0;
				_ky = 0.0;
				_kx = 0.0;
				beta_ref = 0.0;
			};
			double beta_ref;
		protected:
			double _kpar,_kx,_ky;
			double alpha_ref;

	};
}
#endif // DISPRELN_H
