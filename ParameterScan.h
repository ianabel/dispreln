#include "DispReln.h"

namespace DispReln {
template<typename T> struct ky_scanner 
{
	T impl;
	ky_scanner( ky_scanner<T> const& other ) : impl( other.impl ) {};
	ky_scanner( T const & in ) : impl( in ) {};
	Complex operator()( Complex xi ) { return impl( xi );};
	void set_param( Real ky )
	{
		impl.set_ky( ky );
	}
};

template<typename T> struct kx_scanner 
{
	T impl;
	kx_scanner( kx_scanner<T> const& other ) : impl( other.impl ) {};
	kx_scanner( T const & in ) : impl( in ) {};
	Complex operator()( Complex xi ) { return impl( xi );};
	void set_param( Real kx )
	{
		impl.set_kx( kx );
	}
};

template<typename T> struct kpar_scanner 
{
	T impl;
	kpar_scanner( kpar_scanner<T> const& other ) : impl( other.impl ) {};
	kpar_scanner( T const & in ) : impl( in ) {};
	Complex operator()( Complex xi ) { return impl( xi );};
	void set_param( Real kpar )
	{
		impl.set_kpar( kpar );
	}
};

template<typename T> struct eta_i_scanner 
{
	T impl;
	eta_i_scanner( eta_i_scanner<T> const& other ) : impl( other.impl ) {};
	eta_i_scanner( T const & in ) : impl( in ) {};
	Complex operator()( Complex xi ) { return impl( xi );};
	void set_param( Real eta_i )
	{
		impl.SpeciesList.front.s.tprim = eta_i*impl.SpeciesList.front.s.tprim;
	}
};
}
