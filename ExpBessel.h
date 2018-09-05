

namespace DispReln {
	template<unsigned int> double ExpBessel( double );

	


	// Using 8.451 5 of Gradshteyn and Ryzhik (8th Ed) and multiplying by exp(-z), we have the asymptotic series
	// Exp(-z)*I_n(z) ~= (2 pi z)^{-1/2} * Sum_k ( Gamma[ n + k + 1/2 ]*(-1)^k / (2 z)^k k! Gamma[ n - k + 1/2 ] )
	// where we will only evaluate this for z > 0 and sufficiently large z that the second term in the asymptotic expansion is
	// beyond the resolution of our floating-point type.
	//
	// Factoring out 512 from z, we have 
	//
	// Exp(-z)*I_n(z) ~= sqrt(x/b) * Sum_k ( a_k x^k )
	//
	// with x = 512.0/z
	// a_k = Gamma[ n + k + 1/2 ]*(-1)^k / (1024)^k k! Gamma[ n - k + 1/2 ]
	// b = sqrt( 1024 pi )
	//
	// Tabulating a_k(0) and a_k(1) with Mathematica, one obtains
	static const double a0[] = {1.,0.000244140625,2.682209014892578125e-7,5.4569682106375694275e-10,1.6320278461989801144e-12,6.4547976338924506479e-15,3.1780212959838319026e-17};
	static const double a1[] = {1.,-0.000732421875,-4.470348358154296875e-7,-7.6397554948925971985e-10,-2.0983215165415458614e-12,-7.8891971080907730141e-15,-3.7558433497990740667e-17};
	static const double b = 56.718523228976512874;

	// Hence
	template<> inline double ExpBessel<0>( double z ) {
		if ( z < 512.0 )
			return std::exp( -z )*std::cyl_bessel_i( 0.0, z );
		double x = 512.0/z;

		double EI = ( a0[ 0 ] + x*( a0[ 1 ] + x*( a0[ 2 ] + x*( a0[ 3 ] + x*( a0[ 4 ] + x * ( a0[ 5 ] + x*a0[ 6 ] ) ) ) ) ) );

		return EI*( std::sqrt( x )/b );
	}

	template<> inline double ExpBessel<1>( double z ) {
		if ( z < 512.0 )
			return std::exp( -z )*std::cyl_bessel_i( 1.0, z );
		double x = 512.0/z;

		double EI = ( a1[ 0 ] + x*( a1[ 1 ] + x*( a1[ 2 ] + x*( a1[ 3 ] + x*( a1[ 4 ] + x * ( a1[ 5 ] + x*a1[ 6 ] ) ) ) ) ) );

		return EI*( std::sqrt( x )/b );
	}

	// We don't need this, but if someone else ever does, they can write it :)
	template<unsigned int i> inline double ExpBessel( double z ) {
		if ( z < 512.0 )
			return std::exp( -z )*std::cyl_bessel_i( 0.0, z );
		else
			throw std::logic_error( "Unimplemented." );
	}

}

