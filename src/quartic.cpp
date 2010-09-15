// found at:   http://www.maths.tcd.ie/~dwmalone/p/561-99/quartic.c
// similar to: http://chem.skku.ac.kr/~wkpark/tutor/chem/pdb2pov/povray31/source/polysolv.c

#include <stdio.h>
#include <math.h>
#include <iostream>

using std::cout;
using std::endl;

int solve_linear(const double *poly,double *roots);
int solve_quadratic(const double *poly,double *roots);
int solve_cubic(const double *poly,double *roots);
int solve_quartic(const double *poly,double *roots);


#ifndef cbrt
inline double cbrt(double val) {
	return ((val < 0.0) ? -pow(-val, 1.0/3.0) : pow(val, 1.0/3.0));
} // should be in math.h?
#endif


#define SIN60 0.86602540378443864675
#define TOO_CLOSE 0.00001
#define SAMESIGN(a,b) ( ((a)>=0) ^ ((b)<0) )

int solve_linear(const double *poly,double *roots)
{
	if( poly[1] != 0.0 ) {
		roots[0] = - poly[0] / poly[1] ;
		return 1;
	}
	return 0;
}

int solve_quadratic(const double *poly,double *roots)
{
	double t, b,c;
	if( poly[2] == 0.0 ) return solve_linear(poly,roots);
	b = poly[1] / poly[2] / 2.0;
	c = poly[0] / poly[2];
	t = b*b - c;
	if( t < 0.0 ) return 0;

	if( fabs(t) < (TOO_CLOSE*TOO_CLOSE) )  {
		roots[0] = -b;
		return 1;
	}
	t = sqrt(t);
	roots[0] = -b - t;
	roots[1] = -b + t;
	return 2;
}

int solve_cubic(const double *poly,double *roots)
{
	double a0,a1,a2,q,r,d;
	if( poly[3] == 0.0 ) return solve_quadratic(poly,roots);
	a0 = poly[0] / poly[3];
	a1 = poly[1] / poly[3];
	a2 = poly[2] / (poly[3] * 3.0);
	q = a1 / 3.0 - a2*a2;
	r = ( a1*a2 - a0 ) * 0.5  - a2*a2*a2 ;
	
	// This is a simpler case.
	if( fabs(r) < TOO_CLOSE ) {
		if( q >= 0.0 ) {
			roots[0] = -a2;
			return 1;
		}
		else {
			double t = sqrt( -3.0 * q );
			roots[0] = -a2 - t;
			roots[1] = -a2;
			roots[2] = -a2 + t;
			return 3;
		}
	}

	// Now we get on with normal cases.
	d = q*q*q + r*r;

	// One double real root and one normal.
	if( fabs(d) < ( TOO_CLOSE * TOO_CLOSE ) ) {
		r = cbrt(r);
		roots[0] = 2.0*r - a2;
		roots[1] = -r - a2;
		return 2;
	}

	// One real and two complex.
	if( d > 0.0 ) {
		double s1,s2;
		d = sqrt(d);
		s1 = cbrt(r + d);
		s2 = cbrt(r - d); // d > r?
		roots[0] = s1 + s2 - a2;
		return 1;
	}

	// Three real, surprisingly this is the most difficult one.
	d = sqrt(-d);
	q = atan2(d,r) / 3.0;
	r = cbrt(_hypot(d,r));
	{
		double s,t;
		s =  (2.0 * r) * cos(q);
		t = fabs(2.0 * r * sin(q)) * SIN60 ;
		roots[0] =  s + -a2 ;
		roots[1] = -( 0.5 * s ) -a2 - t;
		roots[2] = -( 0.5 * s ) -a2 + t;
		return 3;
	}
}

int solve_quartic(const double *poly,double *roots)
{
	double a3,a2,a1,a0;
	double resolvant[4];
	double resolvant_roots[3];
	int resolvant_num_roots;
	int i;
	if( poly[4] == 0.0 ) return solve_cubic(poly,roots);
	a0 = poly[0] / poly[4];
	a1 = poly[1] / poly[4];
	a2 = poly[2] / poly[4];
	a3 = poly[3] / poly[4];
	resolvant[0] = -( a1*a1 + a0 * a3 * a3  - 4.0 * a0 * a2 );
	resolvant[1] = a1*a3 - 4.0 * a0;
	resolvant[2] = -a2;
	resolvant[3] = 1.0;
	resolvant_num_roots = solve_cubic(resolvant,resolvant_roots);

	for( i = 0 ; i < resolvant_num_roots ; i++ ) {
		double u,s,t;
		u = resolvant_roots[i];
		s = a3*a3 / 4.0 + u - a2;
		t = u * u / 4.0 - a0;

		if( s >= 0.0 && t >= 0.0 ) {
			double quad[3];
			int root_count;
			s = sqrt(s);
			t = sqrt(t);
			if( !SAMESIGN(s*t,a3 * u / 2.0 - a1) ) t = -t;
			quad[0] = u / 2.0 + t;
			quad[1] = a3 / 2.0 + s;
			quad[2] = 1.0;
			root_count = solve_quadratic(quad,roots);
			quad[0] = u / 2.0 - t;
			quad[1] = a3 / 2.0 - s;
			root_count += solve_quadratic(quad,roots+root_count);
			return root_count;
		}
	}
	return 0;
}


