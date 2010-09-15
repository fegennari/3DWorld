// found at: http://research.microsoft.com/%7Eawf/graphics/ray-torus.html

#include "3DWorld.h" // for PI


int solve_quartic(const double *poly,double *roots); // quartic.cpp

bool line_intersect_torus(double ax, double ay, double az,
						  double bx, double by, double bz,
						  double R, double r, double vlength, float &t)
//              
// Intersect ray
//      [ax]     [bx]              
//	    [ay] + t [by]
//	    [az]     [bz]
// with torus at origin, major radius R, minor radius r
//
{
	// This struct is syntactic sugar so that a.b,
	// (the dot product) looks just right (:-)
	struct { double a,b; } a;
	a.b = ax*bx + ay*by + az*bz;
	a.a = ax*ax + ay*ay + az*az;

	// Set up quartic in t:
	//
	//  4     3     2
	// t + A t + B t + C t + D = 0
	//
	double R2 = R*R;
	double K = a.a - r*r - R2;
	double A = 4*a.b;
	double B = 2*(2*a.b*a.b + K + 2*R2*bz*bz);
	double C = 4*(K*a.b + 2*R2*az*bz);
	double D = K*K + 4*R2*(az*az - r*r);

	// Solve quartic...
	double roots[4], poly[5] = {D, C, B, A, 1.0};
	int nroots = solve_quartic(poly,roots);
	t = 1.0;
	bool found(0);

	while(nroots--) {
		float const t0((float)roots[nroots]/vlength);
		//double const x(ax + t*bx), y(ay + t*by), l(R*(PI_TWO - atan2(y,x)));
		//if (l <= vlength && l >= 0)

		if (t0 <= t && t0 >= 0.0) {
			t     = t0;
			found = 1;
		}
	}
	return found;
}


