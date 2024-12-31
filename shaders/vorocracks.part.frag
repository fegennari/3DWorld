// variant of Vorocracks: https://shadertoy.com/view/lsVyRy
// integrated with cracks here: https://www.shadertoy.com/view/Xd3fRN

float RATIO = 1.,              // stone length/width ratio
      CRACK_depth = 3.,
      CRACK_zebra_scale = 1.,  // fractal shape of the fault zebra
      CRACK_zebra_amp = .67,
      CRACK_profile = 1.,      // fault vertical shape  1.  .2 
      CRACK_slope = 100.,       //                      10.  1.4
      CRACK_width = .0;

// === Voronoi =====================================================

#define hash22(p)  fract( 18.5453 * sin( p * mat2(127.1,311.7,269.5,183.3)) )

// --- Voronoi distance to borders. inspired by https://www.shadertoy.com/view/ldl3W8
vec3 voronoiB( vec2 u )  // returns len + id
{
    vec2 iu = floor(u), C, P;
	float m = 1e9,d;
    for( int k=0; k < 9; k++ ) {
        vec2  p = iu + vec2(k%3-1,k/3-1),
              o = hash22(p),
      	      r = p - u + o;
		d = dot(r,r);
        if( d < m ) m = d, C = p-iu, P = r;
    }
    m = 1e9;
    
    for( int k=0; k < 25; k++ ) {
        vec2 p = iu+C + vec2(k%5-2,k/5-2),
		     o = hash22(p),
             r = p-u + o;

        if( dot(P-r,P-r)>1e-5 )
        m = min( m, .5*dot( (P+r), normalize(r-P) ) );
    }
    return vec3( m, P+u );
}

// === pseudo Perlin noise =============================================
#define rot(a) mat2(cos(a),-sin(a),sin(a),cos(a))
    
// --- 2D
#define hash21(p) fract(sin(dot(p,vec2(127.1,311.7)))*43758.5453123)

float noise2_(vec2 p) {
    vec2 i = floor(p);
    vec2 f = fract(p); f = f*f*(3.-2.*f); // smoothstep
    return mix( mix(hash21(i+vec2(0,0)),hash21(i+vec2(1,0)),f.x),
                mix(hash21(i+vec2(0,1)),hash21(i+vec2(1,1)),f.x), f.y);
}

#define noise22(p) vec2(noise2_(p),noise2_(p+17.7))

vec2 fbm22(vec2 p) {
    vec2 v = vec2(0);
    float a = .5;
    mat2 R = rot(.37);

    for (int i = 0; i < 6; i++, p*=2.,a/=2.) 
        p *= R,
        v += a * noise22(p);
    return v;
}
    
// ======================================================

float get_crack_weight(in vec2 tc) {
    vec2 U = tc;
    vec3 H0;
    float weight = 0.0;

    for(float i=0.; i<CRACK_depth ; i++) {
        vec2 V =  U / vec2(RATIO,1),                  // voronoi cell shape
             D = CRACK_zebra_amp * fbm22(U/CRACK_zebra_scale) * CRACK_zebra_scale;
        vec3  H = voronoiB( V + D ); if (i==0.) H0=H;
        float d = H.x;                                // distance to cracks
        d = min( 1., CRACK_slope * pow(max(0.,d-CRACK_width),CRACK_profile) );
  
        weight += (1.-d) / exp2(i);
        U *= 1.5 * rot(.37);
    }
    return 1.0 - weight;
}