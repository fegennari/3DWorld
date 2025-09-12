
vec3 colorize(float val) {
	float a = 5*val, b = 7*val, c = 11*val;
	return vec3((a - int(a)), (b - int(b)), (c - int(c)));
}
vec2 conjugate   (vec2 v) {return vec2((v.x*v.x - v.y*v.y), -2.0*v.x*v.y);}
vec2 complex_mult(vec2 v) {return vec2((v.x*v.x - v.y*v.y), (v.x*v.y + v.y*v.x));}

// https://www.shadertoy.com/view/wdtXDs
vec4 gen_ZzArt(vec2 pos, float seed) {
	vec4 f = vec4(10.0*seed);
	vec4 a=pos.xyxy;
	a.xy*=vec2(-8.434,-5.070);
	a.xy+=vec2(1.605,3.840);
	a.wz*=vec2(-8.434,-5.070);
	a.wz+=vec2(1.605,3.840);
	vec4 b=a;
	b.zywx*=atan(a.xxzy);
	a.xwzy+=(b.wyzw+ f.x);
	b.zwxy+=atan(a.zzxz);
	b.wyxz+=sin(a.wzzy);
	a.wxyz=(b.xwyz);
	a.xzwy*=sin(a.wzzx);
	a.zxyw=tan(a.wxxx+ f.y);
	a.xwzy/=abs(a.xwzz);
	// Convert to HSL with Smooth HSV by iq
	a.x=a.x*-0.635+0.101;
	a.y=a.y*0.292;
	b=clamp(abs(mod(a.x*6.+vec4(0,4,2,1),6.)-3.)-1.,0.,1.);
	b=b*b*(3.-2.*b);
	return a.z*mix(vec4(1),b,a.y);
}

// https://www.shadertoy.com/view/DlVczV
vec4 vortex(vec2 p, float seed) {
	seed = 4.0 + 16.0*seed;
	vec2 v=vec2(1.0);
	vec4 O=vec4(0.0);
    //Loop through arcs (i=radius, P=pi, l=length)
    for(float i=.2,l; i<1.;
    //Pick color for each arc
    O+=(cos(i*5.+vec4(0,1,2,3))+1.)*
    //Shade and attenuate light
    (1.+v.y/(l=length(v)+.003))/l)
        //Compute polar coordinate position
        v=vec2(mod(atan(p.y,p.x)+i+i*seed,6.28)-3.14,1)*length(p)-i,
        //Clamp to light length
        v.x-=clamp(v.x+=i,-i,i),
        //Iterate radius
        i+=.05;
    //Tanh tonemap: shadertoy.com/view/ms3BD7
    return tanh(O/1e2);
}

// https://www.shadertoy.com/view/WtjyzR
#define NUM_LAYERS 16.
#define ITER 23
vec4 tex(vec3 p, float seed) {
    float t = seed+78.;
    vec4 o = vec4(p.xyz,3.*sin(t*.1));
    vec4 dec = vec4 (1.,.9,.1,.15) + vec4(.06*cos(t*.1),0,0,.14*cos(t*.23));
    for (int i=0 ; i++ < ITER;) o.xzyw = abs(o/dot(o,o)- dec);
    return o;
}
vec4 colorful(vec2 pos, float seed) {
	seed *= 10.0;
    vec3 col = vec3(0);
    float t= seed* .3;

	for(float i=0.; i<=1.; i+=1./NUM_LAYERS) {
        float d = fract(i+t); // depth
        float s = mix(5.,.5,d); // scale
        float f = d * smoothstep(1.,.9,d); //fade
        col+= tex(vec3(pos*s,i*4.), seed).xyz*f;
    }
    col/=NUM_LAYERS;
    col*=vec3(2,1.,2.);
    col=pow(col,vec3(.5 ));
    return vec4(col,1.0);
}

vec4 gen_abstract_art(vec2 tc, vec3 seed) {
	// seed.r selects the fractal/mode, seed.g selects the location, and seed.b selects the zoom level
	int mode = int(5.99*seed.r);
	vec2 pos = 2.0*tc - vec2(1.0); // [-1.0, 1.0]
	if (mode == 3) {return gen_ZzArt(pos, seed.g);}
	if (mode == 4) {return vortex   (pos, seed.g);}
	if (mode == 5) {return colorful (pos, seed.g);}
	int cix  = min(4, int(4.0*seed.g));
	pos.y    = -pos.y; // invert Y
	vec2 c; // center of window

	if (mode == 0) { // madelbrot
		vec2 centers[4] = {vec2(-1.41117, 0.0), vec2(-1.3953, 0.018367), vec2(-1.23102, 0.167902), vec2(0.437802, 0.340398)};
		c = centers[cix];
	}
	else if (mode == 1) { // tricorn
		vec2 centers[4] = {vec2(-1.41975, 0.0), vec2(0.7093, 1.22854), vec2(0.742104, 1.2854), vec2(-1.36524, -0.0388551)};
		c = centers[cix];
	}
	else { // burning ship
		vec2 centers[4] = {vec2(-1.57535, -0.00717238), vec2(0.807586, -1.40662), vec2(-0.676457, -1.10419), vec2(0.810333, -1.40356)};
		c = centers[cix];
	}
	c       += (0.001 + 0.009*seed.b)*pos;
	vec2 z   = vec2(0.0, 0.0);
	uint val = 0;

	for (; val < 200; ++val) {
		if (dot(z, z) > 4.0) break;
		if      (mode == 0) {z = complex_mult    (z)  + c;} // madelbrot
		else if (mode == 1) {z = conjugate       (z)  + c;} // tricorn
		else                {z = complex_mult(abs(z)) + c;} // burning ship
	}
	float val2 = (float(val) - log2(log2(dot(z, z))) + 1.0)/200.0; // from http://www.iquilezles.org/www/articles/mset_smooth/mset_smooth.htm
	//float val2 = val/200.0;
	return vec4(colorize(val2), 1.0);
}

