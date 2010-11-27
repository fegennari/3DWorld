uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
varying vec3 eye, vpos;

#define TEST_CLIP_T(reg, va, vb, vd, vc) {float t = ((va) - (vb))/(vd); if ((vc) > 0.0) \
                                          {if (t > tmin) tmin = t;} else {if (t < tmax) tmax = t;}}

void main()
{
	if (use_texgen) {
		setup_texgen(0);
	}
	else {
		gl_TexCoord[0] = gl_MultiTexCoord0;
	}	
	gl_Position = ftransform();
	gl_FrontColor = gl_Color;
	
	if (!smoke_enabled) {
		eye  = vec3(0,0,0);
		vpos = vec3(0,0,0);
		return;
	}
	
	// clip to scene bounds
	vec3 v1 = gl_Vertex.xyz;
	vec3 v2 = (inverse(gl_ModelViewMatrix) * vec4(0.0, 0.0, 0.0, 1.0)).xyz; // world space
	float tmin = 0.0, tmax = 1.0;
	vec3 dv = v2 - v1;
	TEST_CLIP_T(0x01, smoke_bb[0], v1.x, dv.x,  dv.x); // -x plane
	TEST_CLIP_T(0x02, smoke_bb[1], v1.x, dv.x, -dv.x); // +x plane
	TEST_CLIP_T(0x04, smoke_bb[2], v1.y, dv.y,  dv.y); // -y plane
	TEST_CLIP_T(0x08, smoke_bb[3], v1.y, dv.y, -dv.y); // +y plane
	TEST_CLIP_T(0x10, smoke_bb[4], v1.z, dv.z,  dv.z); // -z plane
	TEST_CLIP_T(0x20, smoke_bb[5], v1.z, dv.z, -dv.z); // +z plane
	
	if (tmin >= tmax) { // clipped away
		eye  = v1;
		vpos = v1;
	}
	else {
		eye  = v1 + dv*tmax;
		vpos = v1 + dv*tmin;
	}
} 
