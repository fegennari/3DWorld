uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
uniform float step_delta;
attribute float shadowed;
varying vec3 eye, vpos, spos;

const int MAX_LIGHTS = 10;
attribute float num_lights;
//uniform float data[8*MAX_LIGHTS];

/*float get_intensity_at(in vec3 pos, in int off) {
	float radius = data[off+3];
	if (radius == 0.0) return data[off+7]; // no falloff
	if (abs(pos.z - data[off+2]) > radius) return 0.0; // fast test
	vec3 center = vec3(data[off+0], data[off+1], data[off+2]);
	float dist = length(pos - center);
	if (dist > radius) return 0.0;
	float rscale = (radius - dist)/radius;
	return rscale*rscale*data[off+7]; // quadratic 1/r^2 attenuation
}*/

attribute vec4 dlp0, dlc0;
attribute vec4 dlp1, dlc1;
attribute vec4 dlp2, dlc2;
attribute vec4 dlp3, dlc3;
attribute vec4 dlp4, dlc4;

vec4 add_color(in vec3 pos, in vec4 lp, in vec4 lc) {
	float weight;
	float radius = 0.5*lp.w;
	
	if (radius == 0.0) {
		weight = 1.0; // no falloff
	}
	else {
		vec3 center = lp.xyz;
		float dist = length(pos - center);
		
		if (dist > radius) {
			weight = 0.0;
		}
		else {
			float rscale = (radius - dist)/radius;
			weight = rscale*rscale; // quadratic 1/r^2 attenuation
		}
	}
	return vec4(lc.rgb, 0.0)*weight*lc.a;
}

#define ADD_LIGHT(en, N) if (en && (s&(1<<N)) == 0) color.rgb += add_light_comp(normal, N).rgb;

void main()
{
	if (use_texgen) {
		setup_texgen(0);
	}
	else {
		gl_TexCoord[0] = gl_MultiTexCoord0;
	}	
	gl_Position = ftransform();
	vec3 normal = normalize(gl_NormalMatrix * gl_Normal); // eye space
	vec4 color  = gl_Color;
	int s = int(round(shadowed));
	ADD_LIGHT(enable_light0, 0);
	ADD_LIGHT(enable_light1, 1);
	
	int nl = int(round(num_lights));
	/*for (int i = 0; i < nl; ++i) {
		int off = 8*i;
		float cscale = get_intensity_at(gl_Vertex.xyz, off);
		color += vec4(data[off+4], data[off+5], data[off+6], 0.0)*cscale;
	}*/
	if (nl > 0) color += add_color(gl_Vertex.xyz, dlp0, dlc0);
	if (nl > 1) color += add_color(gl_Vertex.xyz, dlp1, dlc1);
	if (nl > 2) color += add_color(gl_Vertex.xyz, dlp2, dlc2);
	if (nl > 3) color += add_color(gl_Vertex.xyz, dlp3, dlc3);
	if (nl > 4) color += add_color(gl_Vertex.xyz, dlp4, dlc4);
	gl_FrontColor = clamp(color, 0.0, 1.0);
	
	spos = gl_Vertex.xyz + (0.25*step_delta)*normalize(gl_Normal); // move slightly away from the vertex
	
	if (!smoke_enabled) {
		eye  = vec3(0,0,0);
		vpos = vec3(0,0,0);
		return;
	}
	vec3 v2 = (inverse(gl_ModelViewMatrix) * vec4(0.0, 0.0, 0.0, 1.0)).xyz; // world space
	pt_pair res = clip_line(gl_Vertex.xyz, v2, smoke_bb);
	eye  = res.v1;
	vpos = res.v2;
} 
