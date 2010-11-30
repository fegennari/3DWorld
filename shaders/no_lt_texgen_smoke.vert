uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
varying vec3 eye, vpos, orig_vertex;

const int MAX_LIGHTS = 256;
// store light_source as: center.xyz, radius, color.rgba
uniform int num_lights = 0;
uniform sampler2D dlights_tex;
uniform float x_scene_size, y_scene_size, czmin, czmax; // scene bounds (world space)

const float CTHRESH = 0.02;
const float t_scale = 1.0/MAX_LIGHTS;

vec4 get_tex_val(int i, int which) {
	return texture2D(dlights_tex, vec2(0.5*which, i*t_scale));
}

float get_intensity_at(in vec3 pos, in vec3 off, in vec3 scale, in int i) {
	vec4 pos_r = get_tex_val(i,0);
	float radius = pos_r.w;
	if (radius == 0.0) return get_tex_val(i,1).a; // no falloff
	vec3 center = pos_r.xyz*scale + off;
	if (abs(pos.z - center.z) > radius) return 0.0; // fast test
	float dist = length(pos - center);
	if (dist > radius) return 0.0;
	float rscale = (radius - dist)/radius;
	return rscale*rscale*get_tex_val(i,1).a; // quadratic 1/r^2 attenuation
}

void main()
{
	if (use_texgen) {
		setup_texgen(0);
	}
	else {
		gl_TexCoord[0] = gl_MultiTexCoord0;
	}
	orig_vertex = gl_Vertex.xyz;
	gl_Position = ftransform();
	//gl_FrontColor = gl_Color;
	vec4 color = gl_Color;
	
	/*if (num_lights > 0) {
		vec3 llc = vec3(-x_scene_size, -y_scene_size, czmin);
		vec3 urc = vec3( x_scene_size,  y_scene_size, czmax);
		vec3 scale = (urc - llc);
		
		for (int i = 0; i < num_lights; ++i) {
			float cscale = get_intensity_at(gl_Vertex.xyz, llc, scale, i);
			if (cscale < CTHRESH) continue;
			color += vec4(get_tex_val(i,1).rgb, 0.0)*cscale;
		}
	}*/
	gl_FrontColor = color;
	
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
