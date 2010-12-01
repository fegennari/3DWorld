uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
varying vec3 eye, vpos;

#extension GL_ARB_uniform_buffer_object : enable

// store light_source as: center.xyz, radius, color.rgba
//                               012  3       4567
const int MAX_LIGHTS = 10;
uniform int num_lights = 0;

layout(std140) uniform uniform_block {
  float data[8*MAX_LIGHTS];
};

const float CTHRESH = 0.02;

float get_intensity_at(in vec3 pos, in int off) {
	float radius = data[off+3];
	if (radius == 0.0) return data[off+7]; // no falloff
	if (abs(pos.z - data[off+2]) > radius) return 0.0; // fast test
	vec3 center = vec3(data[off+0], data[off+1], data[off+2]);
	float dist = length(pos - center);
	if (dist > radius) return 0.0;
	float rscale = (radius - dist)/radius;
	return rscale*rscale*data[off+7]; // quadratic 1/r^2 attenuation
}

void main()
{
	if (use_texgen) {
		setup_texgen(0);
	}
	else {
		gl_TexCoord[0] = gl_MultiTexCoord0;
	}	
	gl_Position = ftransform();
	//gl_FrontColor = gl_Color;
	vec4 color = gl_Color;
	
	for (int i = 0; i < num_lights; ++i) {
		int off = 8*i;
		float cscale = get_intensity_at(gl_Vertex.xyz, off);
		if (cscale < CTHRESH) continue;
		color += vec4(data[off+4], data[off+5], data[off+6], 0.0)*cscale;
	}
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
