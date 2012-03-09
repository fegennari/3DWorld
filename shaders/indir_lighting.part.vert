uniform float half_dxy;
uniform float indir_vert_offset = 0.25;
uniform vec3 scene_llc, scene_scale; // scene bounds (world space)
varying vec3 spos; // world space

void setup_indir_lighting(in vec3 normal) {
	if (indir_lighting) {
		spos = gl_Vertex.xyz + (indir_vert_offset*half_dxy)*normal; // move slightly away from the vertex
		spos = clamp((spos - scene_llc)/scene_scale, 0.0, 1.0); // should be in [0.0, 1.0] range
	}
} 
