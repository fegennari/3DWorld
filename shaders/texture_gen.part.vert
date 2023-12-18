uniform float tex_coord_weight = 0.0;
uniform vec4 texgen_s = vec4(1.0, 0.0, 0.0, 0.0);
uniform vec4 texgen_t = vec4(0.0, 1.0, 0.0, 0.0);
uniform vec3 texgen_origin = vec3(0.0);
out vec2 tc;

void setup_texgen_st() {
	vec4 vertex = fg_Vertex - vec4(texgen_origin, 0.0);
	tc = vec2(dot(vertex, texgen_s), dot(vertex, texgen_t));
}

// this variant works better for cylinders as the X and Y values won't cancel at 45 degree edges
void setup_texgen_st_no_xy_cancel() {
	vec4 vertex = fg_Vertex - vec4(texgen_origin, 0.0);
	vertex.x =  abs(vertex.x);
	vertex.y = -abs(vertex.y);
	tc = vec2(dot(vertex, texgen_s), dot(vertex, texgen_t));
}

uniform int tc_start_ix = 0;

void set_tc0_from_vert_id() { // 0,0 1,0 1,1 0,1
	int tix = (gl_VertexID + tc_start_ix) & 3;
	if      (tix == 0) {tc = vec2(0,0);}
	else if (tix == 1) {tc = vec2(1,0);}
	else if (tix == 2) {tc = vec2(1,1);}
	else               {tc = vec2(0,1);}
}

void set_bent_quad_tc0_from_vert_id() {
	int tix = (gl_VertexID + tc_start_ix) & 7;
	if      (tix == 0) {tc = vec2(.5,1);}
	else if (tix == 1) {tc = vec2( 0,1);}
	else if (tix == 2) {tc = vec2( 0,0);}
	else if (tix == 3) {tc = vec2(.5,0);}
	else if (tix == 4) {tc = vec2( 1,1);}
	else if (tix == 5) {tc = vec2(.5,1);}
	else if (tix == 6) {tc = vec2(.5,0);}
	else               {tc = vec2( 1,0);}
}

void set_tc0_blend_from_tc_vert_id() {
#ifdef ENABLE_TEX_COORD_WEIGHT
	if (tex_coord_weight == 0.0) {set_tc0_from_vert_id();} else {tc = tex_coord_weight*fg_TexCoord;}
#else
	set_tc0_from_vert_id();
#endif
}
