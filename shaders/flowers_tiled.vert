uniform vec2 xlate = vec2(0.0);
uniform float x1, y1, dx_inv, dy_inv;
uniform vec4 clip_box1, clip_box2; // {x1 y1 x2 y2} => {x y z w}
uniform sampler2D height_tex;

out vec4 vertex, epos;
out vec3 eye_norm;

void main() {

	set_tc0_from_vert_id();
	vec3 gwdelta = get_grass_wind_delta(fg_Vertex.xyz, 0.5); // wind_amt=0.5
	eye_norm     = fg_NormalMatrix * normalize(normalize(fg_Normal) + gwdelta/height); // eye space, height comes from wind.part
	vertex       = fg_Vertex + vec4(gwdelta, 0.0);
	vertex.z    += texture(height_tex, vec2((fg_Vertex.x - x1)*dx_inv, (fg_Vertex.y - y1)*dy_inv)).r;
#ifdef ENABLE_VERTEX_CLIP
	float vx = vertex.x + xlate.x;
	float vy = vertex.y + xlate.y;
	if ((vx > clip_box1.x && vy > clip_box1.y && vx < clip_box1.z && vy < clip_box1.w) ||
	    (vx > clip_box2.x && vy > clip_box2.y && vx < clip_box2.z && vy < clip_box2.w)) {vertex.z -= height;}
#endif
	epos         = fg_ModelViewMatrix * (vertex + vec4(xlate, 0.0, 0.0));
	gl_Position  = fg_ProjectionMatrix * epos;
	fg_Color_vf  = fg_Color;
}
