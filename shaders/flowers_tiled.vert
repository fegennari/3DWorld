uniform vec2 xlate = vec2(0.0);
uniform float x1, y1, dx_inv, dy_inv;
uniform vec4 clip_box1, clip_box2; // {x1 y1 x2 y2} => {x y z w}
uniform sampler2D height_tex;
//uniform float height = 1.0; // from wind.part

out vec4 vertex, epos;
out vec3 eye_norm;

void main() {

	set_tc0_from_vert_id();
	vec3 gwdelta = get_grass_wind_delta(fg_Vertex.xyz, 0.5); // wind_amt=0.5
	eye_norm     = fg_NormalMatrix * normalize(normalize(fg_Normal) + gwdelta/height); // eye space
	vertex       = fg_Vertex + vec4(gwdelta, 0.0);
	vertex.z    += texture(height_tex, vec2((fg_Vertex.x - x1)*dx_inv, (fg_Vertex.y - y1)*dy_inv)).r;
#ifdef ENABLE_VERTEX_CLIP
	vec2 v = fg_Vertex.xy + xlate.xy; // Note: using original vertex rather than wind modified 'vertex' avoids flowers stretching in the wind when near a building
	if ((v.x > clip_box1.x && v.y > clip_box1.y && v.x < clip_box1.z && v.y < clip_box1.w) ||
	    (v.x > clip_box2.x && v.y > clip_box2.y && v.x < clip_box2.z && v.y < clip_box2.w)) {vertex.z -= height;}
#endif
	epos         = fg_ModelViewMatrix * (vertex + vec4(xlate, 0.0, 0.0));
	gl_Position  = fg_ProjectionMatrix * epos;
	fg_Color_vf  = fg_Color;
}
