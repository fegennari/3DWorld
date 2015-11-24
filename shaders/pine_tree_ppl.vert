#ifdef ENABLE_INSTANCING
in vec3 xlate;
uniform float vertex_scale = 1.0;
#endif

out float world_space_zval;
out vec4 epos;
out vec3 normal; // eye space

void main() {
	set_tc0_from_vert_id();
	vec4 vertex = fg_Vertex;
#ifdef ENABLE_INSTANCING
	vertex.xyz *= vertex_scale;
	vertex.xyz += xlate;
#endif
	add_leaf_wind(vertex);
	world_space_zval = vertex.z;
	epos             = fg_ModelViewMatrix  * vertex;
	gl_Position      = fg_ProjectionMatrix * epos;
	gl_FogFragCoord  = length(epos.xyz); // set standard fog coord
	normal           = normalize(fg_NormalMatrix * fg_Normal);
	normal          *= -sign(dot(normal, epos.xyz)); // two-sided lighting
	fg_Color_vf      = fg_Color;
}
