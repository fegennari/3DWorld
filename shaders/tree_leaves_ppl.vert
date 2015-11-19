uniform float normal_scale = 1.0;
uniform float water_depth = 0.0;
uniform vec3 world_space_offset = vec3(0.0);

out vec4 epos;
out vec3 normal; // eye space
out vec3 ws_pos;
out vec3 ws_normal;

void main() {
	set_tc0_from_vert_id();
	vec4 vpos   = fg_Vertex;
	add_leaf_wind(vpos);
	gl_Position = fg_ModelViewProjectionMatrix * vpos;
	fg_Color_vf = fg_Color;

	ws_pos    = fg_Vertex.xyz + world_space_offset;
	epos      = fg_ModelViewMatrix * fg_Vertex;
	ws_normal = normalize(fg_Normal);
	normal    = normalize(fg_NormalMatrix * ws_normal * normal_scale); // eye space
	if (dot(normal, epos.xyz) > 0.0) {normal = -normal;} // facing away from the eye, so reverse (could use faceforward())

	if (underwater) {
		vec3 eye        = fg_ModelViewMatrixInverse[3].xyz;
		gl_FogFragCoord = length(mix(eye, fg_Vertex.xyz, min(1.0, water_depth/max(0.0001, (fg_Vertex.z - eye.z)))) - eye);
	}
	else {
		gl_FogFragCoord = length(epos.xyz);
	}
} 
