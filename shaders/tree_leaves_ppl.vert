uniform float normal_scale = 1.0;
uniform float water_depth  = 0.0;

out vec4 epos;
out vec3 normal; // eye space
out vec3 ws_pos, ws_normal;

void main() {
	set_tc0_blend_from_tc_vert_id();
	set_ws_pos_and_normal(ws_pos, ws_normal); // from world_space_offset_rot.part.vert
	vec4 vpos   = fg_Vertex;
	add_leaf_wind(vpos);
	gl_Position = fg_ProjectionMatrix * (fg_ModelViewMatrix * vpos); // Note: faster than using fg_ModelViewProjectionMatrix (avoids CPU mult+upload)
	fg_Color_vf = fg_Color;
	epos        = fg_ModelViewMatrix * fg_Vertex;
	normal      = normalize(fg_NormalMatrix * fg_Normal * normal_scale); // eye space

	if (underwater) {
		vec3 eye        = fg_ModelViewMatrixInverse[3].xyz;
		gl_FogFragCoord = length(mix(eye, fg_Vertex.xyz, min(1.0, water_depth/max(0.0001, (fg_Vertex.z - eye.z)))) - eye);
	}
	else {
		gl_FogFragCoord = length(epos.xyz);
	}
} 
