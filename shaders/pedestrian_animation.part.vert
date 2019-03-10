const float PI = 3.14159;

uniform float animation_time     = 0.0;
uniform float animation_scale    = 1.0;
uniform float model_delta_height = 0.0;
uniform int animation_id = 0;

mat4 do_rotation(vec3 axis, float angle) {
	vec3 a = normalize(axis);
	float s = sin(angle);
	float c = cos(angle);
	float oc = 1.0 - c;
	float sx = s * a.x;
	float sy = s * a.y;
	float sz = s * a.z;
	float ocx = oc * a.x;
	float ocy = oc * a.y;
	float ocz = oc * a.z;
	float ocxx = ocx * a.x;
	float ocxy = ocx * a.y;
	float ocxz = ocx * a.z;
	float ocyy = ocy * a.y;
	float ocyz = ocy * a.z;
	float oczz = ocz * a.z;
	return mat4(vec4(ocxx+c, ocxy-sz, ocxz+sy, 0.0), vec4(ocxy+sz, ocyy+c, ocyz-sx, 0.0), vec4(ocxz-sy, ocyz+sx, oczz+c, 0.0), vec4(0.0, 0.0, 0.0, 1.0));
}

void apply_vertex_animation(inout vec4 vertex, inout vec3 normal, in vec2 tc) {
	if (animation_id == 0 || animation_time == 0.0) return; // animation disabled
	float anim_scale = 0.01*animation_scale;
	float anim_val   = 150.0*animation_time;

	if (animation_id == 1) { // bunny hop
		vertex.y += 1.0*anim_scale*abs(sin(anim_val));
	}
	else if (animation_id == 2) { // flip
		float v = fract(0.2*anim_val);
		if (v < 0.5) {
			float rot_height = 2.25*anim_scale + model_delta_height;
			vertex.y += 0.8*anim_scale*sin(2.0*PI*v);
			vertex.y -= rot_height; // rotate about this point
			mat4 m    = do_rotation(vec3(1,0,0), 2.0*PI*smoothstep(0.0, 1.0, 2.0*v));
			vertex   *= m;
			normal    = normalize(normal * mat3(m)); // is this okay, or do we need inverse transform?
			vertex.y += rot_height;
		}
	}
	else if (animation_id == 3) { // twirl
		float v = fract(0.2*anim_val);
		if (v < 0.5) {
			vertex.y += 0.5*anim_scale*sin(2.0*PI*v);
			mat4 m    = do_rotation(vec3(0,1,0), 2.0*PI*smoothstep(0.0, 1.0, 2.0*v));
			vertex   *= m;
			normal    = normalize(normal * mat3(m));
		}
	}
	else if (animation_id == 4) { // march
		if (vertex.y < 0.7*anim_scale + model_delta_height) {
			float v   = sin(2.5*anim_val);
			vertex.y += 0.4*anim_scale*max(((vertex.x > 0.0) ? v : -v), 0.0);
		}
	}
	else if (animation_id == 5 || animation_id == 6) { // walk
		float rot_height = 1.3*anim_scale + model_delta_height*((animation_id == 5) ? -1.0 : 1.0);
		if (vertex.y < rot_height) {
			float v   = sin(2.5*anim_val);
			vertex.y -= rot_height; // rotate about this point
			mat4 m    = do_rotation(vec3(((vertex.x > 0.0) ? 1.0 : -1.0),0,0), ((animation_id == 5) ? 1.0 : 0.3)*v);
			normal    = normalize(normal * mat3(m));
			vertex   *= m;
			vertex.y += rot_height;
		}
	}
}
