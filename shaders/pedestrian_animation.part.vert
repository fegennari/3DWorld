const float PI = 3.14159;

uniform float animation_time     = 0.0;
uniform float animation_scale    = 1.0;
uniform float model_delta_height = 0.0;
uniform int animation_id = 0;

mat3 do_rotation(vec3 a, float angle) { // Note: axis does not need to be normalized
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
	return mat3(vec3(ocxx+c, ocxy-sz, ocxz+sy), vec3(ocxy+sz, ocyy+c, ocyz-sx), vec3(ocxz-sy, ocyz+sx, oczz+c));
}

void rotate_about(inout vec3 vertex, inout vec3 normal, in float yval, in vec3 axis, in float angle) {
	vertex.y -= yval; // rotate about this point
	mat3 m    = do_rotation(axis, angle);
	vertex   *= m;
	normal   *= m; // rotate only, no scale - no inverse transpose needed
	vertex.y += yval;
}

void apply_vertex_animation(inout vec4 vertex, inout vec3 normal) {
	if (animation_id == 0 || animation_time == 0.0) return; // animation disabled
	float anim_scale = 0.01*animation_scale;
	float anim_val   = 150.0*animation_time;

	if (animation_id == 2) { // bunny hop
		vertex.y += 1.0*anim_scale*abs(sin(anim_val));
	}
	else if (animation_id == 3) { // flip
		float v = fract(0.2*anim_val);
		if (v < 0.5) {
			float rot_height = 2.25*anim_scale + model_delta_height;
			vertex.y += 0.8*anim_scale*sin(2.0*PI*v);
			rotate_about(vertex.xyz, normal, rot_height, vec3(1,0,0), 2.0*PI*smoothstep(0.0, 1.0, 2.0*v));
		}
	}
	else if (animation_id == 4) { // twirl
		float v = fract(0.2*anim_val);
		if (v < 0.5) {
			vertex.y += 0.5*anim_scale*sin(2.0*PI*v);
			rotate_about(vertex.xyz, normal, 0.0, vec3(0,1,0), 2.0*PI*smoothstep(0.0, 1.0, 2.0*v));
		}
	}
	else if (animation_id == 5) { // march
		if (vertex.y < 0.7*anim_scale + model_delta_height) {
			float v   = sin(2.5*anim_val);
			vertex.y += 0.4*anim_scale*max(((vertex.x > 0.0) ? v : -v), 0.0);
		}
	}
	else if (animation_id == 1 || animation_id == 6) { // walk
		float hip_height = 1.25*anim_scale + model_delta_height*((animation_id == 6) ? -1.0 : 1.0);
		float hip_scale  = 1.0 - clamp((10.0*(vertex.y - hip_height) + 0.5), 0.0, 1.0);
		
		if (hip_scale > 0.0) {
			float time = 2.0*anim_val;
			float v    = hip_scale * sin(time);
			vec3 axis  = vec3(((vertex.x > 0.0) ? 1.0 : -1.0),0,0);
			
			if (animation_id == 1) { // normal walk
				float dv = cos(time); // derivative of v

				if (dv*axis.x < 0.0) { // moving forward, bend the knee
					float knee_height = 0.55*hip_height;
					if (vertex.y < knee_height) {rotate_about(vertex.xyz, normal, knee_height, axis, -0.8*dv);} // knee joint
				}
				if (v*axis.x > 0.0) {v *= 0.5;} // leg in back, half the angle
				rotate_about(vertex.xyz, normal, hip_height, axis, 0.5*v); // hip joint
			}
			else {
				rotate_about(vertex.xyz, normal, hip_height, axis, v); // alien
			}
		}
	}
	else if (animation_id == 7) { // rats
		// y = up/down, x = left/right, z = front/back; x and z are centered around 0, y is about [0, height]
		float height = 14.0*model_delta_height;

		if (vertex.y < 0.14*height && vertex.z > -1.2*height) { // legs
			float lr_sign = ((vertex.x < 0.0) ? -1.0 : 1.0); // left/right
			float fb_sign = ((vertex.z < 0.0) ? -1.0 : 1.0); // front back
			vertex.y += max(0.0, 0.05*anim_scale*abs(sin(5.0*anim_val + 0.25*PI*lr_sign*fb_sign)));
		}
	}
}
