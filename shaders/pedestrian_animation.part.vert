const float PI = 3.14159;

uniform float animation_time     = 0.0;
uniform float animation_scale    = 1.0;
uniform float model_delta_height = 0.0;
uniform float rotate_amount      = 1.0; // for swings
uniform int   animation_id       = 0;
// 0=none, 1=walk, 2=bunny hop, 3=flip, 4=twirl, 5=march, 6=alien walk, 7=rat, 8=spider, 9=bones animation, 10=helicopter rotate, 11=swingset

#ifdef USE_BONE_ANIMATIONS
layout(location = 4) in uvec4 bone_ids;
layout(location = 5) in vec4 bone_weights;
const int MAX_MODEL_BONES = 200; // must agree with the value used in model3d::setup_bone_transforms()
uniform mat4 bones[MAX_MODEL_BONES];
#endif // USE_BONE_ANIMATIONS

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

// Note: this is no longer just for pedestrians, it now has cases for rats and spiders
void apply_vertex_animation(inout vec4 vertex, inout vec3 normal, in vec2 tc) {

#ifdef USE_BONE_ANIMATIONS
	if (animation_id == 9) { // bone animation
		mat4 bone_transform = bones[bone_ids[0]] * bone_weights[0];
		bone_transform     += bones[bone_ids[1]] * bone_weights[1];
		bone_transform     += bones[bone_ids[2]] * bone_weights[2];
		bone_transform     += bones[bone_ids[3]] * bone_weights[3];
		vertex = bone_transform * vertex;
		normal = inverse(transpose(mat3(bone_transform))) * normal;
		return;
	}
#endif // USE_BONE_ANIMATIONS
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
		float time       = 2.0*anim_val;
		
		if (hip_scale > 0.0) {
			float v   = hip_scale * sin(time);
			vec3 axis = vec3(((vertex.x > 0.0) ? 1.0 : -1.0),0,0);
			
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
		vertex.y += 0.02*anim_scale*sin(2.0*time); // small amount of bob up and down when walking
	}
	else if (animation_id == 7) { // rats
		// y = up/down, x = left/right, z = front/back; x and z are centered around 0, y is about [0, height]
		float height  = 14.0*model_delta_height;
		float leg_top = 0.4*height;

		if (vertex.y < leg_top && vertex.z > -1.2*height && abs(vertex.x) > ((vertex.z < 0.0) ? 0.25 : 0.15)*height) { // legs
			float move_amt  = 0.012*anim_scale*(leg_top - vertex.y)/model_delta_height;
			float lr_sign   = ((vertex.x < 0.0) ? -1.0 : 1.0); // left/right
			float fb_sign   = ((vertex.z < 0.0) ? -1.0 : 1.0); // front back
			float cycle_pos = 0.06*anim_val + 0.5*PI*lr_sign*fb_sign;
			vertex.y += max(0.0, 0.5*move_amt*(cos(cycle_pos) - 0.05)); // up and down; always positive
			vertex.z += move_amt*sin(cycle_pos); // forward and backward
		}
	}
	else if (animation_id == 8) { // spiders
		// tc.x is the position of the leg from front to back: {0.0, 0.25, 0.5, 0.75}
		// tc.y is the joint index: negative for left, positive for right; 0.0 for body joint, 0.3 for knee, 0.7 for ankle, 1.0 for foot
		// first and thrid legs move together; second and fourth legs move together
		float leg_dist  = abs(tc.y);
		float lr_off    = ((tc.y < 0.0) ? 0.0 : 1.0);
		float cycle_pos = 0.04*anim_val + 2.0*PI*(2.0*tc.x + 0.5*lr_off);
		float up_amt    = max(sin(-cycle_pos), 0.0);
		vertex.x += 0.25*animation_scale*leg_dist*cos(cycle_pos); // move forward and backward
		vertex.z += ((leg_dist > 0.1) ? 1.2*animation_scale*(leg_dist - 0.5)*up_amt : 0.0); // knee moves down, angle and foot move up, body joint stays in place
		vertex.y += 0.50*animation_scale*tc.y*up_amt;
	}
	else if (animation_id == 10) { // helicopter rotate
		if (vertex.y > 0.95*anim_scale) {
			mat3 m      = do_rotation(vec3(0.0, 1.0, 0.0), anim_val);
			vertex.xyz *= m;
			normal     *= m; // rotate only, no scale - no inverse transpose needed
		}
	}
	else if (animation_id == 11) { // swingset
		// y = up/down, x = left/right, z = front/back
		float pivot_y = 9.12*model_delta_height;

		if (vertex.y < pivot_y && abs(vertex.x) < 5.0*model_delta_height) { // swings part
			float angle = 0.5*rotate_amount*((vertex.x < 0.0) ? -1.0 : 1.0)*sin(anim_val); // each side is 180 degrees out of phase
			rotate_about(vertex.xyz, normal, pivot_y, vec3(1.0, 0.0, 0.0), angle);
		}
	}
	// else error/skip
}
