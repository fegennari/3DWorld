#ifdef USE_BUMP_MAP
attribute vec3 tangent;

varying vec3 tangent_v, binorm_v;
varying vec3 ts_pos; // tangent space pos

void setup_tbn(in vec3 eye_norm, in vec3 eye_pos)
{
	// Building the matrix Eye Space -> Tangent Space
	tangent_v = normalize(gl_NormalMatrix * tangent);
	binorm_v  = cross(eye_norm, tangent_v);

	ts_pos.x  = dot(eye_pos, tangent_v);
	ts_pos.y  = dot(eye_pos, binorm_v);
	ts_pos.z  = dot(eye_pos, eye_norm);
}
#endif
