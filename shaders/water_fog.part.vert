uniform float water_plane_z;

void set_fog_coord(in vec4 vertex, in vec3 eye)
{
	// clip the line to the water plane if the eye is above the water
	// could use different terms/fog scaling/color for inside/outside water?
	float t   = min(1.0, (eye.z - water_plane_z)/max(0.0001, (eye.z - vertex.z)));
	vec3 vert = mix(eye, vertex.xyz, ((t < 0.0) ? 1.0 : t));
	gl_FogFragCoord = length(vert - eye);
}
