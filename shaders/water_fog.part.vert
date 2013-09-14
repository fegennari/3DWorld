uniform float water_plane_z;

void set_fog_coord()
{
	// clip the line to the water plane if the eye is above the water
	// could use different terms/fog scaling/color for inside/outside water?
	vec4 eye = gl_ModelViewMatrixInverse[3]; // world space
	float t = min(1.0, (eye.z - water_plane_z)/max(0.0001, (eye.z - gl_Vertex.z)));
	vec4 clipped_vert = mix(eye, gl_Vertex, ((t < 0.0) ? 1.0 : t));
	gl_FogFragCoord = length((gl_ModelViewMatrix * clipped_vert).xyz);
} 
