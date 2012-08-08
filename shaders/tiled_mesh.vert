uniform float water_plane_z;
varying vec3 vertex, epos;

void main()
{
	setup_texgen(0);
	setup_texgen(1);
	vertex = gl_Vertex.xyz;
	epos   = (gl_ModelViewMatrix * gl_Vertex).xyz;
	gl_Position = ftransform();

	// calculate fog coord
	// clip the line to the water plane if the eye is above the water
	// could use different terms/fog scaling/color for inside/outside water?
	vec4 eye = gl_ModelViewMatrixInverse[3]; // world space
	float t = min(1.0, (eye.z - water_plane_z)/max(0.0, (eye.z - gl_Vertex.z)));
	vec4 clipped_vert = mix(eye, gl_Vertex, ((t < 0.0) ? 1.0 : t));
	gl_FogFragCoord = length((gl_ModelViewMatrix * clipped_vert).xyz);
} 
