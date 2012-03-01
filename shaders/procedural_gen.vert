uniform float half_dxy;
uniform float indir_vert_offset = 0.25;

varying vec3 eye, vpos, spos, normal; // world space

void main()
{
	// *** texture stuff here ***
	gl_Position   = ftransform();
	gl_FrontColor = gl_Color;
	normal   = normalize(gl_Normal);
	eye_norm = normalize(gl_NormalMatrix * normal);
	epos     = gl_ModelViewMatrix * gl_Vertex;
	vpos     = gl_Vertex.xyz;
	spos     = gl_Vertex.xyz + (indir_vert_offset*half_dxy)*normal; // move slightly away from the vertex
	eye      = gl_ModelViewMatrixInverse[3].xyz; // world space
	set_fog(); // set standard fog coord
} 
