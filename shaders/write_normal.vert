varying vec3 normal;

void main()
{
	gl_Position = fg_ftransform();
	normal      = fg_Normal; // world space (not normalized)
} 
