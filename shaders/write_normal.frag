in vec3 normal;

void main()
{
	fg_FragColor = vec4(0.5*(normal + 1.0), 1.0); // Note: normal not normalized
}
