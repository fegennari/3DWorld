// simple toon fragment shader
// www.lighthouse3d.com

varying vec3 normal, lightDir;

void main()
{
	float intensity;
	vec3 n;
	vec4 color;

	n = normalize(normal);
	intensity = max(dot(lightDir,n),0.0); 

	if (intensity > 0.98)
		color = vec4(0.8,0.8,0.8,1.0);
	else if (intensity > 0.5)
		color = vec4(0.4,0.4,0.8,1.0);	
	else if (intensity > 0.25)
		color = vec4(0.2,0.2,0.4,1.0);
	else
		color = vec4(0.1,0.1,0.1,1.0);		
		
	gl_FragColor = color;
}
