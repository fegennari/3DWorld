varying vec4 pos;
uniform float far_clip = 100.0; // far clipping plane

void main()
{
	gl_FragColor = vec4(pos.z/far_clip, 0, 0, 1); // write depth (zval) into red channel
}
