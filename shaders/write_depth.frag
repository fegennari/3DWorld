uniform float far_clip = 100.0; // far clipping plane
in vec4 pos;

void main()
{
	fg_FragColor = vec4(pos.z/far_clip, 0, 0, 1); // write depth (zval) into red channel
}
