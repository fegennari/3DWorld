// input: point.xyz, height, color
// output: 30 pine tree quads
uniform float mesh_scale = 1.0;
varying out vec3 normal;

const int N_PT_LEVELS = 6;
const int N_PT_RINGS  = 5;
const float rd = 0.5;
const float PI = 3.14159;

void main()
{
	gl_FrontColor = gl_FrontColorIn[0]; // all colors are the same
	vec4 pos      = gl_PositionIn[0]; // center of the tree
	float height  = gl_TexCoordIn[0][7].s;
	float SQRT2   = sqrt(2.0);
	float theta0  = (int(1.0E6*height)%360)*(PI/180.0); // pseudo-random rotation

	for (int j = 0; j < N_PT_LEVELS; ++j) {
		float sz = 0.5*(height + 0.03/mesh_scale)*((N_PT_LEVELS - j - 0.4)/float(N_PT_LEVELS));
		float z = (j + 1.8)*height/(N_PT_LEVELS + 2.8) - rd*sz;

		for (int k = 0; k < N_PT_RINGS; ++k) {
			float theta = 2.0*PI*(3.3*j + k/float(N_PT_RINGS)) + theta0;
			float tv = theta - 0.25*PI;
			float s = sin(tv)/SQRT2;
			float c = cos(tv)/SQRT2;
			vec4 pts[4];

			for (int i = 0; i < 4; ++i) {
				float tv2 = theta - 0.5*PI*i;
				pts[i] = pos + vec4(sz*(sin(tv2)+s), sz*(cos(tv2)+c), rd*sz*(-1.0+2.0*(i>>1))+z, 0.0);
			}
			normal = normalize(gl_NormalMatrix * cross((pts[1] - pts[0]).xyz, (pts[2] - pts[1]).xyz));
			gl_Position = gl_ModelViewProjectionMatrix * pts[0]; gl_TexCoord[0].st = vec2(0.0, 1.0); EmitVertex();
			gl_Position = gl_ModelViewProjectionMatrix * pts[1]; gl_TexCoord[0].st = vec2(0.0, 0.0); EmitVertex();
			gl_Position = gl_ModelViewProjectionMatrix * pts[2]; gl_TexCoord[0].st = vec2(1.0, 1.0); EmitVertex();
			gl_Position = gl_ModelViewProjectionMatrix * pts[3]; gl_TexCoord[0].st = vec2(1.0, 0.0); EmitVertex();
			EndPrimitive();
		}
	}
}