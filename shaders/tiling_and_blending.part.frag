// from Deliot2019
// https://eheitzresearch.wordpress.com/738-2/

// Textures
uniform sampler2D Tinput; // Gaussian input T(I)
uniform sampler2D invT; // Inverse histogram transformation T^{-1}

// Decorrelated color space vectors and origin
uniform vec3 _colorSpaceVector1;
uniform vec3 _colorSpaceVector2;
uniform vec3 _colorSpaceVector3;
uniform vec3 _colorSpaceOrigin;

vec3 ReturnToOriginalColorSpace(vec3 color) {
	vec3 result = _colorSpaceOrigin + _colorSpaceVector1 * color.r + _colorSpaceVector2 * color.g + _colorSpaceVector3 * color.b;
	return result;
}

// Compute local triangle barycentric coordinates and vertex IDs
void TriangleGrid(vec2 uv, out float w1, out float w2, out float w3, out ivec2 vertex1, out ivec2 vertex2, out ivec2 vertex3) {
	// Scaling of the input
	uv *= 3.464; // 2 * sqrt(3)

	// Skew input space into simplex triangle grid
	const mat2 gridToSkewedGrid = mat2(1.0, 0.0, -0.57735027, 1.15470054);
	vec2 skewedCoord = gridToSkewedGrid * uv;

	// Compute local triangle vertex IDs and local barycentric coordinates
	ivec2 baseId = ivec2(floor(skewedCoord));
	vec3 temp = vec3(fract(skewedCoord), 0);
	temp.z = 1.0 - temp.x - temp.y;
	if (temp.z > 0.0) {
		w1 = temp.z;
		w2 = temp.y;
		w3 = temp.x;
		vertex1 = baseId;
		vertex2 = baseId + ivec2(0, 1);
		vertex3 = baseId + ivec2(1, 0);
	}
	else {
		w1 = -temp.z;
		w2 = 1.0 - temp.y;
		w3 = 1.0 - temp.x;
		vertex1 = baseId + ivec2(1, 1);
		vertex2 = baseId + ivec2(1, 0);
		vertex3 = baseId + ivec2(0, 1);
	}
}

vec2 hash(vec2 p) {
	return fract(sin((p) * mat2(127.1, 311.7, 269.5, 183.3) )*43758.5453);
}

// By-Example procedural noise at uv
vec3 ByExampleProceduralNoise(vec2 uv) {
	// Get triangle info
	float w1, w2, w3;
	ivec2 vertex1, vertex2, vertex3;
	TriangleGrid(uv, w1, w2, w3, vertex1, vertex2, vertex3);
		
	// Assign random offset to each triangle vertex
	vec2 uv1 = uv + hash(vertex1);
	vec2 uv2 = uv + hash(vertex2);
	vec2 uv3 = uv + hash(vertex3);

	// Precompute UV derivatives 
	vec2 duvdx = dFdx(uv);
	vec2 duvdy = dFdy(uv);

	// Fetch Gaussian input
	vec3 G1 = textureGrad(Tinput, uv1, duvdx, duvdy).rgb;
	vec3 G2 = textureGrad(Tinput, uv2, duvdx, duvdy).rgb;
	vec3 G3 = textureGrad(Tinput, uv3, duvdx, duvdy).rgb;

	// Variance-preserving blending
	vec3 G = w1*G1 + w2*G2 + w3*G3;
	G = G - vec3(0.5);
	G = G * inversesqrt(w1*w1 + w2*w2 + w3*w3);
	//if(_DXTScalers.x >= 0.0) G = G * _DXTScalers; // Only with DXT compression (Section 1.6)
	G = G + vec3(0.5);

	// Compute LOD level to fetch the prefiltered look-up table invT
	float LOD = textureQueryLod(Tinput, uv).y / float(textureSize(invT, 0).y);

	// Fetch prefiltered LUT (T^{-1})
	vec3 color;
	color.r	= texture(invT, vec2(G.r, LOD)).r;
	color.g	= texture(invT, vec2(G.g, LOD)).g;
	color.b	= texture(invT, vec2(G.b, LOD)).b;
	
	// Original color space
	color = ReturnToOriginalColorSpace(color);

	return color;
}

