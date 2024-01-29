uniform float time = 0.0;

in vec2 tc;

// https://www.shadertoy.com/view/sdKyWc
float rand(vec2 co){
    return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}

float q(float a, vec2 b, vec2 c, vec2 xy) {
	vec2 qq = - (xy - b) * (xy - b)/(4.0 * c * c);
	return a * exp(qq.x + qq.y);
}

vec3 colormap(float q, float depth) {
    return 
        vec3(1.0, depth * 0.5 + 0.5, depth * 0.3 + 0.5)
      + vec3(0.2, 0.7, 0.5) * smoothstep(vec3(0.4, 0.69, 0.7), vec3(0.69, 1.6, 5.5), vec3(q*depth))
      + vec3(0.6, -0.6, -0.9) * smoothstep(0.69, 0.7, q);
}

float blob(float a, vec2 center, float size, vec2 xy, vec2 range, vec2 speed) {
    vec2 b = center + range * cos( rand(center) + time / speed);
    vec2 c = mix(vec2(2.0 * size, 0.0), vec2(0.0, 2.0 * size), (0.5 - 0.08 * sin(2.0 * time/(speed.y))));
	return q(a, b, c, xy);
}

void main() {
    vec2 uv = tc - vec2(0.5);
    float fq = 0.0;
    for (int i = 0; i < 25; i += 1) {
        fq += blob(
            mix(0.9, 1.0 ,rand(vec2(i, 1))),
            vec2(mix(-0.5, 0.5, rand(vec2(i, 2))), 0.0),
            mix(0.03, 0.07, rand(vec2(i, 3))),
            uv,
            vec2(mix(0.0, 0.6, rand(vec2(i, 4))), 0.5),
            vec2(mix(10.0, 22.0, rand(vec2(i, 5))), mix(5.0, 7.0, rand(vec2(i, 6))))
        );
    }
    fg_FragColor = gl_Color * vec4(colormap(fq, 1. - tc.y), 1.0);
}
