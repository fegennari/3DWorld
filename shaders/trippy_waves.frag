uniform sampler2D frame_buffer_tex;
uniform float time      = 0.0;
uniform float intensity = 1.0;

in vec2 tc;

// https://www.shadertoy.com/view/Xt2BDG
mat2 rotate(float angle) {
    float c = cos(angle);
    float s = sin(angle);
    return mat2(c, -s, s, c);
}
void main() {
    vec2 uv = 2.0*tc - 1.0;
    vec3 result = vec3(0,0,0);
    float t = 1.;
    float vtime = 0.05*time;
    float offset = -5. * vtime;
    float base = 100. * length(uv);
    float d = sin(-vtime + 15. * length(uv));
    d *= d * d;
    mat2 rot = rotate(5. * length(uv));
    uv += .5;
    uv = abs(rot * uv);
    
    for (int p = 0; p < 3; p++) {
        result[p] = sin(offset + t * base) - cos(20. * uv.x) - cos(20. * uv.y);
        t += 0.05;
    }
    result *= result;
    result = 1. - result;
    vec3 color   = texture(frame_buffer_tex, tc).rgb;
    fg_FragColor = vec4((color + 0.5*intensity*result*d), 1.0);
}
