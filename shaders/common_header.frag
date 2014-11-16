layout(location = 0) out vec4 fg_FragColor;
in vec4 fg_Color_vf;
#define gl_Color fg_Color_vf

in float fg_FogFragCoord; // not always used
#define gl_FogFragCoord fg_FogFragCoord
