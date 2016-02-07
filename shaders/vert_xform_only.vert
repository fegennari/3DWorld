void main() {
#ifdef POS_FROM_EPOS_MULT // needed for model3d z-prepass
	gl_Position = fg_ProjectionMatrix * (fg_ModelViewMatrix * fg_Vertex);
#else
	gl_Position = fg_ftransform();
#endif
	fg_Color_vf = fg_Color;
}
