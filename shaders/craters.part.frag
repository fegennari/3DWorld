// add craters by modifying the normal
// assumes the object center is at (0,0,0) in world space
// normal is in eye space
// vertex is in world space (or could be the normal)
void adjust_normal_for_craters(inout vec3 norm, in vec3 vertex) {
	
	float v0 = 1.0; // using a variable here is slow
	vec3 dir = normalize(vertex); // world space normal

	for (int i = 0; i < 50; ++i) { // Note: inefficient, but fast enough for a single object render
		vec3 center = rand_vec3(v0);
		vec3 dir2   = dir - normalize(center);
		float dist  = length(dir2);
		float rad1  = 0.07*(0.25 + 0.75*rand_01(v0+3.0));
		float rad2  = 1.5*rad1;
		v0         += 4.0;
		
		if (dist < rad2) { // at crater (parabola)
			vec3 cnorm = normalize(fg_NormalMatrix * dir2/dist);
			float cwt;

			if (dist < rad1) { // inside crater
				cwt  = 0.75*dist/rad1; // higher power?
				cnorm = -cnorm;
			}
			else { // on rim of crater
				cwt  = 0.5*sqrt(1.0 - (dist - rad1)/(rad2 - rad1));
			}
			norm = normalize((1.0 - cwt)*norm + cwt*cnorm);
		}
	}
}

