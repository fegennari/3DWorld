const int SINES_PER_FREQ  = 12;
const int MAX_FREQ_BINS   = 5;
const int TOT_NUM_SINES   = SINES_PER_FREQ*MAX_FREQ_BINS;
const int NUM_SINE_PARAMS = 2*3+1; // 2*NUM_DIMENSIONS+1
const int SINE_DATA_SIZE  = NUM_SINE_PARAMS*TOT_NUM_SINES;
uniform float rdata[SINE_DATA_SIZE]; // too large for a uniform?

float procedural_eval(in vec3 pos) { // equivalent to noise_gen_3d::get_val()
	float val = 0.0;
	
	for (int k = 0; k < TOT_NUM_SINES; ++k) {
		int index2 = NUM_SINE_PARAMS*k;
		float x = sin(rdata[index2+1]*pos.x + rdata[index2+2]);
		float y = sin(rdata[index2+3]*pos.y + rdata[index2+4]);
		float z = sin(rdata[index2+5]*pos.z + rdata[index2+6]);
		val += rdata[index2]*x*y*z;
	}
	return val;
}
