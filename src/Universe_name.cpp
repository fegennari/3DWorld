// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 4/27/06
#include "universe.h"


string const v_com = "a e i o ";
string const v_str = "u ai io ";
string const v_beg = "au ea ei eo eu ou ya ye yo ";
string const v_mid = "ao au ea ee ei eo eu ia ie oa oi oo ou ue ";
string const v_end = "ay ee ey ia ie oo oy ue ";
string const c_com = "l n r s t ";
string const c_str = "b b c c d d f f g g h h j k m m p p v w w x ch ch sp st st th th ";
string const c_beg = "bl br cl cr dr fl fr gl gr ph pl pr sk sh sh tr tr wh qu str thr ";
string const c_mid = "bb dd ff gg pp rr rr tt tt rc ";
string const c_end = "ck ck ct gh ld ld ll ln ln lm lp lt mp nc nd nd ng nk rk rs rs rt rt ss ss ion gth nch ";

string const v_start  = v_com + v_com + v_com + v_com + v_com + v_str + v_beg;
string const v_middle = v_com + v_com + v_com + v_com + v_com + v_str + v_mid;
string const v_ending = v_com + v_com + v_com + v_com + v_com + v_str + v_end;
string const c_start  = c_com + c_com + c_com + c_com + c_str + c_beg;
string const c_middle = c_com + c_com + c_com + c_com + c_str + c_beg + c_mid + c_end;
string const c_ending = c_com + c_com + c_com + c_com + c_str + c_end;

vector<string> n_start[2], n_middle[2], n_ending[2];



void parse_str_list(string const &str, vector<string> &vs) { // str must end with a space

	string cur;

	for (unsigned i = 0; i < str.size(); ++i) {
		if (str[i] == ' ') {
			assert(!cur.empty());
			vs.push_back(cur);
			cur.clear();
		}
		else {
			cur.push_back(str[i]);
		}
	}
	assert(cur.empty());
}


void parse_str_tables() {

	static bool parsed(0);
	if (parsed) return;
	parsed = 1;
	parse_str_list(v_start,  n_start [0]);
	parse_str_list(v_middle, n_middle[0]);
	parse_str_list(v_ending, n_ending[0]);
	parse_str_list(c_start,  n_start [1]);
	parse_str_list(c_middle, n_middle[1]);
	parse_str_list(c_ending, n_ending[1]);
}


void named_obj::gen_name(s_object const &sobj) {

	parse_str_tables();
	bool at_end(0), vc((rand2() % 100) < 90); // 0 = vowel, 1 = consonant
	name = "";

	for (unsigned i = 0; !at_end; ++i, vc ^= 1) {
		at_end = (i >= 5 || name.size() >= 8 || (i > 1 && (rand2() % 100) < int((vc ? 10 : 5) + 16*i)));
		vector<string> const &str((i == 0) ? n_start[vc] : (at_end ? n_ending[vc] : n_middle[vc]));
		assert(!str.empty());
		name += str[rand2() % str.size()];
	}
	assert(!name.empty());
	name[0] += ('A' - 'a'); // make uppercase
	lookup_given_name(sobj); // already named, overwrite the old value (but need to preserve random number generator state)
	//cout << name << "  ";
}




