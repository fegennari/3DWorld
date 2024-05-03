// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 4/27/06
#include "universe.h"

class name_gen_t {
	vector<string> n_start[2], n_middle[2], n_ending[2];

	static void parse_str_list(string const &str, vector<string> &vs) { // str must end with a space
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
public:
	void init(string const &v_start, string const &v_middle, string const &v_ending, string const &c_start, string const &c_middle, string const &c_ending) {
		parse_str_list(v_start,  n_start [0]);
		parse_str_list(v_middle, n_middle[0]);
		parse_str_list(v_ending, n_ending[0]);
		parse_str_list(c_start,  n_start [1]);
		parse_str_list(c_middle, n_middle[1]);
		parse_str_list(c_ending, n_ending[1]);
	}
	string gen_name(rand_gen_t &rgen) {
		bool at_end(0), vc((rgen.rand() % 100) < 90); // 0 = vowel, 1 = consonant
		string name;

		for (unsigned i = 0; !at_end; ++i, vc ^= 1) {
			at_end = (i >= 5 || name.size() >= 8 || (i > 1 && (rgen.rand() % 100) < int((vc ? 10 : 5) + 16*i)));
			vector<string> const &str((i == 0) ? n_start[vc] : (at_end ? n_ending[vc] : n_middle[vc]));
			assert(!str.empty());
			name += str[rgen.rand() % str.size()];
		}
		assert(!name.empty());
		name[0] += ('A' - 'a'); // make uppercase
		//cout << name << endl; // TESTING
		return name;
	}
	bool valid() const {return !n_start[0].empty();}
};
name_gen_t name_gen_universe, name_gen_city;

void parse_universe_name_str_tables() {
	if (name_gen_universe.valid()) return; // already parsed

	{ // Suggested phoneme distribution from GPT4
		// Vowels
		string const v_com = "a e i o u ae ai ao au ea eo eu ia ie io oa oi ou ua ui ue ";
		string const v_str = "ii uu yu yi aia oio uai eie ";
		string const v_beg = "aii eau ieo oau yae yea yio ";
		string const v_mid = "aia aio aie aea iea ioa iou oai oua oue ";
		string const v_end = "aie aye oie uie uai yai yio yea ";
		// Consonants
		string const c_com = "l n r s t z v b d g h k m p w x ";
		string const c_str = "ph th ch sh zh gh bh dh jh kh lh mh rh wh ";
		string const c_beg = "bl br cl cr dr fl fr gl gr pl pr tr sl sr sk st sp sh sw thr fr ";
		string const c_mid = "bb dd ff gg ll nn mm pp rr ss tt zz xh chl chr phl phr thr shl shw ";
		string const c_end = "ct ck nd ng nk nt rt rk rn rm rp rb rd rg rtst sh ";

		string const v_start  = v_com + v_com + v_str + v_beg + v_beg;
		string const v_middle = v_com + v_com + v_str + v_mid + v_mid;
		string const v_ending = v_com + v_com + v_str + v_end + v_end;
		string const c_start  = c_com + c_com + c_str + c_beg + c_beg;
		string const c_middle = c_com + c_com + c_str + c_beg + c_mid  + c_mid + c_end;
		string const c_ending = c_com + c_com + c_str + c_end + c_end;
		name_gen_universe.init(v_start, v_middle, v_ending, c_start, c_middle, c_ending);
	}
	{ // original strings
		// Vowels
		string const v_com = "a e i o ";
		string const v_str = "u ai io ";
		string const v_beg = "au ea ei eo eu ou ya ye yo ";
		string const v_mid = "ao au ea ee ei eo eu ia ie oa oi oo ou ue ";
		string const v_end = "ay ee ey ia ie oo oy ue ion ";
		// Consonants
		string const c_com = "l n r s t ";
		string const c_str = "b b c c d d f f g g h h j k m m p p v w w x ch ch sp st st th th ";
		string const c_beg = "bl br cl cr dr fl fr gl gr ph pl pr sk sh sh tr tr wh q str thr ";
		string const c_mid = "bb dd ff gg pp rr rr tt tt rc ";
		string const c_end = "ck ck ct gh ld ld ll ln ln lm lp lt mp nc nd nd ng nk rk rs rs rt rt ss ss gth nch ";

		string const v_start  = v_com + v_com + v_com + v_com + v_com + v_str + v_beg;
		string const v_middle = v_com + v_com + v_com + v_com + v_com + v_str + v_mid;
		string const v_ending = v_com + v_com + v_com + v_com + v_com + v_str + v_end;
		string const c_start  = c_com + c_com + c_com + c_com + c_str + c_beg;
		string const c_middle = c_com + c_com + c_com + c_com + c_str + c_beg + c_mid + c_end;
		string const c_ending = c_com + c_com + c_com + c_com + c_str + c_end;
		name_gen_city.init(v_start, v_middle, v_ending, c_start, c_middle, c_ending);
	}
}

string gen_random_name(rand_gen_t &rgen, unsigned min_len, bool for_universe) {
	parse_universe_name_str_tables();
	string name;

	for (unsigned n = 0; n < 100; ++n) {
		name = (for_universe ? name_gen_universe : name_gen_city).gen_name(rgen);
		if (name.size() >= min_len) break;
	}
	return name;
}

extern rand_gen_t global_rand_gen;

void named_obj::gen_name(s_object const &sobj) {
	name = gen_random_name(global_rand_gen, 0, 1); // min_len=0, for_universe=1
	lookup_given_name(sobj); // already named, overwrite the old value (but need to preserve random number generator state)
	//cout << name << "  ";
}




