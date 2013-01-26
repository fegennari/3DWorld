// 3D World - Universe Stellar Object Modification Management
// by Frank Gennari
// 11/18/06
#include "universe.h"
#include "asteroid.h" // required for vector<usateroid_field>::_Tidy()?
#include <iostream>
#include <fstream>

using namespace std;


string const property_tag("BEGIN_PROPERTY"), end_tag("END");
string const def_modmap_fn("universe.modmap");


modmap modmaps[N_UMODS];


bool import_default_modmap() {

	if (!import_modmap(def_modmap_fn)) {
		cerr << "Could not import modmap file '" << def_modmap_fn << "'." << endl;
		return 0;
	}
	return 1;
}


bool import_modmap(string const &filename) {

	unsigned num;
	string str;
	ifstream in(filename.c_str());

	if (!in.good()) {
		cerr << "Failed to open modmap file '" << filename << "' for reading" << endl;
		return 0;
	}
	if (!(in >> num) || num != N_UMODS) {
		cerr << "Failed to read modmap file '" << filename << "' header" << endl;
		return 0;
	}
	for (unsigned i = 0; i < N_UMODS; ++i) {
		if (!in.good() || !(in >> str >> num) || str != property_tag) {
			cerr << "Failed to read modmap file '" << filename << "' header for modmap " << i << endl;
			return 0;
		}
		modmap_val_t val;
		s_object sobj;
		modmaps[i].clear();
		
		for (unsigned j = 0; j < num; ++j) {
			if (!sobj.read(in) || !(in >> val)) {
				cerr << "Failed to read modmap file '" << filename << "' entry " << modmaps[i].size() << " for modmap " << i << endl;
				return 0;
			}
			modmaps[i][sobj] = val;
		}
	}
	return (in.good() && (in >> str) && str == end_tag);
}


bool export_modmap(string const &filename) {

	ofstream out(filename.c_str());
	if (!out.good() || !(out << N_UMODS << endl)) return 0;

	for (unsigned i = 0; i < N_UMODS; ++i) {
		if (!out.good() || !(out << property_tag << " " << modmaps[i].size() << endl)) return 0;
		
		for (modmap::const_iterator it = modmaps[i].begin(); it != modmaps[i].end(); ++it) {
			if (!it->first.write(out) || !(out << " " << it->second << endl)) return 0;
		}
	}
	return (out.good() && (out << end_tag << endl));
}


bool s_object::is_destroyed() const {

	return (modmaps[MOD_DESTROYED].find(*this) != modmaps[MOD_DESTROYED].end());
}


void s_object::register_destroyed_sobj() const {

	if (type != UTYPE_NONE) modmaps[MOD_DESTROYED][get_shifted_sobj(*this)] = "1";
}


unsigned s_object::get_owner() const {

	modmap::const_iterator it(modmaps[MOD_OWNER].find(*this));
	if (it == modmaps[MOD_OWNER].end() || it->second.empty()) return NO_OWNER;
	return unsigned(it->second[0] - '0');
}


void s_object::set_owner(unsigned owner) const {

	if (owner == NO_OWNER) {
		modmaps[MOD_OWNER].erase(get_shifted_sobj(*this)); // should be OK even if doesn't exist (but should exist)
		return;
	}
	string &str(modmaps[MOD_OWNER][get_shifted_sobj(*this)]);
	str.clear(); // remove previous owner (if any)
	str.push_back(char(owner) + '0'); // add new owner
}


bool named_obj::rename(s_object const &sobj, string const &name_) {

	name = name_;
	modmaps[MOD_NAME][sobj] = name;
	return 1;
}


bool named_obj::lookup_given_name(s_object const &sobj) {

	modmap::const_iterator it(modmaps[MOD_NAME].find(sobj));
	if (it == modmaps[MOD_NAME].end()) return 0;
	name = it->second;
	return 1;
}





