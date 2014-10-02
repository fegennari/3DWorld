// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 9/3/03

#include "spillover.h"

using std::cout;
using std::endl;


void spillover::init(unsigned max_index) {

	clear();
	data.resize(max_index);
	seen.resize(max_index, 0);
}

void spillover::insert(unsigned index1, unsigned index2) { // insert index2 into index1 (source, dest)
	
	assert(index1 < data.size() && index2 < data.size());
	assert(index1 != index2);
	data[index1].insert(index2);
}

void spillover::remove(unsigned index1, unsigned index2) { // remove index2 from index1
	
	assert(index1 < data.size() && index2 < data.size());
	assert(index1 != index2);
	data[index1].erase(index2);
}

void spillover::remove_all_i(unsigned index1) { // remove outgoing edges

	assert(index1 < data.size());
	data[index1].clear();
}

void spillover::remove_connected(unsigned index1) { // remove incoming edges

	assert(index1 < data.size());
	set<unsigned> const sdata(data[index1]); // have to copy the set so that we can modify the original

	for (set<unsigned>::const_iterator j = sdata.begin(); j != sdata.end(); ++j) {
		if (member(*j, index1)) remove(index1, *j);
	}
}

bool spillover::member(unsigned index1, unsigned index2) const { // index2 is a member of index1

	assert(index1 < data.size() && index2 < data.size());
	assert(index1 != index2);
	return (data[index1].find(index2) != data[index1].end());
}

bool spillover::member_recur(unsigned index1, unsigned index2, unsigned level) const { // index2 is a member of index1

	assert(seen.size() == data.size());
	assert(index1 < data.size() && index2 < data.size());
	assert(index1 != index2);
	if (level == 0) {++cur_seen_ix;} // invalidate seen values
	seen[index1] = cur_seen_ix;
	if (data[index1].find(index2) != data[index1].end()) return 1;

	for (set<unsigned>::const_iterator i = data[index1].begin(); i != data[index1].end(); ++i) {
		if (seen[*i] == cur_seen_ix)           continue; // already seen
		if (member_recur(*i, index2, level+1)) return 1;
	}
	return 0;
}

bool spillover::member2way(unsigned index1, unsigned index2) const {
	
	return (member_recur(index1, index2) && member_recur(index2, index1));
}

void spillover::get_fanout(unsigned index1, set<unsigned> &fanout) const { // recursive

	assert(index1 < data.size());
	fanout.insert(index1);
	
	for (set<unsigned>::const_iterator i = data[index1].begin(); i != data[index1].end(); ++i) {
		if (fanout.find(*i) == fanout.end()) get_fanout(*i, fanout);
	}
}

void spillover::get_connected_components(unsigned index1, vector<unsigned> &cc) const {

	cc.resize(0);
	assert(index1 < data.size());
	if (data[index1].empty()) return;
	set<unsigned> fanout;
	get_fanout(index1, fanout);

	for (set<unsigned>::const_iterator i = fanout.begin(); i != fanout.end(); ++i) {
		if (*i != index1 && member_recur(*i, index1)) cc.push_back(*i);
	}
}




