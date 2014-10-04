// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 9/3/03

#include "spillover.h"

using std::cout;
using std::endl;


void spillover::init(unsigned max_index) {
	clear();
	data.resize(max_index);
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

	for (auto j = sdata.begin(); j != sdata.end(); ++j) {
		if (member(*j, index1)) {remove(index1, *j);}
	}
}

bool spillover::member(unsigned index1, unsigned index2) const { // index2 is a member of index1
	assert(index1 < data.size() && index2 < data.size());
	assert(index1 != index2);
	return (data[index1].find(index2) != data[index1].end());
}

bool spillover::member_deep(unsigned index1, unsigned index2) {
	assert(index2 < data.size());
	++cur_seen_ix; // invalidate seen values
	return member_recur(index1, index2);
}

bool spillover::member_recur(unsigned index1, unsigned index2, bool use_cache, vector<unsigned char> *used) { // index2 is a member of index1

	assert(index1 != index2);
	data[index1].seen = cur_seen_ix;
	if (used != nullptr && (*used)[index1]) return 0; // already used
	if (use_cache && data[index1].unconnected == cur_connected) return 0;
	if (use_cache && data[index1].connected   == cur_connected) return 1;

	for (auto i = data[index1].begin(); i != data[index1].end(); ++i) {
		assert(*i < data.size());
		if (*i == index2) return 1; // found it
		if (data[*i].seen == cur_seen_ix) continue; // already seen
		
		if (member_recur(*i, index2, use_cache, used)) {
			if (use_cache) {data[*i].connected = cur_connected;}
			return 1;
		}
	}
	return 0;
}

bool spillover::member2way(unsigned index1, unsigned index2) {
	return (member_deep(index1, index2) && member_deep(index2, index1));
}

void spillover::get_fanout(unsigned index1, vector<unsigned> &fanout, vector<unsigned char> *used) {
	
	for (auto i = data[index1].begin(); i != data[index1].end(); ++i) {
		assert(*i < data.size());
		if (used != nullptr && (*used)[*i]) continue; // already used
		if (data[*i].seen == cur_seen_ix)   continue; // already seen
		data[*i].seen = cur_seen_ix;
		fanout.push_back(*i);
		get_fanout(*i, fanout, used);
	}
}

void spillover::get_connected_components(unsigned index1, vector<unsigned> &cc, vector<unsigned char> *used) {

	cc.resize(0);
	assert(index1 < data.size());
	if (data[index1].empty()) return;
	++cur_seen_ix;   // invalidate seen values
	++cur_connected; // invalidate connected/unconnected
	data[index1].seen = cur_seen_ix;
	fanout.resize(0);
	get_fanout(index1, fanout, used); // excludes index1

	for (auto i = fanout.begin(); i != fanout.end(); ++i) {
		assert(*i != index1);
		++cur_seen_ix; // invalidate seen values
		bool const ret(member_recur(*i, index1, 1, used));
		(ret ? data[*i].connected : data[*i].unconnected) = cur_connected;
		if (ret) {cc.push_back(*i);}
	}
}




