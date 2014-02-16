// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 9/3/03

#ifndef _SPILLOVER_H_
#define _SPILLOVER_H_

#include "3DWorld.h" // need iterater #defs


class spillover {

public:
	void clear() {data.resize(0);}
	void init(unsigned max_index);
	void insert(unsigned index1, unsigned index2);
	void remove(unsigned index1, unsigned index2);
	void remove_all_i(unsigned index1);
	void remove_connected(unsigned index1);
	bool member(unsigned index1, unsigned index2) const;
	bool member_recur(unsigned index1, unsigned index2, unsigned level=0) const;
	bool member2way(unsigned index1, unsigned index2) const;
	void get_fanout(unsigned index1, set<unsigned> &fanout) const;
	void get_connected_components(unsigned index1, vector<unsigned> &cc) const;

private:
	vector<set<unsigned> > data;
	mutable set<unsigned> seen;
};


#endif

