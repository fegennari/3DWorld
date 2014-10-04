// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 9/3/03

#ifndef _SPILLOVER_H_
#define _SPILLOVER_H_

#include "3DWorld.h" // need iterator #defs


class spillover {

public:
	spillover() : cur_seen_ix(1), cur_connected(1) {}
	void clear() {data.clear(); cur_seen_ix = cur_connected = 1;}
	void init(unsigned max_index);
	void insert(unsigned index1, unsigned index2);
	void remove(unsigned index1, unsigned index2);
	void remove_all_i(unsigned index1);
	void remove_connected(unsigned index1);
	bool member(unsigned index1, unsigned index2) const;
	bool member_deep(unsigned index1, unsigned index2);
	bool member_recur(unsigned index1, unsigned index2, bool use_cache=0, vector<unsigned char> *used=nullptr);
	bool member2way(unsigned index1, unsigned index2);
	void get_fanout(unsigned index1, vector<unsigned> &fanout, vector<unsigned char> *used);
	void get_connected_components(unsigned index1, vector<unsigned> &cc, vector<unsigned char> *used=nullptr);

private:
	struct graph_node : public set<unsigned> {
		unsigned seen, connected, unconnected;
		graph_node() : seen(0), connected(0), unconnected(0) {}
	};
	vector<graph_node> data;
	vector<unsigned> fanout;
	unsigned cur_seen_ix, cur_connected;
};


#endif

