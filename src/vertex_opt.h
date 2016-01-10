// 3D World - GPU Indexed Vertex Optimization Algorithms Header
// by Frank Gennari
// 6/5/13

#ifndef _VERT_OPT_H_
#define _VERT_OPT_H_

#include "3DWorld.h"


class vert_optimizer {

	vector<unsigned> &indices;
	unsigned num_verts, npts_per_prim;

	struct vbuf_entry_t {
		unsigned ix, pos;
		vbuf_entry_t() : ix((unsigned)-1), pos(0) {}
	};

	template<unsigned N> struct vert_block_t {
		unsigned v[N];

		unsigned min_ix() const {return min(min(v[0], v[1]), ((N == 3) ? v[2] : min(v[2], v[3])));}
		bool operator<(vert_block_t<N> const &b) const {return (min_ix() < b.min_ix());}
		
		static void sort_by_min_ix(vector<unsigned> &ixs) {
			vert_block_t<N> *start((vert_block_t<N> *)(&ixs.front()));
			vert_block_t<N> *end  ((vert_block_t<N> *)(&ixs.front() + ixs.size()));
			std::sort(start, end);
		}
	};

	float calc_acmr() const;

public:
	vert_optimizer(vector<unsigned> &indices_, unsigned num_verts_, unsigned npts_per_prim_) :
	  indices(indices_), num_verts(num_verts_), npts_per_prim(npts_per_prim_) {}
	void run(bool full_opt, bool verbose);
};

#endif // _VERT_OPT_H_
