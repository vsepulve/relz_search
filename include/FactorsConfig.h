#ifndef _FACTORS_CONFIG_H
#define _FACTORS_CONFIG_H

#include <set>
#include <map>
#include <list>
#include <vector>
#include <unordered_map>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rmq_support.hpp>
#include <sdsl/inv_perm_support.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bp_support_g.hpp>
#include <sdsl/bp_support_gg.hpp>

using namespace sdsl;
using namespace std;

// Type for bitarray S of the relz index, of size (reference_text + Z)
//typedef rrr_vector<127> bits_s_type;
typedef bit_vector bits_s_type;

// Type for bitarray B of the relz index, of size (collection_text)
//typedef rrr_vector<127> bits_b_type;
typedef sd_vector<> bits_b_type;

// Type for the fm-index
// Default Values: wt_huff<>, 32 (Sample density for suffix array), 64 (Sample density for inverse suffix array)
//typedef csa_wt<wt_huff<>, 32, 64, sa_order_sa_sampling<>, isa_sampling<>> fm_index_type;
typedef csa_wt<wt_huff<>, 16, 32, sa_order_sa_sampling<>, isa_sampling<>> fm_index_type;

// Type for the rmq structure
// typedef rmq_succinct_sct<false, bp_support_sada<256,32,rank_support_v5<> > > rmq_type;
// typedef rmq_succinct_sct<false, bp_support_sada<128,32,rank_support_v5<> > > rmq_type;
typedef rmq_succinct_sct<false, bp_support_gg<> > rmq_type;

// Type for the Wavelet Tree structure used to combine ranges
//typedef wt_int<rrr_vector<127>> wt_type;
typedef wt_int<rrr_vector<63>> wt_type;



#endif //_FACTORS_CONFIG_H





