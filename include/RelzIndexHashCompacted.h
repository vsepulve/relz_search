#ifndef _RELZ_INDEX_HASH_COMPACTED_H
#define _RELZ_INDEX_HASH_COMPACTED_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <algorithm>
#include <vector>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rmq_support.hpp>
#include <sdsl/inv_perm_support.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "FactorsConfig.h"
#include "FactorsIterator.h"
#include "FactorsFastIteratorComparator.h"
#include "NanoTimer.h"

#include "BitsUtils.h"
#include "KarpRabin.h"
#include "KarpRabinFactorsSuffixes.h"
//#include "KarpRabinFactorsSuffixesv1.h"
#include "KarpRabinFactorsSuffixesv2.h"
#include "HashTrie.h"
#include "HashTriev2.h"

using namespace sdsl;
using namespace std;

class RelzIndexHashCompacted {

private: 
	
	unsigned int len_text;
	unsigned int n_factors;
	unsigned int len_ref;
	char *ref_text;
	
	fm_index_type fm_index;
	
	rmq_type rmq;
	
	bits_s_type bits_s;
	bits_s_type::select_1_type select1_s;
	bits_s_type::select_0_type select0_s;
	
	bits_b_type bits_b;
	bits_b_type::select_1_type select1_b;
	bits_b_type::select_0_type select0_b;
	
	wt_type wt;
	
	int_vector<> pi;
	int_vector<> pi_inv;
	
	int_vector<> arr_x;
	int_vector<> arr_y;
	
	// Prueba de aceleracion de recursive_rmq almacenando los datos de los factores descomprimidos
	static const bool precompute_rmq = false;
	vector<unsigned int> arr_tu;
	vector<unsigned int> arr_pu;
	vector<unsigned int> arr_lu;
	
	// Simplifacion del calculo de posicion de inicio de factor en coleccion 
	vector<unsigned int> factors_start;
	
	void recursive_rmq(unsigned int ini, unsigned int fin, unsigned int min_pos, unsigned int occ_ref, vector<unsigned int> &results);
	
public: 
	RelzIndexHashCompacted();
	RelzIndexHashCompacted(KarpRabin *_karp_rabin);
	RelzIndexHashCompacted(vector<pair<unsigned int, unsigned int> > &factors, char *full_text, unsigned int _len_text, const char *_ref_text, unsigned int _len_ref, KarpRabin *_karp_rabin);
	~RelzIndexHashCompacted();
	
	// void find(const string &pattern, vector<unsigned int> &results);
	void findTimes(const string &pattern, vector<unsigned int> &results);
	
	unsigned long long querytime_p1;
	unsigned long long querytime_p2;
	unsigned long long querytime_p3x;
	unsigned long long querytime_p3y;
	unsigned long long querytime_p4;
	
	unsigned int occs_a;
	unsigned int occs_b;
	unsigned int occs_c;
	
	KarpRabinFactorsSuffixes *getKRFactors(){
		return kr_factors;
	}
	
	KarpRabin *karp_rabin;
	vector<unsigned long long> arr_kr_ref;
	vector<unsigned long long> arr_kr_s;
	KarpRabinFactorsSuffixes *kr_factors;
	
//	HashTrie tree_y;
//	HashTrieRev tree_x;
	
	HashTriev2 tree_y;
	HashTriev2Rev tree_x;
	
	void printSize();
	
	void save(const string &file_base);
	
	void load(const string &file_base, KarpRabin *_karp_rabin);
	
	
	
};







#endif //_RELZ_INDEX_HASH_H
