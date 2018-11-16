#ifndef _RELZ_INDEX_HASH_H
#define _RELZ_INDEX_HASH_H

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
#include "HashTrie.h"
#include "KarpRabinFactorsSuffixes.h"

using namespace sdsl;
using namespace std;

class RelzIndexHash {

private: 
	
	unsigned int len_text;
	const char *ref_text;
	unsigned int len_ref;
	unsigned int n_factors;
	
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
	static const bool precompute_rmq = true;
	vector<unsigned int> arr_tu;
	vector<unsigned int> arr_pu;
	vector<unsigned int> arr_lu;
	
	// Solo para pruebas
//	bit_vector arr_s;
//	bit_vector arr_b;
//	int_vector<> ez;
//	int_vector<> values_wt;
	// Simplifacion del calculo de posicion de inicio de factor en coleccion 
	vector<unsigned int> factors_start;
	
	// Cache de iteradores
	unordered_map<unsigned int, FactorsIterator> mapa_iterators;
	unordered_map<unsigned int, FactorsIteratorReverse> mapa_iterators_rev;
	
	// Falta la estructura para agregar efectivamente los resultados, quiza un vector de posiciones
	void recursive_rmq(unsigned int ini, unsigned int fin, unsigned int min_pos, unsigned int occ_ref, vector<unsigned int> &results);
	
	char getChar(unsigned int factor, unsigned int pos);
	
	char getCharRev(unsigned int factor, unsigned int pos);
	
public: 
	RelzIndexHash();
	RelzIndexHash(KarpRabin *_karp_rabin);
	RelzIndexHash(vector<pair<unsigned int, unsigned int> > &factors, char *full_text, unsigned int _len_text, const char *_ref_text, unsigned int _len_ref, KarpRabin *_karp_rabin, const char *index_base_file);
	~RelzIndexHash();
	
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
	HashTrie tree_y;
	HashTrieRev tree_x;
	KarpRabinFactorsSuffixes *kr_factors;
	
	
	void printSize();
	
	
	
	
};







#endif //_RELZ_INDEX_HASH_H