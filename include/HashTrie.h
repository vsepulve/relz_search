#ifndef _HASH_TRIE_v3_H
#define _HASH_TRIE_v3_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <iostream>
#include <sstream>
#include <fstream>

#include <algorithm>
#include <vector>
#include <unordered_map>
#include <map>
#include <set>
#include <memory>

#include <sparsehash/dense_hash_map>
#include <sparsehash/sparse_hash_map>

#include "NanoTimer.h"
#include "KarpRabin.h"
#include "KarpRabinFactorsSuffixes.h"
#include "FactorsIteratorCompacted.h"

#define NOT_FOUND 0xffffffff

// Maximum length for hashes to traverse global table
// it must be in the form 00..01..11 (2^M-1), here we are using 8 bits
#define MAX_LEN_HASH 0xff

using namespace std;

using google::dense_hash_map;
using google::sparse_hash_map;

class HashTrieNode{

public: 

	// Uncompressed node's data
	unsigned int len;
	unsigned int min;
	unsigned int hash;
	char first;
	
	unsigned int max;
	
	// Structure for chidls, sorted by first char
	vector<HashTrieNode> childs_vector;
	
	HashTrieNode();
	~HashTrieNode();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *arr_y, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len);
	
	void buildRev(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *arr_x, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len);
	
	unsigned int totalChilds(unsigned int &max_len, unsigned int &max_childs, unsigned int &max_height, unsigned int height);
	
/*	void compactData(unsigned int &next_pos, int_vector<> &positions_childs, int_vector<> &n_childs, int_vector<> &len_childs, int_vector<> &min_childs, int_vector<> &hash_childs, int_vector<> &first_childs);*/
	void compactData(unsigned int &next_pos, int_vector<> &positions_childs, int_vector<> &n_childs, int_vector<> &len_childs, int_vector<> &min_childs, int_vector<> &hash_childs, int_vector<> &first_childs, int_vector<> &arr_max_childs, int_vector<> &len_path, unsigned int ini_path, int_vector<> &len_hash);
	
};

class HashTrie{

private: 
	KarpRabin *karp_rabin;
	int_vector<> *arr_factors;
	unsigned int len_text;
	
	CompactedText *compacted_text;
	
	// map<unsigned int, unsigned int> global_hash;
	// unordered_map<unsigned int, unsigned int> global_hash;
	sparse_hash_map<unsigned int, unsigned int> global_hash;
	// dense_hash_map<unsigned int, unsigned int> global_hash;
	
	// Para esta busqueda necesito max_child explicito tambien
	int_vector<> arr_max_childs;
	// Con esto se podra calcular len_childs en recorrido, reemplazandolo
	int_vector<> len_path;
	// Arreglo con el largo de los hash de prefijos
	// No estoy seguro si es posible evitar este valor explicitamente (dado que no estoy navegando por el arbol)
	int_vector<> len_hash;
	
	// Select structures for S and B to process factors data (pos, len)
	bits_s_type::select_1_type *select1_s;
	bits_b_type::select_1_type *select1_b;
	bits_b_type::select_0_type *select0_b;
	
	// Permutation to turn positional factors to its id
	int_vector<> *pi_inv;
	
	// Datos serializados de los nodos
	int_vector<> positions_childs;
	int_vector<> n_childs;
	int_vector<> len_childs;
	int_vector<> min_childs;
	int_vector<> hash_childs;
	int_vector<> first_childs;
	
	void compactData(HashTrieNode &root_node);
	
	// Clone of HashTrieNode::findChild, but using compacted arrays and node_pos as the position for the current node
	// Returns NOT_FOUND in case of not found
	unsigned int findChildInternal(unsigned int node_pos, char c);

	// Clone of HashTrieNode::getRange, but using compacted arrays and node_pos as the position for the current node
	pair<unsigned int, unsigned int> getRangeInternal(unsigned int node_pos, vector<unsigned long long> &kr_pat_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, unsigned int cur_max, const string &pattern);
	pair<unsigned int, unsigned int> getRangeInternalNoHash(unsigned int node_pos, vector<unsigned long long> &kr_pat_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, unsigned int cur_max, const string &pattern);
	pair<unsigned int, unsigned int> getRangeTableInternal(vector<unsigned long long> &kr_pat_vector, unsigned int pos, KarpRabin *karp_rabin, const string &pattern);
	
	// Clone of HashTrieNode::getRange, but using compacted arrays and node_pos as the position for the current node
	// Version for Reverse factors
	pair<unsigned int, unsigned int> getRangeRevInternal(unsigned int node_pos, vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, unsigned int cur_max, const string &pattern_rev);
	pair<unsigned int, unsigned int> getRangeRevInternalNoHash(unsigned int node_pos, vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, unsigned int cur_max, const string &pattern_rev);
	pair<unsigned int, unsigned int> getRangeTableRevInternal(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, KarpRabin *karp_rabin, const string &pattern_rev);
	
	// Clone of HashTrieNode::print, but using compacted arrays and node_pos as the position for the current node
	void printInternal(unsigned int node_pos, unsigned int level);
	
	void prepareHashMap(unsigned int node_pos, unsigned int path_len, unordered_map<unsigned int, pair<unsigned int, unsigned int>> &marked_hash);
	void prepareHashMapRev(unsigned int node_pos, unsigned int path_len, unordered_map<unsigned int, pair<unsigned int, unsigned int>> &marked_hash);
	
public: 
	
	HashTrie();
	HashTrie(const char *full_text, unsigned int _len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_factors, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *kr_factors, bool reverse = false);
	virtual ~HashTrie();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_factors, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *kr_factors, bool reverse = false);
	
	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, const string &pattern, bool use_hash = true);
	
	pair<unsigned int, unsigned int> getRangeTable(vector<unsigned long long> &kr_pat_vector, unsigned int pos, const string &pattern);
	
	pair<unsigned int, unsigned int> getRangeRev(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, const string &pattern_rev, bool use_hash = true);
	
	pair<unsigned int, unsigned int> getRangeTableRev(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, const string &pattern_rev);
	
	void print();
	
	void printSize();
	
	void save(const string &file);
	
	void load(unsigned int _len_text, KarpRabin *_karp_rabin, CompactedText *_compacted_text, bits_s_type::select_1_type *_select1_s, bits_b_type::select_1_type *_select1_b, bits_b_type::select_0_type *_select0_b, int_vector<> *_pi_inv, int_vector<> *_arr_factors, const string &file);
	
	unsigned int totalChilds(){
		return positions_childs.size();
	}
	
	unsigned int getSizeBytes();
	
	unsigned long long hash_nano;
	
	static unsigned int codeChar(char c);
	
	static constexpr char decodeChar[] = "ACGT";
	
	unsigned int n_collisions;
	unsigned long long nano_search_global;
	unsigned long long nano_search_collisions;
	
};





#endif //_HASH_TRIE_H





