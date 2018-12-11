#ifndef _HASH_TRIE_v2_H
#define _HASH_TRIE_v2_H

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

#include "NanoTimer.h"
#include "KarpRabin.h"
#include "KarpRabinFactorsSuffixes.h"
#include "KarpRabinFactorsSuffixesv2.h"

#define NOT_FOUND 0xffffffff

using namespace std;

class HashTriev3Node{
	
	// Returns NOT_FOUND in case of not found
	unsigned int findChild(char c);

public: 

	// Uncompressed node's data
	unsigned int len;
	unsigned int min;
	unsigned int hash;
	char first;
	
	// Structure for chidls, sorted by first char
	vector<HashTriev3Node> childs_vector;
	
	HashTriev3Node();
	~HashTriev3Node();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *arr_y, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len);
	
	void print(unsigned int level);
	
	unsigned int totalChilds(unsigned int &max_len, unsigned int &max_childs, unsigned int &max_height, unsigned int height);
	
	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, int_vector<> *arr_y, unsigned int cur_max, const string &pattern, unsigned long long *hash_nano);
	
	void save(fstream &writer);
	
	void load(fstream &reader);
	
};

class HashTriev3{

private: 
	KarpRabin *karp_rabin;
	KarpRabinFactorsSuffixes *kr_factors;
	int_vector<> *arr_y;
	
	// Datos serializados de los nodos
	int_vector<> positions_childs;
	int_vector<> n_childs;
	int_vector<> len_childs;
	int_vector<> min_childs;
	int_vector<> hash_childs;
	vector<char> first_childs;
	
public: 
	HashTriev3Node root;
	
	HashTriev3();
	HashTriev3(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	virtual ~HashTriev3();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	
	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, const string &pattern);
	
	void print();
	
	void printSize();
	
	void save(const string &file);
	
	void load(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors, int_vector<> *arr_y, const string &file);
	
	unsigned long long hash_nano;
	
};

class HashTriev3RevNode{
	
	// Returns NOT_FOUND in case of not found
	unsigned int findChild(char c);

public: 

	// Uncompressed node's data
	unsigned int len;
	unsigned int min;
//	unsigned int max;
	unsigned int hash;
	// Position of the factor in the collection
//	unsigned int min_factor_pos;
	char first;
	
	// Structure for chidls, indexed by first char
//	unordered_map<char, shared_ptr<HashTriev3RevNode>> childs;

	// Structure for chidls, sorted by first char
	vector<HashTriev3RevNode> childs_vector;
	
	HashTriev3RevNode();
	~HashTriev3RevNode();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *arr_x, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len);
	
	void print(unsigned int level);
	
	unsigned int totalChilds(unsigned int &max_len, unsigned int &max_childs, unsigned int &max_height, unsigned int height);
	
//	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, int_vector<> *arr_x, const string &pattern_rev);
	
	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, int_vector<> *arr_x, unsigned int cur_max, const string &pattern_rev);
	
	void save(fstream &writer);
	
	void load(fstream &reader);
	
};

class HashTriev3Rev{

private: 
	KarpRabin *karp_rabin;
	KarpRabinFactorsSuffixes *kr_factors;
	int_vector<> *arr_x;
	
public: 
	HashTriev3RevNode root;
	
	HashTriev3Rev();
	HashTriev3Rev(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	virtual ~HashTriev3Rev();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	
	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, const string &pattern_rev);
	
	void print();
	
	void printSize();
	
	void save(const string &file);
	
	void load(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors, int_vector<> *_arr_x, const string &file);
	
};






#endif //_HASH_TRIE_H





