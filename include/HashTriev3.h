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

#include "NanoTimer.h"
#include "KarpRabin.h"
#include "KarpRabinFactorsSuffixes.h"
#include "KarpRabinFactorsSuffixesv2.h"

#define NOT_FOUND 0xffffffff

using namespace std;

class HashTriev3Node{

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
	
	unsigned int totalChilds(unsigned int &max_len, unsigned int &max_childs, unsigned int &max_height, unsigned int height);
	
	void compactData(unsigned int &next_pos, int_vector<> &positions_childs, int_vector<> &n_childs, int_vector<> &len_childs, int_vector<> &min_childs, int_vector<> &hash_childs, vector<char> &first_childs);
	
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
	
	void compactData(HashTriev3Node &root_node);
	
	// Clone of HashTriev3Node::findChild, but using compacted arrays and node_pos as the position for the current node
	// Returns NOT_FOUND in case of not found
	unsigned int findChildInternal(unsigned int node_pos, char c);

	// Clone of HashTriev3Node::getRange, but using compacted arrays and node_pos as the position for the current node
	pair<unsigned int, unsigned int> getRangeInternal(unsigned int node_pos, vector<unsigned long long> &kr_pat_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, int_vector<> *arr_y, unsigned int cur_max, const string &pattern);
	
	// Clone of HashTriev3Node::print, but using compacted arrays and node_pos as the position for the current node
	void printInternal(unsigned int node_pos, unsigned int level);
	
public: 
	
	HashTriev3();
	HashTriev3(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	virtual ~HashTriev3();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	
	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, const string &pattern);
	
	void print();
	
	void printSize();
	
	void save(const string &file);
	
	void load(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors, int_vector<> *arr_y, const string &file);
	
	unsigned int totalChilds(){
		return positions_childs.size();
	}
	
	unsigned long long hash_nano;
	
};

class HashTriev3RevNode{
	
	// Returns NOT_FOUND in case of not found
	unsigned int findChild(char c);

public: 

	// Uncompressed node's data
	unsigned int len;
	unsigned int min;
	unsigned int hash;
	char first;
	
	// Structure for chidls, sorted by first char
	vector<HashTriev3RevNode> childs_vector;
	
	HashTriev3RevNode();
	~HashTriev3RevNode();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *arr_x, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len);
	
	void print(unsigned int level);
	
	unsigned int totalChilds(unsigned int &max_len, unsigned int &max_childs, unsigned int &max_height, unsigned int height);
	
	void compactData(unsigned int &next_pos, int_vector<> &positions_childs, int_vector<> &n_childs, int_vector<> &len_childs, int_vector<> &min_childs, int_vector<> &hash_childs, vector<char> &first_childs);
	
//	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, int_vector<> *arr_x, const string &pattern_rev);
	
//	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, int_vector<> *arr_x, unsigned int cur_max, const string &pattern_rev);
//	
//	void save(fstream &writer);
//	
//	void load(fstream &reader);
	
};

class HashTriev3Rev{

private: 
	KarpRabin *karp_rabin;
	KarpRabinFactorsSuffixes *kr_factors;
	int_vector<> *arr_x;
	
	// Datos serializados de los nodos
	int_vector<> positions_childs;
	int_vector<> n_childs;
	int_vector<> len_childs;
	int_vector<> min_childs;
	int_vector<> hash_childs;
	vector<char> first_childs;
	
	void compactData(HashTriev3RevNode &root_node);
	
	// Clone of HashTriev3Node::findChild, but using compacted arrays and node_pos as the position for the current node
	// Returns NOT_FOUND in case of not found
	unsigned int findChildInternal(unsigned int node_pos, char c);

	// Clone of HashTriev3Node::getRange, but using compacted arrays and node_pos as the position for the current node
	pair<unsigned int, unsigned int> getRangeInternal(unsigned int node_pos, vector<unsigned long long> &kr_pat_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, int_vector<> *arr_y, unsigned int cur_max, const string &pattern);
	
	// Clone of HashTriev3Node::print, but using compacted arrays and node_pos as the position for the current node
	void printInternal(unsigned int node_pos, unsigned int level);
	
public: 
//	HashTriev3RevNode root;
	
	HashTriev3Rev();
	HashTriev3Rev(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	virtual ~HashTriev3Rev();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	
	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, const string &pattern_rev);
	
	void print();
	
	void printSize();
	
	void save(const string &file);
	
	void load(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors, int_vector<> *_arr_x, const string &file);
	
	unsigned int totalChilds(){
		return positions_childs.size();
	}
	
};






#endif //_HASH_TRIE_H





