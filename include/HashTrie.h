#ifndef _HASH_TRIE_CHAR_H
#define _HASH_TRIE_CHAR_H

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

class HashTrieNode{
	
	// Returns 0xffffffff in case of not found
	unsigned int findChild(char c);

public: 

	// Uncompressed node's data
	unsigned int len;
	unsigned int min;
	unsigned int max;
	unsigned int hash;
	// Position of the factor in the collection
	unsigned int min_factor_pos;
	char first;
	
	// Structure for chidls, indexed by first char
	map<char, shared_ptr<HashTrieNode>> childs;
	
//	vector<pair<char, HashTrieNode>> childs_vector;
	vector<HashTrieNode> childs_vector;
	
	HashTrieNode();
	~HashTrieNode();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *arr_y, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len);
	
	void print(unsigned int level);
	
	unsigned int totalChilds(unsigned int &max_len, unsigned int &max_childs, unsigned int &max_height, unsigned int height);
	
	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, int_vector<> *arr_y, const string &pattern, unsigned long long *hash_nano);
	
//	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, int_vector<> *arr_y, unsigned int cur_max, const string &pattern, unsigned long long *hash_nano);
	
	void save(fstream &writer);
	
	void load(fstream &reader);
	
};

class HashTrie{

private: 
	KarpRabin *karp_rabin;
	KarpRabinFactorsSuffixes *kr_factors;
	int_vector<> *arr_y;
	
public: 
	HashTrieNode root;
	
	HashTrie();
	HashTrie(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	virtual ~HashTrie();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	
	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, const string &pattern);
	
	void print();
	
	void printSize();
	
	void save(const string &file);
	
	void load(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors, int_vector<> *arr_y, const string &file);
	
	unsigned long long hash_nano;
	
};

class HashTrieRevNode{

public: 

	// Uncompressed node's data
	unsigned int len;
	unsigned int min;
	unsigned int max;
	unsigned int hash;
	// Position of the factor in the collection
	unsigned int min_factor_pos;
	
	// To simplify tge code, I add the text to evaluate prefix hash
	// Obviously this text can be accesed from the reference at querytime (its only 1 factor)
//	string text;
	
	// Structure for chidls, indexed by first char
	unordered_map<char, shared_ptr<HashTrieRevNode>> childs;
	
	HashTrieRevNode();
	~HashTrieRevNode();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *arr_x, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len);
	
	void print(unsigned int level);
	
	unsigned int totalChilds(unsigned int &max_len, unsigned int &max_childs, unsigned int &max_height, unsigned int height);
	
	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, int_vector<> *arr_x, const string &pattern_rev);
	
	void save(fstream &writer);
	
	void load(fstream &reader);
	
};

class HashTrieRev{

private: 
	KarpRabin *karp_rabin;
	KarpRabinFactorsSuffixes *kr_factors;
	int_vector<> *arr_x;
	
public: 
	HashTrieRevNode root;
	
	HashTrieRev();
	HashTrieRev(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	virtual ~HashTrieRev();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	
	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, const string &pattern_rev);
	
	void print();
	
	void printSize();
	
	void save(const string &file);
	
	void load(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors, int_vector<> *_arr_x, const string &file);
	
};






#endif //_HASH_TRIE_H





