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

using namespace std;

class HashTrieCharNode{

public: 

	// Uncompressed node's data
	unsigned int len;
	unsigned int min;
	unsigned int max;
	unsigned int hash;
	// Position of the factor in the collection
	unsigned int min_factor_pos;
	
	// Structure for chidls, indexed by first char
	map<char, shared_ptr<HashTrieCharNode>> childs;
	
	HashTrieCharNode();
	~HashTrieCharNode();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len);
	
	void print(unsigned int level);
	
	unsigned int totalChilds(unsigned int &max_len, unsigned int &max_childs, unsigned int &max_height, unsigned int height);
	
	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, const string &pattern, unsigned long long *hash_nano);
	
	void save(fstream &writer);
	
	void load(fstream &reader, unsigned int processed, vector<unsigned int> &factors_start, const char *full_text);
	
};

class HashTrieChar{

private: 
	KarpRabin *karp_rabin;
	KarpRabinFactorsSuffixes *kr_factors;
	
public: 
	HashTrieCharNode root;
	
	HashTrieChar();
	HashTrieChar(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	virtual ~HashTrieChar();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	
	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, const string &pattern);
	
	void print();
	
	void printSize();
	
	void save(const string &file);
	
	void load(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors, const string &file, vector<unsigned int> &factors_start, const char *full_text);
	
	unsigned long long hash_nano;
	
};

class HashTrieCharRevNode{

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
	string text;
	
	// Structure for chidls, indexed by first char
	unordered_map<char, shared_ptr<HashTrieCharRevNode>> childs;
	
	HashTrieCharRevNode();
	~HashTrieCharRevNode();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_x, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len);
	
	void print(unsigned int level);
	
	unsigned int totalChilds(unsigned int &max_len, unsigned int &max_childs, unsigned int &max_height, unsigned int height);
	
	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, const string &pattern_rev);
	
	void save(fstream &writer);
	
	void load(fstream &reader);
	
};

class HashTrieCharRev{

private: 
	KarpRabin *karp_rabin;
	KarpRabinFactorsSuffixes *kr_factors;
	
public: 
	HashTrieCharRevNode root;
	
	HashTrieCharRev();
	HashTrieCharRev(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	virtual ~HashTrieCharRev();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	
	pair<unsigned int, unsigned int> getRange(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, const string &pattern_rev);
	
	void print();
	
	void printSize();
	
	void save(const string &file);
	
	void load(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors, const string &file);
	
};






#endif //_HASH_TRIE_H





