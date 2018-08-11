#ifndef _HASH_TRIE_H
#define _HASH_TRIE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <iostream>
#include <sstream>
#include <fstream>

#include <vector>
#include <unordered_map>
#include <map>
#include <set>
#include <memory>

#include "KarpRabin.h"
#include "KarpRabinFactorsSuffixes.h"

using namespace std;

//class HashTrieNode;

class HashTrieNode{

public: 

	// Datos descomprimidos del nodo
//	unsigned long long hash;
	unsigned int len;
	unsigned int min;
	unsigned int max;
//	unsigned int factor;
	
	// Estructura para los hijos
//	map<unsigned long long, HashTrieNode> childs;
	unordered_map<unsigned int, shared_ptr<HashTrieNode>> childs;
	vector<unsigned int> childs_lenghts;
	
	unsigned int getMaxComp(const char *str_1, unsigned int len_1, const char *str_2, unsigned int len_2);
	
	HashTrieNode();
	~HashTrieNode();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, bool _factor_only, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len);
	
	void print(unsigned int level);
	
	pair<unsigned int, unsigned int> getRange(const char *pattern, unsigned int pat_len, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors);
};

class HashTrie{

private: 
	KarpRabin *karp_rabin;
	KarpRabinFactorsSuffixes *kr_factors;
	HashTrieNode root;
	bool factor_only;
	
public: 
	
	HashTrie();
	HashTrie(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, bool _factor_only, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	virtual ~HashTrie();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, bool _factor_only, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	
	// Notar que este metodo despues debe usar la estructura de hash de prefijos del patron completo para acelerar sus hash
	pair<unsigned int, unsigned int> getRange(const string &pattern);
	
	void print();
	
};







#endif //_HASH_TRIE_H





