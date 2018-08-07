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
#include <memory>

#include "KarpRabin.h"

using namespace std;

//class HashTrieNode;

class HashTrieNode{

public: 

	// Datos descomprimidos del nodo
	unsigned long long hash;
	unsigned int len;
	unsigned int min;
	unsigned int max;
	unsigned int factor;
	
	// Estructura para los hijos
//	map<unsigned long long, HashTrieNode> childs;
	unordered_map<unsigned int, shared_ptr<HashTrieNode>> childs;
	
	unsigned int getMaxComp(const char *str_1, unsigned int len_1, const char *str_2, unsigned int len_2);
	
	HashTrieNode();
	~HashTrieNode();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, bool _factor_only, KarpRabin *_karp_rabin, unsigned int min, unsigned int max, unsigned int processed_len);
	
	void print(unsigned int level);
	
};

class HashTrie{

private: 
	KarpRabin *karp_rabin;
	HashTrieNode root;
	bool factor_only;
	
public: 
	
	HashTrie();
	HashTrie(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, bool _factor_only, KarpRabin *_karp_rabin);
	virtual ~HashTrie();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, bool _factor_only, KarpRabin *_karp_rabin);
	
	void print();
	
};







#endif //_HASH_TRIE_H





