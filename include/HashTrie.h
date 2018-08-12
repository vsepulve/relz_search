#ifndef _HASH_TRIE_H
#define _HASH_TRIE_H

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

#include "KarpRabin.h"
#include "KarpRabinFactorsSuffixes.h"

using namespace std;

class HashTrieNode{

public: 

	// Datos descomprimidos del nodo
	unsigned int len;
	unsigned int min;
	unsigned int max;
	// Posicion del factor de referencia en la coleccion (no en el arreglo lexicografico)
	unsigned int min_factor_pos;
	
	// Estructura para los hijos
	vector<unsigned int> childs_lenghts;
	unordered_map<unsigned long long, shared_ptr<HashTrieNode>> childs;
	
	// Estructra secundadaria para los hijos
	// Arreglo de hijos ordenados por largo del texto para acelerar la busqueda
	vector<pair<unsigned int, vector<shared_ptr<HashTrieNode>>>> childs_pairs;
	
	unsigned int getMaxComp(const char *str_1, unsigned int len_1, const char *str_2, unsigned int len_2);
	
	HashTrieNode();
	~HashTrieNode();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len);
	
	void print(unsigned int level);
	
	pair<unsigned int, unsigned int> getRange(const char *pattern, unsigned int pat_len, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors);
	
	void save(fstream &writer);
	
	void load(fstream &reader);
	
	void prepareChilds();
	
};

class HashTrie{

private: 
	KarpRabin *karp_rabin;
	KarpRabinFactorsSuffixes *kr_factors;
	HashTrieNode root;
	
public: 
	
	HashTrie();
	HashTrie(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	HashTrie(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	virtual ~HashTrie();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_y, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	
	// Notar que este metodo despues debe usar la estructura de hash de prefijos del patron completo para acelerar sus hash
	pair<unsigned int, unsigned int> getRange(const string &pattern);
	
	void print();
	
	void save(const string &file);
	
	void load(const string &file);
	
	void prepareChilds();
	
};

class HashTrieRevNode{

public: 

	// Datos descomprimidos del nodo
	unsigned int len;
	unsigned int min;
	unsigned int max;
	// Posicion del factor de referencia en la coleccion (no en el arreglo lexicografico)
	unsigned int min_factor_pos;
	
	// Estructura para los hijos
	vector<unsigned int> childs_lenghts;
	unordered_map<unsigned long long, shared_ptr<HashTrieRevNode>> childs;
	
	// Estructra secundadaria para los hijos
	// Arreglo de hijos ordenados por largo del texto para acelerar la busqueda
	vector<pair<unsigned int, vector<shared_ptr<HashTrieRevNode>>>> childs_pairs;
	
	unsigned int getMaxComp(const char *str_1, unsigned int len_1, const char *str_2, unsigned int len_2);
	
	HashTrieRevNode();
	~HashTrieRevNode();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_x, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len);
	
	void print(unsigned int level);
	
	pair<unsigned int, unsigned int> getRange(const char *pattern, unsigned int pat_len, unsigned int processed, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors);
	
	void save(fstream &writer);
	
	void load(fstream &reader);
	
	void prepareChilds();
	
};

class HashTrieRev{

private: 
	KarpRabin *karp_rabin;
	KarpRabinFactorsSuffixes *kr_factors;
	HashTrieRevNode root;
	
public: 
	
	HashTrieRev();
	HashTrieRev(KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	HashTrieRev(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	virtual ~HashTrieRev();
	
	void build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, vector<unsigned int> &arr_x, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *_kr_factors);
	
	// Notar que este metodo despues debe usar la estructura de hash de prefijos del patron completo para acelerar sus hash
	pair<unsigned int, unsigned int> getRange(const string &pattern);
	
	void print();
	
	void save(const string &file);
	
	void load(const string &file);
	
	void prepareChilds();
	
};






#endif //_HASH_TRIE_H





