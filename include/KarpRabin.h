#ifndef _KARP_RABIN_H
#define _KARP_RABIN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <iostream>
#include <sstream>
#include <fstream>

#include <vector>

using namespace std;

class KarpRabin{

protected: 
	unsigned int voc_bits;
	unsigned int kr_mod;
	unsigned int table_size;
	unsigned long long *pow_table;
	
	// Direct version (*x, y times)
	unsigned long long ullpow(unsigned long long x, unsigned int y);
	
	// Base 2 bits version, recursive to avoid overflow
	unsigned long long ullpow2_rec(unsigned int bits, unsigned int y);
	
	// Version using precomputed table
	unsigned long long ullpow2_table(unsigned int bits, unsigned int y);
	
	// Returns the start from the first factor greater than start, from the sorted array
	unsigned int nextFactor(unsigned int start, vector<unsigned int> &factors_start);
	
public: 
	
	KarpRabin();
	KarpRabin(unsigned int _voc_bits, unsigned int _kr_mod, unsigned int _table_size = 1000000);
	virtual ~KarpRabin();
	
	// Evaluate the full hash in str.length() operations
	unsigned long long hash(const string &str);
	
	// Evaluate the full hash in str.length() operations
	unsigned long long hash(const char *str, unsigned long long str_len);
	
	// This version uses, if str_len is too big, the hash from the factors
	// Note that this version is only available during the building process
	unsigned long long hash(const char *full_text, unsigned int start, unsigned int str_len, vector<unsigned int> &factors_start);
	
	// Evaluate the hash of the concatenation in constant time
	unsigned long long concat(unsigned long long kr1, unsigned long long kr2, unsigned int len2);
	
	// Evaluate the hash of the subtract in constant time
	unsigned long long subtract_prefix(unsigned long long kr12, unsigned long long kr1, unsigned int len2);
	
};







#endif //_TEXT_FILTER_H





