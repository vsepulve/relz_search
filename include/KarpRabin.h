#ifndef _KARP_RABIN_H
#define _KARP_RABIN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

class KarpRabin{

protected: 
	unsigned int voc_bits;
	unsigned int kr_mod;
	
public: 
	
	// Direct version (*x, y times)
	unsigned long long ullpow(unsigned long long x, unsigned int y);
	
	// Base 2 bits version, recursive to avoid overflow
	unsigned long long ullpow2_rec(unsigned int bits, unsigned int y);
	
	KarpRabin();
	KarpRabin(unsigned int _voc_bits, unsigned int _kr_mod);
	virtual ~KarpRabin();
	
	// Evaluate the full hash in str.length() operations
	unsigned long long hash(const string &str);
	
	// Evaluate the hash of the concatenation in constant time
	unsigned long long concat(unsigned long long kr1, unsigned long long kr2, unsigned int len2);
	
};







#endif //_TEXT_FILTER_H





