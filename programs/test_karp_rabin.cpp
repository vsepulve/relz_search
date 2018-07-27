#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>

#include <algorithm>
#include <vector>

//#include <sdsl/suffix_arrays.hpp>
//#include <sdsl/bit_vectors.hpp>
//#include <sdsl/rmq_support.hpp>
//#include <sdsl/inv_perm_support.hpp>
//#include <sdsl/wavelet_trees.hpp>

//#include "ReferenceIndexBasic.h"
//#include "CompressorSingleBuffer.h"
//#include "CoderBlocksRelz.h"
//#include "DecoderBlocksRelz.h"
//#include "TextFilterFull.h"

//#include "FactorsIterator.h"
//#include "FactorsIteratorComparator.h"

//using namespace sdsl;
using namespace std;

#define MOD 787

unsigned long long ullpow(int x, int y){
	unsigned long long ret = 1;
	for(int i = 0; i < y; i++){
		ret = (ret*x) % MOD;
	}
	return ret;
}

//unsigned long long ullpow2(int bits, int y){
//	return ((1ull << ((bits*y) % 64) ) % MOD);
//}

unsigned long long ullpow2_rec(unsigned int bits, unsigned int y){
	if( bits * y < 64 ){
		return ((1ull << (bits*y) ) % MOD);
	}
	else{
//		cout << "Rec (" << (bits * y) << " => " << (y/2) << " | " << (y/2 + (y & 0x1)) << ")\n";
		return ( ullpow2_rec(bits, y/2) * ullpow2_rec(bits, y/2 + (y & 0x1)) ) % MOD;
	}
}

unsigned long long rolling_hash(string &str){
	unsigned long long hash = 0;
	size_t str_len = str.length();
	for(unsigned int i = 0, k = str_len-1; i < str_len; i++, k--) {
//		hash = hash + (str[i] * ullpow(256, k)) % MOD;
		hash = hash + (str[i] * ullpow2_rec(8, k)) % MOD;
		hash = hash % MOD;
	}
	return hash;
}

unsigned long long rolling_hash_concat(unsigned long long kb1, unsigned long long kb2, unsigned int len2){
//	return (kb2 + (kb1 * ullpow(256, len2)) % MOD) % MOD;
	return (kb2 + (kb1 * ullpow2_rec(8, len2)) % MOD) % MOD;
}

//unsigned long long karp_rabin_direct(string &s, unsigned int param_p = 787){
//	unsigned long long res = 0;
//	// 8 para alfabeto de 256 = 2^8
//	// Notar que quizas sea razonable usar base 2 para alfabeto 4 = 2^2
//	unsigned long long bits = 8;
//	unsigned long long b = 1;
//	for( unsigned int i = 0; i < s.length(); ++i ){
//		unsigned long long part = (unsigned long long)( (unsigned char)s[i] ) * b;
//		part %= param_p;
//		res += part;
//		b <<= bits;
//		if( b == 0 ){
//			b = 1;
//		}
//	}
//	
//	return res;
//}

//unsigned long long karp_rabin_concat(unsigned long long kb1, unsigned long long kb2, unsigned int len1, unsigned int param_p = 787){
//	unsigned long long bits = 8;
//	unsigned long long b = 1;
//	b <<= ((bits * len1) % 64);
//	return kb1 + (kb2 * b) % param_p;
//}

int main() {
	
	cout << "Inicio\n";
	
	for( unsigned int i = 0; i < 100; ++i ){
	cout << "ullpow(256, " << i << "): " << ullpow(256, i) << " / " << ullpow2_rec(8, i) << "\n";
	}
	
	string s;
	
	s = "A";
	cout << "karp_rabin_direct(\"" << s << "\"): " << rolling_hash(s) << "\n";
	
	s = "AB";
	cout << "karp_rabin_direct(\"" << s << "\"): " << rolling_hash(s) << "\n";
	
	s = "ABC";
	cout << "karp_rabin_direct(\"" << s << "\"): " << rolling_hash(s) << "\n";
	
	s = "LA";
	cout << "karp_rabin_direct(\"" << s << "\"): " << rolling_hash(s) << "\n";
	
	s = "LALA";
	cout << "karp_rabin_direct(\"" << s << "\"): " << rolling_hash(s) << "\n";
	
	string s1, s2;
	unsigned long long kb1, kb2, kb_real, kb_concat;
	
	s1 = "LA";
	s2 = "LA";
	s = s1;
	s += s2;
	kb1 = rolling_hash(s1);
	kb2 = rolling_hash(s2);
	kb_real = rolling_hash(s);
	kb_concat = rolling_hash_concat( kb1, kb2, s2.length() );
	cout << "karp_rabin_concat(\"" << s1 << "\", \"" << s2 << "\"): " <<  kb_concat << " / " << kb_real << " (kb1: " << kb1 << ", kb2: " << kb2 << ")\n";
	
	s1 = "LA";
	s2 = "LANDA";
	s = s1;
	s += s2;
	kb1 = rolling_hash(s1);
	kb2 = rolling_hash(s2);
	kb_real = rolling_hash(s);
	kb_concat = rolling_hash_concat( kb1, kb2, s2.length() );
	cout << "karp_rabin_concat(\"" << s1 << "\", \"" << s2 << "\"): " <<  kb_concat << " / " << kb_real << " (kb1: " << kb1 << ", kb2: " << kb2 << ")\n";
	
	s1 = "ALABA";
	s2 = "LAALABARDA";
	s = s1;
	s += s2;
	kb1 = rolling_hash(s1);
	kb2 = rolling_hash(s2);
	kb_real = rolling_hash(s);
	kb_concat = rolling_hash_concat( kb1, kb2, s2.length() );
	cout << "karp_rabin_concat(\"" << s1 << "\", \"" << s2 << "\"): " <<  kb_concat << " / " << kb_real << " (kb1: " << kb1 << ", kb2: " << kb2 << ")\n";
	
	s1 = "asdafasfd";
	s2 = "ldfdfkjvndkfjnvdkjfnv";
	s = s1;
	s += s2;
	kb1 = rolling_hash(s1);
	kb2 = rolling_hash(s2);
	kb_real = rolling_hash(s);
	kb_concat = rolling_hash_concat( kb1, kb2, s2.length() );
	cout << "karp_rabin_concat(\"" << s1 << "\", \"" << s2 << "\"): " <<  kb_concat << " / " << kb_real << " (kb1: " << kb1 << ", kb2: " << kb2 << ")\n";
	
	s1 = "asdafasfdldfdfkjvndkfjnvdkjfnvldfdfkjvndkfjnvdkjfnv";
	s2 = "ldfdfkjvndkfjnvdkjfnvldfdfkjvndkfjnvdkjfnvldfdfkjvndkfjnvdkjfnvldfdfkjvndkfjnvdkjfnv";
	s = s1;
	s += s2;
	kb1 = rolling_hash(s1);
	kb2 = rolling_hash(s2);
	kb_real = rolling_hash(s);
	kb_concat = rolling_hash_concat( kb1, kb2, s2.length() );
	cout << "karp_rabin_concat(\"" << s1 << "\", \"" << s2 << "\"): " <<  kb_concat << " / " << kb_real << " (kb1: " << kb1 << ", kb2: " << kb2 << ")\n";

	cout << "Fin\n";
}




















