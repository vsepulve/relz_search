#include "KarpRabin.h"

KarpRabin::KarpRabin(){
	voc_bits = 8;
	kr_mod = 787;
}

KarpRabin::KarpRabin(unsigned int _voc_bits, unsigned int _kr_mod){
	voc_bits = _voc_bits;
	kr_mod = _kr_mod;
}

KarpRabin::~KarpRabin(){
}
	
// Direct version (*x, y times)
unsigned long long KarpRabin::ullpow(unsigned long long x, unsigned int y){
	unsigned long long ret = 1;
	for(unsigned int i = 0; i < y; i++){
		ret = (ret*x) % kr_mod;
	}
	return ret;
}

// Base 2 bits version (2^bits^y), recursive to avoid overflow
unsigned long long KarpRabin::ullpow2_rec(unsigned int bits, unsigned int y){
	if( bits * y < 64 ){
		return ((1ull << (bits*y) ) % kr_mod);
	}
	else{
//		cout << "KarpRabin::ullpow2_rec (" << (bits * y) << " => " << (y/2) << " | " << (y/2 + (y & 0x1)) << ")\n";
		return ( ullpow2_rec(bits, y>>1) * ullpow2_rec(bits, (y>>1) + (y & 0x1)) ) % kr_mod;
	}
}

// Evaluate the full hash in str.length() operations
unsigned long long KarpRabin::hash(const string &str){
	unsigned long long ret = 0;
	size_t str_len = str.length();
	for(unsigned int i = 0, k = str_len-1; i < str_len; i++, k--) {
//		ret = ret + (str[i] * ullpow(1<<voc_bits, k)) % kr_mod;
		ret = ret + (str[i] * ullpow2_rec(voc_bits, k)) % kr_mod;
		ret = ret % kr_mod;
	}
	return ret;
}

// Evaluate the hash of the concatenation in constant time
unsigned long long KarpRabin::concat(unsigned long long kr1, unsigned long long kr2, unsigned int len2){
//	return (kb2 + (kb1 * ullpow(1<<voc_bits, len2)) % kr_mod) % kr_mod;
	return (kr2 + (kr1 * ullpow2_rec(voc_bits, len2)) % kr_mod) % kr_mod;
}


















