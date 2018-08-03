#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>

#include <algorithm>
#include <vector>

using namespace std;

#include "KarpRabin.h"

int main() {
	
	cout << "Inicio\n";
	
	unsigned int bits = 8;
	unsigned int mod = 787;
//	unsigned int bits = 2;
//	unsigned int mod = 787;
//	unsigned int bits = 2;
//	unsigned int mod = 223;
	KarpRabin karp_rabin(bits, mod);
	
//	unsigned int base = 1<<bits;
//	for( unsigned int i = 0; i < 100; ++i ){
//		cout << "ullpow(" << base << ", " << i << "): " << karp_rabin.ullpow(base, i) << " / " << karp_rabin.ullpow2_rec(bits, i) << "\n";
//		if( karp_rabin.ullpow(base, i) != karp_rabin.ullpow2_rec(bits, i) ){
//			cout << "Error!\n";
//			return 0;
//		}
//	}
	
	string s;
	
	s = "A";
	cout << "karp_rabin_direct(\"" << s << "\"): " << karp_rabin.hash(s) << "\n";
	
	s = "AB";
	cout << "karp_rabin_direct(\"" << s << "\"): " << karp_rabin.hash(s) << "\n";
	
	s = "ABC";
	cout << "karp_rabin_direct(\"" << s << "\"): " << karp_rabin.hash(s) << "\n";
	
	s = "LA";
	cout << "karp_rabin_direct(\"" << s << "\"): " << karp_rabin.hash(s) << "\n";
	
	s = "LALA";
	cout << "karp_rabin_direct(\"" << s << "\"): " << karp_rabin.hash(s) << "\n";
	
	string s1, s2;
	unsigned long long kr1, kr2, kr_real, kr_concat;
	
	s1 = "LA";
	s2 = "LA";
	s = s1;
	s += s2;
	kr1 = karp_rabin.hash(s1);
	kr2 = karp_rabin.hash(s2);
	kr_real = karp_rabin.hash(s);
	kr_concat = karp_rabin.concat( kr1, kr2, s2.length() );
	cout << "karp_rabin_concat(\"" << s1 << "\", \"" << s2 << "\"): " <<  kr_concat << " / " << kr_real << " (kr1: " << kr1 << ", kr2: " << kr2 << ")\n";
	
	s1 = "LA";
	s2 = "LANDA";
	s = s1;
	s += s2;
	kr1 = karp_rabin.hash(s1);
	kr2 = karp_rabin.hash(s2);
	kr_real = karp_rabin.hash(s);
	kr_concat = karp_rabin.concat( kr1, kr2, s2.length() );
	cout << "karp_rabin_concat(\"" << s1 << "\", \"" << s2 << "\"): " <<  kr_concat << " / " << kr_real << " (kr1: " << kr1 << ", kr2: " << kr2 << ")\n";
	
	s1 = "ALABA";
	s2 = "LAALABARDA";
	s = s1;
	s += s2;
	kr1 = karp_rabin.hash(s1);
	kr2 = karp_rabin.hash(s2);
	kr_real = karp_rabin.hash(s);
	kr_concat = karp_rabin.concat( kr1, kr2, s2.length() );
	cout << "karp_rabin_concat(\"" << s1 << "\", \"" << s2 << "\"): " <<  kr_concat << " / " << kr_real << " (kr1: " << kr1 << ", kr2: " << kr2 << ")\n";
	
	s1 = "asdafasfd";
	s2 = "ldfdfkjvndkfjnvdkjfnv";
	s = s1;
	s += s2;
	kr1 = karp_rabin.hash(s1);
	kr2 = karp_rabin.hash(s2);
	kr_real = karp_rabin.hash(s);
	kr_concat = karp_rabin.concat( kr1, kr2, s2.length() );
	cout << "karp_rabin_concat(\"" << s1 << "\", \"" << s2 << "\"): " <<  kr_concat << " / " << kr_real << " (kr1: " << kr1 << ", kr2: " << kr2 << ")\n";
	
	s1 = "asdafasfdldfdfkjvndkfjnvdkjfnvldfdfkjvndkfjnvdkjfnv";
	s2 = "ldfdfkjvndkfjnvdkjfnvldfdfkjvndkfjnvdkjfnvldfdfkjvndkfjnvdkjfnvldfdfkjvndkfjnvdkjfnv";
	s = s1;
	s += s2;
	kr1 = karp_rabin.hash(s1);
	kr2 = karp_rabin.hash(s2);
	kr_real = karp_rabin.hash(s);
	kr_concat = karp_rabin.concat( kr1, kr2, s2.length() );
	cout << "karp_rabin_concat(\"" << s1 << "\", \"" << s2 << "\"): " <<  kr_concat << " / " << kr_real << " (kr1: " << kr1 << ", kr2: " << kr2 << ")\n";

	cout << "Fin\n";
}




















