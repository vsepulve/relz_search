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
//	unsigned int mod = 15485863;
//	unsigned int bits = 2;
	unsigned int mod = 787;
//	unsigned int bits = 2;
//	unsigned int mod = 223;
	KarpRabin karp_rabin(bits, mod);
	
//	unsigned int base = 1<<bits;
//	for( unsigned int i = 0; i < 10000; ++i ){
//		cout << "ullpow(" << base << ", " << i << "): " << karp_rabin.ullpow(base, i) << " / " << karp_rabin.ullpow2_rec(bits, i) << " / " << karp_rabin.ullpow2_table(bits, i) << "\n";
//		if( karp_rabin.ullpow(base, i) != karp_rabin.ullpow2_rec(bits, i) ){
//			cout << "Error!\n";
//			return 0;
//		}
//		if( karp_rabin.ullpow(base, i) != karp_rabin.ullpow2_table(bits, i) ){
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

	
	s1 = "ALA";
	s2 = "BARDA";
	string s12 = s1;
	s12 += s2;
	unsigned long long kr1_real = karp_rabin.hash(s1);
	cout << "kr1_real: " << kr1_real << " (len1 " << s1.length() << ")\n";
	unsigned long long kr2_real = karp_rabin.hash(s2);
	cout << "kr2_real: " << kr2_real << " (len2 " << s2.length() << ")\n";
	unsigned long long kr12_real = karp_rabin.hash(s12);
	cout << "kr12_real: " << kr12_real << "\n";
//	unsigned long long kr_sub2 = karp_rabin.subtract( kr12_real, kr2_real, s1.length() );
//	cout << "kr12 - kr2: " << kr_sub2 << " / " << kr1_real << "\n";
	unsigned long long kr_sub1 = karp_rabin.subtract_prefix( kr12_real, kr1_real, s2.length() );
	cout << "kr12 - kr1: " << kr_sub1 << " / " << kr2_real << "\n";
	
	

	cout << "Fin\n";
}




















