#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>

#include <algorithm>
#include <vector>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rmq_support.hpp>
#include <sdsl/inv_perm_support.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bp_support_g.hpp>
#include <sdsl/bp_support_gg.hpp>

using namespace std;
using namespace sdsl;

#include "KarpRabin.h"
#include "FactorsConfig.h"
#include "KarpRabinFactorsSuffixesv2.h"

int main() {
	
	cout << "Inicio\n";
	
	
	unsigned int bits = 8;
	unsigned int mod = 15485863;
	KarpRabin karp_rabin(bits, mod, 1000000);
//	KarpRabin karp_rabin(bits, mod, 1100000000);
	
	string file_base = "index_yeast_compacted_test1/index_hash_compacted";
	
	string index_basic_file = file_base + ".base";
	fstream reader(index_basic_file, fstream::in);
	// Version of the index
	unsigned char version = 0;
	reader.read((char*)&version, 1);
	if( version != 3 ){
		cout << "Wrong Version\n";
		return 0;
	}
	unsigned int len_text = 0;
	reader.read((char*)&len_text, sizeof(int));
	unsigned int n_factors = 0;
	reader.read((char*)&n_factors, sizeof(int));
	unsigned int len_ref = 0;
	reader.read((char*)&len_ref, sizeof(int));
	cout << "len_text: " << len_text << ", n_factors: " << n_factors << ", len_ref: " << len_ref << "\n";
	char *ref_text = new char[len_ref + 1];
	reader.read((char*)ref_text, len_ref);
	ref_text[len_ref] = 0;
	reader.close();
	
	bits_s_type bits_s;
	bits_s_type::select_1_type select1_s;
	string bits_s_file = file_base + ".arrs";
	load_from_file(bits_s, bits_s_file);
	select1_s = bits_s_type::select_1_type(&bits_s);
	
	bits_b_type bits_b;
	bits_b_type::select_1_type select1_b;
	bits_b_type::select_0_type select0_b;
	string bits_b_file = file_base + ".arrb";
	load_from_file(bits_b, bits_b_file);
	select1_b = bits_b_type::select_1_type(&bits_b);
	select0_b = bits_b_type::select_0_type(&bits_b);
	
	int_vector<> pi_inv;
	string pi_inv_file = file_base + ".pi_inv";
	load_from_file(pi_inv, pi_inv_file);
	
	string krs_file = file_base + ".krsuffixes";
	
	KarpRabinFactorsSuffixesv2 kr_factors(krs_file, &karp_rabin, ref_text, &select1_s, &select1_b, &select0_b, &pi_inv);
	
	unsigned int hash = kr_factors.hash(1000, 1000000, 1000000);
	cout << "hash: " << hash << "\n";
	
	
	/*
	unsigned int n_tests = 1000;
	unsigned int min_len = 500000;
	
	vector<unsigned int> hash_table;
	for(unsigned int i = 0; i < n_tests; ++i){
		unsigned int hash = kr_factors.hash(1000, 1000000, min_len * i);
		hash_table.push_back(hash);
	}
	
	NanoTimer timer;
	for(unsigned int i = 0; i < n_tests; ++i){
		unsigned int hash = kr_factors.hash(1000, 1000000, min_len * i);
		if( hash != hash_table[i] ){
			cerr << "Error - kr_factors.hash(1000, 1000000, " << (min_len * i) << "): " << hash << " != " << hash_table[i] << "\n";
			exit(0);
		}
	}
	cout << "Time 1: " << timer.getMilisec() << "\n";
	timer.reset();
	
	for(unsigned int i = 0; i < n_tests; ++i){
		unsigned int hash = kr_factors.hashBin(1000, 1000000, min_len * i);
		if( hash != hash_table[i] ){
			cerr << "Error - kr_factors.hashBin(1000, 1000000, " << (min_len * i) << "): " << hash << " != " << hash_table[i] << "\n";
			exit(0);
		}
	}
	cout << "Time 2: " << timer.getMilisec() << "\n";
	timer.reset();
	
	for(unsigned int i = 0; i < n_tests; ++i){
		unsigned int hash = kr_factors.hashFast(1000, 1000000, min_len * i);
		if( hash != hash_table[i] ){
			cerr << "Error - kr_factors.hashFast(1000, 1000000, " << (min_len * i) << "): " << hash << " != " << hash_table[i] << "\n";
			exit(0);
		}
	}
	cout << "Time 3: " << timer.getMilisec() << "\n";
	timer.reset();
	*/
	
	
	
//	cout << "kr_factors.hash: " << kr_factors.hash(1000, 1000000, 500000000) << "\n";
//	cout << "kr_factors.hashBin: " << kr_factors.hashBin(1000, 1000000, 500000000) << "\n";
//	cout << "kr_factors.hashFast: " << kr_factors.hashFast(1000, 1000000, 500000000) << "\n";
//	
//	cout << "kr_factors.hash: " << kr_factors.hash(6647, 1000000, 500000000) << "\n";
//	cout << "kr_factors.hashBin: " << kr_factors.hashBin(6647, 1000000, 500000000) << "\n";
//	cout << "kr_factors.hashFast: " << kr_factors.hashFast(6647, 1000000, 500000000) << "\n";
//	
//	cout << "Times: " << kr_factors.nano1 << ", " << kr_factors.nano2 << ", " << kr_factors.nano3 << ", " << kr_factors.nano4 << "\n";
	
	/*
	string input_file = "index_yeast_compacted_test1/index_hash_compacted.tree_y.first_old";
	string output_file = "index_yeast_compacted_test1/index_hash_compacted.tree_y.first";
	
	// Test de conversion de archivo first, de char a int compactado
	fstream reader(input_file, fstream::in);
	
	reader.seekg(0, reader.end);
	unsigned int length = reader.tellg();
	reader.seekg(0, reader.beg);
	
	vector<char> first_chars;
	for(unsigned int i = 0; i < length; ++i){
		char c;
		reader.read(&c, 1);
		first_chars.push_back(c);
	}
	reader.close();
	
	for(unsigned int i = 0; i < 10; ++i){
		cout << "first_chars[" << i << "]: " << first_chars[i] << "\n";
	}
	for(unsigned int i = length-10; i < length; ++i){
		cout << "first_chars[" << i << "]: " << first_chars[i] << "\n";
	}
	
	int_vector<> first_compacetd(length);
	for(unsigned int i = 0; i < length; ++i){
		first_compacetd[i] = 0;
		if( first_chars[i] == 'C' ){
			first_compacetd[i] = 1;
		}
		else if( first_chars[i] == 'G' ){
			first_compacetd[i] = 2;
		}
		else if( first_chars[i] == 'T' ){
			first_compacetd[i] = 3;
		}
	}
	sdsl::util::bit_compress(first_compacetd);
	
	cout << "bits/node hash_childs: " << (8.0*size_in_bytes(first_compacetd)/length) << "\n";
	
	store_to_file(first_compacetd, output_file);
	*/
	
	/*
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
	*/
	

	cout << "Fin\n";
}




















