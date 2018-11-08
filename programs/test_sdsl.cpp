#include <stdio.h>
#include <stdlib.h>
#include <random>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rmq_support.hpp>
#include <sdsl/inv_perm_support.hpp>
#include <sdsl/wavelet_trees.hpp>

using namespace sdsl;
using namespace std;

int main() {
//	csa_wt<> fm_index;
//	construct_im(fm_index, "mississippi!", 1);
//	std::cout << "'si' occurs " << count(fm_index,"si") << " times.\n";
//	store_to_file(fm_index,"fm_index-file.sdsl");
//	std::ofstream out("fm_index-file.sdsl.html");
//	write_structure<HTML_FORMAT>(fm_index,out);
	
	string ref = "alabardas";
	string text = "alba";
	text += "balas";
	text += "lalabas";
	string pat = "ba";
	unsigned int len_r = 9;
	unsigned int z = 8;
	
	bit_vector arr_s(len_r + z, 0);
	arr_s[0] = 1;
	arr_s[2] = 1;
	arr_s[3] = 1;
	arr_s[4] = 1;
	arr_s[7] = 1;
	arr_s[8] = 1;
	arr_s[14] = 1;
	arr_s[15] = 1;
	
	rrr_vector<127> rrr_s(arr_s);
	rrr_vector<127>::select_1_type select1(&rrr_s);
	cout << "Posicion de primer 1: " << select1(1) << "\n";
	cout << "Posicion de tercer 1: " << select1(3) << "\n";
	cout << "Posicion de quinto 1: " << select1(5) << "\n";
	
	// Notar que la posicion del select DEBE empezar desde 1, no desde 0
	// De este modo, hay que sumar 1 a las posiciones de la ref para buscar en S
	rrr_vector<127>::select_0_type select0(&rrr_s);
	cout << "Posicion de primer 0: " << select0(1) << "\n";
	cout << "Posicion de tercer 0: " << select0(3) << "\n";
	cout << "Posicion de quinto 0: " << select0(5) << "\n";
	cout << "Posicion de 0th 0: " << select0(0) << "\n";
	
	// En esta primera prueba, almaceno la permutacion de z descomprimida
//	vector<unsigned int> pi(8);
	int_vector<> pi(8);
	pi[0] = 0;
	pi[1] = 3;
	pi[2] = 5;
	pi[3] = 6;
	pi[4] = 1;
	pi[5] = 2;
	pi[6] = 4;
	pi[7] = 7;
	inv_perm_support<> perm(&pi);
	cout << "pi[0]: " << pi[0] << ", perm[0]: " << perm[0] << "\n";
	cout << "pi[2]: " << pi[2] << ", perm[2]: " << perm[2] << "\n";
	cout << "pi[5]: " << pi[5] << ", perm[5]: " << perm[5] << "\n";
	
	// Para esta prueba tambien guardo ez descomprimido
//	vector<unsigned int> ez(8);
	int_vector<> ez(8);
	ez[0] = 1;
	ez[1] = 2;
	ez[2] = 2;
	ez[3] = 4;
	ez[4] = 4;
	ez[5] = 4;
	ez[6] = 8;
	ez[7] = 8;
	// rmq_succinct_sct<> rmq(&ez);
	rmq_succinct_sct<false, bp_support_sada<256,32,rank_support_v5<> > > rmq(&ez);
	// rmq_maximum_sct<> rmq(&ez);
	
	csa_wt<> fm_index;
	// Construccion con datos en memoria, en un string
	construct_im(fm_index, ref, 1);
	// Construccion con datos en un archivo
//	construct(fm_index, file, 1);
	
	cout << "Texto de ref: \"" << ref << "\"\n";
	
	string query = "ba";
	size_t m = query.size();
	size_t occs = sdsl::count(fm_index, query.begin(), query.end());
	cout << "# occs de \"" << query << "\": " << occs << "\n";
	if( occs > 0 ){
		auto locations = locate(fm_index, query.begin(), query.begin()+m);
		sort(locations.begin(), locations.end());
		for( unsigned int i = 0; i < occs; ++i ){
			unsigned int occ_i = locations[i];
			cout << "occ[" << i << "]: " << occ_i << " (" << ref.substr(occ_i, m) << ")\n";
			// Comprobar los factores que cuben esta ocurrencia (el string ref[occ_i, occ_i + m - 1])
			unsigned int select = select0(occ_i + 1);
			unsigned int pos_ez = select - 1 - occ_i;
			cout << "select: " << select << " => pos_ez: " << pos_ez << "\n";
			
			// Ahora la busqueda (recursiva) en el rmq (entre 0 y pos_ez)
			unsigned int max = rmq(0, pos_ez);
			cout << "max: " << max << "\n";

		}
	}
	
//	std::cout << "'si' occurs " << count(fm_index,"si") << " times.\n";
//	store_to_file(fm_index,"fm_index-file.sdsl");
//	std::ofstream out("fm_index-file.sdsl.html");
//	write_structure<HTML_FORMAT>(fm_index,out);
	
	
	cout << "Prueba de WT\n";
	
	wt_int<rrr_vector<63>> wt;
	construct_im(wt, int_vector<> {0, 6, 3, 5, 7, 1, 2, 4});
	// La consulta de rango recibe [min_x, max_x, min_y, max_y] inclusivos (es decir, horizontal primero, vertical despues)
	// Los resultados son pares [x, y] (horizontal primero, vertical despues)
	// Los contadores de valores y posiciones empiezan en 0
	cout << "Buscando en [1, 2] x [1, 5]:\n";
	auto res = wt.range_search_2d(1, 2, 1, 5);
	for (auto point : res.second){
		cout << "(" << point.first << ", " << point.second << ")\n";
	}
	
	cout << "Buscando en [1, 3] x [1, 5]:\n";
	res = wt.range_search_2d(1, 3, 1, 5);
	for (auto point : res.second){
		cout << "(" << point.first << ", " << point.second << ")\n";
	}
	
	
	unsigned int n_ints = 100000;
	unsigned int max_int = 100000;
	cout << "n_ints: " << n_ints << ", max_int: " << max_int << "\n";
	
	int_vector<32> arr_32(n_ints);
	cout << "size(arr_32): " << size_in_bytes(arr_32) << " bytes (" << ((double)(size_in_bytes(arr_32))*8.0/n_ints) << " bits/int)\n";
	
	int_vector<30> arr_30(n_ints);
	cout << "size(arr_30): " << size_in_bytes(arr_30) << " bytes (" << ((double)(size_in_bytes(arr_30))*8.0/n_ints) << " bits/int)\n";
	
	int_vector<28> arr_28(n_ints);
	cout << "size(arr_28): " << size_in_bytes(arr_28) << " bytes (" << ((double)(size_in_bytes(arr_28))*8.0/n_ints) << " bits/int)\n";
	
	int_vector<15> arr_15(n_ints);
	cout << "size(arr_15): " << size_in_bytes(arr_15) << " bytes (" << ((double)(size_in_bytes(arr_15))*8.0/n_ints) << " bits/int)\n";

	
	int_vector<> arr(n_ints);
	
	random_device seed;
	mt19937 generator(seed());
	uniform_int_distribution<> dist(0, max_int - 1);
	
	for(unsigned int i = 0; i < n_ints; ++i){
		arr[i] = dist(generator);
	}
	
	cout << "size(arr): " << size_in_bytes(arr) << " bytes (" << ((double)(size_in_bytes(arr))*8.0/n_ints) << " bits/int)\n";
	for(unsigned int i = 0; i < ( (n_ints<10)?n_ints:10 ); ++i){
		cout << "arr[" << i << "]: " << arr[i] << "\n";
	}
	
	sdsl::util::bit_compress(arr);
	cout << "size(arr): " << size_in_bytes(arr) << " bytes (" << ((double)(size_in_bytes(arr))*8.0/n_ints) << " bits/int)\n";
	for(unsigned int i = 0; i < ( (n_ints<10)?n_ints:10 ); ++i){
		cout << "arr[" << i << "]: " << arr[i] << "\n";
	}

}




















