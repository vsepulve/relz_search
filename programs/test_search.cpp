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

#include "ReferenceIndexBasic.h"
#include "CompressorSingleBuffer.h"
#include "CoderBlocksRelz.h"
#include "DecoderBlocksRelz.h"
#include "TextFilterFull.h"

using namespace sdsl;
using namespace std;

void recursive_rmq(unsigned int ini, unsigned int fin, unsigned int crit, rmq_succinct_sct<false, bp_support_sada<256,32,rank_support_v5<> > > &rmq, int_vector<> &ez, inv_perm_support<> &perm){
	cout << " -> recursive_rmq(" << ini << ", " << fin << ")\n";
	
	unsigned int pos_max = rmq(ini, fin);
	
	cout << " -> max pos Ez: " << pos_max << " (Ez: " << ez[pos_max] << ", factor: " << perm[pos_max] << ")\n";
	if( ez[pos_max] < crit ){
		cout << "Omitiendo\n";
		return;
	}
	else{
		cout << "Agregando\n";
	}
	
	if( (pos_max > 0) && (ini < pos_max) ){
		recursive_rmq(ini, pos_max-1, crit, rmq, ez, perm);
	}
	if( pos_max < fin ){
		recursive_rmq(pos_max+1, fin, crit, rmq, ez, perm);
	}
}

int main() {
	
	string ref = "ALABARDAS";
	string text = "ALBA";
	text += "BALAS";
	text += "LALABAS";
	string input = "test_text.txt";
	string output = "test_output.relz";
	string serialized_reference = "test_ref.bin";
	
	vector<pair<unsigned int, unsigned int> > factors;
	
	unsigned int len_ref = ref.length();
	unsigned int z = 0;
	
	// Construir SA referencia
	// Esto eventualmente se convertira en un FM-index (csa_wt)
	ReferenceIndex *reference = new ReferenceIndexBasic(ref.c_str(), 1);
	reference->save(serialized_reference.c_str());
	
	// Preparar Compresor
	TextFilter *filter = new TextFilterFull();
	CompressorSingleBuffer compressor(
		output.c_str(), 
		new CoderBlocksRelz(reference), 
		new DecoderBlocksRelz(reference->getText()), 
		filter
		);
	
	compressor.compress(input.c_str(), 1, 1000000, 0, &factors);
	
	cout << "Factors: \n";
	// Factores en version ini, fin (absoluto) y ordenados por ini
//	vector<pair<unsigned int, unsigned int> > factors_sort;
	vector<pair<unsigned int, pair<unsigned int, unsigned int> > > factors_sort;
	unsigned cur_pos = 0;
	for( pair<unsigned int, unsigned int> factor : factors ){
		cout << "(" << factor.first << ", " << factor.second << ", " << cur_pos << ")\n";
		factors_sort.push_back( 
			pair<unsigned int, pair<unsigned int, unsigned int> >(
				factor.first, pair<unsigned int, unsigned int>(factor.first + factor.second - 1, cur_pos++)
				)
			);
	}
	sort(factors_sort.begin(), factors_sort.end());
	cout << "Factors Sorted: \n";
	for( pair<unsigned int, pair<unsigned int, unsigned int> > factor : factors_sort ){
		cout << "(" << factor.first << ", " << factor.second.first << ", " << factor.second.second << ")\n";
	}
	z = factors_sort.size();
	
	// Bit vector S
	bit_vector arr_s(len_ref + z, 0);
	unsigned cur_ref = 0;
	cur_pos = 0;
	for( unsigned int i = 0; i < z; ++i ){
		unsigned int ini = factors_sort[i].first;
//		unsigned int fin = factors_sort[i].second.first;
		if( ini == cur_ref ){
			arr_s[cur_pos++] = 1;
		}
		else{
			arr_s[cur_pos++] = 0;
			++cur_ref;
			--i;
		}
	}
	cout << "Bit Array S: \n";
	for( unsigned int i = 0; i < len_ref + z; ++i ){
		cout << "arr_s[" << i << "]: " << arr_s[i] << "\n";
	}
	
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
	
	
	// Permutacion 
	int_vector<> pi(z);
	int_vector<> pi_inv(z);
	for( unsigned int i = 0; i < z; ++i ){
		pi[i] = factors_sort[i].second.second;
		pi_inv[ factors_sort[i].second.second ] = i;
	}
	inv_perm_support<> perm_inv(&pi);
	inv_perm_support<> perm(&pi_inv);
	cout << "Permutation: \n";
	for( unsigned int i = 0; i < z; ++i ){
		cout << "pi[" << i << "]: " << pi[i] << ", perm[" << i << "]: " << perm[i] << ", perm_inv[" << i << "]: " << perm_inv[i] << "\n";
	}
	
	// Posiciones finales Ez
//	vector<unsigned int> ez(8);
	int_vector<> ez(z);
	for( unsigned int i = 0; i < z; ++i ){
		ez[i] = factors_sort[i].second.first;
	}
	// rmq_succinct_sct<> rmq(&ez);
	rmq_succinct_sct<false, bp_support_sada<256,32,rank_support_v5<> > > rmq(&ez);
	// rmq_maximum_sct<> rmq(&ez);
	
	csa_wt<> fm_index;
	// Construccion con datos en memoria, en un string
	construct_im(fm_index, ref, 1);
	// Construccion con datos en un archivo
//	construct(fm_index, file, 1);
	
	cout << "Texto de ref: \"" << ref << "\"\n";
	
//	string query = "LA";
	string query = "BA";
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
//			unsigned int pos_max = rmq(0, pos_ez);
//			cout << "max pos Ez: " << pos_max << " (Ez: " << ez[pos_max] << ", factor: " << perm[pos_max] << ")\n";
			recursive_rmq(0, pos_ez, (occ_i + m - 1), rmq, ez, perm);

		}
	}
	
//	std::cout << "'si' occurs " << count(fm_index,"si") << " times.\n";
//	store_to_file(fm_index,"fm_index-file.sdsl");
//	std::ofstream out("fm_index-file.sdsl.html");
//	write_structure<HTML_FORMAT>(fm_index,out);
	
	
	
	
	
	
	
	
	
	
	



	delete reference;




}




















