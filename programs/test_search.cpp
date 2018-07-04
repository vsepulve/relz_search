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

void recursive_rmq_v2(unsigned int ini, unsigned int fin, unsigned int crit, rmq_succinct_sct<false, bp_support_sada<256,32,rank_support_v5<> > > &rmq, rrr_vector<127>::select_1_type &select1_s, rrr_vector<127>::select_1_type &select1_b, rrr_vector<127>::select_0_type &select0_b, inv_perm_support<> &perm){
	cout << " -> recursive_rmq(" << ini << ", " << fin << ")\n";
	
	unsigned int pos_max = rmq(ini, fin);
	
	unsigned int tu = select1_s(pos_max + 1) - pos_max;
	unsigned int pu = select1_b(perm[pos_max] + 1);
	unsigned int lu = select1_b(perm[pos_max] + 2) - pu;
	
//	cout << " -> max pos Ez: " << pos_max << " (Ez: " << ez[pos_max] << ", factor: " << perm[pos_max] << ")\n";
	cout << " -> max pos Ez: " << pos_max << " (tu: " << tu << ", pu: " << pu << ", lu: " << lu << ")\n";
	if( tu + lu < crit ){
		cout << "Omitiendo\n";
		return;
	}
	else{
		cout << "Agregando\n";
	}
	
	if( (pos_max > 0) && (ini < pos_max) ){
		recursive_rmq_v2(ini, pos_max-1, crit, rmq, select1_s, select1_b, select0_b, perm);
	}
	if( pos_max < fin ){
		recursive_rmq_v2(pos_max+1, fin, crit, rmq, select1_s, select1_b, select0_b, perm);
	}
}

class FactorsIterator {

private: 
	
	rrr_vector<127>::select_1_type *select1_s;
	rrr_vector<127>::select_1_type *select1_b;
	rrr_vector<127>::select_0_type *select0_b;
	inv_perm_support<> *perm;
	char *ref_text;
	
	// Posicion actual en la referencia para el factor actual
	unsigned int cur_pos;
	
	// Factor actual, con posiciones para la referencia
	unsigned int start_f;
	unsigned int n_factors;
	unsigned int cur_f;
	unsigned int cur_f_ini;
	unsigned int cur_f_fin;
	
	void loadFactor(unsigned int f){
		cout << "FactorsIterator::loadFactor\n";
		cur_f = f;
		unsigned int tu = select1_s->operator()(cur_f + 1) - cur_f;
		unsigned int pu = select1_b->operator()( (*perm)[cur_f] + 1);
		unsigned int lu = select1_b->operator()( (*perm)[cur_f] + 2) - pu;
		cout << "FactorsIterator::loadFactor - tu: " << tu << ", pu: " << pu << ", lu: " << lu << "\n";
		cur_f_ini = tu;
		cur_f_fin = tu + lu - 1;
		cur_pos = cur_f_ini;
	}
	
public: 
	
	FactorsIterator( unsigned int _start_f, unsigned int _n_factors, 
			rrr_vector<127>::select_1_type *_select1_s, 
			rrr_vector<127>::select_1_type *_select1_b, 
			rrr_vector<127>::select_0_type *_select0_b, 
			inv_perm_support<> *_perm, 
			char *_ref_text ){
		start_f = _start_f;
		n_factors = _n_factors;
		select1_s = _select1_s;
		select1_b = _select1_b;
		select0_b = _select0_b;
		ref_text = _ref_text;
		cur_pos = 0;
		cur_f = 0;
		cur_f_ini = 0;
		cur_f_fin = 0;
		reset();
	}
	
	void reset(){
		cout << "FactorsIterator::reset\n";
		loadFactor(start_f);
	}
	
	char next(){
		cout << "FactorsIterator::next - cur_pos: " << cur_pos << ", cur_f_fin: " << cur_f_fin << ", cur_f: " << cur_f << " / " << n_factors << "\n";
//		return 0;
		char ret = ref_text[cur_pos++];
		if( (cur_pos > cur_f_fin) && cur_f < n_factors ){
			loadFactor(++cur_f);
		}
		return ret;
	}
	
	bool hasNext(){
		cout << "FactorsIterator::hasNext - " << cur_pos << " <= " << cur_f_fin << "?\n";
//		return false;
		if( cur_pos <= cur_f_fin ){
			return true;
		}
		return false;
	}
	
};

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
	unsigned int len_text = text.length();
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
	rrr_vector<127>::select_1_type select1_s(&rrr_s);
	cout << "Posicion de primer 1: " << select1_s(1) << "\n";
	cout << "Posicion de tercer 1: " << select1_s(3) << "\n";
	cout << "Posicion de quinto 1: " << select1_s(5) << "\n";
	
	// Notar que la posicion del select DEBE empezar desde 1, no desde 0
	// De este modo, hay que sumar 1 a las posiciones de la ref para buscar en S
	rrr_vector<127>::select_0_type select0_s(&rrr_s);
	cout << "Posicion de primer 0: " << select0_s(1) << "\n";
	cout << "Posicion de tercer 0: " << select0_s(3) << "\n";
	cout << "Posicion de quinto 0: " << select0_s(5) << "\n";
	cout << "Posicion de 0th 0: " << select0_s(0) << "\n";
	
	
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
	
	// Bit vector B (inicio de las frases en texto)
	bit_vector arr_b(len_text, 0);
	unsigned int pos_text = 0;
	for( unsigned int i = 0; i < z; ++i ){
//		unsigned int ini = factors[i].first;
		unsigned int len = factors[i].second;
		arr_b[ pos_text ] = 1;
		pos_text += len;
	}
	cout << "Bit Vector B: \n";
	for( unsigned int i = 0; i < len_text; ++i ){
		cout << "B[" << i << "]: " << arr_b[i] << "\n";
	}
	rrr_vector<127> rrr_b(arr_b);
	rrr_vector<127>::select_1_type select1_b(&rrr_b);
	rrr_vector<127>::select_0_type select0_b(&rrr_b);
	
	// Construccion del FM Index (aunque use el SA original para la compresion, esto es para la busqueda)
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
			unsigned int select = select0_s(occ_i + 1);
			unsigned int pos_ez = select - 1 - occ_i;
			cout << "select: " << select << " => pos_ez: " << pos_ez << "\n";
			
			// Ahora la busqueda (recursiva) en el rmq (entre 0 y pos_ez)
//			unsigned int pos_max = rmq(0, pos_ez);
//			cout << "max pos Ez: " << pos_max << " (Ez: " << ez[pos_max] << ", factor: " << perm[pos_max] << ")\n";
			cout << "----- Search V1 -----\n";
			recursive_rmq(0, pos_ez, (occ_i + m - 1), rmq, ez, perm);
			cout << "----- Search V2 -----\n";
			recursive_rmq_v2(0, pos_ez, (occ_i + m), rmq, select1_s, select1_b, select0_b, perm);
			cout << "----- -----\n";
			
		}
	}
	
	// Indice secundario
	cout << "Preparando estructuras de indice secundario\n";
	
	// Creo que seria ideal preparar iteradores de factor directo y reverso, que accedan a la referencia
	// Esas estructuras deberian poder usar las estructuras comprimidas para evaluar la info de los factores
	// Cada iterador internamente puede mantener los valores actuales dada un cur_pos actual
	// Basta con que los iteradores retornen el char de cierta pos FactorsIterator::get(unsigned int pos) o "char FactorsIterator::next()"
	// A parte del next, necesitaria una forma de controlar el final del iterator, quizas "bool FactorsIterator::hasNext()"
	
	
	



	delete reference;




}




















