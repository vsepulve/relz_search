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

#include "ReferenceIndexBasic.h"
#include "CompressorSingleBuffer.h"
#include "CoderBlocksRelz.h"
#include "DecoderBlocksRelz.h"
#include "TextFilterFull.h"

#include "FactorsIndex.h"
#include "FactorsIterator.h"
#include "FactorsIteratorComparator.h"

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[]){

	if(argc != 4){
		cout<<"\nModo de Uso: relz_search serialized_ref sequence_file output_relz\n";
		return 0;
	}
	
	const char *ref_file = argv[1];
	const char *input = argv[2];
	const char *output = argv[3];
	
	ReferenceIndex *reference = new ReferenceIndexBasic();
	reference->load(ref_file);
	
	// Preparar Compresor
	TextFilter *filter = new TextFilterFull();
	CompressorSingleBuffer compressor(
		output, 
		new CoderBlocksRelz(reference), 
		new DecoderBlocksRelz(reference->getText()), 
		filter
	);
	
	vector<pair<unsigned int, unsigned int> > factors;
	compressor.compress(input, 1, 1000000, 0, &factors);
	
	unsigned int len_ref = reference->getLength();
	unsigned int len_text = compressor.getTextSize();
//	unsigned int z = 0;
	
	const char *ref = reference->getText();
	
	/*
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
	*/
	
	 FactorsIndex index(factors, len_text, ref, len_ref);
//	 index.test("AB");
	 index.find("AB");
	
	
	
	
	
	
	/*
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
			// cout << "occ[" << i << "]: " << occ_i << " (" << ref.substr(occ_i, m) << ")\n";
			cout << "occ[" << i << "]: " << occ_i << " \n";
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
	
	
	cout << "Probando Iterador\n";
	
	FactorsIterator it( 2, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	cout << "-----\n";
	
	while( it.hasNext() ){
		char c = it.next();
		cout << "it.next(): " << c << " (text_pos " << it.position() << " / " << it.length() << ")\n";
		cout << "-----\n";
	}
	
	cout << "Fin prueba iterador\n";
	cout << "-----\n";
	
	cout << "Probando Iterador Reverso\n";
	
	FactorsIteratorReverse it_rev( 2, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	cout << "-----\n";
	
	while( it_rev.hasNext() ){
		char c = it_rev.next();
		cout << "it_rev.next(): " << c << " (text_pos " << it_rev.position() << " / " << it_rev.length() << ")\n";
		cout << "-----\n";
	}
	
	cout << "Fin prueba iterador\n";
	cout << "-----\n";

	
	cout << "Probando Iterador Reverso -1\n";
	it_rev = FactorsIteratorReverse( -1, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	while( it_rev.hasNext() ){
		char c = it_rev.next();
		cout << "it_rev.next(): " << c << " (text_pos " << it_rev.position() << " / " << it_rev.length() << ")\n";
		cout << "-----\n";
	}

	cout << "Probando Iterador Reverso 0\n";
	it_rev = FactorsIteratorReverse( 0, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	while( it_rev.hasNext() ){
		char c = it_rev.next();
		cout << "it_rev.next(): " << c << " (text_pos " << it_rev.position() << " / " << it_rev.length() << ")\n";
		cout << "-----\n";
	}

	cout << "Probando Iterador Reverso 1\n";
	it_rev = FactorsIteratorReverse( 1, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	while( it_rev.hasNext() ){
		char c = it_rev.next();
		cout << "it_rev.next(): " << c << " (text_pos " << it_rev.position() << " / " << it_rev.length() << ")\n";
		cout << "-----\n";
	}

	cout << "Probando Iterador Reverso 5\n";
	it_rev = FactorsIteratorReverse( 5, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	while( it_rev.hasNext() ){
		char c = it_rev.next();
		cout << "it_rev.next(): " << c << " (text_pos " << it_rev.position() << " / " << it_rev.length() << ")\n";
		cout << "-----\n";
	}

	cout << "Probando Iterador Reverso 7\n";
	it_rev = FactorsIteratorReverse( 7, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	while( it_rev.hasNext() ){
		char c = it_rev.next();
		cout << "it_rev.next(): " << c << " (text_pos " << it_rev.position() << " / " << it_rev.length() << ")\n";
		cout << "-----\n";
	}
	
	cout << "-----\n";
		
	cout << "Probando Iterador Reverso 8\n";
	it_rev = FactorsIteratorReverse( 8, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	while( it_rev.hasNext() ){
		char c = it_rev.next();
		cout << "it_rev.next(): " << c << " (text_pos " << it_rev.position() << " / " << it_rev.length() << ")\n";
		cout << "-----\n";
	}
	cout << "Probando reset\n";
	it_rev.reset();
	while( it_rev.hasNext() ){
		char c = it_rev.next();
		cout << "it_rev.next(): " << c << " (text_pos " << it_rev.position() << " / " << it_rev.length() << ")\n";
		cout << "-----\n";
	}

	cout << "-----\n";

	cout << "Probando Iterador Reverso 5\n";
	it_rev = FactorsIteratorReverse( 5, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	while( it_rev.hasNext() ){
		char c = it_rev.next();
		cout << "it_rev.next(): " << c << " (text_pos " << it_rev.position() << " / " << it_rev.length() << ")\n";
		cout << "-----\n";
	}
	cout << "Probando reset\n";
	it_rev.reset();
	while( it_rev.hasNext() ){
		char c = it_rev.next();
		cout << "it_rev.next(): " << c << " (text_pos " << it_rev.position() << " / " << it_rev.length() << ")\n";
		cout << "-----\n";
	}
	
	cout << "-----\n";
	
	cout << "Probando acceso posicional\n";
	unsigned int f = 5;
	for( unsigned int i = 0; i < 10; ++i ){
		char c = getChar( f, i, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
		if( c == 0){
			break;
		}
		cout << "factor_" << f << "[" << i << "]: " << c << " \n";
	}
	cout << "-----\n";
	f = 7;
	for( unsigned int i = 0; i < 10; ++i ){
		char c = getChar( f, i, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
		if( c == 0){
			break;
		}
		cout << "factor_" << f << "[" << i << "]: " << c << " \n";
	}
	cout << "-----\n";
	f = 3;
	for( unsigned int i = 0; i < 10; ++i ){
		char c = getChar( f, i, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
		if( c == 0){
			break;
		}
		cout << "factor_" << f << "[" << i << "]: " << c << " \n";
	}
	cout << "-----\n";
	
	cout << "Probando acceso posicional Reverso\n";
	f = 3;
	for( unsigned int i = 0; i < 10; ++i ){
		char c = getCharRev( f-1, i, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
		if( c == 0){
			break;
		}
		cout << "factor_rev_" << f << "[" << i << "]: " << c << " \n";
	}
	cout << "-----\n";
	f = 5;
	for( unsigned int i = 0; i < 10; ++i ){
		char c = getCharRev( f-1, i, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
		if( c == 0){
			break;
		}
		cout << "factor_rev_" << f << "[" << i << "]: " << c << " \n";
	}
	cout << "-----\n";
	cout << "factor_rev_" << f << "[" << 3 << "]: " << getCharRev( f-1, 3, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text) << " \n";
	cout << "factor_rev_" << f << "[" << 5 << "]: " << getCharRev( f-1, 5, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text) << " \n";
	cout << "factor_rev_" << f << "[" << 2 << "]: " << getCharRev( f-1, 2, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text) << " \n";
	cout << "factor_rev_" << f << "[" << 0 << "]: " << getCharRev( f-1, 0, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text) << " \n";
	cout << "-----\n";
	
	cout << "Preparando arr X\n";
	vector<unsigned int> arr_x(z);
	for( unsigned int i = 0; i < z; ++i ){
		arr_x[i] = i;
	}
	FactorsIteratorReverseComparator comp_rev(z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	stable_sort(arr_x.begin(), arr_x.end(), comp_rev);
	int_vector<> pre_x_inv(z);
	for( unsigned int i = 0; i < z; ++i ){
		pre_x_inv[ arr_x[i] ] = i;
	}
	inv_perm_support<> perm_x(&pre_x_inv);
	for( unsigned int i = 0; i < z; ++i ){
		pre_x_inv[ arr_x[i] ] = i;
		cout << " arr_x[" << i << "]: " << arr_x[i] << " (perm_x[" << i << "]: " << perm_x[i] << ") \n";
	}
	cout << "-----\n";
	
	cout << "Preparando arr Y\n";
	vector<unsigned int> arr_y(z);
	for( unsigned int i = 0; i < z; ++i ){
		arr_y[i] = i;
	}
	FactorsIteratorComparator comp(z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	stable_sort(arr_y.begin(), arr_y.end(), comp);
	int_vector<> pre_y(z);
	int_vector<> pre_y_inv(z);
	for( unsigned int i = 0; i < z; ++i ){
		pre_y[i] = arr_y[i];
		pre_y_inv[ arr_y[i] ] = i;
	}
	inv_perm_support<> perm_y_inv(&pre_y);
	inv_perm_support<> perm_y(&pre_y_inv);
	for( unsigned int i = 0; i < z; ++i ){
		cout << " arr_y[" << i << "]: " << arr_y[i] << " (perm_y[" << i << "]: " << perm_y[i] << ", perm_y_inv[" << i << "]: " << perm_y_inv[i] << ")\n";
	}
	cout << "-----\n";
	
	cout << "Preparando WT\n";
	int_vector<> values_wt(z);
	for( unsigned int i = 0; i < z; ++i ){
		values_wt[i] = perm_y_inv[ arr_x[ i ] ];
		cout << " values_wt[" << i << "]: " << values_wt[i] << " \n";
	}
	
	wt_int<rrr_vector<63>> wt;
	construct_im(wt, values_wt);
	
//	cout << "Buscando en [1, 2] x [1, 5]:\n";
//	auto res = wt.range_search_2d(1, 2, 1, 5);
//	for (auto point : res.second){
//		cout << "(" << point.first << ", " << point.second << ")\n";
//	}
	
	cout << "Buscando en [1, 5] x [1, 2]:\n";
	auto res = wt.range_search_2d(1, 5, 1, 2);
	for (auto point : res.second){
		cout << "(" << point.first << ", " << point.second << ") => factor " << perm_y[point.second] << "\n";
		// Aqui tengo el id posicional del factor
		// Puedo sacar sus datos con las formulas para tu, pu, lu
		// Es necesario leer tambien el factor anterior tambien, pero la posicion se tiene
		// Notar que las posiciones de corte en m las conozco porque itero por ella (en m^2)
	}
	
	cout << "Realizando busquedas reales en el WT\n";
	// El codigo de la busqueda de rangos deberia estar basado en el codigo de reference
	
	cout << "Prueba de patron \"A\"\n";
	getRangeY("A", perm_y, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	cout << "-----\n";
	
	cout << "Prueba de patron \"B\"\n";
	getRangeY("B", perm_y, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	cout << "-----\n";
	
	cout << "Prueba de patron \"BA\"\n";
	getRangeY("BA", perm_y, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	cout << "-----\n";
	
	cout << "Prueba de patron \"BAL\"\n";
	getRangeY("BAL", perm_y, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	cout << "-----\n";
	
	cout << "Prueba de patron \"BALB\"\n";
	getRangeY("BALB", perm_y, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	cout << "-----\n";
	
	cout << "Prueba de patron \"Z\"\n";
	getRangeY("Z", perm_y, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	cout << "-----\n";
	
	cout << "Prueba de patron \"0\"\n";
	getRangeY("0", perm_y, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	cout << "-----\n";
	
	cout << "Prueba de patron \"\"\n";
	getRangeY("", perm_y, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	cout << "-----\n";
	
	
	string pattern = "AB";
	for(unsigned int i = 1; i < pattern.length(); ++i){
		string p1 = pattern.substr(0, i);
		string p2 = pattern.substr(i, pattern.length() - i);
		cout << "Corte de \"" << pattern << "\": (" << p1 << "| " << p2 << ")\n";
		pair<unsigned int, unsigned int> r1 = getRangeX(p1.c_str(), perm_x, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
		pair<unsigned int, unsigned int> r2 = getRangeY(p2.c_str(), perm_y, z, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
		
		if( r1.second == (unsigned int)(-1) || r1.second < r1.first
			|| r2.second == (unsigned int)(-1) || r2.second < r1.first ){
			cout << "Rangos Invalidos, omitiendo...\n";
			continue;
		}
		
		cout << "Buscando en [" << r1.first << ", " << r1.second << "] x [" << r2.first << ", " << r2.second << "]:\n";
		auto res = wt.range_search_2d(r1.first, r1.second, r2.first, r2.second);
		for (auto point : res.second){
			cout << "(" << point.first << ", " << point.second << ") => factor " << perm_y[point.second] << "\n";
		}
		
	}
	*/
	
	delete reference;




}




















