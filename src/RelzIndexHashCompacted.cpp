#include "RelzIndexHashCompacted.h"

RelzIndexHashCompacted::RelzIndexHashCompacted(){
	len_text = 0;
	len_ref = 0;
	n_factors = 0;
	karp_rabin = NULL;
	ref_text = NULL;
}

RelzIndexHashCompacted::RelzIndexHashCompacted(KarpRabin *_karp_rabin){
	len_text = 0;
	len_ref = 0;
	n_factors = 0;
	karp_rabin = _karp_rabin;
	ref_text = NULL;
}

RelzIndexHashCompacted::RelzIndexHashCompacted(vector<pair<unsigned int, unsigned int> > &factors, char *full_text, unsigned int _len_text, const char *_ref_text, unsigned int _len_ref, KarpRabin *_karp_rabin){
	
	len_text = _len_text;
	len_ref = _len_ref;
	n_factors = factors.size();
	karp_rabin = _karp_rabin;
	
	ref_text = new char[_len_ref + 1];
	memcpy(ref_text, _ref_text, _len_ref);
	ref_text[_len_ref] = 0;
	
	NanoTimer timer;
	
	cout << "RelzIndexHashCompacted - Inicio (factors: " << factors.size() << ", len_text: " << len_text << ", len_ref: " << len_ref << ")\n";
//	cout << "RelzIndexHashCompacted - Ref: " << ref_text << "\n";
//	cout << "RelzIndexHashCompacted - Text: " << full_text << "\n";
	
	cout << "RelzIndexHashCompacted - Preparing Factors\n";
	// Factores en version ini, fin (absoluto) y ordenados por ini
	vector<pair<unsigned int, pair<unsigned int, unsigned int> > > factors_sort;
//	vector<unsigned int> factors_start;
	unsigned int cur_start = 0;
	unsigned int cur_pos = 0;
	for( pair<unsigned int, unsigned int> factor : factors ){
//		cout << "(" << factor.first << ", " << factor.second << ", " << cur_pos << ") - cur_start: " << cur_start << "\n";
		factors_sort.push_back( 
			pair<unsigned int, pair<unsigned int, unsigned int> >(
				factor.first, pair<unsigned int, unsigned int>(factor.first + factor.second - 1, cur_pos++)
				)
			);
		factors_start.push_back(cur_start);
		cur_start += factor.second;
	}
	sort(factors_sort.begin(), factors_sort.end());
	cout << "RelzIndexHashCompacted - Factors Sorted prepared in " << timer.getMilisec() << "\n";
	timer.reset();
//	cout << "Factors Sorted: \n";
//	for( unsigned int i = 0; i < factors_sort.size(); ++i ){
//		cout << "factors_sort[" << i << "]: (" << factors_sort[i].first << ", " << factors_sort[i].second.first << ", " << factors_sort[i].second.second << ")\n";
//	}
	
	// Bit vector S
	cout << "RelzIndexHashCompacted - Preparing Vector S\n";
	bit_vector arr_s(len_ref + n_factors, 0);
	unsigned cur_ref = 0;
	cur_pos = 0;
	for( unsigned int i = 0; i < n_factors; ++i ){
		unsigned int ini = factors_sort[i].first;
		if( ini == cur_ref ){
			arr_s[cur_pos++] = 1;
		}
		else{
			arr_s[cur_pos++] = 0;
			++cur_ref;
			--i;
		}
	}
	cout << "RelzIndexHashCompacted - Vector S prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	bits_s = bits_s_type(arr_s);
	select1_s = bits_s_type::select_1_type(&bits_s);
	select0_s = bits_s_type::select_0_type(&bits_s);
	
	// Permutacion 
	cout << "RelzIndexHashCompacted - Preparing Permutation PI\n";
	pi = int_vector<>(n_factors);
	pi_inv = int_vector<>(n_factors);
	for( unsigned int i = 0; i < n_factors; ++i ){
		pi[i] = factors_sort[i].second.second;
		pi_inv[ factors_sort[i].second.second ] = i;
	}
	
	cout << "RelzIndexHashCompacted - PI prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	// Posiciones finales Ez
	int_vector<> ez(n_factors);
	for( unsigned int i = 0; i < n_factors; ++i ){
		ez[i] = factors_sort[i].second.first;
	}
	rmq = rmq_type(&ez);
	
	// Bit vector B (inicio de las frases en texto)
	cout << "RelzIndexHashCompacted - Preparing Vector B\n";
	bit_vector arr_b(len_text, 0);
	unsigned int pos_text = 0;
	for( unsigned int i = 0; i < n_factors; ++i ){
		unsigned int len = factors[i].second;
		arr_b[ pos_text ] = 1;
		pos_text += len;
	}
	cout << "RelzIndexHashCompacted - Vector B prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	bits_b = bits_b_type(arr_b);
	select1_b = bits_b_type::select_1_type(&bits_b);
	select0_b = bits_b_type::select_0_type(&bits_b);
	
	// Construccion del FM Index (aunque use el SA original para la compresion, esto es para la busqueda)
	// Construccion con datos en memoria, en un string
	cout << "RelzIndexHashCompacted - Preparing fm_index\n";
	construct_im(fm_index, _ref_text, 1);
	cout << "RelzIndexHashCompacted - fm_index prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	// Preparacion de permutaciones X e Y
	cout << "RelzIndexHashCompacted - Preparing arr X\n";
	vector<unsigned int> arr_x_original(n_factors);
	for( unsigned int i = 0; i < n_factors; ++i ){
		arr_x_original[i] = i;
	}
	FactorsFastIteratorReverseComparator comp_rev(&factors_start, full_text, len_text);
	stable_sort(arr_x_original.begin(), arr_x_original.end(), comp_rev);
	arr_x = int_vector<>(n_factors);
	for( unsigned int i = 0; i < n_factors; ++i ){
		arr_x[i] = arr_x_original[i];
	}
	
	for( unsigned int i = 0; i < n_factors; ++i ){
		cout << " arr_x[" << i << "]: " << arr_x[i] << " -> ";
		char c = 0;
		FactorsIteratorReverse it(arr_x[i] - 1, n_factors, &select1_s, &select1_b, &select0_b, &pi_inv, _ref_text, &fm_index, len_text);
		for(unsigned int k = 0; k < 20 && it.hasNext() && (c = it.next()) != 0; ++k ) 
			cout << c;
		cout << "\n";
	}
	cout << "-----\n";
	
	cout << "RelzIndexHashCompacted - Preparing arr Y\n";
	vector<unsigned int> arr_y_original(n_factors);
	for( unsigned int i = 0; i < n_factors; ++i ){
		arr_y_original[i] = i;
	}
	
	FactorsFastIteratorComparator comp(&factors_start, full_text, len_text);
	stable_sort(arr_y_original.begin(), arr_y_original.end(), comp);
	
	arr_y = int_vector<>(n_factors);
	int_vector<> arr_y_inv(n_factors);
	for( unsigned int i = 0; i < n_factors; ++i ){
		arr_y[i] = arr_y_original[i];
		arr_y_inv[ arr_y_original[i] ] = i;
	}
	
	for( unsigned int i = 0; i < n_factors; ++i ){
		cout << " arr_y[" << i << "]: " << arr_y[i] << " -> ";
		char c = 0;
		FactorsIterator it(arr_y[i], n_factors, &select1_s, &select1_b, &select0_b, &pi_inv, _ref_text, &fm_index, len_text);
		for(unsigned int k = 0; k < 20 && it.hasNext() && (c = it.next()) != 0; ++k ) 
			cout << c;
		cout << " (" << it.length() << ")\n";
	}
	cout << "-----\n";

	cout << "RelzIndexHashCompacted - X & Y prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	cout << "RelzIndexHashCompacted - Preparing WT\n";
	int_vector<> values_wt(n_factors);
	for( unsigned int i = 0; i < n_factors; ++i ){
		values_wt[i] = arr_y_inv[ arr_x[ i ] ];
	}
	
	construct_im(wt, values_wt);
	cout << "RelzIndexHashCompacted - WT prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	// Prueba de aceleracion de recursive_rmq almacenando los datos de los factores descomprimidos
	if(precompute_rmq){
		for( unsigned int i = 0; i < n_factors; ++i ){
			unsigned int tu = select1_s(i + 1) - i;
			unsigned int pu = select1_b(pi[i] + 1);
			unsigned int lu = select1_b(pi[i] + 2) - pu;
			arr_tu.push_back(tu);
			arr_pu.push_back(pu);
			arr_lu.push_back(lu);
		}
	}
	
	cout << "RelzIndexHashCompacted - Testing factors_start\n";
	for(unsigned int i = 0; i < n_factors; ++i){
		if(i < 10 || i > n_factors-10 ){
			cout << "factors_start[" << i << "]: " << factors_start[i] << " (" << select1_b(i+1) << ")\n";
		}
		if( factors_start[i] != select1_b(i+1) ){
			cerr << "Error\n";
			exit(0);
		}
	}
	
	cout << "RelzIndexHashCompacted - Preparing KarpRobin Structures\n";
	
	// En principio se necesita un arreglo de un hash por cada log(n) caracteres de la referencia
	// Ademas necesito la firma de cada sufijo de factor quizas?
	
	BitsUtils bits_utils;
	unsigned int log_n = bits_utils.n_bits(len_ref);
	cout << "RelzIndexHashCompacted - Adding hash for Reference prefixes\n";
	cout << "log_n: " << log_n << " de " << len_ref << "\n";
	arr_kr_ref.push_back( karp_rabin->hash(ref_text, log_n) );
	unsigned int processed_text = log_n;
	while( processed_text < len_ref ){
		unsigned int word_len = log_n;
		if( len_ref - processed_text < word_len ){
			word_len = len_ref - processed_text;
		}
		unsigned long long kr1 = arr_kr_ref.back();
		unsigned long long kr2 = karp_rabin->hash(ref_text + processed_text, word_len);
		unsigned long long kr_total = karp_rabin->concat(kr1, kr2, word_len);
//		unsigned long long kr_test = karp_rabin->hash(ref_text, processed_text + word_len);
//		cout << "Agregando " << kr_total << " / " << kr_test << " (kr1: " << kr1 << ", " << kr2 << ", word_len: " << word_len << ")\n";
		arr_kr_ref.push_back(kr_total);
		processed_text += word_len;
	}
	
	// factors_start almacena la posicion de inicio de cada frase (=> los caracteres ANTERIORES forman el prefijo)
	// Por eso agrego un 0 para la primera frase (prefijo nulo)
	cout << "RelzIndexHashCompacted - Adding hash for Frases prefixes\n";
	arr_kr_s.push_back(0);
	for(unsigned int i = 1; i < factors_start.size(); ++i){
		unsigned int word_len = factors_start[i] - factors_start[i-1];
		unsigned long long kr1 = arr_kr_s.back();
		unsigned long long kr2 = karp_rabin->hash(full_text + factors_start[i-1], word_len);
		unsigned long long kr_total = karp_rabin->concat(kr1, kr2, word_len);
		
//		string s(full_text, factors_start[i]);
//		cout << "RelzIndexHashCompacted - Adding hash de \"" << s << "\" " << kr_total << " / " << karp_rabin->hash(s) << "\n";
		
		arr_kr_s.push_back(kr_total);
	}
	// Agrego tambien el hash de la coleccion completa como un factor ficticio con id = n_factors
	// Podria necesitarlo despues para simplificar los calculos
	unsigned int word_len = len_text - factors_start.back();
	unsigned long long kr1 = arr_kr_s.back();
	unsigned long long kr2 = karp_rabin->hash(full_text + factors_start.back(), word_len);
	unsigned long long kr_total = karp_rabin->concat(kr1, kr2, word_len);
	arr_kr_s.push_back(kr_total);
	
	cout << "RelzIndexHashCompacted - KarpRobin prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	// Original
//	kr_factors = new KarpRabinFactorsSuffixesv1(n_factors, &arr_kr_s, karp_rabin, ref_text, &pi_inv, &arr_tu, &arr_pu, &arr_lu, &factors_start);
	
	// New
	kr_factors = new KarpRabinFactorsSuffixesv2(n_factors, &arr_kr_s, karp_rabin, ref_text, &select1_s, &select1_b, &select0_b, &pi_inv, &factors_start);
	
	// Para esta fase, en CONSTRUCCION usare datos descomprimidos para simplificarlo
	// Obviamente esto es olo para construccion y los datos usados no se almacenan, solo los datos de los nodos
	// Tambien, si hay un archivo para tree_y almacenado, lo cargo en lugar de construirlo
	
	cout << "RelzIndexHashCompacted - Building Tree Y\n";
	timer.reset();
//	tree_y.build(full_text, len_text, factors_start, arr_y_original, karp_rabin, kr_factors);
	tree_y.build(full_text, len_text, factors_start, &arr_y, karp_rabin, kr_factors);
	cout << "RelzIndexHashCompacted - Tree Y finished in (" << timer.getMilisec() << " ms)\n";
	tree_y.print();
	
	cout << "RelzIndexHashCompacted - Building Tree X\n";
	timer.reset();
//	tree_x.build(full_text, len_text, factors_start, arr_x_original, karp_rabin, kr_factors);
	tree_x.build(full_text, len_text, factors_start, &arr_x, karp_rabin, kr_factors);
	cout << "RelzIndexHashCompacted - Tree X finished in (" << timer.getMilisec() << " ms)\n";
	tree_x.print();
	
	cout << "RelzIndexHashCompacted - Trees prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	// Compactacion de arreglos descomprimidos
	sdsl::util::bit_compress(pi);
	sdsl::util::bit_compress(pi_inv);
	sdsl::util::bit_compress(arr_x);
	sdsl::util::bit_compress(arr_y);
	
	cout << "RelzIndexHashCompacted - End\n";
	
}

RelzIndexHashCompacted::~RelzIndexHashCompacted(){
	delete [] ref_text;
	ref_text = NULL;
}

void RelzIndexHashCompacted::printSize(){
	double total_bytes = 0;
	
	// texto descomprimido
//	unsigned int len_ref = ref_text->length();
	total_bytes += len_ref;
	cout << "RelzIndexReference::printSize - Reference Text: " << (8.0*len_ref/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(fm_index);
	cout << "RelzIndexHashCompacted::printSize - fm_index: " << (8.0*size_in_bytes(fm_index)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(rmq);
	cout << "RelzIndexHashCompacted::printSize - rmq: " << (8.0*size_in_bytes(rmq)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(bits_s);
	cout << "RelzIndexHashCompacted::printSize - bits_s: " << (8.0*size_in_bytes(bits_s)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(bits_b);
	cout << "RelzIndexHashCompacted::printSize - bits_b: " << (8.0*size_in_bytes(bits_b)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(pi);
	cout << "RelzIndexHashCompacted::printSize - pi: " << (8.0*size_in_bytes(pi)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(pi_inv);
	cout << "RelzIndexHashCompacted::printSize - pi_inv: " << (8.0*size_in_bytes(pi_inv)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(arr_x);
	cout << "RelzIndexHashCompacted::printSize - arr_x: " << (8.0*size_in_bytes(arr_x)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(arr_y);
	cout << "RelzIndexHashCompacted::printSize - arr_y: " << (8.0*size_in_bytes(arr_y)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(wt);
	cout << "RelzIndexHashCompacted::printSize - wt: " << (8.0*size_in_bytes(wt)/len_text) << " bps\n";
	
//	tree_y
	unsigned int max_len = 0;
	unsigned int max_childs = 0;
	unsigned int max_height = 0;
	unsigned int total_childs = 0;
	total_childs = tree_y.root.totalChilds(max_len, max_childs, max_height, 0);
	total_bytes += ((double)total_childs)*12.625;
//	cout << "RelzIndexHashCompacted::printSize - tree_y: " << (((double)total_childs)*12.625/(1024*1024)) << " MB\n";
	cout << "RelzIndexHashCompacted::printSize - tree_y: " << (8.0*((double)total_childs)*12.625/(len_text)) << " bps (total_childs: " << total_childs << ", max_len: " << max_len << ", max_childs: " << max_childs << ", max_height: " << max_height << ")\n";
	
//	tree_x
	max_len = 0;
	max_childs = 0;
	max_height = 0;
	total_childs = 0;
	total_childs = tree_x.root.totalChilds(max_len, max_childs, max_height, 0);
	total_bytes += ((double)total_childs)*12.625;
//	cout << "RelzIndexHashCompacted::printSize - tree_x: " << (((double)total_childs)*12.625/(1024*1024)) << " MB\n";
	cout << "RelzIndexHashCompacted::printSize - tree_x: " << (8.0*((double)total_childs)*12.625/(len_text)) << " bps (total_childs: " << total_childs << ", max_len: " << max_len << ", max_childs: " << max_childs << ", max_height: " << max_height << ")\n";
	
	cout << "RelzIndexHashCompacted::printSize - Total " << total_bytes << " (" << (total_bytes/(1024*1024)) << " MB)\n";
	
}

void RelzIndexHashCompacted::findTimes(const string &pattern, vector<unsigned int> &results){
	
	cout << "RelzIndexHashCompacted::findTimes - Start (\"" << pattern << "\")\n";
	NanoTimer timer;
	
	cout << "RelzIndexHashCompacted::findTimes - Section A, reference\n";
	
	size_t m = pattern.size();
	size_t occs = sdsl::count(fm_index, pattern.begin(), pattern.end());
	occs_a += occs;
	vector<int_vector<64>> arr_locations;
	if( occs > 0 ){
		arr_locations.push_back(locate(fm_index, pattern.begin(), pattern.begin()+m));
		sort(arr_locations.back().begin(), arr_locations.back().end());
	}
	querytime_p1 += timer.getNanosec();
	timer.reset();
	
	for( int_vector<64> locations : arr_locations ){
		for( unsigned int i = 0; i < occs; ++i ){
			unsigned int occ_i = locations[i];
			// Comprobar los factores que cuben esta ocurrencia (el string ref[occ_i, occ_i + m - 1])
			unsigned int select = select0_s(occ_i + 1);
			unsigned int pos_ez = select - 1 - occ_i;
			// Now the recursive search in rmq (0 - pos_ez)
			if( occ_i >= select ){
				continue;
			}
			recursive_rmq(0, pos_ez, (occ_i + m), occ_i, results);
		}
	}
	querytime_p2 += timer.getNanosec();
	
	cout << "RelzIndexHashCompacted::findTimes - Section B, ranges\n";
	
	vector<unsigned long long> kr_pat_vector;
	vector<unsigned long long> kr_pat_rev_vector;
	karp_rabin->hashPrefixes(pattern, kr_pat_vector);
	karp_rabin->hashPrefixesRev(pattern, kr_pat_rev_vector);
	
	string pattern_rev;
	for(unsigned int i = 0; i < pattern.length(); ++i){
		pattern_rev += pattern[pattern.length() - 1 - i];
	}
//	cout << "-----  pattern: " << pattern << " -----\n";
//	cout << "-----  pattern_rev: " << pattern_rev << " -----\n";
	
	for(unsigned int i = 1; i < pattern.length(); ++i){
		timer.reset();
		
		cout << "-----  tree_x.getRange -----\n";
		pair<unsigned int, unsigned int> r1 = tree_x.getRange(kr_pat_rev_vector, i, pattern_rev);
		querytime_p3x += timer.getNanosec();
		timer.reset();
		cout << "-----\n";
		
		if( r1.first == (unsigned int)(-1) || r1.second == (unsigned int)(-1) || r1.second < r1.first ){
			continue;
		}
		
		cout << "-----  tree_y.getRange -----\n";
		pair<unsigned int, unsigned int> r2 = tree_y.getRange(kr_pat_vector, i, pattern);
		querytime_p3y += timer.getNanosec();
		timer.reset();
		cout << "-----\n";
		
		if( r2.first == (unsigned int)(-1) || r2.second == (unsigned int)(-1) || r2.second < r2.first ){
			continue;
		}
		
		cout << "RelzIndexHashCompacted::findTimes - Searching in [" << r1.first << ", " << r1.second << "] x [" << r2.first << ", " << r2.second << "]:\n";
		auto res = wt.range_search_2d(r1.first, r1.second, r2.first, r2.second);
		cout << "RelzIndexHashCompacted::findTimes - Adding " << res.second.size() << " points\n";
		for (auto point : res.second){
			unsigned int f = arr_y[point.second];
			unsigned int pu = select1_b(f + 1);
			
			// Verificacion
			bool omit = false;
			
			FactorsIteratorReverse it_x(f-1, n_factors, &select1_s, &select1_b, &select0_b, &pi_inv, ref_text, &fm_index, len_text);
//			cout << "text_x: ";
			for(unsigned int pos = 0; pos < i; ++pos){
				char c = it_x.next();
//				cout << c << "|" << pattern[ i - 1 - pos ] << " ";
				if( c != pattern[ i - 1 - pos ] ){
					omit = true;
					break;
				}
			}
//			cout << "\n";
			
			FactorsIterator it_y(f, n_factors, &select1_s, &select1_b, &select0_b, &pi_inv, ref_text, &fm_index, len_text);
//			cout << "text_y: ";
			for(unsigned int pos = 0; pos < pattern.length()-i; ++pos){
				char c = it_y.next();
//				cout << c << "|" << pattern[ i + pos ] << " ";
				if( c != pattern[ i + pos ] ){
					omit = true;
					break;
				}
			}
//			cout << "\n";
			
			if( omit ){
//				cout << "RelzIndexHashCompacted::findTimes - Omiting bad result.\n";
				continue;
			}
			
			results.push_back(pu - i);
			++occs_c;
		}
		querytime_p4 += timer.getNanosec();
		
	}
	
	cout << "RelzIndexHashCompacted::findTimes - End\n";
	
}

void RelzIndexHashCompacted::recursive_rmq(unsigned int ini, unsigned int fin, unsigned int min_pos, unsigned int occ_ref, vector<unsigned int> &results){
//	cout << "RelzIndexHashCompacted::recursive_rmq - " << ini << ", " << fin << "\n";
	
	unsigned int pos_max = rmq(ini, fin);
	
//	cout << "FactorsIndex::recursive_rmq - Computing factor\n";
	assert(pos_max < n_factors);
	unsigned int tu = 0;
	unsigned int pu = 0;
	unsigned int lu = 0;
	if( precompute_rmq ){
		tu = arr_tu[pos_max];
		pu = arr_pu[pos_max];
		lu = arr_lu[pos_max];
	}
	else{
		tu = select1_s(pos_max + 1) - pos_max;
		pu = select1_b(pi[pos_max] + 1);
		lu = select1_b(pi[pos_max] + 2) - pu;
	}
	
//	cout << "RelzIndexHashCompacted::recursive_rmq - tu: " << tu << ", pu: " << pu << ", lu: " << lu << "\n";
	
	if( tu + lu < min_pos ){
		return;
	}
	else{
//		cout << " -> Adding " << (pu + (occ_ref - tu)) << "\n";
		results.push_back(pu + (occ_ref - tu));
		++occs_b;
	}
	
	if( (pos_max > 0) && (ini < pos_max) ){
		recursive_rmq(ini, pos_max-1, min_pos, occ_ref, results);
	}
	if( pos_max < fin ){
		recursive_rmq(pos_max+1, fin, min_pos, occ_ref, results);
	}
}

void RelzIndexHashCompacted::save(const string &file_base){
	cout << "RelzIndexHashCompacted::save - Start (base \"" << file_base << "\")\n";
	
	// Base
	string index_basic_file = file_base + ".base";
	fstream writer(index_basic_file, fstream::out | fstream::trunc);
	// Version of the index
	unsigned char version = 3;
	writer.write((char*)&version, 1);
	// len_text
	writer.write((char*)&len_text, sizeof(int));
	// n_factors
	writer.write((char*)&n_factors, sizeof(int));
	// len_ref
	writer.write((char*)&len_ref, sizeof(int));
	// Reference Text
	writer.write((char*)ref_text, len_ref);
	// Close Base
	writer.close();
	
	// fm_index
	string fm_index_file = file_base + ".fm";
	store_to_file(fm_index, fm_index_file);
	
	// rmq
	string rmq_file = file_base + ".rmq";
	store_to_file(rmq, rmq_file);
	
	// bits_s
	string bits_s_file = file_base + ".arrs";
	store_to_file(bits_s, bits_s_file);
	
	// bits_b
	string bits_b_file = file_base + ".arrb";
	store_to_file(bits_b, bits_b_file);
	
	// wt
	string wt_file = file_base + ".wt";
	store_to_file(wt, wt_file);
	
	// pi
	string pi_file = file_base + ".pi";
	store_to_file(pi, pi_file);
	
	// pi_inv
	string pi_inv_file = file_base + ".pi_inv";
	store_to_file(pi_inv, pi_inv_file);
	
	// arr_x
	string x_file = file_base + ".x";
	store_to_file(arr_x, x_file);
	
	// arr_y
	string y_file = file_base + ".y";
	store_to_file(arr_y, y_file);
	
	
	// KarpRabinFactorsSuffixes
	string krs_file = file_base + ".krsuffixes";
	kr_factors->save(krs_file);
	
	// tree_x
	string tree_x_file = file_base + ".tree.x";
	tree_x.save(tree_x_file);
	
	// tree_y
	string tree_y_file = file_base + ".tree.y";
	tree_y.save(tree_y_file);
	
	
	
	cout << "RelzIndexHashCompacted::save - End\n";
}

void RelzIndexHashCompacted::load(const string &file_base, KarpRabin *_karp_rabin){
	cout << "RelzIndexHashCompacted::load - Start (base \"" << file_base << "\")\n";
	
	// Base
	string index_basic_file = file_base + ".base";
	fstream reader(index_basic_file, fstream::in);
	// Version of the index
	unsigned char version = 0;
	reader.read((char*)&version, 1);
	if( version != 3 ){
		cout << "RelzIndexHashCompacted::load - Wrong Version\n";
		return;
	}
	// len_text
	reader.read((char*)&len_text, sizeof(int));
	// n_factors
	reader.read((char*)&n_factors, sizeof(int));
	// len_ref
	reader.read((char*)&len_ref, sizeof(int));
	// Reference Text
	if( ref_text != NULL ){
		delete [] ref_text;
		ref_text = NULL;
	}
	ref_text = new char[len_ref + 1];
	reader.read((char*)ref_text, len_ref);
	ref_text[len_ref] = 0;
	// Close Base
	reader.close();
	
	// fm_index
	cout << "RelzIndexHashCompacted::load - fm_index\n";
	string fm_index_file = file_base + ".fm";
	load_from_file(fm_index, fm_index_file);
	
	// rmq
	cout << "RelzIndexHashCompacted::load - rmq\n";
	string rmq_file = file_base + ".rmq";
	load_from_file(rmq, rmq_file);
	
	// bits_s
	cout << "RelzIndexHashCompacted::load - bits_s\n";
	string bits_s_file = file_base + ".arrs";
	load_from_file(bits_s, bits_s_file);
	select1_s = bits_s_type::select_1_type(&bits_s);
	select0_s = bits_s_type::select_0_type(&bits_s);
	
	// bits_b
	cout << "RelzIndexHashCompacted::load - bits_b\n";
	string bits_b_file = file_base + ".arrb";
	load_from_file(bits_b, bits_b_file);
	select1_b = bits_b_type::select_1_type(&bits_b);
	select0_b = bits_b_type::select_0_type(&bits_b);
	
	// pi
	string pi_file = file_base + ".pi";
	load_from_file(pi, pi_file);
	
	// pi_inv
	string pi_inv_file = file_base + ".pi_inv";
	load_from_file(pi_inv, pi_inv_file);
	
	// arr_x
	string x_file = file_base + ".x";
	load_from_file(arr_x, x_file);
	
	// arr_y
	string y_file = file_base + ".y";
	load_from_file(arr_y, y_file);
	
	// wt
	cout << "RelzIndexHashCompacted::load - wt\n";
	string wt_file = file_base + ".wt";
	load_from_file(wt, wt_file);
	
	if(precompute_rmq){
		for( unsigned int i = 0; i < n_factors; ++i ){
			unsigned int tu = select1_s(i + 1) - i;
			unsigned int pu = select1_b(pi[i] + 1);
			unsigned int lu = select1_b(pi[i] + 2) - pu;
			arr_tu.push_back(tu);
			arr_pu.push_back(pu);
			arr_lu.push_back(lu);
		}
	}
	
	// KarpRabin (from parameter)
	karp_rabin = _karp_rabin;
	
	// KarpRabinFactorsSuffixes																				
	string krs_file = file_base + ".krsuffixes";
	// Instead of a load, we use a builder here to load, to pass the more specific parameters
	// kr_factors->load(krs_file);
	kr_factors = new KarpRabinFactorsSuffixesv2(krs_file, karp_rabin, ref_text, &select1_s, &select1_b, &select0_b, &pi_inv);
	
	// tree_x
	string tree_x_file = file_base + ".tree.x";
	tree_x.load(karp_rabin, kr_factors, &arr_x, tree_x_file);
	
	// tree_y
	string tree_y_file = file_base + ".tree.y";
	tree_y.load(karp_rabin, kr_factors, &arr_y, tree_y_file);
	
	
	
	cout << "RelzIndexHashCompacted::load - End\n";
	
}


