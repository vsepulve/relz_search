#include "RelzIndexHash.h"

RelzIndexHash::RelzIndexHash(){
	len_text = 0;
	ref_text = NULL;
	len_ref = 0;
	n_factors = 0;
	karp_rabin = NULL;
}

//RelzIndexHash::RelzIndexHash(vector<pair<unsigned int, unsigned int> > &factors, char *full_text, unsigned int _len_text, const char *_ref_text, unsigned int _len_ref, KarpRabin *_karp_rabin, const char *kr_frases_file, bool load_kr_frases){

RelzIndexHash::RelzIndexHash(vector<pair<unsigned int, unsigned int> > &factors, char *full_text, unsigned int _len_text, const char *_ref_text, unsigned int _len_ref, KarpRabin *_karp_rabin, const char *index_base_file){
	
	len_text = _len_text;
	ref_text = _ref_text;
	len_ref = _len_ref;
	n_factors = factors.size();
	karp_rabin = _karp_rabin;
	
	NanoTimer timer;
	
	cout << "RelzIndexHash - Inicio (factors: " << factors.size() << ", len_text: " << len_text << ", len_ref: " << len_ref << ")\n";
//	cout << "RelzIndexHash - Ref: " << ref_text << "\n";
//	cout << "RelzIndexHash - Text: " << full_text << "\n";
	
	cout << "RelzIndexHash - Preparing Factors\n";
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
	cout << "RelzIndexHash - Factors Sorted prepared in " << timer.getMilisec() << "\n";
	timer.reset();
//	cout << "Factors Sorted: \n";
//	for( unsigned int i = 0; i < factors_sort.size(); ++i ){
//		cout << "factors_sort[" << i << "]: (" << factors_sort[i].first << ", " << factors_sort[i].second.first << ", " << factors_sort[i].second.second << ")\n";
//	}
	
	// Bit vector S
	cout << "RelzIndexHash - Preparing Vector S\n";
	bit_vector _arr_s(len_ref + n_factors, 0);
	arr_s = _arr_s;
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
	cout << "RelzIndexHash - Vector S prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	bits_s = bits_s_type(arr_s);
	select1_s = bits_s_type::select_1_type(&bits_s);
	select0_s = bits_s_type::select_0_type(&bits_s);
	
	// Permutacion 
	cout << "RelzIndexHash - Preparing Permutation PI\n";
	pi = int_vector<>(n_factors);
	pi_inv = int_vector<>(n_factors);
	for( unsigned int i = 0; i < n_factors; ++i ){
		pi[i] = factors_sort[i].second.second;
		pi_inv[ factors_sort[i].second.second ] = i;
	}
	
	cout << "RelzIndexHash - PI prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	// Posiciones finales Ez
	int_vector<> _ez(n_factors);
	ez = _ez;
	for( unsigned int i = 0; i < n_factors; ++i ){
		ez[i] = factors_sort[i].second.first;
	}
	rmq = rmq_type(&ez);
	
	// Bit vector B (inicio de las frases en texto)
	cout << "RelzIndexHash - Preparing Vector B\n";
	bit_vector _arr_b(len_text, 0);
	arr_b = _arr_b;
	unsigned int pos_text = 0;
	for( unsigned int i = 0; i < n_factors; ++i ){
		unsigned int len = factors[i].second;
		arr_b[ pos_text ] = 1;
		pos_text += len;
	}
	cout << "RelzIndexHash - Vector B prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	bits_b = bits_b_type(arr_b);
	select1_b = bits_b_type::select_1_type(&bits_b);
	select0_b = bits_b_type::select_0_type(&bits_b);
	
	// Construccion del FM Index (aunque use el SA original para la compresion, esto es para la busqueda)
	// Construccion con datos en memoria, en un string
	cout << "RelzIndexHash - Preparing fm_index\n";
	construct_im(fm_index, ref_text, 1);
	cout << "RelzIndexHash - fm_index prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	// Preparacion de permutaciones X e Y
	cout << "RelzIndexHash - Preparing arr X\n";
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
		for(unsigned int k = 0; k < 20 && (c = getCharRev(arr_x[i] - 1, k)) != 0; ++k ) 
			cout << c;
		cout << "\n";
	}
	cout << "-----\n";
	
	cout << "RelzIndexHash - Preparing arr Y\n";
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
		for(unsigned int k = 0; k < 20 && (c = getChar(arr_y[i], k)) != 0; ++k ) 
			cout << c;
		cout << " (" << (len_text - factors_start[arr_y[i]]) << ")\n";
	}
	cout << "-----\n";

	cout << "RelzIndexHash - X & Y prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	cout << "RelzIndexHash - Preparing WT\n";
	int_vector<> _values_wt(n_factors);
	values_wt = _values_wt;
	for( unsigned int i = 0; i < n_factors; ++i ){
		values_wt[i] = arr_y_inv[ arr_x[ i ] ];
	}
	
	construct_im(wt, values_wt);
	cout << "RelzIndexHash - WT prepared in " << timer.getMilisec() << "\n";
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
	
	cout << "RelzIndexHash - Preparing KarpRobin Structures\n";
	
	// En principio se necesita un arreglo de un hash por cada log(n) caracteres de la referencia
	// Ademas necesito la firma de cada sufijo de factor quizas?
	
	BitsUtils bits_utils;
	unsigned int log_n = bits_utils.n_bits(len_ref);
	cout << "RelzIndexHash - Adding hash for Reference prefixes\n";
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
//		if( arr_kr_ref.size() >= 10 ){
//			break;
//		}
	}
	
	// Notar que el arreglo almacena los has de los PREFIJOS de cada frase
	// factors_start almacena la posicion de inicio de cada frase (=> los caracteres ANTERIORES forman el prefijo)
	// Por eso agrego un 0 para la primera frase (prefijo nulo)
	cout << "RelzIndexHash - Adding hash for Frases prefixes\n";
	arr_kr_s.push_back(0);
	for(unsigned int i = 1; i < factors_start.size(); ++i){
		unsigned int word_len = factors_start[i] - factors_start[i-1];
		unsigned long long kr1 = arr_kr_s.back();
		unsigned long long kr2 = karp_rabin->hash(full_text + factors_start[i-1], word_len);
		unsigned long long kr_total = karp_rabin->concat(kr1, kr2, word_len);
		
//		string s(full_text, factors_start[i]);
//		cout << "RelzIndexHash - Adding hash de \"" << s << "\" " << kr_total << " / " << karp_rabin->hash(s) << "\n";
		
		arr_kr_s.push_back(kr_total);
	}
	// Agrego tambien el hash de la coleccion completa como un factor ficticio con id = n_factors
	// Podria necesitarlo despues para simplificar los calculos
	unsigned int word_len = len_text - factors_start.back();
	unsigned long long kr1 = arr_kr_s.back();
	unsigned long long kr2 = karp_rabin->hash(full_text + factors_start.back(), word_len);
	unsigned long long kr_total = karp_rabin->concat(kr1, kr2, word_len);
	arr_kr_s.push_back(kr_total);
	
	cout << "RelzIndexHash - KarpRobin prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	kr_factors = new KarpRabinFactorsSuffixes(n_factors, &arr_kr_s, karp_rabin, ref_text, &pi_inv, &arr_tu, &arr_pu, &arr_lu, &factors_start);
	
	kr_factors->hash(1, 3, 9);
	cout << "-----\n";
	kr_factors->hashFast(1, 3, 9);
	cout << "-----\n";
	const char *test_str = "ALASLALAB";
	cout << "RelzIndexHash - Testing KarpRabinFactorsSuffixes, hash: " << karp_rabin->hash(test_str, strlen(test_str)) << "\n";
	
	
	string index_y(index_base_file, strlen(index_base_file));
	index_y += ".index.y";
	string index_x(index_base_file, strlen(index_base_file));
	index_x += ".index.x";
	cout << "RelzIndexHash - Preparing Trees (files " << index_y << " and " << index_x << ")\n";
	
	timer.reset();
	// Para esta fase, en CONSTRUCCION usare datos descomprimidos para simplificarlo
	// Obviamente esto es olo para construccion y los datos usados no se almacenan, solo los datos de los nodos
	// Tambien, si hay un archivo para tree_y almacenado, lo cargo en lugar de construirlo
	bool load_y = false;
	fstream reader(index_y, fstream::in);
	if( reader.good() ){
		reader.seekg(0, reader.end);
		unsigned long long file_size = reader.tellg();
		reader.close();
		if( file_size > 4 ){
			load_y = true;
		}
	}
	load_y = false;
	if( load_y ){
		cout << "RelzIndexHash - Loading Tree Y\n";
//		tree_y.load(karp_rabin, kr_factors, index_y);
		tree_y.load(karp_rabin, kr_factors, index_y, factors_start, full_text);
	}
	else{
		cout << "RelzIndexHash - Building Tree Y\n";
		tree_y.build(full_text, len_text, factors_start, arr_y_original, karp_rabin, kr_factors);
		tree_y.save(index_y);
	}
	cout << "RelzIndexHash - Tree Y finished in (" << timer.getMilisec() << " ms)\n";
	tree_y.print();
//	tree_y.prepareChilds();
//	tree_y.printSize();
	
	cout << "RelzIndexHash - Building Tree X\n";
	timer.reset();
	tree_x.build(full_text, len_text, factors_start, arr_x_original, karp_rabin, kr_factors);
	cout << "RelzIndexHash - Tree X finished in (" << timer.getMilisec() << " ms)\n";
	tree_x.print();
//	tree_x.save(index_x);
//	tree_x.prepareChilds();
//	tree_x.printSize();
	
	vector<unsigned long long> pat_vector;
	karp_rabin->hashPrefixes("ALABARDA", pat_vector);
	
	cout << "RelzIndexHash - Trees prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	cout << "RelzIndexHash - End\n";
	
}

RelzIndexHash::~RelzIndexHash(){
	
}

void RelzIndexHash::printSize(){
	double total_bytes = 0;
	
	// texto descomprimido
//	if( ! omit_text ){
//		total_bytes += len_ref;
//		cout << "RelzIndexHash::printSize - Reference Text: " << (len_ref/(1024*1024)) << " MB\n";
//	}
	
//	if( precompute_rmq ){
//		// 3 integers => 12 bytes per factor
//		total_bytes += n_factors * 12;
//		cout << "RelzIndexHash::printSize - Factors: " << (n_factors*12/(1024*1024)) << " MB\n";
//	}
	
	total_bytes += size_in_bytes(fm_index);
	cout << "RelzIndexHash::printSize - fm_index: " << (8.0*size_in_bytes(fm_index)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(rmq);
	cout << "RelzIndexHash::printSize - rmq: " << (8.0*size_in_bytes(rmq)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(bits_s);
	cout << "RelzIndexHash::printSize - bits_s: " << (8.0*size_in_bytes(bits_s)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(bits_b);
	cout << "RelzIndexHash::printSize - bits_b: " << (8.0*size_in_bytes(bits_b)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(pi);
	cout << "RelzIndexHash::printSize - pi: " << (8.0*size_in_bytes(pi)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(pi_inv);
	cout << "RelzIndexHash::printSize - pi_inv: " << (8.0*size_in_bytes(pi_inv)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(arr_x);
	cout << "RelzIndexHash::printSize - arr_x: " << (8.0*size_in_bytes(arr_x)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(arr_y);
	cout << "RelzIndexHash::printSize - arr_y: " << (8.0*size_in_bytes(arr_y)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(wt);
	cout << "RelzIndexHash::printSize - wt: " << (8.0*size_in_bytes(wt)/len_text) << " bps\n";
	
//	tree_y
	unsigned int max_len = 0;
	unsigned int max_childs = 0;
	unsigned int max_height = 0;
	unsigned int total_childs = 0;
	total_childs = tree_y.root.totalChilds(max_len, max_childs, max_height, 0);
	total_bytes += ((double)total_childs)*12.625;
//	cout << "RelzIndexHash::printSize - tree_y: " << (((double)total_childs)*12.625/(1024*1024)) << " MB\n";
	cout << "RelzIndexHash::printSize - tree_y: " << (8.0*((double)total_childs)*12.625/(len_text)) << " bps (total_childs: " << total_childs << ", max_len: " << max_len << ", max_childs: " << max_childs << ", max_height: " << max_height << ")\n";
	
//	tree_x
	max_len = 0;
	max_childs = 0;
	max_height = 0;
	total_childs = 0;
	total_childs = tree_x.root.totalChilds(max_len, max_childs, max_height, 0);
	total_bytes += ((double)total_childs)*12.625;
//	cout << "RelzIndexHash::printSize - tree_x: " << (((double)total_childs)*12.625/(1024*1024)) << " MB\n";
	cout << "RelzIndexHash::printSize - tree_x: " << (8.0*((double)total_childs)*12.625/(len_text)) << " bps (total_childs: " << total_childs << ", max_len: " << max_len << ", max_childs: " << max_childs << ", max_height: " << max_height << ")\n";
	
	cout << "RelzIndexHash::printSize - Total " << total_bytes << " (" << (total_bytes/(1024*1024)) << " MB)\n";
	
}

void RelzIndexHash::findTimes(const string &pattern, vector<unsigned int> &results){
	
	cout << "RelzIndexHash::findTimes - Start (\"" << pattern << "\")\n";
	NanoTimer timer;
	
	cout << "RelzIndexHash::findTimes - Section A, reference\n";
	
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
	
	cout << "RelzIndexHash::findTimes - Section B, ranges\n";
	
	vector<unsigned long long> kr_pat_vector;
	vector<unsigned long long> kr_pat_rev_vector;
	karp_rabin->hashPrefixes(pattern, kr_pat_vector);
	karp_rabin->hashPrefixesRev(pattern, kr_pat_rev_vector);
	
	string pattern_rev;
	for(unsigned int i = 0; i < pattern.length(); ++i){
		pattern_rev += pattern[pattern.length() - 1 - i];
	}
	cout << "-----  pattern: " << pattern << " -----\n";
	cout << "-----  pattern_rev: " << pattern_rev << " -----\n";
	
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
		
		cout << "RelzIndexHash::findTimes - Searching in [" << r1.first << ", " << r1.second << "] x [" << r2.first << ", " << r2.second << "]:\n";
		auto res = wt.range_search_2d(r1.first, r1.second, r2.first, r2.second);
		cout << "RelzIndexHash::findTimes - Adding " << res.second.size() << " points\n";
		for (auto point : res.second){
			unsigned int f = arr_y[point.second];
			unsigned int pu = select1_b(f + 1);
			
			// Verificacion
			bool omit = false;
			
			FactorsIteratorReverse it_x(f-1, n_factors, &select1_s, &select1_b, &select0_b, &pi, &pi_inv, ref_text, &fm_index, len_text);
			cout << "text_x: ";
			for(unsigned int pos = 0; pos < i; ++pos){
				char c = it_x.next();
				cout << c << "|" << pattern[ i - 1 - pos ] << " ";
				if( c != pattern[ i - 1 - pos ] ){
					omit = true;
					break;
				}
			}
			cout << "\n";
			
			FactorsIterator it_y(f, n_factors, &select1_s, &select1_b, &select0_b, &pi, &pi_inv, ref_text, &fm_index, len_text);
			cout << "text_y: ";
			for(unsigned int pos = 0; pos < pattern.length()-i; ++pos){
				char c = it_y.next();
				cout << c << "|" << pattern[ i + pos ] << " ";
				if( c != pattern[ i + pos ] ){
					omit = true;
					break;
				}
			}
			cout << "\n";
			
			if( omit ){
				cout << "RelzIndexHash::findTimes - Omiting bad result.\n";
				continue;
			}
			
			results.push_back(pu - i);
			++occs_c;
		}
		querytime_p4 += timer.getNanosec();
		
	}
	
	cout << "RelzIndexHash::findTimes - End\n";
	
}

/*
void RelzIndexHash::find(const string &pattern, vector<unsigned int> &results){

//	cout << "RelzIndexHash::find - Start\n";
	
//	cout << "RelzIndexHash::find - Section A, reference\n";
	
	size_t m = pattern.size();
	size_t occs = sdsl::count(fm_index, pattern.begin(), pattern.end());
	if( occs > 0 ){
		auto locations = locate(fm_index, pattern.begin(), pattern.begin()+m);
		sort(locations.begin(), locations.end());
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
	
//	cout << "RelzIndexHash::find - Section B, ranges\n";
	for(unsigned int i = 1; i < pattern.length(); ++i){
		string p1 = pattern.substr(0, i);
		string p1_rev = "";
		for( unsigned int k = 0; k < p1.length(); ++k ){
			p1_rev += p1[ p1.length() - 1 - k ];
		}
		string p2 = pattern.substr(i, pattern.length() - i);
		pair<unsigned int, unsigned int> r1 = getRangeX(p1_rev.c_str());
		pair<unsigned int, unsigned int> r2 = getRangeY(p2.c_str());
		
		if( r1.second == (unsigned int)(-1) || r1.second < r1.first
			|| r2.second == (unsigned int)(-1) || r2.second < r2.first ){
//			cout << "RelzIndexHash::find - Invalid ranges, omitting...\n";
			continue;
		}
		
//		cout << "RelzIndexHash::find - Searching in [" << r1.first << ", " << r1.second << "] x [" << r2.first << ", " << r2.second << "]:\n";
		auto res = wt.range_search_2d(r1.first, r1.second, r2.first, r2.second);
		for (auto point : res.second){
			unsigned int f = arr_y[point.second];
			unsigned int pu = select1_b(f + 1);
//			cout << " -> Adding " << (pu - p1.length()) << "\n";
			results.push_back(pu - p1.length());
		}
		
	}
//	cout << "RelzIndexHash::find - End\n";
	
}
*/

void RelzIndexHash::recursive_rmq(unsigned int ini, unsigned int fin, unsigned int min_pos, unsigned int occ_ref, vector<unsigned int> &results){
//	cout << "RelzIndexHash::recursive_rmq - " << ini << ", " << fin << "\n";
	
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
	
//	cout << "RelzIndexHash::recursive_rmq - tu: " << tu << ", pu: " << pu << ", lu: " << lu << "\n";
	
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

char RelzIndexHash::getChar(unsigned int factor, unsigned int pos){
	
	// Iterators cache
	if( mapa_iterators.find(factor) == mapa_iterators.end() ){
		mapa_iterators[factor] = FactorsIterator(factor, n_factors, &select1_s, &select1_b, &select0_b, &pi, &pi_inv, ref_text, &fm_index, len_text);
	}
	FactorsIterator it = mapa_iterators[factor];
	if( pos >= it.length() ){
		return 0;
	}
	if( it.position() > pos ){
		it.reset();
	}
	
	char c = 0;
	while( it.position() <= pos ){
		c = it.next();
	}
	return c;
}

char RelzIndexHash::getCharRev(unsigned int factor, unsigned int pos){
	
	if( factor == (unsigned int)(-1) ){
		return 0;
	}
	
	// Iterators cache
	if( mapa_iterators_rev.find(factor) == mapa_iterators_rev.end() ){
		mapa_iterators_rev[factor] = FactorsIteratorReverse(factor, n_factors, &select1_s, &select1_b, &select0_b, &pi, &pi_inv, ref_text, &fm_index, len_text);
	}
	FactorsIteratorReverse it = mapa_iterators_rev[factor];
	if( pos >= it.length() ){
		return 0;
	}
	if( it.position() > pos ){
		it.reset();
	}
	char c = 0;
	while( it.position() <= pos ){
		c = it.next();
	}
	return c;
}


