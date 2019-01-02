#include "RelzIndex.h"

RelzIndex::RelzIndex(){
	len_text = 0;
	ref_text = NULL;
	len_ref = 0;
	n_factors = 0;
	omit_text = false;
}

RelzIndex::RelzIndex(vector<pair<unsigned int, unsigned int> > &factors, char *full_text, unsigned int _len_text, const char *_ref_text, unsigned int _len_ref, bool _omit_text){
	
	len_text = _len_text;
	len_ref = _len_ref;
	n_factors = factors.size();
	omit_text = _omit_text;
	if( omit_text ){
		ref_text = (char*)_ref_text;
	}
	else{
		ref_text = new char[_len_ref + 1];
		memcpy(ref_text, _ref_text, _len_ref);
		ref_text[_len_ref] = 0;
	}
	
	NanoTimer timer;
	
	cout << "RelzIndex - Inicio (factors: " << factors.size() << ", len_text: " << len_text << ", len_ref: " << len_ref << ")\n";
	
	cout << "RelzIndex - Preparing Factors\n";
	// Factores en version ini, fin (absoluto) y ordenados por ini
	vector<pair<unsigned int, pair<unsigned int, unsigned int> > > factors_sort;
	vector<unsigned int> factors_start;
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
	cout << "RelzIndex - Factors Sorted prepared in " << timer.getMilisec() << "\n";
	timer.reset();
//	cout << "Factors Sorted: \n";
//	for( pair<unsigned int, pair<unsigned int, unsigned int> > factor : factors_sort ){
//		cout << "(" << factor.first << ", " << factor.second.first << ", " << factor.second.second << ")\n";
//	}
	
	// Bit vector S
	cout << "RelzIndex - Preparing Vector S\n";
	bit_vector arr_s = bit_vector(len_ref + n_factors, 0);
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
	cout << "RelzIndex - Vector S prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	bits_s = bits_s_type(arr_s);
	select1_s = bits_s_type::select_1_type(&bits_s);
	select0_s = bits_s_type::select_0_type(&bits_s);
	
	// Permutacion 
	cout << "RelzIndex - Preparing Permutation PI\n";
	pi = int_vector<>(n_factors);
	pi_inv = int_vector<>(n_factors);
	for( unsigned int i = 0; i < n_factors; ++i ){
		pi[i] = factors_sort[i].second.second;
		pi_inv[ factors_sort[i].second.second ] = i;
	}
	
	cout << "RelzIndex - PI prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	// Posiciones finales Ez
	int_vector<> ez = int_vector<>(n_factors);
	for( unsigned int i = 0; i < n_factors; ++i ){
		ez[i] = factors_sort[i].second.first;
	}
	rmq = rmq_type(&ez);
	
	// Bit vector B (inicio de las frases en texto)
	cout << "RelzIndex - Preparing Vector B\n";
	bit_vector arr_b(len_text, 0);
	unsigned int pos_text = 0;
	for( unsigned int i = 0; i < n_factors; ++i ){
		unsigned int len = factors[i].second;
		arr_b[ pos_text ] = 1;
		pos_text += len;
	}
	cout << "RelzIndex - Vector B prepared in " << timer.getMilisec() << "\n";
	timer.reset();

	bits_b = bits_b_type(arr_b);
	select1_b = bits_b_type::select_1_type(&bits_b);
	select0_b = bits_b_type::select_0_type(&bits_b);
	
	// Construccion del FM Index (aunque use el SA original para la compresion, esto es para la busqueda)
	// Construccion con datos en memoria, en un string
	cout << "RelzIndex - Preparing fm_index\n";
	construct_im(fm_index, _ref_text, 1);
	cout << "RelzIndex - fm_index prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	// Preparacion de permutaciones X e Y
	cout << "RelzIndex - Preparing arr X\n";
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
	
//	for( unsigned int i = 0; i < n_factors; ++i ){
//		cout << " arr_x[" << i << "]: " << arr_x[i] << " -> ";
//		char c = 0;
//		FactorsIteratorReverse it(arr_x[i] - 1, n_factors, &select1_s, &select1_b, &select0_b, &pi_inv, (omit_text)?NULL:ref_text, &fm_index, len_text);
//		for(unsigned int k = 0; k < 20 && it.hasNext() && (c = it.next()) != 0; ++k ) 
//			cout << c;
//		cout << "\n";
//	}
//	cout << "-----\n";
	
	cout << "RelzIndex - Preparing arr Y\n";
	
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
	
//	for( unsigned int i = 0; i < n_factors; ++i ){
//		cout << " arr_y[" << i << "]: " << arr_y[i] << " -> ";
//		char c = 0;
//		FactorsIterator it(arr_y[i], n_factors, &select1_s, &select1_b, &select0_b, &pi_inv, (omit_text)?NULL:ref_text, &fm_index, len_text);
//		for(unsigned int k = 0; k < 20 && it.hasNext() && (c = it.next()) != 0; ++k ) 
//			cout << c;
//		cout << " (" << it.length() << ")\n";
//	}
//	cout << "-----\n";
	
	cout << "RelzIndex - X & Y prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	cout << "RelzIndex - Preparing WT\n";
	int_vector<> values_wt(n_factors);
	for( unsigned int i = 0; i < n_factors; ++i ){
		values_wt[i] = arr_y_inv[ arr_x[ i ] ];
	}
	construct_im(wt, values_wt);
	cout << "RelzIndex - WT prepared in " << timer.getMilisec() << "\n";
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
	
	// Compactacion de arreglos descomprimidos
	sdsl::util::bit_compress(pi);
	sdsl::util::bit_compress(pi_inv);
	sdsl::util::bit_compress(arr_x);
	sdsl::util::bit_compress(arr_y);
	
	cout << "RelzIndex - End\n";
	
}

RelzIndex::~RelzIndex(){
	if( ! omit_text && ref_text != NULL ){
		delete [] ref_text;
		ref_text = NULL;
	}
}

void RelzIndex::findTimes(const string &pattern, vector<unsigned int> &results){

//	cout << "RelzIndex::findTimes - Start (\"" << pattern << "\")\n";
	NanoTimer timer;
	
//	cout << "RelzIndex::findTimes - Section A, reference\n";
	
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
	
//	cout << "RelzIndex::findTimes - Section B, ranges\n";
	for(unsigned int i = 1; i < pattern.length(); ++i){
		timer.reset();
		string p1 = pattern.substr(0, i);
		string p1_rev = "";
		for( unsigned int k = 0; k < p1.length(); ++k ){
			p1_rev += p1[ p1.length() - 1 - k ];
		}
		string p2 = pattern.substr(i, pattern.length() - i);
		
		// Rango X
		pair<unsigned int, unsigned int> r1 = getRangeX(p1_rev.c_str());
		querytime_p3 += timer.getNanosec();
		timer.reset();
		
		if( r1.first == (unsigned int)(-1) || r1.second == (unsigned int)(-1) || r1.second < r1.first ){
//			cout << "RelzIndex::findTimes - getRangeX invalido, omitiendo\n";
			continue;
		}
		
		// Rango Y
		pair<unsigned int, unsigned int> r2 = getRangeY(p2.c_str());
		querytime_p3 += timer.getNanosec();
		timer.reset();
		
		if( r2.first == (unsigned int)(-1) || r2.second == (unsigned int)(-1) || r2.second < r2.first ){
//			cout << "RelzIndex::findTimes - getRangeY invalido, omitiendo\n";
			continue;
		}
		
//		cout << "RelzIndex::findTimes - Searching in [" << r1.first << ", " << r1.second << "] x [" << r2.first << ", " << r2.second << "]:\n";
		auto res = wt.range_search_2d(r1.first, r1.second, r2.first, r2.second);
		for (auto point : res.second){
			unsigned int f = arr_y[point.second];
			unsigned int pu = select1_b(f + 1);
//			cout << "RelzIndex::findTimes - Adding pos " << (pu - p1.length()) << "\n";
			results.push_back(pu - p1.length());
			++occs_c;
		}
		
		querytime_p4 += timer.getNanosec();
	}
//	cout << "RelzIndex::findTimes - End\n";
	
}

void RelzIndex::recursive_rmq(unsigned int ini, unsigned int fin, unsigned int min_pos, unsigned int occ_ref, vector<unsigned int> &results){
//	cout << "RelzIndex::recursive_rmq - " << ini << ", " << fin << "\n";
	
	unsigned int pos_max = rmq(ini, fin);
	
//	cout << "RelzIndex::recursive_rmq - Computing factor\n";
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
	
//	cout << "RelzIndex::recursive_rmq - tu: " << tu << ", pu: " << pu << ", lu: " << lu << "\n";
	
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

template <typename ItereatorType>
bool RelzIndex::factorLess(unsigned int factor, const char *pattern, unsigned int len, bool equal){
	if( factor == (unsigned int)(-1) ){
		return true;
	}
	ItereatorType it(factor, n_factors, &select1_s, &select1_b, &select0_b, &pi_inv, (omit_text)?NULL:ref_text, &fm_index, len_text);
	if( it.length() == 0 ){
		return true;
	}
	it.setMaxLength(len);
	
	char c1 = it.next();
	unsigned int pos = 0;
	char c2 = pattern[pos++];
	bool next = true;
//	cout << "RelzIndex::factorLess - " << c1 << " vs " << c2 << "\n";
	while( next && (c1 == c2) ){
		if( it.hasNext() ){
			c1 = it.next();
		}
		else{
			c1 = 0;
			next = false;
		}
		if(pos < len){
			c2 = pattern[pos++];
		}
		else{
			c2 = 0;
			next = false;
		}
//		cout << "RelzIndex::factorLess - " << c1 << " vs " << c2 << "\n";
	}
	if( equal && (c2 == 0) ){
//		cout << "RelzIndex::factorLess - res " << (c1 <= c2) << "\n";
		return true;
	}
	else{
//		cout << "RelzIndex::factorLess - res " << (c1 < c2) << "\n";
		return (c1 < c2);
	}
}

// Notar que, a diferencia de la busqueda en referencia, esta debe ser completa
// Es decir, solo importa el rango que contiene al patron completo
pair<unsigned int, unsigned int> RelzIndex::getRangeY(const char *pattern){
	
	// Version de Revision completa
	
	unsigned int izq = -1;
	unsigned int der = -1;
	unsigned int pat_len = strlen(pattern);
	unsigned int l, h, m, fm;
	
	// Busqueda binaria del lado izquierdo
//	cout << "getRangeY - Busqueda Izquierda\n";
	l = 0;
	h = n_factors-1;
//	cout << "getRangeY - l: " << l << ", h: " << h << "\n";
	while(l < h){
		m = l + ((h-l)>>1);
		fm = arr_y[m];
//		cout << "getRangeY - FactorRev " << fm << " < pattern?\n";
		if( factorLess<FactorsIterator>(fm, pattern, pat_len) ){
//			cout << "getRangeY - caso 1: l = " << (m+1) << "\n";
			l = m+1;
		}
		else{
//			cout << "getRangeY - caso 2: h = " << m << "\n";
			h = m;
		}
	}
//	cout << "getRangeY - end h: " << h << "\n";
	izq = h;
	fm = arr_y[izq];
	if( factorLess<FactorsIterator>(fm, pattern, pat_len) ){
		++izq;
	}
//	cout << "getRangeY - izq: " << izq << "\n";
	
	// Busqueda binaria del lado derecho
//	cout << "getRangeY - Busqueda derecha\n";
	l = izq;
	h = n_factors-1;
//	cout << "getRangeY - l: " << l << ", h: " << h << "\n";
	while(l < h){
		m = l + ((h-l)>>1);
		fm = arr_y[m];
//		cout << "getRangeY - FactorRev " << fm << " < pattern?\n";
		if( factorLess<FactorsIterator>(fm, pattern, pat_len, true) ){
//			cout << "getRangeY - caso 1: l = " << (m+1) << "\n";
			l = m+1;
		}
		else{
//			cout << "getRangeY - caso 2: h = " << m << "\n";
			h = m;
		}
	}
//	cout << "getRangeY - end h: " << h << "\n";
	der = h;
	fm = arr_y[der];
	if( !factorLess<FactorsIterator>(fm, pattern, pat_len, true) ){
		--der;
	}
//	cout << "getRangeY - der: " << der << "\n";
	
//	cout << "getRangeY - result: (" << izq << ", " << der << ")\n";
	return pair<unsigned int, unsigned int>(izq, der);
}

pair<unsigned int, unsigned int> RelzIndex::getRangeX(const char *pattern){
	
	// Version de Revision completa
	
	unsigned int izq = -1;
	unsigned int der = -1;
	unsigned int pat_len = strlen(pattern);
	unsigned int l, h, m, fm;
	
	// Busqueda binaria del lado izquierdo
//	cout << "getRangeX - Busqueda Izquierda\n";
	l = 0;
	h = n_factors-1;
//	cout << "getRangeX - l: " << l << ", h: " << h << "\n";
	while(l < h){
		m = l + ((h-l)>>1);
		fm = arr_x[m];
//		cout << "getRangeX - FactorRev " << fm << " < pattern?\n";
		if( factorLess<FactorsIteratorReverse>(fm-1, pattern, pat_len) ){
//			cout << "getRangeX - caso 1: l = " << (m+1) << "\n";
			l = m+1;
		}
		else{
//			cout << "getRangeX - caso 2: h = " << m << "\n";
			h = m;
		}
	}
//	cout << "getRangeX - end h: " << h << "\n";
	izq = h;
	fm = arr_x[izq];
	if( factorLess<FactorsIteratorReverse>(fm-1, pattern, pat_len) ){
		++izq;
	}
//	cout << "getRangeX - izq: " << izq << "\n";
	
	// Busqueda binaria del lado derecho
//	cout << "getRangeX - Busqueda derecha\n";
	l = izq;
	h = n_factors-1;
//	cout << "getRangeX - l: " << l << ", h: " << h << "\n";
	while(l < h){
		m = l + ((h-l)>>1);
		fm = arr_x[m];
//		cout << "getRangeX - FactorRev " << fm << " < pattern?\n";
		if( factorLess<FactorsIteratorReverse>(fm-1, pattern, pat_len, true) ){
//			cout << "getRangeX - caso 1: l = " << (m+1) << "\n";
			l = m+1;
		}
		else{
//			cout << "getRangeX - caso 2: h = " << m << "\n";
			h = m;
		}
	}
//	cout << "getRangeX - end h: " << h << "\n";
	der = h;
	fm = arr_x[der];
	if( (der > 0) && !factorLess<FactorsIteratorReverse>(fm-1, pattern, pat_len, true) ){
		--der;
	}
//	cout << "getRangeX - der: " << der << "\n";
	
//	cout << "getRangeX - result: (" << izq << ", " << der << ")\n";
	return pair<unsigned int, unsigned int>(izq, der);
}

void RelzIndex::printSize(){
	double total_bytes = 0;
	
	// texto descomprimido
	if( ! omit_text ){
		total_bytes += len_ref;
		cout << "RelzIndex::printSize - Reference Text: " << (8.0*len_ref/len_text) << " bps\n";
	}
	
	total_bytes += size_in_bytes(fm_index);
	cout << "RelzIndex::printSize - fm_index: " << (8.0*size_in_bytes(fm_index)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(rmq);
	cout << "RelzIndex::printSize - rmq: " << (8.0*size_in_bytes(rmq)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(bits_s);
	cout << "RelzIndex::printSize - bits_s: " << (8.0*size_in_bytes(bits_s)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(bits_b);
	cout << "RelzIndex::printSize - bits_b: " << (8.0*size_in_bytes(bits_b)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(pi);
	cout << "RelzIndex::printSize - pi: " << (8.0*size_in_bytes(pi)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(pi_inv);
	cout << "RelzIndex::printSize - pi_inv: " << (8.0*size_in_bytes(pi_inv)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(arr_x);
	cout << "RelzIndex::printSize - arr_x: " << (8.0*size_in_bytes(arr_x)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(arr_y);
	cout << "RelzIndex::printSize - arr_y: " << (8.0*size_in_bytes(arr_y)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(wt);
	cout << "RelzIndex::printSize - wt: " << (8.0*size_in_bytes(wt)/len_text) << " bps\n";
	
	cout << "RelzIndex::printSize - Total " << total_bytes/(1024*1024) << " MB (" << (8.0*total_bytes/len_text) << " bps)\n";
	
}

void RelzIndex::save(const string &file_base){
	cout << "RelzIndex::save - Start (base \"" << file_base << "\")\n";
	
	// Base
	string index_basic_file = file_base + ".base";
	fstream writer(index_basic_file, fstream::out | fstream::trunc);
	// Version of the index
	unsigned char version = 1;
	writer.write((char*)&version, 1);
	// len_text
	writer.write((char*)&len_text, sizeof(int));
	// n_factors
	writer.write((char*)&n_factors, sizeof(int));
	// omit_text
	writer.write((char*)&omit_text, 1);
	// len_ref
	writer.write((char*)&len_ref, sizeof(int));
	// Reference Text
	if( ! omit_text ){
		writer.write((char*)ref_text, len_ref);
	}
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
	
	// wt
	string wt_file = file_base + ".wt";
	store_to_file(wt, wt_file);
	
	cout << "RelzIndex::save - End\n";
}

void RelzIndex::load(const string &file_base){
	cout << "RelzIndex::load - Start (base \"" << file_base << "\")\n";
	
	// Base
	string index_basic_file = file_base + ".base";
	fstream reader(index_basic_file, fstream::in);
	// Version of the index
	unsigned char version = 0;
	reader.read((char*)&version, 1);
	if( version != 1 ){
		cout << "RelzIndex::load - Wrong Version\n";
		return;
	}
	// len_text
	reader.read((char*)&len_text, sizeof(int));
	// n_factors
	reader.read((char*)&n_factors, sizeof(int));
	// omit_text
	reader.read((char*)&omit_text, 1);
	// len_ref
	reader.read((char*)&len_ref, sizeof(int));
	// Reference Text
	if( ref_text != NULL ){
		delete [] ref_text;
		ref_text = NULL;
	}
	if( ! omit_text ){
		ref_text = new char[len_ref + 1];
		reader.read((char*)ref_text, len_ref);
		ref_text[len_ref] = 0;
	}
	// Close Base
	reader.close();
	
	// fm_index
	cout << "RelzIndex::load - fm_index\n";
	string fm_index_file = file_base + ".fm";
	load_from_file(fm_index, fm_index_file);
	
	// rmq
	cout << "RelzIndex::load - rmq\n";
	string rmq_file = file_base + ".rmq";
	load_from_file(rmq, rmq_file);
	
	// bits_s
	cout << "RelzIndex::load - bits_s\n";
	string bits_s_file = file_base + ".arrs";
	load_from_file(bits_s, bits_s_file);
	select1_s = bits_s_type::select_1_type(&bits_s);
	select0_s = bits_s_type::select_0_type(&bits_s);
	
	// bits_b
	cout << "RelzIndex::load - bits_b\n";
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
	cout << "RelzIndex::load - wt\n";
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
	
	cout << "RelzIndex::load - End\n";
	
}




















