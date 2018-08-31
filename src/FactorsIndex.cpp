#include "FactorsIndex.h"

FactorsIndex::FactorsIndex(){
	len_text = 0;
	ref_text = NULL;
	len_ref = 0;
	n_factors = 0;
	omit_text = false;
	acelerar_rmq = false;
}

FactorsIndex::FactorsIndex(vector<pair<unsigned int, unsigned int> > &factors, char *full_text, unsigned int _len_text, const char *_ref_text, unsigned int _len_ref, bool _omit_text){
	
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
	acelerar_rmq = false;
	
	NanoTimer timer;
	
	cout << "FactorsIndex - Inicio (factors: " << factors.size() << ", len_text: " << len_text << ", len_ref: " << len_ref << ")\n";
	
	cout << "FactorsIndex - Preparing Factors\n";
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
	cout << "FactorsIndex - Factors Sorted prepared in " << timer.getMilisec() << "\n";
	timer.reset();
//	cout << "Factors Sorted: \n";
//	for( pair<unsigned int, pair<unsigned int, unsigned int> > factor : factors_sort ){
//		cout << "(" << factor.first << ", " << factor.second.first << ", " << factor.second.second << ")\n";
//	}
	
	// Bit vector S
	cout << "FactorsIndex - Preparing Vector S\n";
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
	cout << "FactorsIndex - Vector S prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	bits_s = bits_s_type(arr_s);
	select1_s = bits_s_type::select_1_type(&bits_s);
	select0_s = bits_s_type::select_0_type(&bits_s);
	
	// Permutacion 
	cout << "FactorsIndex - Preparing Permutation PI\n";
	int_vector<> _pi(n_factors);
	int_vector<> _pi_inv(n_factors);
	pi = _pi;
	pi_inv = _pi_inv;	
	for( unsigned int i = 0; i < n_factors; ++i ){
		pi[i] = factors_sort[i].second.second;
		pi_inv[ factors_sort[i].second.second ] = i;
	}
	inv_perm_support<> _perm_inv(&pi);
	inv_perm_support<> _perm(&pi_inv);
	perm_inv = _perm_inv;
	perm = _perm;
	
//	cout << "FactorsIndex - Testing Permutations\n";
//	for( unsigned int i = 0; i < n_factors; ++i ){
//		cout << "perm: " << perm[i] << ", perm_inv: " << perm_inv[i] << ", pi: " << pi[i] << ", pi_inv: " << pi_inv[i] << "\n";
//	}
	
	cout << "FactorsIndex - PI prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	// Posiciones finales Ez
	int_vector<> _ez(n_factors);
	ez = _ez;
	for( unsigned int i = 0; i < n_factors; ++i ){
		ez[i] = factors_sort[i].second.first;
	}
//	rmq = rmq_succinct_sct<false, bp_support_sada<256,32,rank_support_v5<> > >(&ez);
	rmq = rmq_type(&ez);
	
	// Bit vector B (inicio de las frases en texto)
	cout << "FactorsIndex - Preparing Vector B\n";
	bit_vector _arr_b(len_text, 0);
	arr_b = _arr_b;
	unsigned int pos_text = 0;
	for( unsigned int i = 0; i < n_factors; ++i ){
		unsigned int len = factors[i].second;
		arr_b[ pos_text ] = 1;
		pos_text += len;
	}
	cout << "FactorsIndex - Vector B prepared in " << timer.getMilisec() << "\n";
	timer.reset();

	bits_b = bits_b_type(arr_b);
	select1_b = bits_b_type::select_1_type(&bits_b);
	select0_b = bits_b_type::select_0_type(&bits_b);
	
	// Construccion del FM Index (aunque use el SA original para la compresion, esto es para la busqueda)
	// Construccion con datos en memoria, en un string
	cout << "FactorsIndex - Preparing fm_index\n";
	construct_im(fm_index, _ref_text, 1);
	cout << "FactorsIndex - fm_index prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	// Preparacion de permutaciones X e Y
	cout << "FactorsIndex - Preparing arr X\n";
	vector<unsigned int> arr_x(n_factors);
	for( unsigned int i = 0; i < n_factors; ++i ){
		arr_x[i] = i;
	}
//	FactorsIteratorReverseComparator comp_rev(n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref_text, &fm_index, len_text);
	FactorsFastIteratorReverseComparator comp_rev(&factors_start, full_text, len_text);
	cout << "FactorsIndex - stable_sort...\n";
	stable_sort(arr_x.begin(), arr_x.end(), comp_rev);
	cout << "FactorsIndex - Ok\n";
	int_vector<> _pre_x_inv(n_factors);
	pre_x_inv = _pre_x_inv;
	for( unsigned int i = 0; i < n_factors; ++i ){
		pre_x_inv[ arr_x[i] ] = i;
	}
	inv_perm_support<> _perm_x(&pre_x_inv);
	perm_x = _perm_x;
	for( unsigned int i = 0; i < n_factors; ++i ){
		cout << " arr_x[" << i << "]: " << arr_x[i] << " -> ";
		char c = 0;
		for(unsigned int k = 0; k < 10 && (c = getCharRev(arr_x[i] - 1, k)) != 0; ++k ) 
			cout << c;
		cout << "\n";
	}
	cout << "-----\n";
	
	cout << "FactorsIndex - Preparing arr Y\n";
	vector<unsigned int> arr_y(n_factors);
	for( unsigned int i = 0; i < n_factors; ++i ){
		arr_y[i] = i;
	}
//	FactorsIteratorComparator comp(n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref_text, &fm_index, len_text);
	FactorsFastIteratorComparator comp(&factors_start, full_text, len_text);
	stable_sort(arr_y.begin(), arr_y.end(), comp);
	int_vector<> _pre_y(n_factors);
	int_vector<> _pre_y_inv(n_factors);
	pre_y = _pre_y;
	pre_y_inv = _pre_y_inv;
	for( unsigned int i = 0; i < n_factors; ++i ){
		pre_y[i] = arr_y[i];
		pre_y_inv[ arr_y[i] ] = i;
	}
	inv_perm_support<> _perm_y(&pre_y_inv);
	inv_perm_support<> _perm_y_inv(&pre_y);
	perm_y = _perm_y;
	perm_y_inv = _perm_y_inv;
	for( unsigned int i = 0; i < n_factors; ++i ){
		cout << " arr_y[" << i << "]: " << perm_y[i] << " -> ";
		char c = 0;
		for(unsigned int k = 0; k < 10 && (c = getChar(perm_y[i], k)) != 0; ++k ) 
			cout << c;
		cout << "\n";
	}
	cout << "-----\n";
	cout << "FactorsIndex - X & Y prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	cout << "FactorsIndex - Preparing WT\n";
	int_vector<> _values_wt(n_factors);
	values_wt = _values_wt;
	for( unsigned int i = 0; i < n_factors; ++i ){
		values_wt[i] = perm_y_inv[ arr_x[ i ] ];
	}
	
	construct_im(wt, values_wt);
	cout << "FactorsIndex - WT prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	// Prueba de aceleracion de recursive_rmq almacenando los datos de los factores descomprimidos
	acelerar_rmq = true;
	for( unsigned int i = 0; i < n_factors; ++i ){
		unsigned int tu = select1_s(i + 1) - i;
		unsigned int pu = select1_b(perm[i] + 1);
		unsigned int lu = select1_b(perm[i] + 2) - pu;
		arr_tu.push_back(tu);
		arr_pu.push_back(pu);
		arr_lu.push_back(lu);
	}
	
	cout << "FactorsIndex - End\n";
	
}

FactorsIndex::~FactorsIndex(){
	if( ! omit_text && ref_text != NULL ){
		delete [] ref_text;
		ref_text = NULL;
	}
}

void FactorsIndex::findTimes(const string &pattern, vector<unsigned int> &results){

//	cout << "FactorsIndex::findTimes - Start (\"" << pattern << "\")\n";
	NanoTimer timer;
	
//	cout << "FactorsIndex::findTimes - Section A, reference\n";
	
	size_t m = pattern.size();
	size_t occs = sdsl::count(fm_index, pattern.begin(), pattern.end());
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
	
//	cout << "FactorsIndex::findTimes - Section B, ranges\n";
	for(unsigned int i = 1; i < pattern.length(); ++i){
		timer.reset();
		string p1 = pattern.substr(0, i);
		string p1_rev = "";
		for( unsigned int k = 0; k < p1.length(); ++k ){
			p1_rev += p1[ p1.length() - 1 - k ];
		}
		string p2 = pattern.substr(i, pattern.length() - i);
		pair<unsigned int, unsigned int> r1 = getRangeX(p1_rev.c_str());
		pair<unsigned int, unsigned int> r2 = getRangeY(p2.c_str());
		querytime_p3 += timer.getNanosec();
		timer.reset();
		
		if( r1.second == (unsigned int)(-1) || r1.second < r1.first
			|| r2.second == (unsigned int)(-1) || r2.second < r2.first ){
//			cout << "FactorsIndex::findTimes - Invalid ranges, omitting...\n";
			continue;
		}
		
//		cout << "FactorsIndex::findTimes - Searching in [" << r1.first << ", " << r1.second << "] x [" << r2.first << ", " << r2.second << "]:\n";
		auto res = wt.range_search_2d(r1.first, r1.second, r2.first, r2.second);
		for (auto point : res.second){
			unsigned int f = perm_y[point.second];
//			unsigned int cur_perm = perm_inv[f];
//			unsigned int pu = select1_b(perm[cur_perm] + 1);
			unsigned int cur_perm = pi_inv[f];
			unsigned int pu = select1_b(pi[cur_perm] + 1);
			results.push_back(pu - p1.length());
		}
		
		querytime_p4 += timer.getNanosec();
	}
//	cout << "FactorsIndex::findTimes - End\n";
	
}

void FactorsIndex::find(const string &pattern, vector<unsigned int> &results){

//	cout << "FactorsIndex::find - Start (\"" << pattern << "\")\n";
	
//	cout << "FactorsIndex::find - Section A, reference\n";
	
	size_t m = pattern.size();
//	cout << "FactorsIndex::find - 1\n";
	size_t occs = sdsl::count(fm_index, pattern.begin(), pattern.end());
//	cout << "FactorsIndex::find - 2\n";
	if( occs > 0 ){
		auto locations = locate(fm_index, pattern.begin(), pattern.begin()+m);
//		cout << "FactorsIndex::find - 3\n";
		sort(locations.begin(), locations.end());
		for( unsigned int i = 0; i < occs; ++i ){
			unsigned int occ_i = locations[i];
			// Comprobar los factores que cuben esta ocurrencia (el string ref[occ_i, occ_i + m - 1])
			unsigned int select = select0_s(occ_i + 1);
			unsigned int pos_ez = select - 1 - occ_i;
			// Now the recursive search in rmq (0 - pos_ez)
//			cout << "FactorsIndex::find - 4 (occ_i: " << occ_i << ", select: " << select << ", pos_ez: " << pos_ez << ")\n";
			if( occ_i >= select ){
				continue;
			}
			recursive_rmq(0, pos_ez, (occ_i + m), occ_i, results);
//			cout << "FactorsIndex::find - 5\n";
		}
	}
	
//	cout << "FactorsIndex::find - Section B, ranges\n";
	for(unsigned int i = 1; i < pattern.length(); ++i){
		string p1 = pattern.substr(0, i);
		string p1_rev = "";
		for( unsigned int k = 0; k < p1.length(); ++k ){
			p1_rev += p1[ p1.length() - 1 - k ];
		}
		string p2 = pattern.substr(i, pattern.length() - i);
//		cout << "-----  getRangeX (" << p1_rev << ") -----\n";
		pair<unsigned int, unsigned int> r1 = getRangeX(p1_rev.c_str());
//		cout << "-----\n";
		if( r1.first == (unsigned int)(-1) || r1.second == (unsigned int)(-1) || r1.second < r1.first ){
			continue;
		}


//		cout << "-----  getRangeY (" << p1_rev << ") -----\n";
		pair<unsigned int, unsigned int> r2 = getRangeY(p2.c_str());
//		cout << "-----\n";
		
		if( r2.first == (unsigned int)(-1) || r2.second == (unsigned int)(-1) || r2.second < r2.first ){
			continue;
		}
		
//		cout << "FactorsIndex::find - Searching in [" << r1.first << ", " << r1.second << "] x [" << r2.first << ", " << r2.second << "]:\n";
		auto res = wt.range_search_2d(r1.first, r1.second, r2.first, r2.second);
		for (auto point : res.second){
			unsigned int f = perm_y[point.second];
//			unsigned int cur_perm = perm_inv[f];
//			unsigned int pu = select1_b(perm[cur_perm] + 1);
			unsigned int cur_perm = pi_inv[f];
			unsigned int pu = select1_b(pi[cur_perm] + 1);
//			cout << " -> Adding " << (pu - p1.length()) << "\n";
			results.push_back(pu - p1.length());
		}
		
	}
//	cout << "FactorsIndex::find - End\n";
	
}

void FactorsIndex::recursive_rmq(unsigned int ini, unsigned int fin, unsigned int min_pos, unsigned int occ_ref, vector<unsigned int> &results){
//	cout << "FactorsIndex::recursive_rmq - " << ini << ", " << fin << "\n";
	
	unsigned int pos_max = rmq(ini, fin);
	
//	cout << "FactorsIndex::recursive_rmq - Computing factor\n";
//	unsigned int tu = select1_s(pos_max + 1) - pos_max;
//	unsigned int pu = select1_b(perm[pos_max] + 1);
//	unsigned int lu = select1_b(perm[pos_max] + 2) - pu;
	// Prueba de aceleracion
	assert(pos_max < n_factors);
	unsigned int tu = arr_tu[pos_max];
	unsigned int pu = arr_pu[pos_max];
	unsigned int lu = arr_lu[pos_max];
//	cout << "FactorsIndex::recursive_rmq - tu: " << tu << ", pu: " << pu << ", lu: " << lu << "\n";
	
	if( tu + lu < min_pos ){
		return;
	}
	else{
//		cout << " -> Adding " << (pu + (occ_ref - tu)) << "\n";
		results.push_back(pu + (occ_ref - tu));
	}
	
	if( (pos_max > 0) && (ini < pos_max) ){
		recursive_rmq(ini, pos_max-1, min_pos, occ_ref, results);
	}
	if( pos_max < fin ){
		recursive_rmq(pos_max+1, fin, min_pos, occ_ref, results);
	}
}

char FactorsIndex::getChar(unsigned int factor, unsigned int pos, unsigned int max_len){
	
	// Iterators cache
	if( mapa_iterators.find(factor) == mapa_iterators.end() ){
		mapa_iterators[factor] = FactorsIterator(factor, n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, &pi, &pi_inv, (omit_text)?NULL:ref_text, &fm_index, len_text);
	}
	FactorsIterator it = mapa_iterators[factor];
	if( pos >= it.length() ){
		return 0;
	}
	it.setMaxLength(max_len);
	if( it.position() > pos ){
		it.reset();
	}
	
	char c = 0;
	while( it.position() <= pos ){
		c = it.next();
	}
	return c;
}

char FactorsIndex::getCharRev(unsigned int factor, unsigned int pos, unsigned int max_len){
	
	if( factor == (unsigned int)(-1) ){
		return 0;
	}
	
	// Iterators cache
	if( mapa_iterators_rev.find(factor) == mapa_iterators_rev.end() ){
		mapa_iterators_rev[factor] = FactorsIteratorReverse(factor, n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, &pi, &pi_inv, (omit_text)?NULL:ref_text, &fm_index, len_text);
	}
	FactorsIteratorReverse it = mapa_iterators_rev[factor];
	if( pos >= it.length() ){
		return 0;
	}
	it.setMaxLength(max_len);
	if( it.position() > pos ){
		it.reset();
	}
	char c = 0;
	while( it.position() <= pos ){
		c = it.next();
	}
	return c;
}

// Notar que, a diferencia de la busqueda en referencia, esta debe ser completa
// Es decir, solo importa el rango que contiene al patron completo
pair<unsigned int, unsigned int> FactorsIndex::getRangeY(const char *pattern){
	
	unsigned int pat_len = strlen(pattern);
	unsigned int izq = 0;
	unsigned int der = n_factors-1;
	unsigned int cur_pos = 0;
	
//	cout << "getRangeY - Inicio (pat_len: " << pat_len << ", izq: " << izq << ", der: " << der << ")\n";
	
	for( ; cur_pos < pat_len; ++cur_pos ){
//		cout << "getRangeY - cur_pos: " << cur_pos << " (pattern[" << cur_pos << "]: " << pattern[cur_pos] << ")\n";
		
		unsigned int l = izq;
		unsigned int h = der;
		unsigned int m;
		unsigned int fm;
		char c;
		unsigned int text_len;
		
		// Busqueda binaria del lado izquierdo
//		cout << "getRangeY - l: " << l << ", h: " << h << "\n";
		while(l < h){
			m = l + ((h-l)>>1);
			fm = perm_y[m];
			c = getChar(fm, cur_pos, pat_len-cur_pos);
			text_len = mapa_iterators[fm].length();
//			cout << "getRangeY - m: " << m << ", fm: " << fm << ", c: " << c << ", text_len: " << text_len << "\n";
			if( cur_pos > text_len || (unsigned char)(c) < (unsigned char)(pattern[cur_pos]) ){
//				cout << "getRangeY - caso 1: l = " << (m+1) << "\n";
				l = m+1;
			}
			else{
//				cout << "getRangeY - caso 2: h = " << m << "\n";
				h = m;
			}
		}
		izq = h;
		fm = perm_y[izq];
		c = getChar(fm, cur_pos, pat_len-cur_pos);
		text_len = mapa_iterators[fm].length();
		if( (cur_pos < text_len) && (unsigned char)(c) < (unsigned char)(pattern[cur_pos]) ){
			++izq;
		}
//		cout << "getRangeY - izq: " << izq << "\n";
//		cout << "getRangeY - -----\n";
		
		// Busqueda binaria del lado derecho
		l = izq;
		h = der;
//		cout << "getRangeY - l: " << l << ", h: " << h << "\n";
		while(l < h){
			m = l + ((h-l)>>1);
			fm = perm_y[m];
			c = getChar(fm, cur_pos);
			text_len = mapa_iterators[fm].length();
//			cout << "getRangeY - m: " << m << ", fm: " << fm << ", c: " << c << ", text_len: " << text_len << "\n";
			if( cur_pos > text_len || (unsigned char)(c) <= (unsigned char)(pattern[cur_pos]) ){
//				cout << "getRangeY - caso 1: l = " << (m+1) << "\n";
				l = m+1;
			}
			else{
//				cout << "getRangeY - caso 2: h = " << m << "\n";
				h = m;
			}
		}
		der = h;
		fm = perm_y[der];
		c = getChar(fm, cur_pos);
		text_len = mapa_iterators[fm].length();
		if( (der > 0) && (cur_pos < text_len) && (unsigned char)(c) > (unsigned char)(pattern[cur_pos]) ){
			--der;
		}
//		cout << "getRangeY - der: " << der << "\n";
//		cout << "getRangeY - -----\n";
		
	}
	
//	cout << "getRangeY - result: (" << izq << ", " << der << ")\n";
	return pair<unsigned int, unsigned int>(izq, der);
}

pair<unsigned int, unsigned int> FactorsIndex::getRangeX(const char *pattern){
	
	unsigned int pat_len = strlen(pattern);
	unsigned int izq = 0;
	unsigned int der = n_factors-1;
	unsigned int cur_pos = 0;
	
//	cout << "getRangeX - Inicio (pat_len: " << pat_len << ", izq: " << izq << ", der: " << der << ")\n";
	
	for( ; cur_pos < pat_len; ++cur_pos ){
//		cout << "getRangeX - cur_pos: " << cur_pos << " (pattern[" << cur_pos << "]: " << pattern[cur_pos] << ")\n";
		
		unsigned int l = izq;
		unsigned int h = der;
		unsigned int m;
		unsigned int fm;
		char c;
		unsigned int text_len;
		
		// Busqueda binaria del lado izquierdo
//		cout << "getRangeX - l: " << l << ", h: " << h << "\n";
		while(l < h){
			m = l + ((h-l)>>1);
//			cout << "getRangeX - l: perm_x[" << m << "]...\n";
			fm = perm_x[m];
//			cout << "getRangeX - l: getCharRev(" << (fm-1) << ", " << cur_pos << ")\n";
			c = getCharRev(fm-1, cur_pos, pat_len-cur_pos);
			text_len = mapa_iterators_rev[fm-1].length();
//			cout << "getRangeX - m: " << m << ", fm: " << fm << ", c: " << c << ", text_len: " << text_len << "\n";
			if( cur_pos > text_len || (unsigned char)(c) < (unsigned char)(pattern[cur_pos]) ){
//				cout << "getRangeX - caso 1: l = " << (m+1) << "\n";
				l = m+1;
			}
			else{
//				cout << "getRangeX - caso 2: h = " << m << "\n";
				h = m;
			}
		}
//		cout << "getRangeX - end h: " << h << "\n";
		izq = h;
		fm = perm_x[izq];
		c = getCharRev(fm-1, cur_pos, pat_len-cur_pos);
		text_len = mapa_iterators_rev[fm-1].length();
		if( (cur_pos < text_len) && (unsigned char)(c) < (unsigned char)(pattern[cur_pos]) ){
			++izq;
		}
//		cout << "getRangeX - izq: " << izq << "\n";
//		cout << "getRangeX - -----\n";
		
		// Busqueda binaria del lado derecho
		l = izq;
		h = der;
//		cout << "getRangeX - l: " << l << ", h: " << h << "\n";
		while(l < h){
			m = l + ((h-l)>>1);
			fm = perm_x[m];
			c = getCharRev(fm-1, cur_pos);
			text_len = mapa_iterators_rev[fm-1].length();
//			cout << "getRangeX - m: " << m << ", fm: " << fm << ", c: " << c << ", text_len: " << text_len << "\n";
			if( cur_pos > text_len || (unsigned char)(c) <= (unsigned char)(pattern[cur_pos]) ){
//				cout << "getRangeX - caso 1: l = " << (m+1) << "\n";
				l = m+1;
			}
			else{
//				cout << "getRangeX - caso 2: h = " << m << "\n";
				h = m;
			}
		}
		der = h;
		fm = perm_x[der];
		c = getCharRev(fm-1, cur_pos);
		text_len = mapa_iterators_rev[fm-1].length();
		if( (der > 0) && (cur_pos < text_len) && (unsigned char)(c) > (unsigned char)(pattern[cur_pos]) ){
			--der;
		}
//		cout << "getRangeX - der: " << der << "\n";
//		cout << "getRangeX - -----\n";
		
	}
	
//	cout << "getRangeX - result: (" << izq << ", " << der << ")\n";
	return pair<unsigned int, unsigned int>(izq, der);
}

void FactorsIndex::printSize(){
	double total_bytes = 0;
	
	// texto descomprimido
	if( ! omit_text ){
		total_bytes += len_ref;
		cout << "FactorsIndex::printSize - Reference Text: " << (8.0*len_ref/len_text) << " bps\n";
	}
	
//	csa_wt<> fm_index;
	total_bytes += size_in_bytes(fm_index);
	cout << "FactorsIndex::printSize - fm_index: " << (8.0*size_in_bytes(fm_index)/len_text) << " bps\n";
	
//	rmq_succinct_sct<false, bp_support_sada<256,32,rank_support_v5<> > > rmq;
	total_bytes += size_in_bytes(rmq);
	cout << "FactorsIndex::printSize - rmq: " << (8.0*size_in_bytes(rmq)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(bits_s);
	cout << "FactorsIndex::printSize - bits_s: " << (8.0*size_in_bytes(bits_s)/len_text) << " bps\n";
	
//	inv_perm_support<> perm_inv;
	total_bytes += size_in_bytes(perm_inv);
	cout << "FactorsIndex::printSize - perm_inv: " << (8.0*size_in_bytes(perm_inv)/len_text) << " bps\n";
	
//	inv_perm_support<> perm;
	total_bytes += size_in_bytes(perm);
	cout << "FactorsIndex::printSize - perm: " << (8.0*size_in_bytes(perm)/len_text) << " bps\n";
	
	total_bytes += size_in_bytes(bits_b);
	cout << "FactorsIndex::printSize - bits_b: " << (8.0*size_in_bytes(bits_b)/len_text) << " bps\n";
//	
//	inv_perm_support<> perm_x;
	total_bytes += size_in_bytes(perm_x);
	cout << "FactorsIndex::printSize - perm_x: " << (8.0*size_in_bytes(perm_x)/len_text) << " bps\n";
	
//	inv_perm_support<> perm_y;
	total_bytes += size_in_bytes(perm_y);
	cout << "FactorsIndex::printSize - perm_y: " << (8.0*size_in_bytes(perm_y)/len_text) << " bps\n";
	
//	inv_perm_support<> perm_y_inv;
	total_bytes += size_in_bytes(perm_y_inv);
	cout << "FactorsIndex::printSize - perm_y_inv: " << (8.0*size_in_bytes(perm_y_inv)/len_text) << " bps\n";
//	
//	wt_int<rrr_vector<63>> wt;
	total_bytes += size_in_bytes(wt);
	cout << "FactorsIndex::printSize - wt: " << (8.0*size_in_bytes(wt)/len_text) << " bps\n";
	
	cout << "FactorsIndex::printSize - Total " << total_bytes/(1024*1024) << " MB (" << (8.0*total_bytes/len_text) << " bps)\n";
	
}

void FactorsIndex::save(const string &file_base){
	cout << "FactorsIndex::save - Start (base \"" << file_base << "\")\n";
	
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
	
	// perm
	string pi_file = file_base + ".pi";
	store_to_file(perm, pi_file);
	
	// perm_inv
	string pi1_file = file_base + ".pi1";
	store_to_file(perm_inv, pi1_file);
	
	// perm_x
	string x_file = file_base + ".x";
	store_to_file(perm_x, x_file);
	
	// perm_y
	string y_file = file_base + ".y";
	store_to_file(perm_y, y_file);
	
	// perm_y_inv
	string y1_file = file_base + ".y1";
	store_to_file(perm_y_inv, y1_file);
	
	// wt
	string wt_file = file_base + ".wt";
	store_to_file(wt, wt_file);
	
	cout << "FactorsIndex::save - End\n";
}

void FactorsIndex::load(const string &file_base){
	cout << "FactorsIndex::load - Start (base \"" << file_base << "\")\n";
	
	// Base
	string index_basic_file = file_base + ".base";
	fstream reader(index_basic_file, fstream::in);
	// Version of the index
	unsigned char version = 0;
	reader.read((char*)&version, 1);
	if( version != 1 ){
		cout << "FactorsIndex::load - Wrong Version\n";
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
	if( ! omit_text ){
		ref_text = new char[len_ref + 1];
		reader.read((char*)ref_text, len_ref);
		ref_text[len_ref] = 0;
	}
	// Close Base
	reader.close();
	
	// fm_index
	cout << "FactorsIndex::load - fm_index\n";
	string fm_index_file = file_base + ".fm";
	load_from_file(fm_index, fm_index_file);
	
	// rmq
	cout << "FactorsIndex::load - rmq\n";
	string rmq_file = file_base + ".rmq";
	load_from_file(rmq, rmq_file);
	
	// bits_s
	cout << "FactorsIndex::load - bits_s\n";
	string bits_s_file = file_base + ".arrs";
	load_from_file(bits_s, bits_s_file);
	select1_s = bits_s_type::select_1_type(&bits_s);
	select0_s = bits_s_type::select_0_type(&bits_s);
	
	// bits_b
	cout << "FactorsIndex::load - bits_b\n";
	string bits_b_file = file_base + ".arrb";
	load_from_file(bits_b, bits_b_file);
	select1_b = bits_b_type::select_1_type(&bits_b);
	select0_b = bits_b_type::select_0_type(&bits_b);
	
	// perm
	cout << "FactorsIndex::load - perm\n";
	string pi_file = file_base + ".pi";
	load_from_file(perm, pi_file);
	
	// perm_inv
	string pi1_file = file_base + ".pi1";
	load_from_file(perm_inv, pi1_file);
	
	// perm_x
	string x_file = file_base + ".x";
	load_from_file(perm_x, x_file);
	
	// perm_y
	string y_file = file_base + ".y";
	load_from_file(perm_y, y_file);
	
	// perm_y_inv
	string y1_file = file_base + ".y1";
	load_from_file(perm_y_inv, y1_file);
	
	// wt
	cout << "FactorsIndex::load - wt\n";
	string wt_file = file_base + ".wt";
	load_from_file(wt, wt_file);
	
	cout << "FactorsIndex::load - End\n";
	
}




















