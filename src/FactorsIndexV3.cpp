#include "FactorsIndexV3.h"

FactorsIndexV3::FactorsIndexV3(){
	len_text = 0;
	ref_text = NULL;
	len_ref = 0;
	n_factors = 0;
	karp_rabin = NULL;
}

FactorsIndexV3::FactorsIndexV3(vector<pair<unsigned int, unsigned int> > &factors, char *full_text, unsigned int _len_text, const char *_ref_text, unsigned int _len_ref, KarpRabin *_karp_rabin, const char *kr_frases_file, bool load_kr_frases){
	
	len_text = _len_text;
	ref_text = _ref_text;
	len_ref = _len_ref;
	n_factors = factors.size();
	karp_rabin = _karp_rabin;
	NanoTimer timer;
	
	cout << "FactorsIndexV3 - Inicio (factors: " << factors.size() << ", len_text: " << len_text << ", len_ref: " << len_ref << ")\n";
	
	cout << "FactorsIndexV3 - Preparing Factors\n";
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
	cout << "FactorsIndexV3 - Factors Sorted prepared in " << timer.getMilisec() << "\n";
	timer.reset();
//	cout << "Factors Sorted: \n";
//	for( pair<unsigned int, pair<unsigned int, unsigned int> > factor : factors_sort ){
//		cout << "(" << factor.first << ", " << factor.second.first << ", " << factor.second.second << ")\n";
//	}
	
	// Bit vector S
	cout << "FactorsIndexV3 - Preparing Vector S\n";
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
	cout << "FactorsIndexV3 - Vector S prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	rrr_vector<127> _rrr_s(arr_s);
	rrr_s = _rrr_s;
	rrr_vector<127>::select_1_type _select1_s(&rrr_s);
	select1_s = _select1_s;
	
	// Notar que la posicion del select DEBE empezar desde 1, no desde 0
	// De este modo, hay que sumar 1 a las posiciones de la ref para buscar en S
	rrr_vector<127>::select_0_type _select0_s(&rrr_s);
	select0_s = _select0_s;
	
	// Permutacion 
	cout << "FactorsIndexV3 - Preparing Permutation PI\n";
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
	cout << "FactorsIndexV3 - PI prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	// Posiciones finales Ez
	int_vector<> _ez(n_factors);
	ez = _ez;
	for( unsigned int i = 0; i < n_factors; ++i ){
		ez[i] = factors_sort[i].second.first;
	}
	// rmq_succinct_sct<> rmq(&ez);
	rmq_succinct_sct<false, bp_support_sada<256,32,rank_support_v5<> > > _rmq(&ez);
	rmq = _rmq;
	
	// Bit vector B (inicio de las frases en texto)
	cout << "FactorsIndexV3 - Preparing Vector B\n";
	bit_vector _arr_b(len_text, 0);
	arr_b = _arr_b;
	unsigned int pos_text = 0;
	for( unsigned int i = 0; i < n_factors; ++i ){
		unsigned int len = factors[i].second;
		arr_b[ pos_text ] = 1;
		pos_text += len;
	}
	cout << "FactorsIndexV3 - Vector B prepared in " << timer.getMilisec() << "\n";
	timer.reset();

	rrr_vector<127> _rrr_b(arr_b);
	rrr_b = _rrr_b;
	rrr_vector<127>::select_1_type _select1_b(&rrr_b);
	rrr_vector<127>::select_0_type _select0_b(&rrr_b);
	select1_b = _select1_b;
	select0_b = _select0_b;
	
	// Construccion del FM Index (aunque use el SA original para la compresion, esto es para la busqueda)
	// Construccion con datos en memoria, en un string
	cout << "FactorsIndexV3 - Preparing fm_index\n";
	construct_im(fm_index, ref_text, 1);
	cout << "FactorsIndexV3 - fm_index prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	// Preparacion de permutaciones X e Y
	cout << "FactorsIndexV3 - Preparing arr X\n";
	vector<unsigned int> arr_x(n_factors);
	for( unsigned int i = 0; i < n_factors; ++i ){
		arr_x[i] = i;
	}
//	FactorsIteratorReverseComparator comp_rev(n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref_text, &fm_index, len_text);
	FactorsFastIteratorReverseComparator comp_rev(&factors_start, full_text, len_text);
	cout << "FactorsIndexV3 - stable_sort...\n";
	stable_sort(arr_x.begin(), arr_x.end(), comp_rev);
	cout << "FactorsIndexV3 - Ok\n";
	int_vector<> _pre_x_inv(n_factors);
	pre_x_inv = _pre_x_inv;
	for( unsigned int i = 0; i < n_factors; ++i ){
		pre_x_inv[ arr_x[i] ] = i;
	}
	inv_perm_support<> _perm_x(&pre_x_inv);
	perm_x = _perm_x;
//	for( unsigned int i = 0; i < n_factors; ++i ){
//		cout << " arr_x[" << i << "]: " << arr_x[i] << " -> ";
//		char c = 0;
//		for(unsigned int k = 0; k < 10 && (c = getCharRev(arr_x[i] - 1, k)) != 0; ++k ) 
//			cout << c;
//		cout << "\n";
//	}
//	cout << "-----\n";
	
	cout << "FactorsIndexV3 - Preparing arr Y\n";
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
//	for( unsigned int i = 0; i < n_factors; ++i ){
//		cout << " arr_y[" << i << "]: " << perm_y[i] << " -> ";
//		char c = 0;
//		for(unsigned int k = 0; k < 10 && (c = getChar(perm_y[i], k)) != 0; ++k ) 
//			cout << c;
//		cout << "\n";
//	}
//	cout << "-----\n";
	cout << "FactorsIndexV3 - X & Y prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	cout << "FactorsIndexV3 - Preparing WT\n";
	int_vector<> _values_wt(n_factors);
	values_wt = _values_wt;
	for( unsigned int i = 0; i < n_factors; ++i ){
		values_wt[i] = perm_y_inv[ arr_x[ i ] ];
	}
	
	construct_im(wt, values_wt);
	cout << "FactorsIndexV3 - WT prepared in " << timer.getMilisec() << "\n";
	timer.reset();
	
	cout << "FactorsIndexV3 - Preparing KarpRobin Structures\n";
	
	// En principio se necesita un arreglo de un hash por cada log(n) caracteres de la referencia
	// Ademas necesito la firma de cada sufijo de factor quizas?
	
	BitsUtils bits_utils;
	unsigned int log_n = bits_utils.n_bits(len_ref);
	cout << "FactorsIndexV3 - Adding hash for Reference prefixes\n";
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
	cout << "FactorsIndexV3 - Adding hash for Frases prefixes\n";
//	kr_frases_file, bool load_kr_frases
	unsigned int max_line = 1000000;
	char buff[max_line + 1];
	if(load_kr_frases){
		fstream reader(kr_frases_file, fstream::in);
		unsigned int pos = 0;
		unsigned int cur_pos = 0;
		unsigned long long kr = 0;
		while( reader.good() ){
			reader.getline(buff, max_line);
			unsigned int n_read = reader.gcount();
			if( n_read < 1 ){
				continue;
			}
			
			string line(buff);
			stringstream toks(line);
			
			toks >> pos;
			toks >> kr;
			
			if( pos == cur_pos ){
				cout << "Agregando " << kr << " en pos " << pos << "\n";
				arr_kr_s.push_back(kr);
			}
			
			++cur_pos;
		}
		
		reader.close();
	}
	else{
		arr_kr_s.push_back(0);
		for(unsigned int i = 1; i < factors_start.size(); ++i){
			unsigned int word_len = factors_start[i] - factors_start[i-1];
			unsigned long long kr1 = arr_kr_s.back();
			unsigned long long kr2 = karp_rabin->hash(full_text + factors_start[i-1], word_len);
			unsigned long long kr_total = karp_rabin->concat(kr1, kr2, word_len);
			arr_kr_s.push_back(kr_total);
		}
		
		// Si no los cargo, los GUARDO en el archivo
		fstream writer(kr_frases_file, fstream::trunc | fstream::out);
		for(unsigned int i = 0; i < factors_start.size(); ++i){
			sprintf(buff, "%d\t%llu\n", i, arr_kr_s[i]);
			writer << buff;
		}
		writer.close();
		
	}
	
	
	cout << "FactorsIndexV3 - End\n";
	
}

FactorsIndexV3::~FactorsIndexV3(){
	
}

void FactorsIndexV3::find(const string &pattern, vector<unsigned int> &results){

//	cout << "FactorsIndexV3::find - Start\n";
	
//	cout << "FactorsIndexV3::find - Section A, reference\n";
	
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
	
//	cout << "FactorsIndexV3::find - Section B, ranges\n";
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
//			cout << "FactorsIndexV3::find - Invalid ranges, omitting...\n";
			continue;
		}
		
//		cout << "FactorsIndexV3::find - Searching in [" << r1.first << ", " << r1.second << "] x [" << r2.first << ", " << r2.second << "]:\n";
		auto res = wt.range_search_2d(r1.first, r1.second, r2.first, r2.second);
		for (auto point : res.second){
			unsigned int f = perm_y[point.second];
			unsigned int cur_perm = perm_inv[f];
			unsigned int pu = select1_b(perm[cur_perm] + 1);
//			cout << " -> Adding " << (pu - p1.length()) << "\n";
			results.push_back(pu - p1.length());
		}
		
	}
//	cout << "FactorsIndexV3::find - End\n";
	
}

void FactorsIndexV3::recursive_rmq(unsigned int ini, unsigned int fin, unsigned int min_pos, unsigned int occ_ref, vector<unsigned int> &results){
//	cout << "FactorsIndexV3::recursive_rmq - " << ini << ", " << fin << "\n";
	
	unsigned int pos_max = rmq(ini, fin);
	
	unsigned int tu = select1_s(pos_max + 1) - pos_max;
	unsigned int pu = select1_b(perm[pos_max] + 1);
	unsigned int lu = select1_b(perm[pos_max] + 2) - pu;
	
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

char FactorsIndexV3::getChar(unsigned int factor, unsigned int pos){
	
	// Iterators cache
	if( mapa_iterators.find(factor) == mapa_iterators.end() ){
		mapa_iterators[factor] = FactorsIterator(factor, n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref_text, &fm_index, len_text);
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

char FactorsIndexV3::getCharRev(unsigned int factor, unsigned int pos){
	
	if( factor == (unsigned int)(-1) ){
		return 0;
	}
	
	// Iterators cache
	if( mapa_iterators_rev.find(factor) == mapa_iterators_rev.end() ){
		mapa_iterators_rev[factor] = FactorsIteratorReverse(factor, n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref_text, &fm_index, len_text);
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

// Notar que, a diferencia de la busqueda en referencia, esta debe ser completa
// Es decir, solo importa el rango que contiene al patron completo
pair<unsigned int, unsigned int> FactorsIndexV3::getRangeY(const char *pattern){
	
	unsigned int pat_len = strlen(pattern);
	unsigned int izq = 0;
	unsigned int der = n_factors-1;
	
//	cout << "getRangeY - Inicio (pat_len: " << pat_len << ", izq: " << izq << ", der: " << der << ")\n";
	
	for( unsigned int cur_pos = 0; cur_pos < pat_len; ++cur_pos ){
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
			c = getChar(fm, cur_pos);
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
		c = getChar(fm, cur_pos);
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
		if( (cur_pos < text_len) && (unsigned char)(c) > (unsigned char)(pattern[cur_pos]) ){
			--der;
		}
//		cout << "getRangeY - der: " << der << "\n";
//		cout << "getRangeY - -----\n";
		
	}
	
//	cout << "getRangeY - result: (" << izq << ", " << der << ")\n";
	return pair<unsigned int, unsigned int>(izq, der);
}

pair<unsigned int, unsigned int> FactorsIndexV3::getRangeX(const char *pattern){
	
	unsigned int pat_len = strlen(pattern);
	unsigned int izq = 0;
	unsigned int der = n_factors-1;
	
//	cout << "getRangeX - Inicio (pat_len: " << pat_len << ", izq: " << izq << ", der: " << der << ")\n";
	
	for( unsigned int cur_pos = 0; cur_pos < pat_len; ++cur_pos ){
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
			fm = perm_x[m];
			c = getCharRev(fm-1, cur_pos);
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
		izq = h;
		fm = perm_x[izq];
		c = getCharRev(fm-1, cur_pos);
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
		if( (cur_pos < text_len) && (unsigned char)(c) > (unsigned char)(pattern[cur_pos]) ){
			--der;
		}
//		cout << "getRangeX - der: " << der << "\n";
//		cout << "getRangeX - -----\n";
		
	}
	
//	cout << "getRangeX - result: (" << izq << ", " << der << ")\n";
	return pair<unsigned int, unsigned int>(izq, der);
}
