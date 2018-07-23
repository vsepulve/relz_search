#include "FactorsIndex.h"

FactorsIndex::FactorsIndex(){
	len_text = 0;
	ref = NULL;
	len_ref = 0;
	n_factors = 0;
}

FactorsIndex::FactorsIndex(vector<pair<unsigned int, unsigned int> > &factors, unsigned int _len_text, const char *_ref, unsigned int _len_ref){
	
	len_text = _len_text;
	ref = _ref;
	len_ref = _len_ref;
	n_factors = factors.size();
	
	cout << "FactorsIndex - Inicio (factors: " << factors.size() << ", len_text: " << len_text << ", len_ref: " << len_ref << ")\n";
	
//	cout << "Factors: \n";
	// Factores en version ini, fin (absoluto) y ordenados por ini
	vector<pair<unsigned int, pair<unsigned int, unsigned int> > > factors_sort;
	unsigned cur_pos = 0;
	for( pair<unsigned int, unsigned int> factor : factors ){
//		cout << "(" << factor.first << ", " << factor.second << ", " << cur_pos << ")\n";
		factors_sort.push_back( 
			pair<unsigned int, pair<unsigned int, unsigned int> >(
				factor.first, pair<unsigned int, unsigned int>(factor.first + factor.second - 1, cur_pos++)
				)
			);
	}
	sort(factors_sort.begin(), factors_sort.end());
//	cout << "Factors Sorted: \n";
//	for( pair<unsigned int, pair<unsigned int, unsigned int> > factor : factors_sort ){
//		cout << "(" << factor.first << ", " << factor.second.first << ", " << factor.second.second << ")\n";
//	}
	
	// Bit vector S
	bit_vector _arr_s(len_ref + n_factors, 0);
	arr_s = _arr_s;
	unsigned cur_ref = 0;
	cur_pos = 0;
	for( unsigned int i = 0; i < n_factors; ++i ){
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
//	cout << "Bit Array S: \n";
//	for( unsigned int i = 0; i < len_ref + n_factors; ++i ){
//		cout << "arr_s[" << i << "]: " << arr_s[i] << "\n";
//	}
	
	rrr_vector<127> _rrr_s(arr_s);
	rrr_s = _rrr_s;
	rrr_vector<127>::select_1_type _select1_s(&rrr_s);
	select1_s = _select1_s;
//	cout << "Posicion de primer 1: " << select1_s(1) << "\n";
//	cout << "Posicion de tercer 1: " << select1_s(3) << "\n";
//	cout << "Posicion de quinto 1: " << select1_s(5) << "\n";
	
	// Notar que la posicion del select DEBE empezar desde 1, no desde 0
	// De este modo, hay que sumar 1 a las posiciones de la ref para buscar en S
	rrr_vector<127>::select_0_type _select0_s(&rrr_s);
	select0_s = _select0_s;
//	cout << "Posicion de primer 0: " << select0_s(1) << "\n";
//	cout << "Posicion de tercer 0: " << select0_s(3) << "\n";
//	cout << "Posicion de quinto 0: " << select0_s(5) << "\n";
//	cout << "Posicion de 0th 0: " << select0_s(0) << "\n";
	
	// Permutacion 
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
	
//	cout << "Permutation: \n";
//	for( unsigned int i = 0; i < n_factors; ++i ){
//		cout << "pi[" << i << "]: " << pi[i] << ", perm[" << i << "]: " << perm[i] << ", perm_inv[" << i << "]: " << perm_inv[i] << "\n";
//	}
	
	// Posiciones finales Ez
//	vector<unsigned int> ez(8);
	int_vector<> _ez(n_factors);
	ez = _ez;
	for( unsigned int i = 0; i < n_factors; ++i ){
		ez[i] = factors_sort[i].second.first;
	}
	// rmq_succinct_sct<> rmq(&ez);
	rmq_succinct_sct<false, bp_support_sada<256,32,rank_support_v5<> > > _rmq(&ez);
	rmq = _rmq;
	// rmq_maximum_sct<> rmq(&ez);
	
	cout << "Probando RMQ: \n";
	unsigned int r_pos = rmq(0, 5);
	cout << "RMQ(0, 5): " << r_pos << "\n";
	
	// Bit vector B (inicio de las frases en texto)
	bit_vector _arr_b(len_text, 0);
	arr_b = _arr_b;
	unsigned int pos_text = 0;
	for( unsigned int i = 0; i < n_factors; ++i ){
//		unsigned int ini = factors[i].first;
		unsigned int len = factors[i].second;
		arr_b[ pos_text ] = 1;
		pos_text += len;
	}
//	cout << "Bit Vector B: \n";
//	for( unsigned int i = 0; i < len_text; ++i ){
//		cout << "B[" << i << "]: " << arr_b[i] << "\n";
//	}
	rrr_vector<127> _rrr_b(arr_b);
	rrr_b = _rrr_b;
	rrr_vector<127>::select_1_type _select1_b(&rrr_b);
	rrr_vector<127>::select_0_type _select0_b(&rrr_b);
	select1_b = _select1_b;
	select0_b = _select0_b;
	
	// Construccion del FM Index (aunque use el SA original para la compresion, esto es para la busqueda)
//	csa_wt<> fm_index;
	// Construccion con datos en memoria, en un string
	construct_im(fm_index, ref, 1);
	// Construccion con datos en un archivo
//	construct(fm_index, file, 1);
	
	// Preparacion de permutaciones X e Y
	
	cout << "Preparando arr X\n";
	vector<unsigned int> arr_x(n_factors);
	for( unsigned int i = 0; i < n_factors; ++i ){
		arr_x[i] = i;
	}
	FactorsIteratorReverseComparator comp_rev(n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	stable_sort(arr_x.begin(), arr_x.end(), comp_rev);
	int_vector<> _pre_x_inv(n_factors);
	pre_x_inv = _pre_x_inv;
	for( unsigned int i = 0; i < n_factors; ++i ){
		pre_x_inv[ arr_x[i] ] = i;
	}
	inv_perm_support<> _perm_x(&pre_x_inv);
	perm_x = _perm_x;
	for( unsigned int i = 0; i < n_factors; ++i ){
		pre_x_inv[ arr_x[i] ] = i;
//		cout << " arr_x[" << i << "]: " << arr_x[i] << " (perm_x[" << i << "]: " << perm_x[i] << ") \n";
//		cout << " arr_x[" << i << "]: " << arr_x[i] << " -> ";
//		char c = 0;
//		for(unsigned int k = 0; k < 10 && (c = getCharRev(arr_x[i] - 1, k)) != 0; ++k ) 
//			cout << c;
//		cout << "\n";
	}
//	cout << "-----\n";
	
	cout << "Preparando arr Y\n";
	vector<unsigned int> arr_y(n_factors);
	for( unsigned int i = 0; i < n_factors; ++i ){
		arr_y[i] = i;
	}
	FactorsIteratorComparator comp(n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
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
//		cout << " arr_y[" << i << "]: " << arr_y[i] << " (perm_y[" << i << "]: " << perm_y[i] << ", perm_y_inv[" << i << "]: " << perm_y_inv[i] << ")\n";
//		cout << " arr_y[" << i << "]: " << perm_y[i] << " -> ";
//		char c = 0;
//		for(unsigned int k = 0; k < 10 && (c = getChar(perm_y[i], k)) != 0; ++k ) 
//			cout << c;
//		cout << "\n";
	}
//	cout << "-----\n";
	
	cout << "Preparando WT\n";
	int_vector<> _values_wt(n_factors);
	values_wt = _values_wt;
	for( unsigned int i = 0; i < n_factors; ++i ){
		values_wt[i] = perm_y_inv[ arr_x[ i ] ];
//		cout << " values_wt[" << i << "]: " << values_wt[i] << " \n";
	}
	
//	wt_int<rrr_vector<63>> wt;
	construct_im(wt, values_wt);
	
	
	
	cout << "FactorsIndex - End\n";
	
}

FactorsIndex::~FactorsIndex(){
	
}

void FactorsIndex::find(const string &pattern){
	
	cout << "FactorsIndex::find - Start\n";
	
	// Primera parte: Busqueda en referencia (y refinamiento con recursive RMQ)
	cout << "FactorsIndex::find - Section A, reference\n";
	
	size_t m = pattern.size();
	size_t occs = sdsl::count(fm_index, pattern.begin(), pattern.end());
	cout << "FactorsIndex::find - # occs de \"" << pattern << "\": " << occs << "\n";
	if( occs > 0 ){
		auto locations = locate(fm_index, pattern.begin(), pattern.begin()+m);
		sort(locations.begin(), locations.end());
		for( unsigned int i = 0; i < occs; ++i ){
			unsigned int occ_i = locations[i];
			// cout << "occ[" << i << "]: " << occ_i << " (" << ref.substr(occ_i, m) << ")\n";
			cout << "FactorsIndex::find - occ[" << i << "]: " << occ_i << " \n";
			// Comprobar los factores que cuben esta ocurrencia (el string ref[occ_i, occ_i + m - 1])
			unsigned int select = select0_s(occ_i + 1);
			unsigned int pos_ez = select - 1 - occ_i;
			cout << "FactorsIndex::find - select: " << select << " => pos_ez: " << pos_ez << "\n";
			
			// Ahora la busqueda (recursiva) en el rmq (entre 0 y pos_ez)
			recursive_rmq(0, pos_ez, (occ_i + m), occ_i);
			
		}
	}
	
	// Segunda parte: Busqueda de cada par de pattern en rangos X e Y
	cout << "FactorsIndex::find - Section B, ranges\n";
	
	for(unsigned int i = 1; i < pattern.length(); ++i){
		string p1 = pattern.substr(0, i);
		string p1_rev = "";
		for( unsigned int k = 0; k < p1.length(); ++k ){
			p1_rev += p1[ p1.length() - 1 - k ];
		}
		string p2 = pattern.substr(i, pattern.length() - i);
		cout << "FactorsIndex::find - Corte de \"" << pattern << "\": (" << p1 << " -> " << p1_rev << " | " << p2 << ")\n";
		pair<unsigned int, unsigned int> r1 = getRangeX(p1_rev.c_str());
		pair<unsigned int, unsigned int> r2 = getRangeY(p2.c_str());
		
		if( r1.second == (unsigned int)(-1) || r1.second < r1.first
			|| r2.second == (unsigned int)(-1) || r2.second < r2.first ){
			cout << "FactorsIndex::find - Rangos Invalidos, omitiendo...\n";
			continue;
		}
		
		cout << "FactorsIndex::find - Buscando en [" << r1.first << ", " << r1.second << "] x [" << r2.first << ", " << r2.second << "]:\n";
		auto res = wt.range_search_2d(r1.first, r1.second, r2.first, r2.second);
		for (auto point : res.second){
			unsigned int f = perm_y[point.second];
			cout << "FactorsIndex::find - (" << point.first << ", " << point.second << ") => factor " << f << "\n";
			
			unsigned int cur_perm = perm_inv[f];
//			unsigned int tu = select1_s(cur_perm + 1) - cur_perm;
			unsigned int pu = select1_b(perm[cur_perm] + 1);
//			unsigned int lu = select1_b(perm[cur_perm] + 2) - pu;
//			cout << " -> tu: " << tu << ", pu: " << pu << ", lu: " << lu << "\n";
			cout << " -> Agregando " << (pu - p1.length()) << "\n";
			
		}
		
	}
	cout << "FactorsIndex::find - End\n";
	
}


void FactorsIndex::test(const string &pattern){
	
	cout << "Texto de ref: \"" << ref << "\"\n";
	
	cout << "Probando RMQ: \n";
	unsigned int r_pos = rmq(0, 5);
	cout << "RMQ(0, 5): " << r_pos << "\n";
	cout << "Probando perm[4]: " << perm[4] << "\n";
	
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
			cout << "1\n";
			unsigned int select = select0_s(occ_i + 1);
			cout << "2\n";
			unsigned int pos_ez = select - 1 - occ_i;
			cout << "3\n";
			cout << "select: " << select << " => pos_ez: " << pos_ez << "\n";
			
			// Ahora la busqueda (recursiva) en el rmq (entre 0 y pos_ez)
//			unsigned int pos_max = rmq(0, pos_ez);
//			cout << "max pos Ez: " << pos_max << " (Ez: " << ez[pos_max] << ", factor: " << perm[pos_max] << ")\n";
			cout << "----- Search V2 -----\n";
			recursive_rmq(0, pos_ez, (occ_i + m), occ_i);
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
	
	FactorsIterator it( 2, n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	cout << "-----\n";
	
	while( it.hasNext() ){
		char c = it.next();
		cout << "it.next(): " << c << " (text_pos " << it.position() << " / " << it.length() << ")\n";
		cout << "-----\n";
	}
	
	cout << "Fin prueba iterador\n";
	cout << "-----\n";
	
	cout << "Probando Iterador Reverso\n";
	
	FactorsIteratorReverse it_rev( 2, n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	cout << "-----\n";
	
	while( it_rev.hasNext() ){
		char c = it_rev.next();
		cout << "it_rev.next(): " << c << " (text_pos " << it_rev.position() << " / " << it_rev.length() << ")\n";
		cout << "-----\n";
	}
	
	cout << "Fin prueba iterador\n";
	cout << "-----\n";

	
	cout << "Probando Iterador Reverso -1\n";
	it_rev = FactorsIteratorReverse( -1, n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	while( it_rev.hasNext() ){
		char c = it_rev.next();
		cout << "it_rev.next(): " << c << " (text_pos " << it_rev.position() << " / " << it_rev.length() << ")\n";
		cout << "-----\n";
	}

	cout << "Probando Iterador Reverso 0\n";
	it_rev = FactorsIteratorReverse( 0, n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	while( it_rev.hasNext() ){
		char c = it_rev.next();
		cout << "it_rev.next(): " << c << " (text_pos " << it_rev.position() << " / " << it_rev.length() << ")\n";
		cout << "-----\n";
	}

	cout << "Probando Iterador Reverso 1\n";
	it_rev = FactorsIteratorReverse( 1, n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	while( it_rev.hasNext() ){
		char c = it_rev.next();
		cout << "it_rev.next(): " << c << " (text_pos " << it_rev.position() << " / " << it_rev.length() << ")\n";
		cout << "-----\n";
	}

	cout << "Probando Iterador Reverso 5\n";
	it_rev = FactorsIteratorReverse( 5, n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	while( it_rev.hasNext() ){
		char c = it_rev.next();
		cout << "it_rev.next(): " << c << " (text_pos " << it_rev.position() << " / " << it_rev.length() << ")\n";
		cout << "-----\n";
	}

	cout << "Probando Iterador Reverso 7\n";
	it_rev = FactorsIteratorReverse( 7, n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	while( it_rev.hasNext() ){
		char c = it_rev.next();
		cout << "it_rev.next(): " << c << " (text_pos " << it_rev.position() << " / " << it_rev.length() << ")\n";
		cout << "-----\n";
	}
	
	cout << "-----\n";
		
	cout << "Probando Iterador Reverso 8\n";
	it_rev = FactorsIteratorReverse( 8, n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
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
	it_rev = FactorsIteratorReverse( 5, n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
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
		char c = getChar( f, i);
		if( c == 0){
			break;
		}
		cout << "factor_" << f << "[" << i << "]: " << c << " \n";
	}
	cout << "-----\n";
	f = 7;
	for( unsigned int i = 0; i < 10; ++i ){
		char c = getChar( f, i);
		if( c == 0){
			break;
		}
		cout << "factor_" << f << "[" << i << "]: " << c << " \n";
	}
	cout << "-----\n";
	f = 3;
	for( unsigned int i = 0; i < 10; ++i ){
		char c = getChar( f, i);
		if( c == 0){
			break;
		}
		cout << "factor_" << f << "[" << i << "]: " << c << " \n";
	}
	cout << "-----\n";
	
	cout << "Probando acceso posicional Reverso\n";
	f = 3;
	for( unsigned int i = 0; i < 10; ++i ){
		char c = getCharRev( f-1, i);
		if( c == 0){
			break;
		}
		cout << "factor_rev_" << f << "[" << i << "]: " << c << " \n";
	}
	cout << "-----\n";
	f = 5;
	for( unsigned int i = 0; i < 10; ++i ){
		char c = getCharRev( f-1, i);
		if( c == 0){
			break;
		}
		cout << "factor_rev_" << f << "[" << i << "]: " << c << " \n";
	}
	cout << "-----\n";
	cout << "factor_rev_" << f << "[" << 3 << "]: " << getCharRev( f-1, 3) << " \n";
	cout << "factor_rev_" << f << "[" << 5 << "]: " << getCharRev( f-1, 5) << " \n";
	cout << "factor_rev_" << f << "[" << 2 << "]: " << getCharRev( f-1, 2) << " \n";
	cout << "factor_rev_" << f << "[" << 0 << "]: " << getCharRev( f-1, 0) << " \n";
	cout << "-----\n";
	
	
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
	getRangeY("A");
	cout << "-----\n";
	
	cout << "Prueba de patron \"B\"\n";
	getRangeY("B");
	cout << "-----\n";
	
	cout << "Prueba de patron \"BA\"\n";
	getRangeY("BA");
	cout << "-----\n";
	
	cout << "Prueba de patron \"BAL\"\n";
	getRangeY("BAL");
	cout << "-----\n";
	
	cout << "Prueba de patron \"BALB\"\n";
	getRangeY("BALB");
	cout << "-----\n";
	
	cout << "Prueba de patron \"Z\"\n";
	getRangeY("Z");
	cout << "-----\n";
	
	cout << "Prueba de patron \"0\"\n";
	getRangeY("0");
	cout << "-----\n";
	
	cout << "Prueba de patron \"\"\n";
	getRangeY("");
	cout << "-----\n";
	
	
//	string pattern = "AB";
	for(unsigned int i = 1; i < pattern.length(); ++i){
		string p1 = pattern.substr(0, i);
		string p2 = pattern.substr(i, pattern.length() - i);
		cout << "Corte de \"" << pattern << "\": (" << p1 << "| " << p2 << ")\n";
		pair<unsigned int, unsigned int> r1 = getRangeX(p1.c_str());
		pair<unsigned int, unsigned int> r2 = getRangeY(p2.c_str());
		
		if( r1.second == (unsigned int)(-1) || r1.second < r1.first
			|| r2.second == (unsigned int)(-1) || r2.second < r1.first ){
			cout << "Rangos Invalidos, omitiendo...\n";
			continue;
		}
		
		cout << "Buscando en [" << r1.first << ", " << r1.second << "] x [" << r2.first << ", " << r2.second << "]:\n";
		auto res = wt.range_search_2d(r1.first, r1.second, r2.first, r2.second);
		for (auto point : res.second){
			unsigned int f = perm_y[point.second];
			cout << "(" << point.first << ", " << point.second << ") => factor " << f << "\n";
		}
		
	}
	
	
}

void FactorsIndex::recursive_rmq(unsigned int ini, unsigned int fin, unsigned int crit, unsigned int occ_ref){
	cout << "FactorsIndex::recursive_rmq - " << ini << ", " << fin << "\n";
	
	unsigned int pos_max = rmq(ini, fin);
	
	unsigned int tu = select1_s(pos_max + 1) - pos_max;
	unsigned int pu = select1_b(perm[pos_max] + 1);
	unsigned int lu = select1_b(perm[pos_max] + 2) - pu;
	
	cout << "FactorsIndex::recursive_rmq - max pos Ez: " << pos_max << " (tu: " << tu << ", pu: " << pu << ", lu: " << lu << ")\n";
	if( tu + lu < crit ){
		cout << "Omitiendo\n";
		return;
	}
	else{
		cout << " -> Agregando " << (pu + (occ_ref - tu)) << "\n";
	}
	
	if( (pos_max > 0) && (ini < pos_max) ){
		recursive_rmq(ini, pos_max-1, crit, occ_ref);
	}
	if( pos_max < fin ){
		recursive_rmq(pos_max+1, fin, crit, occ_ref);
	}
}

char FactorsIndex::getChar(unsigned int factor, unsigned int pos){
	
	// Segundo enfoque: cache de iteradores completos
	if( mapa_iterators.find(factor) == mapa_iterators.end() ){
		mapa_iterators[factor] = FactorsIterator(factor, n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
	}
	FactorsIterator it = mapa_iterators[factor];
	if( pos >= it.length() ){
		return 0;
	}
	if( it.position() > pos ){
		it.reset();
	}
	
	// Primer enfoque: sin caches 
	// (porque quiza haga un multiiterador con caches internos solo para los valores especificos necesarios)
//	FactorsIterator it( factor, n_factors, select1_s, select1_b, select0_b, perm, perm_inv, ref, len_text );
	char c = 0;
	while( it.position() <= pos ){
		c = it.next();
	}
	return c;
}

char FactorsIndex::getCharRev(unsigned int factor, unsigned int pos){
	
	if( factor == (unsigned int)(-1) ){
		return 0;
	}
	
	// Segundo enfoque: cache de iteradores completos
	if( mapa_iterators_rev.find(factor) == mapa_iterators_rev.end() ){
		mapa_iterators_rev[factor] = FactorsIteratorReverse(factor, n_factors, &select1_s, &select1_b, &select0_b, &perm, &perm_inv, ref, len_text);
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
pair<unsigned int, unsigned int> FactorsIndex::getRangeY(const char *pattern){
	
	unsigned int pat_len = strlen(pattern);
	unsigned int izq = 0;
	unsigned int der = n_factors-1;
	
	cout << "getRangeY - Inicio (pat_len: " << pat_len << ", izq: " << izq << ", der: " << der << ")\n";
	
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
//			if( arr[m] + cur_pos > largo || *(ref + arr[m] + cur_pos) < (unsigned char)(*(text + cur_pos)) ){
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
//			if( arr[m] + cur_pos > largo || *(ref + arr[m] + cur_pos) <= (unsigned char)(*(text + cur_pos)) ){
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
	
	cout << "getRangeY - result: (" << izq << ", " << der << ")\n";
	return pair<unsigned int, unsigned int>(izq, der);
}

pair<unsigned int, unsigned int> FactorsIndex::getRangeX(const char *pattern){
	
	unsigned int pat_len = strlen(pattern);
	unsigned int izq = 0;
	unsigned int der = n_factors-1;
	
	cout << "getRangeX - Inicio (pat_len: " << pat_len << ", izq: " << izq << ", der: " << der << ")\n";
	
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
//			if( arr[m] + cur_pos > largo || *(ref + arr[m] + cur_pos) < (unsigned char)(*(text + cur_pos)) ){
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
//			if( arr[m] + cur_pos > largo || *(ref + arr[m] + cur_pos) <= (unsigned char)(*(text + cur_pos)) ){
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
	
	cout << "getRangeX - result: (" << izq << ", " << der << ")\n";
	return pair<unsigned int, unsigned int>(izq, der);
}

