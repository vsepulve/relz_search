#include "KarpRabinFactorsSuffixesv2.h"

KarpRabinFactorsSuffixesv2::KarpRabinFactorsSuffixesv2(){
	n_factors = 0;
	arr_kr_s = NULL;
	delete_krs = false;
	karp_rabin = NULL;
	ref_text = NULL;
	select1_s = NULL;
	select1_b = NULL;
	select0_b = NULL;
	pi_inv = NULL;
	factors_start = NULL;
	delete_start = false;
	
	max_offset = 0;
	max_length = 0;
	nano1 = 0;
	nano2 = 0;
	nano3 = 0;
	nano4 = 0;
}

KarpRabinFactorsSuffixesv2::KarpRabinFactorsSuffixesv2(unsigned int _n_factors, 
		vector<unsigned long long> *_arr_kr_s, 
		KarpRabin *_karp_rabin, 
		const char *_ref_text, 
		bits_s_type::select_1_type *_select1_s, 
		bits_b_type::select_1_type *_select1_b, 
		bits_b_type::select_0_type *_select0_b, 
		int_vector<> *_pi_inv, 
		vector<unsigned int> *_factors_start){
	n_factors = _n_factors;
	arr_kr_s = _arr_kr_s;
	delete_krs = false;
	karp_rabin = _karp_rabin;
	if( arr_kr_s->size() < n_factors+1 ){
		cout << "KarpRabinFactorsSuffixesv2 - Warning (insufficient prefixes, must include one for the whole collection)\n";
	}
	ref_text = _ref_text;
	select1_s = _select1_s;
	select1_b = _select1_b;
	select0_b = _select0_b;
	pi_inv = _pi_inv;
	factors_start = _factors_start;
	delete_start = false;
	
	max_offset = 0;
	max_length = 0;
	nano1 = 0;
	nano2 = 0;
	nano3 = 0;
	nano4 = 0;
}

KarpRabinFactorsSuffixesv2::KarpRabinFactorsSuffixesv2(const string &file, 
		KarpRabin *_karp_rabin, 
		const char *_ref_text, 
		bits_s_type::select_1_type *_select1_s, 
		bits_b_type::select_1_type *_select1_b, 
		bits_b_type::select_0_type *_select0_b, 
		int_vector<> *_pi_inv){
	// Basic variables loaded from file
	n_factors = 0;
	arr_kr_s = 0;
	delete_krs = false;
	factors_start = NULL;
	delete_start = false;
	load(file);
	// External variables
	karp_rabin = _karp_rabin;
	ref_text = _ref_text;
	select1_s = _select1_s;
	select1_b = _select1_b;
	select0_b = _select0_b;
	pi_inv = _pi_inv;
	
	max_offset = 0;
	max_length = 0;
	nano1 = 0;
	nano2 = 0;
	nano3 = 0;
	nano4 = 0;
}

KarpRabinFactorsSuffixesv2::~KarpRabinFactorsSuffixesv2(){
	n_factors = 0;
	if( delete_krs && arr_kr_s != NULL ){
		arr_kr_s->clear();
		delete arr_kr_s;
	}
	arr_kr_s = NULL;
	if( delete_start && factors_start != NULL){
		factors_start->clear();
		delete factors_start;
	}
	factors_start = NULL;
	karp_rabin = NULL;
	ref_text = NULL;
}

/*
unsigned long long KarpRabinFactorsSuffixesv2::hash(unsigned int factor_ini, unsigned int length){
//	cout << "KarpRabinFactorsSuffixesv2::hash - Start (factor_ini: " << factor_ini << ", length: " << length << ")\n";
	unsigned int cur_len = 0;
	unsigned int factor_cur = factor_ini;
	unsigned long long kr1 = 0;
	unsigned int cur_pi = 0;
	unsigned int tu = 0;
	unsigned int pu = 0;
	unsigned int lu = 0;
	while( cur_len < length ){
		cout << "a\n";
		cur_pi = (*pi_inv)[factor_cur];
		tu = select1_s->operator()(cur_pi + 1) - cur_pi;
		pu = select1_b->operator()(factor_cur + 1);
		lu = select1_b->operator()(factor_cur + 2) - pu;
		cout << "b (tu: " << tu << ", pu: " << pu << ", lu: " << lu << ")\n";
		if( cur_len + lu >= length ){
			break;
		}
		// Agregar factor completo y continuar
		cur_len += lu;
		++factor_cur;
//		cout << "KarpRabinFactorsSuffixesv2::hash - karp_rabin->subtract_prefix (lu: " << lu << ")\n";
		unsigned long long kr_phrase = karp_rabin->subtract_prefix(arr_kr_s->at(factor_cur), arr_kr_s->at(factor_cur-1), lu);
//		cout << "KarpRabinFactorsSuffixesv2::hash - Adding to current hash\n";
		kr1 = karp_rabin->concat(kr1, kr_phrase, lu);
//		cout << "KarpRabinFactorsSuffixesv2::hash - cur kr1: " << kr1 << "\n";
		
	}
	// Calcular el hash de los (length - cur_len) primeros caracteres del factor actual
	// Notar que lo que sigue se puede extraer er KarpRabinReference (desde tu, de largo length - cur_len)
	unsigned long long kr2 = karp_rabin->hash(ref_text + tu, length - cur_len);
//	cout << "KarpRabinFactorsSuffixesv2::hash - kr2: " << kr2 << ", len: " << length - cur_len << "\n";
	
	unsigned long long kr12 = karp_rabin->concat(kr1, kr2, length - cur_len);
//	cout << "KarpRabinFactorsSuffixesv2::hash - End (kr12: " << kr12 << ")\n";
	
	return kr12;
	
}
*/

// Fast version, only valid for offset and length < karp_rabin->getTableSize()
unsigned long long KarpRabinFactorsSuffixesv2::hash(unsigned int factor_ini, unsigned int length){
	cout << "KarpRabinFactorsSuffixesv2::hash - Start (factor_ini: " << factor_ini << ", length: " << length << ", n_factors: " << n_factors << ", n_pi_inv: " << pi_inv->size() << ")\n";
	unsigned long long kr1 = 0;
	
	unsigned int cur_len = 0;
	unsigned int factor_cur = 0;
	unsigned int tu = 0;
	
	// Version con Busqueda Binaria
	unsigned int l = factor_ini;
	// TODO : verificar n_factors - 1
	unsigned int h = n_factors;
	unsigned int factor_cur_bin = factor_ini;
	unsigned int pu_ini = select1_b->operator()(factor_ini + 1);
	
	unsigned int cur_perm_bin = (*pi_inv)[factor_cur_bin];
	unsigned int tu_bin = select1_s->operator()(cur_perm_bin + 1) - cur_perm_bin;
	unsigned int pu_bin = select1_b->operator()(factor_cur_bin + 1);
//	unsigned int lu_bin = select1_b->operator()(factor_cur_bin + 2) - pu_bin;
	
//	cout << "1\n";
	while(l < h){
		factor_cur_bin = l + ((h-l)>>1);
		pu_bin = select1_b->operator()(factor_cur_bin + 1);
		cout << "factor_cur_bin: " << factor_cur_bin << ", factor_ini: " << factor_ini << ", len: " << (pu_bin - pu_ini) << "\n";
		if( (pu_bin - pu_ini) < length ){
			l = factor_cur_bin+1;
		}
		else{
			h = factor_cur_bin;
		}
	}
//	cout << "2 (factor_cur_bin: " << factor_cur_bin << ", h: " << h << ")\n";
//	factor_cur_bin = h;
	if( (factor_cur_bin > factor_ini) && (pu_bin - pu_ini) >= length ){
		--factor_cur_bin;
	}
//	cout << "3 (factor_cur_bin: " << factor_cur_bin << ")\n";
	unsigned int cur_len_bin = (pu_bin - pu_ini);
	cur_perm_bin = (*pi_inv)[factor_cur_bin];
	tu_bin = select1_s->operator()(cur_perm_bin + 1) - cur_perm_bin;
	
//	if( (cur_len_bin != cur_len)
//		|| (factor_cur_bin != factor_cur)
//		|| (tu_bin != tu) ){
//		cout << "KarpRabinFactorsSuffixesv2::hash - Error (cur_len_bin: " << cur_len_bin << " / " << cur_len << ", factor_cur_bin: " << factor_cur_bin << " / " << factor_cur << ", tu_bin: " << tu_bin << " / " << tu << ")\n";
//		exit(0);
//	}
	
	cur_len = cur_len_bin;
	factor_cur = factor_cur_bin;
	tu = tu_bin;
	
	if( factor_cur > factor_ini ){
//		cout << "KarpRabinFactorsSuffixesv2::hash - Preparing kr1 subtract_prefix(" << arr_kr_s->at(factor_cur) << ", " << arr_kr_s->at(factor_ini) << ", " << cur_len << ")\n";
		kr1 = karp_rabin->subtract_prefix(arr_kr_s->at(factor_cur), arr_kr_s->at(factor_ini), cur_len);
	}
//	cout << "KarpRabinFactorsSuffixesv2::hash - kr1: " << kr1 << ", cur_len: " << cur_len << "\n";
	// Calcular el hash de los (length - cur_len) primeros caracteres del factor actual
	// Notar que lo que sigue se puede extraer er KarpRabinReference (desde tu, de largo length - cur_len)
	
	NanoTimer timer;
//	cout << "KarpRabinFactorsSuffixesv2::hashFast - karp_rabin->hash(ref_text + " << tu << ", " << length << " - " << cur_len << ")\n";
	unsigned long long kr2 = karp_rabin->hash(ref_text + tu, length - cur_len);
	kr_nano += timer.getNanosec();
	
//	cout << "KarpRabinFactorsSuffixesv2::hash - kr2: " << kr2 << ", len: " << length - cur_len << "\n";
	unsigned long long kr12 = karp_rabin->concat(kr1, kr2, length - cur_len);
//	cout << "KarpRabinFactorsSuffixesv2::hash - End (kr12: " << kr12 << ")\n";
	return kr12;
}

// Fast version, only valid for offset and length < karp_rabin->getTableSize()
unsigned long long KarpRabinFactorsSuffixesv2::hash(unsigned int factor_ini, unsigned int offset, unsigned int length){
//	cout << "KarpRabinFactorsSuffixesv2::hash - Start (factor_ini: " << factor_ini << ", offset: " << offset << ", length: " << length << ")\n";
	
	if( length == 0 ){
		return 0;
	}
	
	unsigned long long kr_total = hash(factor_ini, offset + length);
//	cout << "KarpRabinFactorsSuffixesv2::hash - kr_total: " << kr_total << "\n";
	unsigned long long kr_prefix = hash(factor_ini, offset);
//	cout << "KarpRabinFactorsSuffixesv2::hash - kr_prefix: " << kr_prefix << "\n";
	unsigned long long kr = karp_rabin->subtract_prefix(kr_total, kr_prefix, length);
//	cout << "KarpRabinFactorsSuffixesv2::hash - kr: " << kr << "\n";
	return kr;
}

/*
unsigned long long KarpRabinFactorsSuffixesv2::hash(unsigned int factor_ini, unsigned int offset, unsigned int length){
	
	// Solo para debug y estadisticas
	if( offset > max_offset ){
		max_offset = offset;
	}
	if( length > max_length ){
		max_length = length;
	}
	
	if( (offset + length) < karp_rabin->getTableSize() ){
//		cout << "KarpRabinFactorsSuffixesv2::hash - hashFast (" << factor_ini << ", " << offset << ", " << length << ")\n";
		return hashFast(factor_ini, offset, length);
	}

//	return hashBin(factor_ini, offset, length);
	
	if( length == 0 ){
		return 0;
	}
	
//	cout << "KarpRabinFactorsSuffixesv2::hash - Start (factor_ini: " << factor_ini << ", offset: " << offset << ", length: " << length << ")\n";
	
	NanoTimer timer;
	
	unsigned int factor_cur = factor_ini;
	
	// Primero omito los factores completos previos a offset
	unsigned int cur_perm = (*pi_inv)[factor_cur];
	unsigned int tu = select1_s->operator()(cur_perm + 1) - cur_perm;
	unsigned int pu = select1_b->operator()(factor_cur + 1);
	unsigned int lu = select1_b->operator()(factor_cur + 2) - pu;
	
	while( lu < offset ){
//		cout << "KarpRabinFactorsSuffixesv2::hash - Ommiting factor " << factor_cur << " (tu: " << tu << ", pu: " << pu << ", lu: " << lu << ")\n";
		offset -= lu;
		++factor_cur;
		
		cur_perm = (*pi_inv)[factor_cur];
		tu = select1_s->operator()(cur_perm + 1) - cur_perm;
		pu = select1_b->operator()(factor_cur + 1);
		lu = select1_b->operator()(factor_cur + 2) - pu;
		
	}
	
	if( lu == offset ){
//		cout << "KarpRabinFactorsSuffixesv2::hash - increasing\n";
		offset = 0;
		++factor_cur;
		
		cur_perm = (*pi_inv)[factor_cur];
		tu = select1_s->operator()(cur_perm + 1) - cur_perm;
		pu = select1_b->operator()(factor_cur + 1);
		lu = select1_b->operator()(factor_cur + 2) - pu;
		
	}
	
	nano1 += timer.getNanosec();
	timer.reset();
	
//	cout << "KarpRabinFactorsSuffixesv2::hash - Preparing to remove first prefix (factor_cur: " << factor_cur << ", offset: " << offset << ", length: " << length << ")\n";
	// Ahora estoy seguro que el offset se aplica a la frase actual
	// Tambien se que offset es menor que lu
	// Calculo el hash de offset y lo resto a la palabra actual
//	cout << "KarpRabinFactorsSuffixesv2::hash - hash 1 (ref_text + " << tu << ", " << offset << ")\n";
	unsigned long long kr_prefix = karp_rabin->hash(ref_text + tu, offset);
//	cout << "KarpRabinFactorsSuffixesv2::hash - kr_prefix: " << kr_prefix << "\n";
	
	// Puede pasar que el largo sea MENOR que el largo de la frase
	// En ese caso, la idea seria omitir el resto
	unsigned long long kr_phrase = 0;
	unsigned long long kr_total = 0;
	if( length + offset <= lu ){
//		cout << "KarpRabinFactorsSuffixesv2::hash - hash 2 (ref_text + " << tu << ", " << length << ")\n";
		kr_phrase = karp_rabin->hash(ref_text + tu, length);
//		cout << "KarpRabinFactorsSuffixesv2::hash - kr_phrase (FINAL) [" << factor_cur << "]: " << kr_phrase << " (len: " << (length + offset) << ", lu: " << lu << ")\n";
		kr_total = karp_rabin->subtract_prefix(kr_phrase, kr_prefix, length + offset);
//		cout << "KarpRabinFactorsSuffixesv2::hash - End (" << kr_total << ")\n";
		return kr_total;
	}
	else{
		kr_phrase = karp_rabin->subtract_prefix(arr_kr_s->at(factor_cur+1), arr_kr_s->at(factor_cur), lu);
//		cout << "KarpRabinFactorsSuffixesv2::hash - kr_phrase [" << factor_cur << "]: " << kr_phrase << "\n";
		kr_total = karp_rabin->subtract_prefix(kr_phrase, kr_prefix, lu - offset);
//		cout << "KarpRabinFactorsSuffixesv2::hash - kr_total: " << kr_total << "\n";
	}
	
	length -= (lu - offset);
	++factor_cur;
	
	nano2 += timer.getNanosec();
	timer.reset();
	
	cur_perm = (*pi_inv)[factor_cur];
	tu = select1_s->operator()(cur_perm + 1) - cur_perm;
	pu = select1_b->operator()(factor_cur + 1);
	lu = select1_b->operator()(factor_cur + 2) - pu;
	
	// Ahora agrego el hash de cada frase completa contenida
	unsigned int sum_len = 0;
	while( lu < length ){
//		cout << "KarpRabinFactorsSuffixesv2::hash - Adding factor " << factor_cur << " (lu: " << lu << ", length: " << length << ")\n";
		
		// Agrego el hash de la frase
		kr_phrase = karp_rabin->subtract_prefix(arr_kr_s->at(factor_cur+1), arr_kr_s->at(factor_cur), lu);
//		cout << "KarpRabinFactorsSuffixesv2::hash - kr_phrase [" << factor_cur << "]: " << kr_phrase << "\n";
		
		kr_total = karp_rabin->concat(kr_total, kr_phrase, lu);
//		cout << "KarpRabinFactorsSuffixesv2::hash - kr_total " << kr_total << "\n";
		
		length -= lu;
		++factor_cur;
		sum_len += lu;
		
		cur_perm = (*pi_inv)[factor_cur];
		tu = select1_s->operator()(cur_perm + 1) - cur_perm;
		pu = select1_b->operator()(factor_cur + 1);
		lu = select1_b->operator()(factor_cur + 2) - pu;
	}
	
//	cout << "KarpRabinFactorsSuffixesv2::hash - sum_len: " << sum_len << "\n";
	
	nano3 += timer.getNanosec();
	timer.reset();
	
//	cout << "KarpRabinFactorsSuffixesv2::hash - Preparing to add last prefix (length: " << length << ")\n";
	// Ahora se que el length que queda es solo un prefijo de la frase actual, de largo menor a lu
//	cout << "KarpRabinFactorsSuffixesv2::hash - hash 3 (ref_text + " << tu << ", " << length << ")\n";
	kr_prefix = karp_rabin->hash(ref_text + tu, length);
//	cout << "KarpRabinFactorsSuffixesv2::hash - kr_prefix: " << kr_prefix << "\n";
	kr_total = karp_rabin->concat(kr_total, kr_prefix, length);
	
	nano4 += timer.getNanosec();
	
//	cout << "KarpRabinFactorsSuffixesv2::hash - End (" << kr_total << ")\n";
	return kr_total;
	
}

unsigned long long KarpRabinFactorsSuffixesv2::hashBin(unsigned int factor_ini, unsigned int offset, unsigned int length){

	if( length == 0 ){
		return 0;
	}
	
	NanoTimer timer;
	
//	cout << "KarpRabinFactorsSuffixesv2::hashBin - Start (factor_ini: " << factor_ini << ", offset: " << offset << ", length: " << length << ")\n";
	
	unsigned int offset_ini = select1_b->operator()(factor_ini + 1);
	
	// Version con Busqueda Binaria
	unsigned int l = factor_ini;
	unsigned int h = pi_inv->size();
	
	unsigned int factor_cur = factor_ini;
	unsigned int cur_perm = (*pi_inv)[factor_cur];
	unsigned int tu = 0;
	unsigned int pu = select1_b->operator()(factor_cur + 1);
	unsigned int lu = 0;
	
	while(l < h){
		factor_cur = l + ((h-l)>>1);
		pu = select1_b->operator()(factor_cur + 1);
//		cout << "KarpRabinFactorsSuffixesv2::hashBin - factor_cur: " << factor_cur << ", pu: " << pu << " / " << (offset + offset_ini) << "\n";
		if( pu < offset + offset_ini ){
			l = factor_cur+1;
		}
		else{
			h = factor_cur;
		}
	}
//	cout << "KarpRabinFactorsSuffixesv2::hashBin - factor_cur: " << factor_cur << " (pu: " << pu << " / " << (offset + offset_ini) << ")\n";
	if( (factor_cur > factor_ini) && pu >= (offset + offset_ini) ){
//		cout << "KarpRabinFactorsSuffixesv2::hashBin - reducing\n";
		--factor_cur;
		pu = select1_b->operator()(factor_cur + 1);
	}
	
//	offset -= (pu - offset_ini);
	cur_perm = (*pi_inv)[factor_cur];
	tu = select1_s->operator()(cur_perm + 1) - cur_perm;
	lu = select1_b->operator()(factor_cur + 2) - pu;
	
	nano1 += timer.getNanosec();
	timer.reset();
	
//	cout << "KarpRabinFactorsSuffixesv2::hashBin - Preparing to remove first prefix (factor_cur: " << factor_cur << ", offset: " << offset << ", length: " << length << ")\n";
	// Ahora estoy seguro que el offset se aplica a la frase actual
	// Tambien se que offset es menor que lu
	// Calculo el hash de offset y lo resto a la palabra actual
//	cout << "KarpRabinFactorsSuffixesv2::hashBin - hash 1 (ref_text + " << tu << ", " << offset + offset_ini - pu << ")\n";
	unsigned long long kr_prefix = karp_rabin->hash(ref_text + tu, offset + offset_ini - pu);
//	cout << "KarpRabinFactorsSuffixesv2::hashBin - kr_prefix: " << kr_prefix << "\n";
	
	// Puede pasar que el largo sea MENOR que el largo de la frase
	// En ese caso, la idea seria omitir el resto
	unsigned long long kr_phrase = 0;
	unsigned long long kr_total = 0;
	if( (length + offset + offset_ini - pu) <= lu ){
//		cout << "KarpRabinFactorsSuffixesv2::hashBin - hash 2 (ref_text + " << tu << ", " << length << ")\n";
		kr_phrase = karp_rabin->hash(ref_text + tu, length);
//		cout << "KarpRabinFactorsSuffixesv2::hashBin - kr_phrase (FINAL) [" << factor_cur << "]: " << kr_phrase << " (len: " << (length + offset) << ", lu: " << lu << ")\n";
		kr_total = karp_rabin->subtract_prefix(kr_phrase, kr_prefix, length + offset + offset_ini - pu);
//		cout << "KarpRabinFactorsSuffixesv2::hashBin - End (" << kr_total << ")\n";
		return kr_total;
	}
	else{
		kr_phrase = karp_rabin->subtract_prefix(arr_kr_s->at(factor_cur+1), arr_kr_s->at(factor_cur), lu);
//		cout << "KarpRabinFactorsSuffixesv2::hashBin - kr_phrase [" << factor_cur << "]: " << kr_phrase << "\n";
		kr_total = karp_rabin->subtract_prefix(kr_phrase, kr_prefix, lu - (offset + offset_ini - pu));
//		cout << "KarpRabinFactorsSuffixesv2::hashBin - kr_total: " << kr_total << "\n";
	}
	
	++factor_cur;
	
	nano2 += timer.getNanosec();
	timer.reset();
	
	// Version con Busqueda Binaria
	l = factor_cur;
	h = pi_inv->size();
	
	factor_ini = factor_cur;
	cur_perm = (*pi_inv)[factor_cur];
	tu = 0;
	pu = select1_b->operator()(factor_cur + 1);
	lu = 0;
	
	while(l < h){
		factor_cur = l + ((h-l)>>1);
		pu = select1_b->operator()(factor_cur + 1);
//		cout << "KarpRabinFactorsSuffixesv2::hashBin - factor_cur: " << factor_cur << ", pu: " << pu << " / " << (offset + offset_ini + length) << "\n";
		if( pu < offset + offset_ini + length ){
			l = factor_cur+1;
		}
		else{
			h = factor_cur;
		}
	}
//	cout << "KarpRabinFactorsSuffixesv2::hashBin - factor_cur: " << factor_cur << " (pu: " << pu << " / " << (offset + offset_ini + length) << ")\n";
	if( pu >= (offset + offset_ini + length) ){
//		cout << "KarpRabinFactorsSuffixesv2::hashBin - reducing\n";
		--factor_cur;
		pu = select1_b->operator()(factor_cur + 1);
	}
	cur_perm = (*pi_inv)[factor_cur];
	tu = select1_s->operator()(cur_perm + 1) - cur_perm;
	lu = select1_b->operator()(factor_cur + 2) - pu;
	
	// Aqui tengo kr_total que es el kr del primer corte
	// Hay que agregar kr_phrase que seria el kr de todas las frases desde factor_ini hasta factor_cur
	// El largo seria pu - pu(factor_ini)
	// Los armo con un substract de kr_full(factor_cur) - kr_full(factor_ini)
	unsigned int sum_len = pu - select1_b->operator()(factor_ini + 1);
//	cout << "KarpRabinFactorsSuffixesv2::hashBin - sum_len: " << sum_len << "\n";
	
	kr_phrase = karp_rabin->subtract_prefix(arr_kr_s->at(factor_cur), arr_kr_s->at(factor_ini), sum_len);
	
//	cout << "KarpRabinFactorsSuffixesv2::hashBin - kr_phrase: " << kr_phrase << "\n";
	
	kr_total = karp_rabin->concat(kr_total, kr_phrase, sum_len);
	
	
	
//	cout << "KarpRabinFactorsSuffixesv2::hashBin - factor_cur: " << factor_cur << ", tu: " << tu << ", pu: " << pu << ", lu: " << lu << "\n";
	
	nano3 += timer.getNanosec();
	timer.reset();
	
//	cout << "KarpRabinFactorsSuffixesv2::hashBin - Preparing to add last prefix (length: " << length << ")\n";
	// Ahora se que el length que queda es solo un prefijo de la frase actual, de largo menor a lu
//	cout << "KarpRabinFactorsSuffixesv2::hashBin - hash 3 (ref_text + " << tu << ", " << length + offset + offset_ini - pu << ")\n";
	kr_prefix = karp_rabin->hash(ref_text + tu, length + offset + offset_ini - pu );
//	cout << "KarpRabinFactorsSuffixesv2::hashBin - kr_prefix: " << kr_prefix << "\n";
	kr_total = karp_rabin->concat(kr_total, kr_prefix, length + offset + offset_ini - pu);
	
	nano4 += timer.getNanosec();
	
//	cout << "KarpRabinFactorsSuffixesv2::hashBin - End (" << kr_total << ")\n";
	return kr_total;
	
}
*/

void KarpRabinFactorsSuffixesv2::save(const string &file){
//	n_factors = 0;
//	arr_kr_s = 0;
//	factors_start = 0;
	
	fstream writer(file, fstream::out | fstream::trunc);
	
	writer.write((char*)&n_factors, sizeof(int));
	
	unsigned int n_arr = arr_kr_s->size();
	writer.write((char*)&n_arr, sizeof(int));
	for(unsigned int i = 0; i < n_arr; ++i){
		unsigned int kr = (unsigned int )(arr_kr_s->at(i));
		writer.write((char*)&kr, sizeof(int));
	}
	
	unsigned long long worst_case = 5 * factors_start->size();
	unsigned char *buff = new unsigned char[ worst_case ];
	unsigned long long cur_byte = 0;
	BitsUtils utils;
	unsigned int last = 0;
	for(unsigned int i = 0; i < factors_start->size(); ++i){
		unsigned int delta = factors_start->at(i) - last;
		last = factors_start->at(i);
		cur_byte += utils.write_varbyte(buff + cur_byte, delta);
	}
	writer.write((char*)&cur_byte, sizeof(long long));
	writer.write((char*)buff, cur_byte);
	
	delete [] buff;
	
	writer.close();
	
//	cout << "KarpRabinFactorsSuffixesv2::save - Verification\n";
//	cout << "n_factors: " << n_factors << " (arr_kr_s: " << arr_kr_s->size() << ", factors_start: " << factors_start->size() << ")\n";
//	for(unsigned int i = 0; i < 10; ++i){
//		cout << "arr_kr_s[" << i << "]: " << arr_kr_s->at(i) << "\n";
//	}
//	for(unsigned int i = arr_kr_s->size() - 10; i < arr_kr_s->size(); ++i){
//		cout << "arr_kr_s[" << i << "]: " << arr_kr_s->at(i) << "\n";
//	}
//	for(unsigned int i = 0; i < 10; ++i){
//		cout << "factors_start[" << i << "]: " << factors_start->at(i) << "\n";
//	}
//	for(unsigned int i = factors_start->size() - 10; i < factors_start->size(); ++i){
//		cout << "factors_start[" << i << "]: " << factors_start->at(i) << "\n";
//	}
//	cout << "KarpRabinFactorsSuffixesv2::save - Verification End\n";
	
	cout << "KarpRabinFactorsSuffixesv2::save - Times: " << nano1 << ", " << nano2 << ", " << nano3 << ", " << nano4 << "\n";
	
}

void KarpRabinFactorsSuffixesv2::load(const string &file){
//	n_factors = 0;
//	arr_kr_s = 0;
//	factors_start = 0;

	fstream reader(file, fstream::in);
	
	n_factors = 0;
	reader.read((char*)&n_factors, sizeof(int));
	
	delete_krs = true;
	arr_kr_s = new vector<unsigned long long>();
	
	unsigned int n_arr = 0;
	reader.read((char*)&n_arr, sizeof(int));
	for(unsigned int i = 0; i < n_arr; ++i){
		unsigned int kr = 0;
		reader.read((char*)&kr, sizeof(int));
		arr_kr_s->push_back(kr);
	}
	
	unsigned long long total_bytes = 0;
	reader.read((char*)&total_bytes, sizeof(long long));
	
	unsigned char *buff = new unsigned char[total_bytes];
	reader.read((char*)buff, total_bytes);
	
	delete_start = true;
	factors_start = new vector<unsigned int>();
	
	BitsUtils utils;
	unsigned long long cur_byte = 0;
	unsigned int sum = 0;
	for(unsigned int i = 0; i < n_factors; ++i){
		unsigned int delta = 0;
		cur_byte += utils.read_varbyte(buff + cur_byte, delta);
		sum += delta;
		factors_start->push_back(sum);
	}
	
	delete [] buff;
	
	reader.close();
	
//	cout << "KarpRabinFactorsSuffixesv2::load - Verification\n";
//	cout << "n_factors: " << n_factors << " (arr_kr_s: " << arr_kr_s->size() << ", factors_start: " << factors_start->size() << ")\n";
//	for(unsigned int i = 0; i < 10; ++i){
//		cout << "arr_kr_s[" << i << "]: " << arr_kr_s->at(i) << "\n";
//	}
//	for(unsigned int i = arr_kr_s->size() - 10; i < arr_kr_s->size(); ++i){
//		cout << "arr_kr_s[" << i << "]: " << arr_kr_s->at(i) << "\n";
//	}
//	for(unsigned int i = 0; i < 10; ++i){
//		cout << "factors_start[" << i << "]: " << factors_start->at(i) << "\n";
//	}
//	for(unsigned int i = factors_start->size() - 10; i < factors_start->size(); ++i){
//		cout << "factors_start[" << i << "]: " << factors_start->at(i) << "\n";
//	}
//	cout << "KarpRabinFactorsSuffixesv2::load - Verification End\n";
	
}














