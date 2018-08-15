#include "KarpRabinFactorsSuffixes.h"

KarpRabinFactorsSuffixes::KarpRabinFactorsSuffixes(){
	n_factors = 0;
	arr_kr_s = NULL;
	karp_rabin = NULL;
	ref_text = NULL;
	perm_inv = NULL;
	arr_tu = NULL;
	arr_pu = NULL;
	arr_lu = NULL;
	factors_start = NULL;
}

KarpRabinFactorsSuffixes::KarpRabinFactorsSuffixes(unsigned int _n_factors, vector<unsigned long long> *_arr_kr_s, KarpRabin *_karp_rabin, const char *_ref_text, inv_perm_support<> *_perm_inv, vector<unsigned int> *_arr_tu, vector<unsigned int> *_arr_pu, vector<unsigned int> *_arr_lu, vector<unsigned int> *_factors_start){
	n_factors = _n_factors;
	arr_kr_s = _arr_kr_s;
	karp_rabin = _karp_rabin;
	if( arr_kr_s->size() < n_factors+1 ){
		cout << "KarpRabinFactorsSuffixes - Warning (insufficient prefixes, must include one for the whole collection)\n";
	}
	ref_text = _ref_text;
	perm_inv = _perm_inv;
	arr_tu = _arr_tu;
	arr_pu = _arr_pu;
	arr_lu = _arr_lu;
	factors_start = _factors_start;
}

KarpRabinFactorsSuffixes::~KarpRabinFactorsSuffixes(){
	n_factors = 0;
	arr_kr_s = NULL;
	karp_rabin = NULL;
	ref_text = NULL;
}

unsigned long long KarpRabinFactorsSuffixes::hash(unsigned int factor_ini, unsigned int length){
//	cout << "KarpRabinFactorsSuffixes::hash - Start (factor_ini: " << factor_ini << ", length: " << length << ")\n";
	unsigned int cur_len = 0;
	unsigned int factor_cur = factor_ini;
	unsigned long long kr1 = 0;
	unsigned int cur_perm = 0;
	unsigned int tu = 0;
	unsigned int lu = 0;
	while( cur_len < length ){
		cout << "a\n";
		cur_perm = (*perm_inv)[factor_cur];
		cout << "b (cur_perm: " << cur_perm << " / " << arr_tu->size() << ")\n";
		tu = arr_tu->at(cur_perm);
		lu = arr_lu->at(cur_perm);
		if( cur_len + lu >= length ){
			break;
		}
		// Agregar factor completo y continuar
		cur_len += lu;
		++factor_cur;
		
//		cout << "KarpRabinFactorsSuffixes::hash - karp_rabin->subtract_prefix (lu: " << lu << ")\n";
		unsigned long long kr_phrase = karp_rabin->subtract_prefix(arr_kr_s->at(factor_cur), arr_kr_s->at(factor_cur-1), lu);
//		cout << "KarpRabinFactorsSuffixes::hash - Adding tu current hash\n";
		kr1 = karp_rabin->concat(kr1, kr_phrase, lu);
//		cout << "KarpRabinFactorsSuffixes::hash - cur kr1: " << kr1 << "\n";
		
	}
	// Calcular el hash de los (length - cur_len) primeros caracteres del factor actual
	// Notar que lo que sigue se puede extraer er KarpRabinReference (desde tu, de largo length - cur_len)
	unsigned long long kr2 = karp_rabin->hash(ref_text + tu, length - cur_len);
//	cout << "KarpRabinFactorsSuffixes::hash - kr2: " << kr2 << ", len: " << length - cur_len << "\n";
	
	unsigned long long kr12 = karp_rabin->concat(kr1, kr2, length - cur_len);
//	cout << "KarpRabinFactorsSuffixes::hash - End (kr12: " << kr12 << ")\n";
	
	return kr12;
	
}


// Fast version, only valid for offset and length < karp_rabin->getTableSize()
unsigned long long KarpRabinFactorsSuffixes::hashFast(unsigned int factor_ini, unsigned int length){
//	cout << "KarpRabinFactorsSuffixes::hashFast - Start (factor_ini: " << factor_ini << ", length: " << length << ")\n";
	unsigned long long kr1 = 0;
	
	/*
	unsigned int cur_len = 0;
	unsigned int factor_cur = factor_ini;
	unsigned int cur_perm = (*perm_inv)[factor_cur];
	unsigned int tu = arr_tu->at(cur_perm);
	unsigned int lu = arr_lu->at(cur_perm);
//	cout << "1...\n";
	while( cur_len + lu < length ){
		// Agregar factor completo y continuar
//		cout << "KarpRabinFactorsSuffixes::hashFast - Adding factor " << factor_cur << " of len " << lu << " (cur_len: " << cur_len << " -> " << (cur_len + lu) << ")\n";
		cur_len += lu;
		++factor_cur;
		cur_perm = (*perm_inv)[factor_cur];
		tu = arr_tu->at(cur_perm);
		lu = arr_lu->at(cur_perm);
	}
	*/
//	cout << "2\n";
	
	// Prueba de Busqueda Binaria
	
	unsigned int l = factor_ini;
	unsigned int h = factors_start->size();
	unsigned int factor_cur_bin = factor_ini;
	unsigned int cur_perm_bin = (*perm_inv)[factor_cur_bin];
	unsigned int tu_bin = arr_tu->at(cur_perm_bin);
	while(l < h){
		factor_cur_bin = l + ((h-l)>>1);
		if( (factors_start->at(factor_cur_bin) - factors_start->at(factor_ini)) < length ){
			l = factor_cur_bin+1;
		}
		else{
			h = factor_cur_bin;
		}
	}
//	factor_cur_bin = h;
	if( (factor_cur_bin > factor_ini) && (factors_start->at(factor_cur_bin) - factors_start->at(factor_ini)) >= length ){
		--factor_cur_bin;
	}
	unsigned int cur_len_bin = (factors_start->at(factor_cur_bin) - factors_start->at(factor_ini));
	cur_perm_bin = (*perm_inv)[factor_cur_bin];
	tu_bin = arr_tu->at(cur_perm_bin);
	
//	cout << "2 bin (factor_cur_bin: " << factor_cur_bin << " / " << factor_cur << ", cur_perm_bin: " << cur_perm_bin << " / " << cur_perm << ", tu_bin: " << tu_bin << " / " << tu << ", cur_len_bin: " << cur_len_bin << " / " << cur_len << ")\n";
	
//	if( (factor_cur_bin != factor_cur)
//		|| (cur_perm_bin != cur_perm)
//		|| (tu_bin != tu) ){
//		cout << "Error en BB (factor_cur_bin: " << factor_cur_bin << " / " << factor_cur << ", cur_perm_bin: " << cur_perm_bin << " / " << cur_perm << ", tu_bin: " << tu_bin << " / " << tu << ")\n";
//		exit(0);
//	}
	
	
	if( factor_cur_bin > factor_ini ){
//		cout << "KarpRabinFactorsSuffixes::hashFast - Preparing kr1 subtract_prefix(" << arr_kr_s->at(factor_cur_bin) << ", " << arr_kr_s->at(factor_ini) << ", " << cur_len_bin << ")\n";
		kr1 = karp_rabin->subtract_prefix(arr_kr_s->at(factor_cur_bin), arr_kr_s->at(factor_ini), cur_len_bin);
	}
//	cout << "KarpRabinFactorsSuffixes::hashFast - kr1: " << kr1 << ", cur_len_bin: " << cur_len_bin << "\n";
	// Calcular el hash de los (length - cur_len) primeros caracteres del factor actual
	// Notar que lo que sigue se puede extraer er KarpRabinReference (desde tu, de largo length - cur_len)
	
	NanoTimer timer;
	unsigned long long kr2 = karp_rabin->hash(ref_text + tu_bin, length - cur_len_bin);
	kr_nano += timer.getNanosec();
	
//	cout << "KarpRabinFactorsSuffixes::hashFast - kr2: " << kr2 << ", len: " << length - cur_len_bin << "\n";
	unsigned long long kr12 = karp_rabin->concat(kr1, kr2, length - cur_len_bin);
//	cout << "KarpRabinFactorsSuffixes::hashFast - End (kr12: " << kr12 << ")\n";
	return kr12;
}

// Fast version, only valid for offset and length < karp_rabin->getTableSize()
unsigned long long KarpRabinFactorsSuffixes::hashFast(unsigned int factor_ini, unsigned int offset, unsigned int length){
//	cout << "KarpRabinFactorsSuffixes::hashFast - Start (factor_ini: " << factor_ini << ", offset: " << offset << ", length: " << length << ")\n";
	
	if( length == 0 ){
		return 0;
	}
	
	unsigned long long kr_total = hashFast(factor_ini, offset + length);
//	cout << "KarpRabinFactorsSuffixes::hashFast - kr_total: " << kr_total << "\n";
	unsigned long long kr_prefix = hashFast(factor_ini, offset);
//	cout << "KarpRabinFactorsSuffixes::hashFast - kr_prefix: " << kr_prefix << "\n";
	unsigned long long kr = karp_rabin->subtract_prefix(kr_total, kr_prefix, length);
//	cout << "KarpRabinFactorsSuffixes::hashFast - kr: " << kr << "\n";
	return kr;
}

unsigned long long KarpRabinFactorsSuffixes::hash(unsigned int factor_ini, unsigned int offset, unsigned int length){
	if( offset < karp_rabin->getTableSize() && length < karp_rabin->getTableSize() ){
		return hashFast(factor_ini, offset, length);
	}

//	cout << "KarpRabinFactorsSuffixes::hash - Start (factor_ini: " << factor_ini << ", offset: " << offset << ", length: " << length << ")\n";
	
	if( length == 0 ){
		return 0;
	}
	
	unsigned int factor_cur = factor_ini;
	
	// Primero omito los factores completos previos a offset
	unsigned int cur_perm = (*perm_inv)[factor_cur];
	unsigned int tu = arr_tu->at(cur_perm);
	unsigned int lu = arr_lu->at(cur_perm);
	while( lu < offset ){
//		cout << "KarpRabinFactorsSuffixes::hash - Ommiting factor " << factor_cur << "\n";
		offset -= lu;
		++factor_cur;
		cur_perm = (*perm_inv)[factor_cur];
		tu = arr_tu->at(cur_perm);
		lu = arr_lu->at(cur_perm);
	}
	
	if( lu == offset ){
		offset = 0;
		++factor_cur;
		cur_perm = (*perm_inv)[factor_cur];
		tu = arr_tu->at(cur_perm);
		lu = arr_lu->at(cur_perm);
	}
	
//	cout << "KarpRabinFactorsSuffixes::hash - Preparing to remove first prefix (offset: " << offset << ", length: " << length << ")\n";
	// Ahora estoy seguro que el offset se aplica a la frase actual
	// Tambien se que offset es menor que lu
	// Calculo el hash de offset y lo resto a la palabra actual
	unsigned long long kr_prefix = karp_rabin->hash(ref_text + tu, offset);
//	cout << "KarpRabinFactorsSuffixes::hash - kr_prefix: " << kr_prefix << "\n";
	
	// Puede pasar que el largo sea MENOR que el largo de la frase
	// En ese caso, la idea seria omitir el resto
	unsigned long long kr_phrase = 0;
	unsigned long long kr_total = 0;
	if( length + offset <= lu ){
		kr_phrase = karp_rabin->hash(ref_text + tu, length);
//		cout << "KarpRabinFactorsSuffixes::hash - kr_phrase (FINAL) [" << factor_cur << "]: " << kr_phrase << " (len: " << (length + offset) << ", lu: " << lu << ")\n";
		kr_total = karp_rabin->subtract_prefix(kr_phrase, kr_prefix, length + offset);
//		cout << "KarpRabinFactorsSuffixes::hash - kr_total (FINAL): " << kr_total << "\n";
		return kr_total;
	}
	else{
		kr_phrase = karp_rabin->subtract_prefix(arr_kr_s->at(factor_cur+1), arr_kr_s->at(factor_cur), lu);
//		cout << "KarpRabinFactorsSuffixes::hash - kr_phrase [" << factor_cur << "]: " << kr_phrase << "\n";
		kr_total = karp_rabin->subtract_prefix(kr_phrase, kr_prefix, lu - offset);
//		cout << "KarpRabinFactorsSuffixes::hash - kr_total: " << kr_total << "\n";
	}
	
	length -= (lu - offset);
	++factor_cur;
	cur_perm = (*perm_inv)[factor_cur];
	tu = arr_tu->at(cur_perm);
	lu = arr_lu->at(cur_perm);
	
	// Ahora agrego el hash de cada frase completa contenida
	while( lu < length ){
//		cout << "KarpRabinFactorsSuffixes::hash - Adding factor " << factor_cur << " (lu: " << lu << ", length: " << length << ")\n";
		
		// Agrego el hash de la frase
		kr_phrase = karp_rabin->subtract_prefix(arr_kr_s->at(factor_cur+1), arr_kr_s->at(factor_cur), lu);
//		cout << "KarpRabinFactorsSuffixes::hash - kr_phrase [" << factor_cur << "]: " << kr_phrase << "\n";
		
		kr_total = karp_rabin->concat(kr_total, kr_phrase, lu);
//		cout << "KarpRabinFactorsSuffixes::hash - kr_total " << kr_total << "\n";
		
		length -= lu;
		++factor_cur;
		cur_perm = (*perm_inv)[factor_cur];
		tu = arr_tu->at(cur_perm);
		lu = arr_lu->at(cur_perm);
	}
	
//	cout << "KarpRabinFactorsSuffixes::hash - Preparing to add last prefix (length: " << length << ")\n";
	// Ahora se que el length que queda es solo un prefijo de la frase actual, de largo menor a lu
	kr_prefix = karp_rabin->hash(ref_text + tu, length);
//	cout << "KarpRabinFactorsSuffixes::hash - kr_prefix: " << kr_prefix << "\n";
	kr_total = karp_rabin->concat(kr_total, kr_prefix, length);
	
//	cout << "KarpRabinFactorsSuffixes::hash - End (kr_total: " << kr_total << ")\n";
	return kr_total;
	
}














