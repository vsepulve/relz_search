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
}

KarpRabinFactorsSuffixes::KarpRabinFactorsSuffixes(unsigned int _n_factors, vector<unsigned long long> *_arr_kr_s, KarpRabin *_karp_rabin, const char *_ref_text, inv_perm_support<> *_perm_inv, vector<unsigned int> *_arr_tu, vector<unsigned int> *_arr_pu, vector<unsigned int> *_arr_lu){
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
//	unsigned int pu = 0;
	unsigned int lu = 0;
	while( cur_len < length ){
		cur_perm = (*perm_inv)[factor_cur];
		tu = arr_tu->at(cur_perm);
//		pu = arr_pu->at(cur_perm);
		lu = arr_lu->at(cur_perm);
		if( cur_len + lu > length ){
			// salir y procesar los caracteres necesarios de este factor
			// Aqui mismo podria procesar de una vez desde factor_ini hasta factor_cur-1 (solo si factor_cur > factor_ini)
			if( factor_cur > factor_ini ){
				kr1 = karp_rabin->subtract_prefix(arr_kr_s->at(factor_cur), arr_kr_s->at(factor_ini), cur_len);
//				cout << "KarpRabinFactorsSuffixes::hash - kr1: " << kr1 << ", cur_len: " << cur_len << "\n";
			}
			break;
		}
		// Agregar factor completo y continuar
		cur_len += lu;
		++factor_cur;
	}
	// Calcular el hash de los (length - cur_len) primeros caracteres del factor actual
	// Notar que lo que sigue se puede extraer er KarpRabinReference (desde tu, de largo length - cur_len)
	unsigned long long kr2 = karp_rabin->hash(ref_text + tu, length - cur_len);
//	cout << "KarpRabinFactorsSuffixes::hash - kr2: " << kr2 << ", len: " << length - cur_len << "\n";
	
	unsigned long long kr12 = karp_rabin->concat(kr1, kr2, length - cur_len);
//	cout << "KarpRabinFactorsSuffixes::hash - End (kr12: " << kr12 << ")\n";
	
	return kr12;
	
}

unsigned long long KarpRabinFactorsSuffixes::hash(unsigned int factor_ini, unsigned int offset, unsigned int length){
//	cout << "KarpRabinFactorsSuffixes::hash - Start (factor_ini: " << factor_ini << ", offset: " << offset << ", length: " << length << ")\n";
	
	// Una opcion es calcular el hash desde factor_ini de largo length
	// Luego el hash desde factor_ini de largo offset
	// Asi se pueden restar para tener el del sujifo
	
	unsigned long long kr_total = hash(factor_ini, length);
//	cout << "KarpRabinFactorsSuffixes::hash - kr_total: " << kr_total << ")\n";
	unsigned long long kr_prefix = hash(factor_ini, offset);
//	cout << "KarpRabinFactorsSuffixes::hash - kr_prefix: " << kr_prefix << ")\n";
	unsigned long long kr = karp_rabin->subtract_prefix(kr_total, kr_prefix, length - offset);
//	cout << "KarpRabinFactorsSuffixes::hash - kr: " << kr << ")\n";
	
	return kr;
}














