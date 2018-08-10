#include "KarpRabin.h"

KarpRabin::KarpRabin(){
	voc_bits = 8;
	kr_mod = 787;
	table_size = 0;
	pow_table = NULL;
}

KarpRabin::KarpRabin(unsigned int _voc_bits, unsigned int _kr_mod, unsigned int _table_size){
	voc_bits = _voc_bits;
	kr_mod = _kr_mod;
	table_size = _table_size + 1;
	pow_table = new unsigned long long[table_size];
	pow_table[0] = 1;
	for(unsigned int i = 1; i < table_size; ++i){
		pow_table[i] = (pow_table[i-1] * (1ull<<voc_bits)) % kr_mod;
	}
}

KarpRabin::~KarpRabin(){
	if( pow_table != NULL ){
		delete [] pow_table;
		pow_table = NULL;
	}
	table_size = 0;
}
	
// Direct version (*x, y times)
unsigned long long KarpRabin::ullpow(unsigned long long x, unsigned int y){
	unsigned long long ret = 1;
	for(unsigned int i = 0; i < y; i++){
		ret = (ret*x) % kr_mod;
	}
	return ret;
}

// Base 2 bits version (2^bits^y), recursive to avoid overflow
unsigned long long KarpRabin::ullpow2_rec(unsigned int bits, unsigned int y){
	if( bits * y < 64 ){
		return ((1ull << (bits*y) ) % kr_mod);
	}
	else{
//		cout << "KarpRabin::ullpow2_rec (" << (bits * y) << " => " << (y/2) << " | " << (y/2 + (y & 0x1)) << ")\n";
		return ( ullpow2_rec(bits, y>>1) * ullpow2_rec(bits, (y>>1) + (y & 0x1)) ) % kr_mod;
	}
}

unsigned long long KarpRabin::ullpow2_table(unsigned int bits, unsigned int y){
	assert(y < table_size);
	assert(pow_table != NULL);
	return pow_table[y];
}


// Evaluate the full hash in str.length() operations
unsigned long long KarpRabin::hash(const string &str){
	unsigned long long ret = 0;
	size_t str_len = str.length();
	for(unsigned int i = 0, k = str_len-1; i < str_len; i++, k--) {
//		ret = ret + (str[i] * ullpow(1<<voc_bits, k)) % kr_mod;
//		ret = ret + (str[i] * ullpow2_rec(voc_bits, k)) % kr_mod;
		ret = ret + ((unsigned long long)(str[i]) * ullpow2_table(voc_bits, k)) % kr_mod;
		ret = ret % kr_mod;
	}
	return ret;
}
unsigned long long KarpRabin::hash(const char *str, unsigned long long str_len){
	
//	string s(str, str_len);
//	cout << "KarpRabin::hash - Start (" << s << ")\n";
	
	unsigned long long ret = 0;
	for(unsigned int i = 0, k = str_len-1; i < str_len; i++, k--) {
//		ret = ret + (str[i] * ullpow(1<<voc_bits, k)) % kr_mod;
//		ret = ret + (str[i] * ullpow2_rec(voc_bits, k)) % kr_mod;
		ret = ret + ((unsigned long long)(str[i]) * ullpow2_table(voc_bits, k)) % kr_mod;
		ret = ret % kr_mod;
	}
	
//	cout << "KarpRabin::hash - End (" << ret << ")\n";
	
	return ret;
}

// Evaluate the hash of the concatenation in constant time
unsigned long long KarpRabin::concat(unsigned long long kr1, unsigned long long kr2, unsigned int len2){
//	return (kb2 + (kb1 * ullpow(1<<voc_bits, len2)) % kr_mod) % kr_mod;
//	return (kr2 + (kr1 * ullpow2_rec(voc_bits, len2)) % kr_mod) % kr_mod;
	return (kr2 + (kr1 * ullpow2_table(voc_bits, len2)) % kr_mod) % kr_mod;
}

// Evaluate the hash of the subtract in constant time
unsigned long long KarpRabin::subtract_prefix(unsigned long long kr12, unsigned long long kr1, unsigned int len2){
//	cout << "KarpRabin::subtract - B^" << len2 << ": " << ullpow2_table(voc_bits, len2) << "\n";
//	cout << "KarpRabin::subtract - kr1 (" << kr1 << ") * B^" << len2 << ": " << ((kr1 * ullpow2_table(voc_bits, len2)) % kr_mod) << "\n";
//	cout << "KarpRabin::subtract - kr12 - kr1 * B^" << len2 << ": " << (kr12 + kr_mod - ((kr1 * ullpow2_table(voc_bits, len2)) % kr_mod)) % kr_mod << "\n";
	return (kr12 + kr_mod - ((kr1 * ullpow2_table(voc_bits, len2)) % kr_mod)) % kr_mod;
}

unsigned int KarpRabin::nextFactor(unsigned int start, vector<unsigned int> &factors_start){
	
//	cout << "KarpRabin::nextFactor - Start (" << start << ")\n";
	
	unsigned int l, h, m, res;
	
	l = 0;
	h = factors_start.size() - 1;
	
	while(l < h){
		m = l + ((h-l)>>1);
//		cout << "KarpRabin::nextFactor - m: " << m << ", factor: " << factors_start[m] << ", start: " << start << "\n";
		if( factors_start[m] <= start ){
//			cout << "KarpRabin::nextFactor - caso 1: l = " << (m+1) << "\n";
			l = m+1;
		}
		else{
//			cout << "KarpRabin::nextFactor - caso 2: h = " << m << "\n";
			h = m;
		}
	}
	res = h;
//	if( (res > 0) && factors_start[m] > start ){
//		--res;
//	}
	
//	cout << "KarpRabin::nextFactor - retornando " << res << " (" << factors_start[res] << " > " << start << " ?)\n";
	return res;
}

unsigned long long KarpRabin::hash(const char *full_text, unsigned int start, unsigned int str_len, vector<unsigned int> &factors_start){
	if( str_len < table_size ){
		return hash(full_text + start, str_len);
	}
	else{
		cout << "KarpRabin::hash - Starting special process using factors (start: " << start << ", str_len: " << str_len << ")\n";
		
		unsigned int next_factor = nextFactor(start, factors_start);
		unsigned int next_start = factors_start[next_factor];
		unsigned long long ret = 0;
		unsigned long long kr1;
		
		unsigned int cur_start = start;
		unsigned int cur_len = next_start - cur_start;
		
//		cout << "KarpRabin::hash - Begining while with cur_start: " << cur_start << ", cur_len: " << cur_len << ", next_start: " << next_start << "\n";
		
		while( cur_start < start + str_len ){
//			cout << "KarpRabin::hash - hash cur_start: " << cur_start << ", cur_len: " << cur_len << "\n";
			kr1 = hash(full_text + cur_start, cur_len);
			// concat
			ret = concat(ret, kr1, cur_len);
			// Preparo siguiente
//			cout << "KarpRabin::hash - Preparing next iteration\n";
			cur_start = next_start;
			++next_factor;
//			cout << "KarpRabin::hash - next_factor: " << next_factor << "/ " << factors_start.size() << "\n";
			if( next_factor < factors_start.size() ){
				next_start = factors_start[next_factor];
			}
			else{
				next_start = start + str_len;
			}
			cur_len = next_start - cur_start;
//			cout << "KarpRabin::hash - Checking cur_len: " << cur_len << ", next_start: " << next_start << ", str_len: " << str_len << ", cur_start: " << cur_start << "\n";
			if( cur_len > str_len - cur_start ){
				cur_len = str_len - cur_start;
			}
//			cout << "KarpRabin::hash - next cur_start: " << cur_start << ", cur_len: " << cur_len << ", next_start: " << next_start << "\n";
		}
		
		return ret;
	}
}
















