#include "HashTrie.h"

constexpr char HashTrie::decodeChar[];

// Node Methods

HashTrieNode::HashTrieNode(){
	len = 0;
	min = 0;
	hash = 0;
//	first = 0;
	first = 'A';
}

HashTrieNode::~HashTrieNode(){
}

void HashTrieNode::compactData(unsigned int &next_pos, int_vector<> &positions_childs, int_vector<> &n_childs, int_vector<> &len_childs, int_vector<> &min_childs, int_vector<> &hash_childs, int_vector<> &first_childs){
	// Guardo la posicion para los hijos de este nodo
	// El valor de next_pos cambiara con cada hijo
	unsigned int cur_pos = next_pos;
	next_pos += childs_vector.size();
//	cout << "HashTrieNode::compactData - childs_vector.size(): " << childs_vector.size() << "\n";
	for(unsigned int i = 0; i < childs_vector.size(); ++i){
		// Data for the child
		positions_childs[cur_pos] = next_pos;
//		cout << "HashTrieNode::compactData - positions_childs[" << cur_pos << "]: " << next_pos << "\n";
		n_childs[cur_pos] = childs_vector[i].childs_vector.size();
		// Omit the length for the leaf nodes
		if( n_childs[cur_pos] == 0 ){
			len_childs[cur_pos] = 0;
		}
		else{
			len_childs[cur_pos] = childs_vector[i].len;
		}
		min_childs[cur_pos] = childs_vector[i].min;
		hash_childs[cur_pos] = childs_vector[i].hash;
		first_childs[cur_pos] = HashTrie::codeChar(childs_vector[i].first);
		// Recursive call
		childs_vector[i].compactData(next_pos, positions_childs, n_childs, len_childs, min_childs, hash_childs, first_childs);
		++cur_pos;
	}
}

// Build for the direct version (Array Y)
void HashTrieNode::build(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *arr_y, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len){
	if( min == max ){
//		cout << "HashTrieNode::build - Leaf\n";
		return;
	}
	if( full_text == NULL || processed_len == len_text ){
//		cout << "HashTrieNode::build - Null text\n";
		return;
	}
	
//	cout << "HashTrieNode::build - Start (full text of " << len_text << ", range[" << min << ", " << max << "], processed_len: " << processed_len << ")\n";
	
	unsigned int min_pos = min;
	unsigned int min_text_start = factors_start[ (*arr_y)[min_pos] ] + processed_len;
	unsigned int min_text_len = 0;
	
	unsigned int cur_pos = 0;
	unsigned int cur_text_start = 0;
	
	unsigned long long hash;
	
	while( min_pos <= max ){
		
		// Primero busco el cur_pos MAYOR que tenga *primer* char igual al de min_pos
		// Despues reviso el largo mayor
		
//		cout << "HashTrieNode::build - Iniciando busqueda binaria en [" << min_pos << ", " << max << "]\n";
		
		unsigned int l = min_pos;
		unsigned int h = max;
		while(l < h){
			cur_pos = l + ((h-l)>>1);
			cur_text_start = factors_start[ (*arr_y)[cur_pos] ] + processed_len;
			if( full_text[cur_text_start] <= full_text[min_text_start] ){
				l = cur_pos+1;
			}
			else{
				h = cur_pos;
			}
		}
		cur_pos = h;
		cur_text_start = factors_start[ (*arr_y)[cur_pos] ] + processed_len;
		if( (cur_pos > min_pos) && full_text[cur_text_start] != full_text[min_text_start] ){
			--cur_pos;
			cur_text_start = factors_start[ (*arr_y)[cur_pos] ] + processed_len;
		}
//		cout << "HashTrieNode::build - cur_pos: " << cur_pos << "\n";
		
		// Si cur_pos == min_pos, el largo es min_text_len
		// Si no, buscar el largo maximo entre min y cur
		
		if( cur_pos > min_pos ){
//			cout << "HashTrieNode::build - Case 1\n";
			// Buscar el mayor
			// Aqui uso el hecho de que full_text termina en '\0'
			min_text_len = 1;
			while( full_text[min_text_start + min_text_len] == full_text[cur_text_start + min_text_len] ){
				++min_text_len;
			}
		}
		else{
//			cout << "HashTrieNode::build - Case 2\n";
			min_text_len = len_text - min_text_start;
		}
//		cout << "HashTrieNode::build - Preparing string s (min_text_start: " << min_text_start << " / " << len_text << ", len: " << min_text_len << ")\n";
		
		// El string es solo para el mensaje de debug
//		string s(full_text + min_text_start, min_text_len);
		char first_char = full_text[min_text_start];
		
//		cout << "HashTrieNode::build - first_char: " << first_char << ")\n";
		
		// Notar que aqui necesito la posicion del factor en la coleccion, no en el arreglo Y
		hash = kr_factors->hash( (*arr_y)[min_pos], processed_len, min_text_len);
//		cout << "HashTrieNode::build - Adding range [" << min_pos << ", " << cur_pos+1 << "), len " << min_text_len << ", hash: " << hash << " / " << karp_rabin->hash(s) << " (\"" << s << "\")\n";
//		cout << "HashTrieNode::build - Adding range [" << min_pos << ", " << cur_pos+1 << "), len " << min_text_len << ", hash: " << hash << "\n";
		// Preparar el hijo, ejecutar la llamda sobre esa instancia
		if( min_text_len == 0 ){
//			cout << "HashTrieNode::build - Omiting child of len 0\n";
		}
		else{
//			cout << "HashTrieNode::build - Preparing Child\n";
			childs_vector.push_back(HashTrieNode());
			childs_vector.back().first = first_char;
			childs_vector.back().len = min_text_len;
			childs_vector.back().min = min_pos;
			childs_vector.back().hash = hash;
//			cout << "HashTrieNode::build - Child Ready, starting recursive build\n";
			childs_vector.back().build(full_text, len_text, factors_start, arr_y, karp_rabin, kr_factors, min_pos, cur_pos, processed_len + min_text_len);
		}
		
//		cout << "HashTrieNode::build - Preparing min_pos = " << cur_pos+1 << " / " << arr_y->size() << "\n";
		min_pos = cur_pos+1;
		if(min_pos >= arr_y->size()){
			break;
		}
		min_text_start = factors_start[ (*arr_y)[min_pos] ] + processed_len;
		
//		cout << "HashTrieNode::build - min_text_start: " << min_text_start << "\n";
	
	}
	
	if( childs_vector.size() > 1 ){
//		cout << "HashTrieNode::build - Sort (" << childs_vector.size() << " childs)\n";
		std::sort(childs_vector.begin(), childs_vector.end(), 
			[](const HashTrieNode &p1, const HashTrieNode &p2) -> bool{return p1.first < p2.first;}
		);
	}
	
//	cout << "----- \n";
	
}

// Build for the reverse version (Array X)
void HashTrieNode::buildRev(const char *full_text, unsigned int len_text, vector<unsigned int> &factors_start, int_vector<> *arr_x, KarpRabin *karp_rabin, KarpRabinFactorsSuffixes *kr_factors, unsigned int min, unsigned int max, unsigned int processed_len){
	if( min == max ){
//		cout << "HashTrieNode::buildRev - Leaf\n";
		return;
	}
	if( full_text == NULL || processed_len == len_text ){
//		cout << "HashTrieNode::buildRev - Null text\n";
		return;
	}
	
//	cout << "HashTrieNode::buildRev - Start (full text of " << len_text << ", range[" << min << ", " << max << "], processed_len: " << processed_len << ")\n";
	
	unsigned int min_pos = min;
	
	unsigned int min_text_len = 0;
	if( (*arr_x)[min_pos] > 0 ){
		min_text_len = factors_start[ (*arr_x)[min_pos] ] - factors_start[ (*arr_x)[min_pos] - 1 ] - processed_len;
	}
	while( min_text_len == 0 && min_pos < max ){
		++min_pos;
		if( (*arr_x)[min_pos] > 0 ){
			min_text_len = factors_start[ (*arr_x)[min_pos] ] - factors_start[ (*arr_x)[min_pos] - 1 ] - processed_len;
		}
	}
	
	unsigned int min_text_start = 0;
	if( (*arr_x)[min_pos] > 0 ){
		min_text_start = factors_start[ (*arr_x)[min_pos] ] - processed_len - 1;
	}
	
	unsigned int cur_pos = 0;
	unsigned int cur_text_start = 0;
	unsigned int cur_text_len = 0;
	unsigned int min_common_text = 0;
	
	unsigned long long hash;
	
	while( min_pos <= max ){
		
		// Primero busco el cur_pos MAYOR que tenga *primer* char igual al de min_pos
		// Despues reviso el largo mayor
		
//		cout << "HashTrieNode::buildRev - Iniciando busqueda binaria en [" << min_pos << ", " << max << "]\n";
		
		unsigned int l = min_pos;
		unsigned int h = max;
		while(l < h){
			cur_pos = l + ((h-l)>>1);
			cur_text_start = factors_start[ (*arr_x)[cur_pos] ] - processed_len - 1;
			if( (*arr_x)[cur_pos] == 0 ){
				cur_text_len = 0;
			}
			else{
				cur_text_len = factors_start[ (*arr_x)[cur_pos] ] - factors_start[ (*arr_x)[cur_pos] - 1 ] - processed_len;
			}
			if( full_text[cur_text_start] <= full_text[min_text_start] ){
				l = cur_pos+1;
			}
			else{
				h = cur_pos;
			}
		}
		cur_pos = h;
		cur_text_start = factors_start[ (*arr_x)[cur_pos] ] - processed_len - 1;
		if( (*arr_x)[cur_pos] == 0 ){
			cur_text_len = 0;
		}
		else{
			cur_text_len = factors_start[ (*arr_x)[cur_pos] ] - factors_start[ (*arr_x)[cur_pos] - 1 ] - processed_len;
		}
		if( (cur_pos > min_pos) && full_text[cur_text_start] != full_text[min_text_start] ){
			--cur_pos;
			cur_text_start = factors_start[ (*arr_x)[cur_pos] ] - processed_len - 1;
			if( (*arr_x)[cur_pos] == 0 ){
				cur_text_len = 0;
			}
			else{
				cur_text_len = factors_start[ (*arr_x)[cur_pos] ] - factors_start[ (*arr_x)[cur_pos] - 1 ] - processed_len;
			}
		}
//		cout << "HashTrieNode::buildRev - cur_pos: " << cur_pos << "\n";
		
		// Si cur_pos == min_pos, el largo es min_text_len
		// Si no, buscar el largo maximo entre min y cur
		
		if( cur_pos > min_pos ){
//			cout << "HashTrieNode::buildRev - Case 1\n";
			// Buscar el mayor
			// Aqui uso el hecho de que full_text termina en '\0'
			min_common_text = 0;
			while( min_common_text < ((min_text_len<cur_text_len)?min_text_len:cur_text_len) 
				&& full_text[min_text_start - min_common_text] == full_text[cur_text_start - min_common_text] 
				){
				++min_common_text;
			}
		}
		else{
//			cout << "HashTrieNode::buildRev - Case 2\n";
			min_common_text = min_text_len;
		}
//		cout << "HashTrieNode::buildRev - min_common_text: " << min_common_text << "\n";
//		cout << "HashTrieNode::buildRev - Preparing string s (min_text_start: " << min_text_start << " / " << len_text << ", len: " << min_common_text << ")\n";
		
		// En el caso de Rev, uso s para el hash directo
		string s;
		for( unsigned int i = 0; i < min_common_text; ++i ){
			s += *(full_text + min_text_start - i);
		}
		hash = karp_rabin->hash(s);
//		cout << "HashTrieNode::buildRev - s: " << s << "\n";
		char first_char = s[0];
//		cout << "HashTrieNode::buildRev - hash: " << hash << "\n";
		
//		cout << "HashTrieNode::buildRev - Adding range [" << min_pos << ", " << cur_pos+1 << "), len " << min_common_text << ", hash: " << hash << "\n";
		// Preparar el hijo, ejecutar la llamda sobre esa instancia
		if( min_common_text == 0 ){
//			cout << "HashTrieNode::buildRev - Omiting child of len 0\n";
		}
		else{
			childs_vector.push_back(HashTrieNode());
			childs_vector.back().first = first_char;
			childs_vector.back().len = min_common_text;
			childs_vector.back().min = min_pos;
			childs_vector.back().hash = hash;
			childs_vector.back().buildRev(full_text, len_text, factors_start, arr_x, karp_rabin, kr_factors, min_pos, cur_pos, processed_len + min_common_text);
		}
		
//		cout << "HashTrieNode::buildRev - Preparing min_pos = " << cur_pos+1 << " / " << arr_x->size() << "\n";
		min_pos = cur_pos+1;
		if(min_pos >= arr_x->size()){
			break;
		}
		min_text_start = 0;
		min_text_len = 0;
		if( (*arr_x)[min_pos] > 0 ){
			min_text_start = factors_start[ (*arr_x)[min_pos] ] - processed_len - 1;
			min_text_len = factors_start[ (*arr_x)[min_pos] ] - factors_start[ (*arr_x)[min_pos] - 1 ] - processed_len;
		}
//		cout << "HashTrieNode::buildRev - min_text_start: " << min_text_start << " de len " << min_text_len << "\n";
	
	}
	
	if( childs_vector.size() > 1 ){
		std::sort(childs_vector.begin(), childs_vector.end(), 
			[](const HashTrieNode &p1, const HashTrieNode &p2) -> bool{return p1.first < p2.first;}
		);
	}
	
//	cout << "----- \n";
	
}

unsigned int HashTrieNode::totalChilds(unsigned int &max_len, unsigned int &max_childs, unsigned int &max_height, unsigned int height){
	unsigned int ret = childs_vector.size();
	if( height > max_height ){
		max_height = height;
	}
	if( len > max_len ){
		max_len = len;
	}
	if( childs_vector.size() > max_childs ){
		max_childs = childs_vector.size();
	}
	for( auto it : childs_vector ){
		ret += it.totalChilds(max_len, max_childs, max_height, 1 + height);
	}
	return ret;
}

// Tree Methods
	
unsigned int HashTrie::codeChar(char c){
	if(c == 'A'){
		return 0;
	}
	else if(c == 'C'){
		return 1;
	}
	else if(c == 'G'){
		return 2;
	}
	else if(c == 'T'){
		return 3;
	}
	else{
		cerr << "HashTrie::codeChar - Unknown char " << c << "\n";
		exit(0);
//		return 0;
	}
}

HashTrie::HashTrie(){
	karp_rabin = NULL;
	arr_factors = NULL;
	compacted_text = NULL;
	select1_s = NULL;
	select1_b = NULL;
	select0_b = NULL;
	pi_inv = NULL;
}

HashTrie::HashTrie(const char *full_text, unsigned int _len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_factors, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *kr_factors, bool reverse){

	len_text = _len_text;
	karp_rabin = _karp_rabin;
	arr_factors = _arr_factors;
	compacted_text = kr_factors->compacted_text;
	select1_s = kr_factors->select1_s;
	select1_b = kr_factors->select1_b;
	select0_b = kr_factors->select0_b;
	pi_inv = kr_factors->pi_inv;
	
	build(full_text, len_text, factors_start, arr_factors, _karp_rabin, kr_factors, reverse);
}

HashTrie::~HashTrie(){
	karp_rabin = NULL;
	arr_factors = NULL;
	compacted_text = NULL;
	select1_s = NULL;
	select1_b = NULL;
	select0_b = NULL;
	pi_inv = NULL;
}

void HashTrie::build(const char *full_text, unsigned int _len_text, vector<unsigned int> &factors_start, int_vector<> *_arr_factors, KarpRabin *_karp_rabin, KarpRabinFactorsSuffixes *kr_factors, bool reverse){

	len_text = _len_text;
	karp_rabin = _karp_rabin;
	arr_factors = _arr_factors;
	compacted_text = kr_factors->compacted_text;
	select1_s = kr_factors->select1_s;
	select1_b = kr_factors->select1_b;
	select0_b = kr_factors->select0_b;
	pi_inv = kr_factors->pi_inv;
	
	cout << "HashTrie::build - Start (full text of " << len_text << ", " << factors_start.size() << " factors)\n";
	HashTrieNode local_root;
	if( reverse ){
		local_root.buildRev(full_text, len_text, factors_start, arr_factors, karp_rabin, kr_factors, 0, factors_start.size()-1, 0);
	}
	else{
		local_root.build(full_text, len_text, factors_start, arr_factors, karp_rabin, kr_factors, 0, factors_start.size()-1, 0);
	}
	compactData(local_root);
	
	unordered_map<unsigned int, pair<unsigned int, unsigned int>> marked_hash;
	if( reverse ){
		prepareHashMapRev(0, 0, marked_hash);
	}
	else{
		prepareHashMap(0, 0, marked_hash);
	}
	
	cout << "HashTrie::build - End\n";
}

void HashTrie::compactData(HashTrieNode &root_node){

	cout << "HashTrie::compactData - Start\n";
	
	// build en un nodo root local
	// Preparacion de arreglos para los datos agrupados
	// Iterar recursivamente por los nodos pasandole los arreglos
	// Desechar el root local

	unsigned int max_len = 0;
	unsigned int max_childs = 0;
	unsigned int max_height = 0;
	// including the root
	unsigned int n_nodes = 1 + root_node.totalChilds(max_len, max_childs, max_height, 0);
	
	cout << "HashTrie::compactData - Total Nodes: " << n_nodes << " (max_len: " << max_len << ", max_childs: " << max_childs << ", max_height: " << max_height << ")\n";
	
	positions_childs.resize(n_nodes);
	n_childs.resize(n_nodes);
	len_childs.resize(n_nodes);
	min_childs.resize(n_nodes);
	hash_childs.resize(n_nodes);
	first_childs.resize(n_nodes);
	
	unsigned int next_pos = 0;
	
	// El llamador guarda la posicion de inicio de la raiz
	// Es decir, cada nodo guarda los datos de sus hijos
	positions_childs[next_pos] = 1;
	cout << "HashTrie::compactData - positions_childs[" << next_pos << "]: 1\n";
	n_childs[next_pos] = root_node.childs_vector.size();
	len_childs[next_pos] = root_node.len;
	min_childs[next_pos] = root_node.min;
	hash_childs[next_pos] = root_node.hash;
//	first_childs[next_pos] = root_node.first;
	first_childs[next_pos] = codeChar(root_node.first);
	++next_pos;
	
	root_node.compactData(next_pos, positions_childs, n_childs, len_childs, min_childs, hash_childs, first_childs);
	
	sdsl::util::bit_compress(positions_childs);
	sdsl::util::bit_compress(n_childs);
	sdsl::util::bit_compress(len_childs);
	sdsl::util::bit_compress(min_childs);
	sdsl::util::bit_compress(hash_childs);
	sdsl::util::bit_compress(first_childs);
	
	float total_bits = 0;
	
	cout << "HashTrie::compactData - bits/node positions_childs: " << (8.0*size_in_bytes(positions_childs)/n_nodes) << "\n";
	cout << "HashTrie::compactData - bits/node n_childs: " << (8.0*size_in_bytes(n_childs)/n_nodes) << "\n";
	cout << "HashTrie::compactData - bits/node len_childs: " << (8.0*size_in_bytes(len_childs)/n_nodes) << "\n";
	cout << "HashTrie::compactData - bits/node min_childs: " << (8.0*size_in_bytes(min_childs)/n_nodes) << "\n";
	cout << "HashTrie::compactData - bits/node hash_childs: " << (8.0*size_in_bytes(hash_childs)/n_nodes) << "\n";
	cout << "HashTrie::compactData - bits/node first_childs: " << (8.0*size_in_bytes(first_childs)/n_nodes) << "\n";
	
	total_bits += (8.0*size_in_bytes(positions_childs)/n_nodes);
	total_bits += (8.0*size_in_bytes(n_childs)/n_nodes);
	total_bits += (8.0*size_in_bytes(len_childs)/n_nodes);
	total_bits += (8.0*size_in_bytes(min_childs)/n_nodes);
	total_bits += (8.0*size_in_bytes(hash_childs)/n_nodes);
//	total_bits += 8;
	total_bits += (8.0*size_in_bytes(first_childs)/n_nodes);
	cout << "HashTrie::compactData - bytes/node: " << total_bits/8.0 << "\n";
	
	// Debug
	for(unsigned int i = 0; i < n_nodes; ++i){
		cout << "HashTrie::compactData - node[" << i << "]: (" << positions_childs[i] << ", " << n_childs[i] << ", " << len_childs[i] << ", " << min_childs[i] << ", " << hash_childs[i] << ", " << decodeChar[ first_childs[i] ] << ")\n";
	}
	
	cout << "HashTrie::compactData - End\n";
}

void HashTrie::prepareHashMapRev(int node_pos, int path_len, unordered_map<unsigned int, pair<unsigned int, unsigned int>> &marked_hash){
	
	for(unsigned int i = 0; i < n_childs[node_pos]; ++i){
		int pos_child_abs = i + positions_childs[node_pos];
		
		int child_len = len_childs[pos_child_abs];
		child_len += path_len;
//		cout << "HashTrie::prepareHashMapRev - child_len: " << child_len << "\n";
		// Si es hoja, ajustar el largo
		unsigned int num_childs = n_childs[pos_child_abs];
		if( num_childs == 0 ){
			unsigned int min_factor_pos = (*arr_factors)[ min_childs[pos_child_abs] ];
			unsigned int pu = select1_b->operator()(min_factor_pos-1 + 1);
			unsigned int lu = select1_b->operator()(min_factor_pos-1 + 2) - pu;
			child_len = lu;
//			cout << "HashTrie::prepareHashMapRev - child_len corrected: " << child_len << "\n";
		}
		
//		int child_len = path_len + len_childs[pos_child_abs];
		char c = decodeChar[ first_childs[pos_child_abs] ];
		cout << "HashTrie::prepareHashMapRev - Child[" << pos_child_abs << "]: " << c << ", child_len " << child_len << "\n";
		prepareHashMapRev(pos_child_abs, child_len, marked_hash);
	}
	
}

void HashTrie::prepareHashMap(int node_pos, int path_len, unordered_map<unsigned int, pair<unsigned int, unsigned int>> &marked_hash){
	
	for(unsigned int i = 0; i < n_childs[node_pos]; ++i){
		int pos_child_abs = i + positions_childs[node_pos];
		
		int child_len = len_childs[pos_child_abs];
		child_len += path_len;
		// Si es hoja, ajustar el largo
		unsigned int num_childs = n_childs[pos_child_abs];
		if( num_childs == 0 ){
			unsigned int min_factor_pos = (*arr_factors)[ min_childs[pos_child_abs] ];
			unsigned int pu = select1_b->operator()(min_factor_pos + 1);
			child_len = len_text - pu;
		}
		
//		int child_len = path_len + len_childs[pos_child_abs];
		char c = decodeChar[ first_childs[pos_child_abs] ];
		cout << "HashTrie::prepareHashMap - Child[" << pos_child_abs << "]: " << c << ", child_len " << child_len << "\n";
		prepareHashMap(pos_child_abs, child_len, marked_hash);
	}
	
}

void HashTrie::print(){
	printInternal(0, 0);
	cout << "Root: " << n_childs[0] << " childs\n";
	cout << "Lengths: ";
	for(unsigned int i = 0; i < n_childs[0]; ++i){
		cout << len_childs[ i + positions_childs[0] ] << " - ";
	}
	cout << "\n";
}

void HashTrie::printInternal(unsigned int node_pos, unsigned int level){
	for(unsigned int i = 0; i < level; ++i){
		cout << "- ";
	}
	cout << "hash: " << hash_childs[node_pos] << ", len: " << len_childs[node_pos] << " , range [" << min_childs[node_pos] << "]\n";
	for(unsigned int i = 0; i < n_childs[node_pos]; ++i){
		cout << decodeChar[ first_childs[ i + positions_childs[node_pos] ] ] << " ";
		printInternal(i + positions_childs[node_pos], level+1);
	}
}

void HashTrie::printSize(){
	unsigned int max_len = 0;
	unsigned int max_childs = 0;
	unsigned int max_height = 0;
	unsigned int total_childs = totalChilds();
	cout << "HashTrie::printSize - totalChilds: " << total_childs << ", max_len: " << max_len << ", max_childs: " << max_childs << ", max_height: " << max_height << "\n";
}

pair<unsigned int, unsigned int> HashTrie::getRange(vector<unsigned long long> &kr_pat_vector, unsigned int pos, const string &pattern, bool use_hash){
	if( use_hash ){
		return getRangeInternal(0, kr_pat_vector, pos, 0, karp_rabin, arr_factors->size() - 1, pattern);
	}
	else{
		return getRangeInternalNoHash(0, kr_pat_vector, pos, 0, karp_rabin, arr_factors->size() - 1, pattern);
	}
}

pair<unsigned int, unsigned int> HashTrie::getRangeRev(vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, const string &pattern_rev, bool use_hash){
	if( use_hash ){
		return getRangeRevInternal(0, kr_pat_rev_vector, pos, 0, karp_rabin, arr_factors->size() - 1, pattern_rev);
	}
	else{
		return getRangeRevInternalNoHash(0, kr_pat_rev_vector, pos, 0, karp_rabin, arr_factors->size() - 1, pattern_rev);
	}
}

unsigned int HashTrie::findChildInternal(unsigned int node_pos, char c){
	for(unsigned int i = 0; i < n_childs[node_pos]; ++i){
		if( c == decodeChar[ first_childs[ i + positions_childs[node_pos] ] ] ){
			return i;
		}
	}
	return NOT_FOUND;
}

pair<unsigned int, unsigned int> HashTrie::getRangeInternal(unsigned int node_pos, vector<unsigned long long> &kr_pat_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, unsigned int cur_max, const string &pattern){
	
	cout << "HashTrie::getRangeInternal - Start (prefixes: " << kr_pat_vector.size() << ", pos: " << pos << ", processed: " << processed << ", node_pos: " << node_pos << ")\n";
	
	unsigned int min = min_childs[node_pos];
	
	if( pos + processed >= kr_pat_vector.size() ){
		cout << "HashTrie::getRangeInternal - [" << min << ", " << cur_max << "]\n";
		return pair<unsigned int, unsigned int>(min, cur_max);
	}
	
	unsigned long long hash_pat = 0;
	unsigned int child_len = 0;
	unsigned int pat_len = kr_pat_vector.size() - pos - processed;
	char first_char_pat = pattern[pos + processed];
//	cout << "HashTrie::getRangeInternal - pat_len: " << pat_len << "\n";
	string pat = pattern.substr(pos + processed, pat_len);
	cout << "HashTrie::getRangeInternal - pat: " << pat << " (first_char_pat: " << first_char_pat << ")\n";
	
	unsigned int pos_child = findChildInternal(node_pos, first_char_pat);
	if( pos_child != NOT_FOUND ){
		unsigned int pos_child_abs = pos_child + positions_childs[node_pos];
		child_len = len_childs[pos_child_abs];
		
		// Si es hoja, ajustar el largo
		unsigned int num_childs = n_childs[pos_child_abs];
		if( num_childs == 0 ){
			unsigned int min_factor_pos = (*arr_factors)[ min_childs[pos_child_abs] ];
			unsigned int pu = select1_b->operator()(min_factor_pos + 1);
			child_len = len_text - pu - processed;
		}
		
		// Ajuste a cur_max
		if( pos_child < n_childs[node_pos] - 1 ){
			cur_max = min_childs[pos_child_abs+1] - 1;
		}
		
		if( child_len <= pat_len ){
//			cout << "HashTrie::getRangeInternal - Case 1, child_len: " << child_len << ", pos_child: " << pos_child << "\n";
			hash_pat = karp_rabin->subtract_prefix(kr_pat_vector[pos + processed + child_len - 1], kr_pat_vector[pos + processed - 1], child_len);
//			string pat_cut = pattern.substr(pos + processed, child_len);
//			cout << "HashTrie::getRangeInternal - pat_cut: " << pat_cut << " hash_pat: " << hash_pat << " / " << karp_rabin->hash(pat_cut) << " (processed: " << processed << ")\n";
			if( hash_pat == hash_childs[pos_child_abs] ){
//				cout << "HashTrie::getRangeInternal - Child found -> [" << min_childs[pos_child_abs] << ", " << cur_max << "]\n";
				return getRangeInternal(pos_child_abs, kr_pat_vector, pos, processed + child_len, karp_rabin, cur_max, pattern);
			}
		}
		else{
//			cout << "HashTrie::getRangeInternal - Case 2, child_len: " << child_len << "\n";
			
			unsigned int min_factor_pos = (*arr_factors)[ min_childs[pos_child_abs] ];
			
			// TEST HASH
			string test_text = "";
			FactorsIteratorCompacted it_y(min_factor_pos, arr_factors->size(), select1_s, select1_b, select0_b, pi_inv, compacted_text, len_text);
			for(unsigned int i = 0; i < processed; ++i){
				it_y.next();
			}
			for(unsigned int i = 0; i < pat_len; ++i){
				test_text += it_y.next();
			}
			unsigned long long hash = karp_rabin->hash(test_text);
//			cout << "HashTrie::getRangeInternal - Test Hash: " << hash_test << " / " << hash << "\n";
//			if( hash_test != hash ){
//				cerr << "HashTrie::getRangeInternal - Error\n";
//				exit(0);
//			}
			
			hash_pat = karp_rabin->subtract_prefix(kr_pat_vector[kr_pat_vector.size() - 1], kr_pat_vector[pos + processed - 1], kr_pat_vector.size() - pos - processed);
//			cout << "HashTrie::getRangeInternal - hash: " << hash << ", hash_pat: " << hash_pat << "\n";
			if( hash == hash_pat ){
				cout << "HashTrie::getRangeInternal - Child found -> [" << min_childs[pos_child_abs] << ", " << cur_max << "]\n";
				return pair<unsigned int, unsigned int>(min_childs[pos_child_abs], cur_max);
			}
		}
	}
	
	cout << "HashTrie::getRangeInternal - Pattern NOT found\n";
	return pair<unsigned int, unsigned int>((unsigned int)(-1), (unsigned int)(-1));
}

pair<unsigned int, unsigned int> HashTrie::getRangeInternalNoHash(unsigned int node_pos, vector<unsigned long long> &kr_pat_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, unsigned int cur_max, const string &pattern){
	
//	cout << "HashTrie::getRangeInternalNoHash - Start (prefixes: " << kr_pat_vector.size() << ", pos: " << pos << ", processed: " << processed << ", node_pos: " << node_pos << ")\n";
	
	unsigned int min = min_childs[node_pos];
	
	if( pos + processed >= kr_pat_vector.size() ){
//		cout << "HashTrie::getRangeInternalNoHash - [" << min << ", " << cur_max << "]\n";
		return pair<unsigned int, unsigned int>(min, cur_max);
	}
	
	unsigned int child_len = 0;
	unsigned int pat_len = kr_pat_vector.size() - pos - processed;
	char first_char_pat = pattern[pos + processed];
//	cout << "HashTrie::getRangeInternalNoHash - pat_len: " << pat_len << "\n";
	string pat = pattern.substr(pos + processed, pat_len);
//	cout << "HashTrie::getRangeInternalNoHash - pat: " << pat << " (first_char_pat: " << first_char_pat << ")\n";
	
	unsigned int pos_child = findChildInternal(node_pos, first_char_pat);
	if( pos_child != NOT_FOUND ){
		unsigned int pos_child_abs = pos_child + positions_childs[node_pos];
		child_len = len_childs[pos_child_abs];
		
		// Si es hoja, ajustar el largo
		unsigned int num_childs = n_childs[pos_child_abs];
		if( num_childs == 0 ){
			unsigned int min_factor_pos = (*arr_factors)[ min_childs[pos_child_abs] ];
			unsigned int pu = select1_b->operator()(min_factor_pos + 1);
			child_len = len_text - pu - processed;
		}
		
		// Ajuste a cur_max
		if( pos_child < n_childs[node_pos] - 1 ){
			cur_max = min_childs[pos_child_abs+1] - 1;
		}
		
		if( child_len <= pat_len ){
//			cout << "HashTrie::getRangeInternalNoHash - Case 1, child_len: " << child_len << ", pos_child: " << pos_child << "\n";
			return getRangeInternalNoHash(pos_child_abs, kr_pat_vector, pos, processed + child_len, karp_rabin, cur_max, pattern);
		}
		else{
//			cout << "HashTrie::getRangeInternalNoHashgetRange - Case 2, child_len: " << child_len << "\n";
			return pair<unsigned int, unsigned int>(min_childs[pos_child_abs], cur_max);
		}
	}
	
//	cout << "HashTrie::getRangeInternalNoHash - Pattern NOT found\n";
	return pair<unsigned int, unsigned int>((unsigned int)(-1), (unsigned int)(-1));
}

pair<unsigned int, unsigned int> HashTrie::getRangeRevInternal(unsigned int node_pos, vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, unsigned int cur_max, const string &pattern_rev){
	
	cout << "HashTrie::getRangeRevInternal - Start (prefixes: " << kr_pat_rev_vector.size() << ", pos: " << pos << ", processed: " << processed << ", pattern_rev: " << pattern_rev << ")\n";
	
	unsigned int min = min_childs[node_pos];
	
	if( processed >= pos ){
		cout << "HashTrie::getRangeRevInternal - [" << min << ", " << cur_max << "]\n";
		return pair<unsigned int, unsigned int>(min, cur_max);
	}
	
	unsigned long long hash_pat = 0;
	unsigned int child_len = 0;
	unsigned int pat_len = pos - processed;
	char first_char_pat = pattern_rev[kr_pat_rev_vector.size() - 1 - pos + processed];
//	cout << "HashTrie::getRangeRevInternal - pat = pattern_rev.substr(" << (kr_pat_rev_vector.size() - 1 - pos + processed) << ", " << pat_len << ");\n";
	string pat = pattern_rev.substr(kr_pat_rev_vector.size() - 1 - pos + processed, pat_len);
	cout << "HashTrie::getRangeRevInternal - pat: " << pat << ", first_char_pat: " << first_char_pat << "\n";
	
	unsigned int pos_child = findChildInternal(node_pos, first_char_pat);
	if( pos_child != NOT_FOUND ){
		unsigned int pos_child_abs = pos_child + positions_childs[node_pos];
		child_len = len_childs[pos_child_abs];
		
		// Si es hoja, ajustar el largo
		unsigned int num_childs = n_childs[pos_child_abs];
		if( num_childs == 0 ){
			unsigned int min_factor_pos = (*arr_factors)[ min_childs[pos_child_abs] ];
			unsigned int pu = select1_b->operator()(min_factor_pos-1 + 1);
			unsigned int lu = select1_b->operator()(min_factor_pos-1 + 2) - pu;
			child_len = lu - processed;
		}
		
		// Ajuste a cur_max
//		cout << "HashTrie::getRangeRevInternal - Adjusting cur_max (pos_child: " << pos_child << " / " << n_childs[node_pos] << ")\n";
		if( pos_child < n_childs[node_pos] - 1 ){
			cur_max = min_childs[pos_child_abs+1] - 1;
//			cout << "HashTrie::getRangeRevInternal - cur_max: " << cur_max << ")\n";
		}
		
		if( child_len <= pat_len ){
//			cout << "HashTrie::getRangeRevInternal - Case 1, child_len: " << child_len << "\n";
			hash_pat = karp_rabin->hash(pattern_rev.c_str() + pattern_rev.length() - pos + processed, child_len);
			string pat_cut = pattern_rev.substr(kr_pat_rev_vector.size() - 1 - pos + processed, child_len);
//			cout << "HashTrie::getRangeRevInternal - pat_cut: " << pat_cut << " hash_pat: " << hash_pat << " / " << karp_rabin->hash(pat_cut) << " (processed: " << processed << ")\n";
			if( hash_pat == hash_childs[pos_child_abs] ){
//				cout << "HashTrie::getRangeRevInternal - Child found -> [" << min_childs[pos_child_abs] << ", " << cur_max << "]\n";
				return getRangeRevInternal(pos_child_abs, kr_pat_rev_vector, pos, processed + child_len, karp_rabin, cur_max, pattern_rev);
			}
		}
		else{
//			cout << "HashTrie::getRangeRevInternal - Case 2, child_len: " << child_len << "\n";
			
			string test_text = "";
			unsigned int min_factor_pos = (*arr_factors)[ min_childs[pos_child_abs] ];
			
			if( min_factor_pos > 0 ){
				
				unsigned int cur_pi = (*pi_inv)[min_factor_pos-1];
				unsigned int tu = select1_s->operator()(cur_pi + 1) - cur_pi;
				unsigned int pu = select1_b->operator()(min_factor_pos-1 + 1);
				unsigned int lu = select1_b->operator()(min_factor_pos-1 + 2) - pu;
				
				if(processed < lu){
					unsigned int len = lu - processed;
					if( child_len < len ){
						len = child_len;
					}
					if( pat_len < len ){
						len = pat_len;
					}
//					cout << "HashTrie::getRange - Adding " << len << " chars\n";
					for(unsigned int i = 0; i < len; ++i){
						test_text += compacted_text->at(tu + lu - processed - i - 1);
					}
				}
			}
			
			unsigned long long hash = karp_rabin->hash(test_text.c_str(), pat_len);
			hash_pat = karp_rabin->hash(pattern_rev.c_str() + pattern_rev.length() - pos + processed, pat_len);
//			cout << "HashTrie::getRange - hash: " << hash << ", hash_pat: " << hash_pat << "\n";
			if( hash == hash_pat ){
				cout << "HashTrie::getRange - Child found -> [" << min_childs[pos_child_abs] << ", " << cur_max << "]\n";
				return pair<unsigned int, unsigned int>(min_childs[pos_child_abs], cur_max);
			}
			
		}
	}
	
	cout << "HashTrie::getRange - Pattern NOT found\n";
	return pair<unsigned int, unsigned int>((unsigned int)(-1), (unsigned int)(-1));
}

pair<unsigned int, unsigned int> HashTrie::getRangeRevInternalNoHash(unsigned int node_pos, vector<unsigned long long> &kr_pat_rev_vector, unsigned int pos, unsigned int processed, KarpRabin *karp_rabin, unsigned int cur_max, const string &pattern_rev){
	
//	cout << "HashTrie::getRangeRevInternalNoHash - Start (prefixes: " << kr_pat_rev_vector.size() << ", pos: " << pos << ", processed: " << processed << ", pattern_rev: " << pattern_rev << ")\n";
	
	unsigned int min = min_childs[node_pos];
	
	if( processed >= pos ){
//		cout << "HashTrie::getRangeRevInternalNoHash - [" << min << ", " << cur_max << "]\n";
		return pair<unsigned int, unsigned int>(min, cur_max);
	}
	
	unsigned int child_len = 0;
	unsigned int pat_len = pos - processed;
	char first_char_pat = pattern_rev[kr_pat_rev_vector.size() - 1 - pos + processed];
//	cout << "HashTrie::getRangeRevInternalNoHash - pat = pattern_rev.substr(" << (kr_pat_rev_vector.size() - 1 - pos + processed) << ", " << pat_len << ");\n";
	string pat = pattern_rev.substr(kr_pat_rev_vector.size() - 1 - pos + processed, pat_len);
//	cout << "HashTrie::getRangeRevInternalNoHash - pat: " << pat << ", first_char_pat: " << first_char_pat << "\n";
	
	unsigned int pos_child = findChildInternal(node_pos, first_char_pat);
	if( pos_child != NOT_FOUND ){
		unsigned int pos_child_abs = pos_child + positions_childs[node_pos];
		child_len = len_childs[pos_child_abs];
		
		// Si es hoja, ajustar el largo
		unsigned int num_childs = n_childs[pos_child_abs];
		if( num_childs == 0 ){
			unsigned int min_factor_pos = (*arr_factors)[ min_childs[pos_child_abs] ];
			unsigned int pu = select1_b->operator()(min_factor_pos-1 + 1);
			unsigned int lu = select1_b->operator()(min_factor_pos-1 + 2) - pu;
			child_len = lu - processed;
		}
		
		// Ajuste a cur_max
//		cout << "HashTrie::getRangeRevInternalNoHash - Adjusting cur_max (pos_child: " << pos_child << " / " << n_childs[node_pos] << ")\n";
		if( pos_child < n_childs[node_pos] - 1 ){
			cur_max = min_childs[pos_child_abs+1] - 1;
//			cout << "HashTrie::getRangeRevInternalNoHash - cur_max: " << cur_max << ")\n";
		}
		
		if( child_len <= pat_len ){
//			cout << "HashTrie::getRangeRevInternalNoHash - Case 1, child_len: " << child_len << "\n";
			return getRangeRevInternalNoHash(pos_child_abs, kr_pat_rev_vector, pos, processed + child_len, karp_rabin, cur_max, pattern_rev);
		}
		else{
//			cout << "HashTrie::getRangeRevInternalNoHash - Case 2, child_len: " << child_len << "\n";
			return pair<unsigned int, unsigned int>(min_childs[pos_child_abs], cur_max);
		}
	}
	
//	cout << "HashTrie::getRangeRevInternalNoHash - Pattern NOT found\n";
	return pair<unsigned int, unsigned int>((unsigned int)(-1), (unsigned int)(-1));
}

void HashTrie::save(const string &file){
	cout << "HashTrie::save - Start (" << file << ")\n";
	
	string pos_file = file + ".pos";
	store_to_file(positions_childs, pos_file);
	
	string nc_file = file + ".nc";
	store_to_file(n_childs, nc_file);
	
	string len_file = file + ".len";
	store_to_file(len_childs, len_file);
	
	string min_file = file + ".min";
	store_to_file(min_childs, min_file);
	
	string hash_file = file + ".hash";
	store_to_file(hash_childs, hash_file);
	
	string first_file = file + ".first";
	store_to_file(first_childs, first_file);
	
	cout << "HashTrie::save - End\n";
}

void HashTrie::load(unsigned int _len_text, KarpRabin *_karp_rabin, CompactedText *_compacted_text, bits_s_type::select_1_type *_select1_s, bits_b_type::select_1_type *_select1_b, bits_b_type::select_0_type *_select0_b, int_vector<> *_pi_inv, int_vector<> *_arr_factors, const string &file){
	cout << "HashTrie::load - Start (" << file << ")\n";
	
	len_text = _len_text;
	karp_rabin = _karp_rabin;
	arr_factors = _arr_factors;
	compacted_text = _compacted_text;
	select1_s = _select1_s;
	select1_b = _select1_b;
	select0_b = _select0_b;
	pi_inv = _pi_inv;
	
	string pos_file = file + ".pos";
	load_from_file(positions_childs, pos_file);
	
	string nc_file = file + ".nc";
	load_from_file(n_childs, nc_file);
	
	string len_file = file + ".len";
	load_from_file(len_childs, len_file);
	
	string min_file = file + ".min";
	load_from_file(min_childs, min_file);
	
	string hash_file = file + ".hash";
	load_from_file(hash_childs, hash_file);
	
	string first_file = file + ".first";
	load_from_file(first_childs, first_file);
	
	cout << "HashTrie::load - End\n";
}

unsigned int HashTrie::getSizeBytes(){
	
	unsigned int total = 0;
	total += size_in_bytes(positions_childs);
	total += size_in_bytes(n_childs);
	total += size_in_bytes(len_childs);
	total += size_in_bytes(min_childs);
	total += size_in_bytes(hash_childs);
//	total += first_childs.size();
	total += size_in_bytes(first_childs);
	
	return total;
}






