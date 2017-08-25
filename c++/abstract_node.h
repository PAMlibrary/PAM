#pragma once
#include "pbbs-include/list_allocator.h"
#include "defs.h"

// Definitions in this file are independent of balance criteria beyond
// maintaining an abstract "rank"
template<class Node>
class avl_tree;

template<class Node>
using tree_interface = avl_tree<Node>;

template<class K, class V, class AugmOp, class Compare>
class node {
public:
    using key_type    = K ;
    using value_type  = V ;
    using key_compare = Compare;
    using entry_type  = std::pair<K, V>;
    using node_type   = node<K, V, AugmOp, Compare>;    
    using allocator   = list_allocator<node_type>;
    using tree_type   = tree_interface<node_type>;
    using aug_type    = typename AugmOp::aug_t;
    using aug_class   = AugmOp;

    node(const entry_type&, node_type* lc, node_type* rc, bool do_update = 1);
    node(const entry_type&, bool initialize = 1);
    node() : ref_cnt(1) {};

    const entry_type get_entry() const { return entry_type(key,value); }
    const K& get_key() const { return key; }
    const V& get_value() const { return value; }
        
    void set_value(const V _value) { value = _value; }
    void set_entry(const entry_type kv) { key = kv.first; value = kv.second;}
	void set_aug(const aug_type ag) { aug_val = ag;}
            
    const aug_type get_aug_val() const { return aug_val; }
    static void* operator new(size_t size) { return allocator::alloc(); }

    inline node_type* copy();
    inline void update();
    static inline node_type* connect_with_update(node_type* t1, node_type* t2, node_type* k) {
		k->lc = t1, k->rc = t2;
		k->update();
		return k;
    }
	
    static inline node_type* connect_without_update(node_type* t1, node_type* t2, node_type* k) {
		k->lc = t1, k->rc = t2;
		return k;
    }
    inline void collect();

    // ordering is designed to save space
    node_type* lc;  // left child
    node_type* rc;  // right child
    K key;
    V value;
    aug_type aug_val; // augmented value
    unsigned char rank; // safe as a height, but not a weight
    tree_size_t node_cnt; // subtree size
    tree_size_t ref_cnt; // reference count
};

template<class K, class V, class AugmOp, class Compare>
inline tree_size_t get_rank(const node<K, V, AugmOp, Compare>* t) {
  return t ? t->rank : 0;
}

template<class K, class V, class AugmOp, class Compare>
typename AugmOp::aug_t get_aug(const node<K, V, AugmOp, Compare>* t) {
  return t ? t->aug_val : AugmOp::get_empty();
}

template<class K, class V, class AugmOp, class Compare>
inline void node<K, V, AugmOp, Compare>::collect() {
    get_key().~key_type();
    get_value().~value_type();
    allocator::free(this);
}

template<class K, class V, class AugmOp, class Compare>
inline node<K, V, AugmOp, Compare>* node<K, V, AugmOp, Compare>::copy() {
    node_type* ret = new node_type(get_entry(), lc, rc, 0);
    ret->rank = rank;
    ret->aug_val = aug_val;
    increase(lc);
    increase(rc);
    return ret;
}

template<class K, class V, class AugmOp, class Compare>
inline void node<K, V, AugmOp, Compare>::update() {
    aug_val = AugmOp::from_entry(get_key(), get_value());
    //if (lc) aug_val = AugmOp::combine(aug_val, lc->aug_val);
	if (lc) aug_val = AugmOp::combine(lc->aug_val, aug_val);// else aug_val = AugmOp::combine(AugmOp::get_empty(), aug_val);
    if (rc) aug_val = AugmOp::combine(aug_val, rc->aug_val);// else aug_val = AugmOp::combine(aug_val, AugmOp::get_empty());

    rank = tree_type::combine_ranks(get_rank(lc), get_rank(rc));
    node_cnt = 1 + get_node_count(lc) + get_node_count(rc);
}

template<class K, class V, class AugmOp, class Compare>
    node<K, V, AugmOp, Compare>::node(const entry_type& kv, node_type* left, node_type* right,
                 bool do_update) {
    set_entry(kv);
    ref_cnt = 1;
    lc = left;
    rc = right;
    if (do_update) update();
}

template<class K, class V, class AugmOp, class Compare>
  node<K, V, AugmOp, Compare>::node(const entry_type& kv, bool initialize) {
    set_entry(kv);
    ref_cnt = 1;
    if (initialize) {
      lc = rc = NULL;
      aug_val = AugmOp::from_entry(get_key(),get_value());
      rank = tree_type::singleton_rank();
      node_cnt = 1;
    }
}

template<class K, class V, class AugmOp, class Compare>
int get_real_height(node<K, V, AugmOp, Compare>* r) {
	if (r == NULL) return 0;
	else return max(get_real_height(r->lc), get_real_height(r->rc))+1;
}


