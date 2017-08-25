#pragma once

#include "common.h"
#include "abstract_node.h"
#include "pbbs-include/binary_search.h"
#include "pbbs-include/sample_sort.h"
#include "pbbs-include/get_time.h"

template<class Node>
struct tree_operations {

  using K           = typename Node::key_type;
  using V           = typename Node::value_type;
  using E           = typename Node::entry_type;
  using tree_type   = typename Node::tree_type;
  using key_compare = typename Node::key_compare;
  using aug_type    = typename Node::aug_type;
  using aug_class   = typename Node::aug_class;
  
  static Node* t_join3(Node* b1, Node* b2, Node* k) {
      return tree_type::t_join(b1, b2, k);
  }

  static bool do_parallel(size_t n, size_t m) {
    if (m > n) std::swap(n,m);
    return (n > 256 && m > 12);
  }

  // fork-join parallel call
  template <class RT, class Lf, class Rf>
  static std::pair<RT,RT> fork(bool do_parallel, Lf left, Rf right) {
    if (do_parallel) {
      RT r = cilk_spawn right();
      RT l = left();
      cilk_sync;
      return std::pair<RT,RT>(l,r);
    } else {
      RT l = left(); 
      RT r = right();
      return std::make_pair(l,r);
    }
  }

  // fork-join parallel call
  template <class Lf, class Rf>
  static void fork_no_result(bool do_parallel, Lf left, Rf right) {
    if (do_parallel) {cilk_spawn right(); left(); cilk_sync; }
    else {left(); right();}
  }

  struct split_info {
      split_info(Node* first, Node* second, bool removed)
          : first(first), second(second), removed(removed) {};

      Node* first;
      Node* second;
      V value;
      bool removed;
  };

  static split_info t_split(Node* bst, const K& e) {
      if (!bst) return split_info(NULL, NULL, false);
      
      Node* lsub = bst->lc, *rsub = bst->rc; 
      const K key = bst->get_key();
      const V value = bst->get_value();

      Node* join = copy_if_needed(bst);
      
      if (comp(key, e)) {
          split_info bstpair = t_split(rsub, e);
          bstpair.first = t_join3(lsub, bstpair.first, join);
          return bstpair;
      }
      else if (comp(e, key)) {
          split_info bstpair = t_split(lsub, e);
          bstpair.second = t_join3(bstpair.second, rsub, join);
          return bstpair;
      } else {
          split_info ret(lsub, rsub, true);
          ret.value = value;
          decrease(join);
          return ret;
      }
  }

  // could be made more efficient, especially for small ranges
  static Node* t_range(Node* b, const K& low, const K& high) {
      split_info left_split = t_split(b, low);
      decrease_recursive(left_split.first);
  
      split_info right_split = t_split(left_split.second, high);
      decrease_recursive(right_split.second);

      Node* ret = right_split.first;
      if (left_split.removed)
            ret = t_insert(ret, E(low,left_split.value));
      if (right_split.removed)
            ret = t_insert(ret, E(high,right_split.value));
  
      return ret;
  }

  static Node* t_range_root(Node* b, const K& key_left, const K& key_right) {
    while (b) {
      if (comp(key_right, b->get_key())) { b = b->lc; continue; } 
      if (comp(b->get_key(), key_left)) { b = b->rc; continue; }
      break;
    }
    return b;
  }

  static Node* t_left(Node* bst, const K& e) {
    if (!bst) return NULL;
    if (comp(e, bst->get_key())) return t_left(bst->lc, e);
    Node* join = new Node(bst->get_entry(), false);
    increase(bst->lc);
    return t_join3(bst->lc, t_left(bst->rc, e), join);
  }

  static Node* t_right(Node* bst, const K& e) {
    if (!bst) return NULL;
    if (comp(bst->get_key(), e)) return t_right(bst->rc, e);
    Node* join = new Node(bst->get_entry(), false);
    increase(bst->rc);
    return t_join3(t_right(bst->lc, e), bst->rc, join);
  }

  static Node* t_range_fast(Node* b, const K& low, const K& high) {
    Node* r = t_range_root(b, low, high);
    if (!r) return NULL;
    Node* join = new Node(r->get_entry(), false);
    return t_join3(t_right(r->lc, low), t_left(r->rc, high), join);
  }

  static Node* t_join2(Node* b1, Node* b2) {
      if (!b1) return b2;
      if (!b2) return b1;
      
      if (b1->rank > b2->rank) {
          Node* join = copy_if_needed(b1);
          return t_join3(join->lc, t_join2(join->rc, b2), join);
       } else {
          Node* join = copy_if_needed(b2);
          return t_join3(t_join2(b1, join->lc), join->rc, join);
      }
  }

  template <class BinaryOp> 
  static Node* t_union(Node* b1, Node* b2, const BinaryOp& op) {
      if (!b1) return b2;
      if (!b2) return b1;

      Node* join = copy_if_needed(b2);

      split_info bsts = t_split(b1, join->get_key());
      if (bsts.removed)
          join->set_value(op(bsts.value, join->get_value()));

      auto P = fork<Node*>(do_parallel(get_node_count(b1), get_node_count(b2)),
        [&] () {return t_union(bsts.first, join->lc, op);},
        [&] () {return t_union(bsts.second, join->rc, op);});
      
      return t_join3(P.first, P.second, join);
  }

  template <class BinaryOp>
  static Node* t_intersect(Node* b1, Node* b2, const BinaryOp& op) {
      if (!b1) { 
          decrease_recursive(b2); 
          return NULL; 
      }
      if (!b2) {
          decrease_recursive(b1);
          return NULL;
      }

      Node* join = copy_if_needed(b2);

      split_info bsts = t_split(b1, join->get_key());

      auto P = fork<Node*>(do_parallel(get_node_count(b1), get_node_count(b2)),
        [&]() {return t_intersect(bsts.first, join->lc, op);},
        [&]() {return t_intersect(bsts.second, join->rc, op);}
      );

      if (bsts.removed) {
          join->set_value(op(bsts.value, join->get_value()));
          return t_join3(P.first, P.second, join);
      } else {
          decrease(join);
          return t_join2(P.first, P.second);
      }
  }

    
  static Node* t_difference(Node* b1, Node* b2) {
    if (!b1) {
      decrease_recursive(b2);
      return NULL;
    }
    if (!b2) return b1;

    Node* join = copy_if_needed(b1);
      
    split_info bsts = t_split(b2, b1->get_key());

    auto P = fork<Node*>(do_parallel(get_node_count(b1), get_node_count(b2)),
      [&]() {return t_difference(join->lc, bsts.first);},
      [&]() {return t_difference(join->rc, bsts.second);});
      
    if (bsts.removed) {
      decrease(join);
      return t_join2(P.first, P.second);
    } else {
      return t_join3(P.first, P.second, join);
    }   
  }

  template<class Func>
  static Node* t_filter(Node* b, const Func& f) {
      if (!b) return NULL;

      Node* join = copy_if_needed(b);

      size_t mn = get_node_count(b);
      auto P = fork<Node*>(mn >= node_limit,
        [&]() {return t_filter(join->lc, f);},
        [&]() {return t_filter(join->rc, f);});

      if (f(b->get_entry())) {
          return t_join3(P.first, P.second, join);
      } else {
          decrease(join);
          return t_join2(P.first, P.second);
      }
  }

  static Node* t_insert(Node* b, const E& e){
      if (!b) return new Node(e);

      Node* tmp = copy_if_needed(b);

      if (comp(b->get_key(), e.first) )
          return t_join3(tmp->lc,t_insert(tmp->rc, e), tmp);
      else if (comp(e.first, b->get_key()) )
          return t_join3(t_insert(tmp->lc, e), tmp->rc, tmp);
      else {
          tmp->set_value(e.second);
          tmp->update();
          return tmp;
      }
  }
  
  template <class Func>
  static Node* t_insert(Node* b, const E& e, const Func& f){
      if (!b) return new Node(e);
      Node* tmp = copy_if_needed(b);

      if (comp(b->get_key(), e.first) )
	return t_join3(tmp->lc,t_insert(tmp->rc, e, f), tmp);
      else if (comp(e.first, b->get_key()) )
	return t_join3(t_insert(tmp->lc, e, f), tmp->rc, tmp);
      else {
          tmp->set_value(f(tmp->get_value(),e.second));
          tmp->update();
          return tmp;
      }
  }
  
  template <class Func>
  static Node* t_insert_lazy(Node* b, const E& e, const Func& f){
      if (!b) return new Node(e);

      Node* tmp = copy_if_needed(b);

      if (comp(b->get_key(), e.first) ) {
		  Node* right = t_insert_lazy(tmp->rc, e);
		  Node* left = tmp->lc;
		  if (Node::tree_type::are_balanced(left, right)) {
			  tmp->set_aug(aug_class::combine(tmp->get_aug_val(), aug_class::from_entry(e.first, e.second)));
			  return Node::connect_without_update(left, right, tmp);
		  } else {
			  return t_join3(left, right, tmp);
		  }
	  }
      else if (comp(e.first, b->get_key()) ) {
		  Node* right = tmp->rc;
		  Node* left = t_insert_lazy(tmp->lc, e);
		  if (Node::tree_type::are_balanced(left, right)) {
				tmp->set_aug(aug_class::combine(tmp->get_aug_val(), aug_class::from_entry(e.first, e.second)));
				return Node::connect_without_update(left, right, tmp);
		  } else {
			return t_join3(left, right, tmp);
		  }
	  }
      else {
          tmp->set_value(f(tmp.get_value(),e.second));
          tmp->update();
          return tmp;
      }
  }
  
  static Node* t_insert_lazy(Node* b, const E& e){
      if (!b) return new Node(e);

      Node* tmp = copy_if_needed(b);

      if (comp(b->get_key(), e.first) ) {
		  Node* right = t_insert_lazy(tmp->rc, e);
		  Node* left = tmp->lc;
		  if (Node::tree_type::are_balanced(left, right)) {
			  aug_type new_aug = aug_class::combine(tmp->get_aug_val(), aug_class::from_entry(e.first, e.second));
			  return Node::connect_without_update(left, right, tmp);
		  } else {
			  return t_join3(left, right, tmp);
		  }
	  }
      else if (comp(e.first, b->get_key()) ) {
		  Node* right = tmp->rc;
		  Node* left = t_insert_lazy(tmp->lc, e);
		  if (Node::tree_type::are_balanced(left, right)) {
				aug_type new_aug = aug_class::combine(tmp->get_aug_val(), aug_class::from_entry(e.first, e.second));
				return Node::connect_without_update(left, right, tmp);
		  } else {
			return t_join3(left, right, tmp);
		  }
	  }
      else {
          tmp->set_value(e.second);
          tmp->update();
          return tmp;
      }
  }

  inline static Node* t_delete(Node* b, const K& k) {
      if (!b) return NULL;
        
      Node* tmp = copy_if_needed(b);

      if (comp(b->get_key(), k)) 
          return t_join3(tmp->lc, t_delete(tmp->rc, k), tmp);
      else if (comp(k, b->get_key())) 
          return t_join3(t_delete(tmp->lc, k), tmp->rc, tmp);
      else {
          Node* r = t_join2(tmp->lc, tmp->rc);
          decrease(tmp);
          return r;
      }
  }

    // Assumes the input is sorted and there are no duplicate keys
  static Node* t_from_sorted_array(E* A, size_t n) {
      if (n <= 0) return NULL;
      if (n == 1) return new Node(A[0]);

      size_t mid = n/2;
      Node* m = new Node(A[mid]);

      auto P = fork<Node*>(n >= node_limit,
      [&]() {return t_from_sorted_array(A, mid);},
      [&]() {return t_from_sorted_array(A+mid+1, n-mid-1);});

      return t_join3(P.first, P.second, m);
  }

  // assumes array A is of lenght n and is sorted with no duplicates
  template <class BinaryOp>
  static Node* t_multi_insert_rec(Node* b, E* A, size_t n,const BinaryOp& op) {
      if (!b) return t_from_sorted_array(A,n);
      if (n == 0) return b;

      Node* join = copy_if_needed(b);

      auto less_first = [] (E a, E b) -> bool {
        return comp(a.first,b.first);};
      V x;
      E b_entry = std::make_pair(b->get_key(),x);
      size_t mid = pbbs::binary_search(sequence<E>(A, n), b_entry, less_first);
      bool dup = (mid < n) && less_first(b_entry,A[mid]);
      if (dup) join->set_value(op(A[mid].second, join->get_value()));
      
      auto P = fork<Node*>(do_parallel(get_node_count(b), n),
        [&] () {return t_multi_insert_rec(join->lc, A, mid, op);},
        [&] () {return t_multi_insert_rec(join->rc, A+mid+dup,
					  n-mid-dup, op);});
      
      return t_join3(P.first, P.second, join);
  }

  template<class NodeType, class Func>
  static void t_forall(Node* b, const Func& f, NodeType*& join_node) {
      if (!b) {
          join_node = NULL;
          return;
      }

      typename NodeType::entry_type entry
          = make_pair(b->get_key(), f(b->get_entry()));

      join_node = new NodeType(entry);

      size_t mn = get_node_count(b);
      auto P = fork<Node*>(mn >= node_limit,
      [&] () {return t_forall(b->lc, f, join_node->lc);},
      [&] () {return t_forall(b->rc, f, join_node->rc);});

      join_node->update();
  }

  static Node* t_find(Node* b, const K& key) {
      while (b) {
          if ( comp(key, b->get_key()) ) b = b->lc;
          else if ( comp(b->get_key(), key) )  b = b->rc;
          else return b;
      }
      return NULL;
  }

  static Node* t_previous(Node* b, const K& key) {
      Node* r = NULL;
      while (b) {
          if ( comp(b->get_key(), key) ) {
              r = b; b = b->rc;
          } else 
              b = b->lc;
      }
      return r;
  }

  static Node* t_next(Node* b, const K& key) {
      Node* r = NULL;
      while (b) {
          if (comp(key, b->get_key()) ) {
              r = b; b = b->lc;
          } else 
              b = b->rc;
      }
      return r;
  }

  static size_t t_rank(Node* b, const K& key) {
      size_t ret = 0;
      while (b) {
          if ( comp(b->get_key(), key) ) {
              ret += 1 + get_node_count(b->lc);
              b = b->rc;
          } else 
              b = b->lc;
      }
      return ret;
  }

  static Node* t_select(Node* b, size_t rank) {    
      size_t lrank = rank;
      while (b) {
          size_t left_size = get_node_count(b->lc);
          if (lrank > left_size) {
              lrank -= left_size + 1;
              b = b->rc;  
          }
          else if (lrank < left_size) 
              b = b->lc;
          else 
              return b;
    }
    return NULL;
  }

  struct range_sum_t {
    aug_type result;
    range_sum_t() : result(aug_class::get_empty()) {}
    void add_entry(K key, V val) {
      result = aug_class::combine(result, aug_class::from_entry(key, val));
    }
    void add_aug_val(aug_type av) {
      result = aug_class::combine(result, av);
    }
  };
    
  template<class aug>
  // the sum right of or at key
  static void t_sum_right(Node* b, const K& key, aug& a) {
      while (b) {
	if (!comp(b->get_key(), key)) {
	  a.add_entry(b->get_key(), b->get_value());
	  if (b->rc) a.add_aug_val((b->rc)->aug_val);
	  b = b->lc;
	} else b = b->rc;
      }
  }

  template<class aug>
  // the sum left of or at key
  static void t_sum_left(Node* b, const K& key, aug& a) {
      while (b) {
	if (!comp(key, b->get_key())) {
	  a.add_entry(b->get_key(), b->get_value());
	  if (b->lc) a.add_aug_val((b->lc)->aug_val);
	  b = b->rc;
	} else b = b->lc;
      }
  }

  template<class aug>
  static void t_sum_range(Node* b, const K& key_left, const K& key_right, aug& a) {
    Node* r = t_range_root(b, key_left, key_right);
    if (r) {
      t_sum_right(r->lc, key_left, a);  // add in left side (right of or at key_left)
      a.add_entry(r->get_key(), r->get_value());   // add in middle
      t_sum_left(r->rc, key_right, a);  // add in right side (left of or at key_right)
    }
  }

  template<typename Func>
  static Node* aug_select(Node* b, const Func& f) {
    if (b == NULL) return NULL;
    if (f(get_aug(b->lc))) {
      if (f(aug_class::from_entry(b->get_key(),b->get_value())))
	return aug_select(b->rc, f);
      return b;
    } return aug_select(b->lc, f);
  }

  // applies the function f(e, i) to every entry and its index i
  // works in parallel in O(log n) depth
  template<typename F>
  static void t_foreach_index(Node* a, size_t start, const F& f) {
    if (!a) return;
    size_t lsize = get_node_count(a->lc);
    fork_no_result(lsize >= node_limit,
      [&] () {t_foreach_index(a->lc, start, f);},
      [&] () {t_foreach_index(a->rc, start + lsize + 1,f);});
    f(a->get_entry(), start+lsize);
  }

  // similar to above but sequential using in-order traversal
  // usefull if using 20th century constructs such as iterators
  template<typename F>
  static void t_foreach_seq(Node* a, const F& f) {
    if (!a) return;
    t_foreach_seq(a->lc, f);
    f(a->get_entry());
    t_foreach_seq(a->rc, f);
  }

  // applies f to all entries at or right of key
  template<typename F>
  static void t_foreach_seq_right(Node* b, const K& key, const F& f) {
    if (b) {
      if (!comp(b->get_key(), key)) {
	t_foreach_seq_right(b->lc, key, f);
	f(b->get_entry());
	t_foreach_seq(b->rc, f);
      } else t_foreach_seq_right(b->rc, key, f);
    }
  }

  // applies f to all entries at or left of key
  template<typename F>
  static void t_foreach_seq_left(Node* b, const K& key, const F& f) {
    if (b) {
      if (!comp(key, b->get_key())) {
	t_foreach_seq(b->lc, f);
	f(b->get_entry());
	t_foreach_seq_left(b->rc, key, f);
      } else t_foreach_seq_left(b->lc, key, f);
    }
  }

  // applies f to all entries between key_left and key_right, inclusive of both
  template<typename F>
  static void t_foreach_seq_range(Node* a, const K& key_left, const K& key_right, 
				  const F& f) {
    Node* r = t_range_root(a, key_left, key_right);
    if (r) {
      t_foreach_seq_right(r->lc, key_left, f);
      f(r->get_entry());
      t_foreach_seq_left(r->rc, key_right, f);
    }
  }

  //  *****************************************************
  //  The following is all for constructing augmented maps from
  //  arrays.  It needs to sort the array and combine duplicates
  //  *****************************************************

  // User supplies a function for reducing all duplicates
  template<class Vin, class Reduce>
  std::pair<E*, size_t> 
  static reduce_duplicates(std::pair<K,Vin>* A, size_t n, Reduce& reduce) {

    // determines the index of start of each block of equal keys
    // and copies values into vals
    bool* is_start = new bool[n];
    Vin* vals = pbbs::new_array_no_init<Vin>(n);
    is_start[0] = 1;
    vals[0] = A[0].second;
    parallel_for(size_t i = 1; i < n; i++) {
      is_start[i] = comp(A[i-1].first, A[i].first);
      new (static_cast<void*>(vals+i)) Vin(A[i].second);
    }

    sequence<tree_size_t> I = 
      pbbs::pack_index<tree_size_t>(sequence<bool>(is_start,n));
    delete[] is_start;

    // combines over each block of equal keys using function reduce
    E* B = pbbs::new_array<E>(I.size());
    parallel_for(size_t i = 0; i < I.size(); i++) {
      size_t start = I[i];
      size_t end = (i==I.size()-1) ? n : I[i+1];
      B[i] = E(A[start].first, reduce(vals+start,vals+end));
    }

    pbbs::delete_array(vals, n);
    return std::pair<E*,size_t>(B,I.size());
  }

  // User supplies a binary associative function for 
  // pairwise combining duplicates
  template<class Bin_Op>
  std::pair<E*, size_t> 
  static combine_duplicates(E* A, size_t n, Bin_Op& bin_op) {

    // determines the index of start of each block of equal keys
    bool* is_start = new bool[n];
    is_start[0] = 1;
    parallel_for(size_t i = 1; i < n; i++)
      is_start[i] = comp(A[i-1].first, A[i].first);
    sequence<tree_size_t> I = 
      pbbs::pack_index<tree_size_t>(sequence<bool>(is_start,n));
    delete[] is_start;

    // combines over each block of equal keys using function reduce
    E* B = pbbs::new_array<E>(I.size());
    parallel_for(size_t i = 0; i < I.size(); i++) {
      size_t start = I[i];
      size_t end = (i==I.size()-1) ? n : I[i+1];
      auto get_val = [&] (size_t i) {return A[i+start].second;};
      B[i] = E(A[start].first, 
	       pbbs::reduce(make_sequence<V>(end-start,get_val),bin_op));
    }
    return std::pair<E*,size_t>(B,I.size());
  }

  template<class Bin_Op>
  static size_t combine_duplicates_seq(E* A, size_t n, Bin_Op& bin_op) {
    size_t j = 0;
    for (size_t i=1; i<n; i++) {
      if (A[j].first < A[i].first)
	A[++j] = A[i];
      else A[j].second = bin_op(A[j].second,A[i].second);
    }
    return j+1;
  }

  // Keeps the first among duplicates
  static std::pair<E*, size_t> remove_duplicates(E* A, size_t n) {

    // determines the index of start of each block of equal keys
    // and copies values into vals
    sequence<bool> Fl(n);
    Fl[0] = true;
    parallel_for (size_t i=1; i < n; i++) 
      Fl[i] = comp(A[i-1].first, A[i].first);

    sequence<E> X = pbbs::pack(sequence<E>(A, n), Fl);

    return std::pair<E*,size_t>(X.as_array(), X.size());
  }

  // Keeps the first among duplicates
  // seqeuntial, so does in place
  static size_t remove_duplicates_seq(E* A, size_t n) {
    size_t  j = 1;
    for (size_t i=1; i<n; i++)
      if (A[j-1].first < A[i].first)
	A[j++] = A[i];
    return  j;
  }

  template<class key_val>
  static void sort_keys(key_val* A, tree_size_t n, 
			bool is_sorted = false, bool sequential = false) {
    auto compare = ([&] (key_val a, key_val b) {
	return comp(a.first, b.first);});
    if (!is_sorted) {
      if (sequential) std::sort(A, A+n, compare);
      else pbbs::sample_sort(A, n, compare);
    }
  }

  template<class Vin, class Reduce>
  static Node* build_reduce(std::pair<K,Vin>* A, size_t n,
			   const Reduce& reduce, bool is_sorted) {
    if (n == 0) return NULL;
    sort_keys(A, n, is_sorted);
    std::pair<E*,size_t> X = reduce_duplicates(A, n, reduce);
    Node* r = t_from_sorted_array(X.first, X.second);
    pbbs::delete_array(X.first, X.second);
    return r;
  }

  static Node* multi_insert(Node* In, E* A, size_t n, 
			    bool is_sorted, bool sequential, bool inplace=false) {
    if (n == 0) return NULL;
    E* B = A;
    if (!inplace) {
      B = pbbs::new_array<E>(n);  
      parallel_for(size_t i =0; i < n; i++) B[i] = A[i];
    }
    sort_keys(B, n, is_sorted);
    Node* r;
    
    if (sequential || n < (1 << 14)) {
      size_t m = remove_duplicates_seq(B,n);
      r = t_multi_insert_rec(In, B, m, get_left<V>());
      //} else if (is_sorted == 2) {
      //return t_multi_insert_rec(In, A, n, get_left<V>());
    } else {      
      std::pair<E*,size_t> X = remove_duplicates(B, n);
      r = t_multi_insert_rec(In, X.first, X.second, 
				   get_left<V>());
      pbbs::delete_array(X.first, X.second);
    }
    if (!inplace) pbbs::delete_array(B,n);
    return r;
  }

  template<class Combine>
  static Node* multi_insert(Node* In, E* A, size_t n, 
			    const Combine& f, 
			    bool is_sorted, bool sequential) {
    if (n == 0) return NULL;
    E* B = pbbs::new_array<E>(n);  
    parallel_for(size_t i =0; i < n; i++) B[i] = A[i];
    Node* r;

    if ((sequential || n < (1 << 14)) && (n < (1 << 20))) {
      sort_keys(B, n, is_sorted, 1);
      size_t m = combine_duplicates_seq(B, n, f);
      r = t_multi_insert_rec(In, B, m, f);
    } else {
      sort_keys(B, n, is_sorted);
      std::pair<E*,size_t> X = combine_duplicates(B, n, f);
      r = t_multi_insert_rec(In, X.first, X.second, f);
      pbbs::delete_array(X.first, X.second);
    }
    pbbs::delete_array(B,n);
    return r;
  }

private:

    static key_compare comp;
};


template<class Node> typename tree_operations<Node>::key_compare 
tree_operations<Node>::comp = tree_operations<Node>::key_compare();
