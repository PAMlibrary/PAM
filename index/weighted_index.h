#pragma once

#include "augmented_map.h"
#include "tree_set.h"

struct strless {
  bool operator()(char* s1, char* s2) {
    return (strcmp(s1,s2) < 0); }
};

auto add = [] (float a, float b) -> float {
    return a + b;};

struct inv_index {
  using word = char*;
  using doc_id = int;
  using weight = float;
 
  struct m_aug {
    using aug_t = weight;
    static aug_t get_empty() {return 0;}
    static aug_t from_entry(doc_id k, weight v) {
      return v;}
    static aug_t combine(aug_t a, aug_t b) {
      return (b > a) ? b : a;}
  };

  using post_elt = pair<doc_id, weight>;
  using post_list = augmented_map<doc_id, weight, m_aug>;
  
  using index_elt = pair<word, post_elt>;
  using index = tree_map<word, post_list, strless>;

  index idx;

  inv_index(index_elt* start, index_elt* end) {
    size_t n = end - start;
    post_list::reserve((size_t) round(.45*n));
    index::reserve(n/300);
    auto reduce = [] (post_elt* s, post_elt* e) {
      return post_list(s,e,add,0,1); };
    idx.build_reduce<post_elt>(start, end, reduce);
  }

  post_list get_list(word w) {
    maybe<post_list> p = idx.find(w);
    if (p) return *p;
    else return post_list();
  }

  post_list And(post_list a, post_list b) {
    return map_intersect(a,b,add);}

  post_list Or(post_list a, post_list b) {
    return map_union(a,b,add);}

  post_list And_Not(post_list a, post_list b) {
    return map_difference(a,b);}

  vector<post_elt> top_k(post_list a, int k) {
    int l = min<int>(k,a.size());
    vector<post_elt> vec(l);
    post_list b = a;
    for (int i=0; i < l; i++) {
      weight m = b.aug_val();
      auto f = [m] (weight v) {return v < m;};
      vec[i] = *b.aug_select(f);
      b.remove(vec[i].first);
    }
    return vec;
  }
};
