// AVL Tree balance operations
#pragma once

#include "common.h"
#include <cstddef> 
#include <utility>
#include "abstract_node.h"

template<class Node>
class avl_tree {
public:

    static tree_size_t combine_ranks(tree_size_t l, tree_size_t r) {
      return std::max(l,r) + 1;
    }

    static tree_size_t singleton_rank() { return 1; }
    
    static Node* t_join(Node* t1, Node* t2, Node* k) {
        if (is_too_heavy(t1, t2)) {
            return right_join(t1, t2, k);
        } else if (is_too_heavy(t2, t1)) {
            return left_join(t1, t2, k);
        } else {
            k->lc = t1, k->rc = t2;
            k->update();
            return k;
        }
    }

    static inline bool is_single_rotation(const Node* t, const bool dir) {
        bool heavier = get_rank(t->lc) > get_rank(t->rc);
        return dir ? heavier : !heavier;
    }

    static inline bool is_too_heavy(const Node* t1, const Node* t2) {
        return get_rank(t1) > get_rank(t2) + 1;
    }
	
    static inline bool are_balanced(const Node* t1, const Node* t2) {
        return (get_rank(t1) <= get_rank(t2) + 1) && (get_rank(t2) <= get_rank(t1) + 1);
    }

    static Node* rebalance(Node* t) {
        Node* ret = t;
        if (is_too_heavy(t->lc, t->rc)) {
            if (is_single_rotation(t->lc, 1)) ret = rotate_right(t);
            else ret = double_rotate_right(t);

        }else if(is_too_heavy(t->rc, t->lc)) {
            if (is_single_rotation(t->rc, 0)) ret = rotate_left(t);
            else ret = double_rotate_left(t);
        }

        return ret;
    }

    static Node* rebalance_right(Node* t) {
        Node* ret = t;
        if (is_too_heavy(t->rc, t->lc)) {
            if (is_single_rotation(t->rc, 0)) ret = rotate_left(t);
            else ret = double_rotate_left(t);
        } else ret->update();

        return ret;
    }

    static Node* rebalance_left(Node* t) {
        Node* ret = t;
        if (is_too_heavy(t->lc, t->rc)) {
            if (is_single_rotation(t->lc, 1)) ret = rotate_right(t);
            else ret = double_rotate_right(t);
        } else ret->update();

        return ret;
    }

    static Node* right_join(Node* t1, Node* t2, Node* k) {
        if (!is_too_heavy(t1, t2))
	       return join_node(t1, t2, k);
      
        Node* ret = copy_if_needed(t1);
        ret->rc = right_join(ret->rc, t2, k);
        return rebalance_right(ret);
    }

    static Node* left_join(Node* t1, Node* t2, Node* k) {
      if (!is_too_heavy(t2, t1))
	       return join_node(t1, t2, k);

        Node* ret = copy_if_needed(t2);
        ret->lc = left_join(t1, ret->lc, k);
        return rebalance_left(ret);
    }
};
