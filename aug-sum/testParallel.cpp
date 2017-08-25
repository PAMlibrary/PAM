#include <cilk/cilk.h>
#include <algorithm>
#include <iostream>
#include <cstring>
#include <vector>
#include <ctime>
#include <random>
#include <cstdio>
#include <set>
#include <map>
#include "augmented_map.h"
#include "pbbs-include/get_time.h"
#include "pbbs-include/sequence_ops.h"
#include "pbbs-include/random_shuffle.h"

using namespace std;

#ifdef LARGE
using key_type = size_t;
#else
using key_type = unsigned int;
#endif

struct aug {
  using aug_t = key_type;
  static aug_t get_empty() { return 0;}
  static aug_t from_entry(key_type k, key_type v) {
    return v;}
  static aug_t combine(aug_t a, aug_t b) {
    return a + b;}
};

//using tmap = augmented_map<key_type, key_type, aug>;
using par = pair<key_type, key_type>;
using tmap = tree_map<key_type, key_type>;

struct mapped {
    key_type k, v;
    mapped(key_type _k, key_type _v) : k(_k), v(_v) {};
    mapped(){};

    bool operator < (const mapped& m) 
        const { return k < m.k; }
    
    bool operator > (const mapped& m) 
        const { return k > m.k; }
    
    bool operator == (const mapped& m) 
        const { return k == m.k; }
};


size_t str_to_int(char* str) {
    return strtol(str, NULL, 10);
}


std::mt19937_64& get_rand_gen() {
    static thread_local std::random_device rd;
    static thread_local std::mt19937_64 generator(rd());
    return generator;
}


par* uniform_input(size_t n, size_t window, bool shuffle = false) {
    par *v = new par[n];

    parallel_for (size_t i = 0; i < n; i++) {
        uniform_int_distribution<> r_keys(1, window);
        key_type k = r_keys(get_rand_gen());
        key_type c = i; //r_keys(get_rand_gen());
        v[i] = make_pair(k, c);
    }

    auto addfirst = [] (par a, par b) -> par {
      return par(a.first+b.first, b.second);};
    auto vv = sequence<par>(v,n);
    pbbs::scan(vv,vv,addfirst,par(0,0),pbbs::fl_scan_inclusive);
    if (shuffle) pbbs::random_shuffle(vv);
    return v;
}

par* uniform_input_unsorted(size_t n, size_t window) {
    par *v = new par[n];

    parallel_for (size_t i = 0; i < n; i++) {
        uniform_int_distribution<> r_keys(1, window);

        key_type k = r_keys(get_rand_gen());
        key_type c = r_keys(get_rand_gen());

        v[i] = make_pair(k, c);
    }

    return v;
}

bool check_union(const tmap& m1, const tmap &m2, const tmap& m3) {
    vector<key_type> e;
    m1.keys(back_inserter(e));
    m2.keys(back_inserter(e));

    sort(e.begin(), e.end());
    e.erase(unique(e.begin(), e.end()), e.end());

    bool ret = (m3.size() == e.size());

    for (auto it : e) 
        ret &= m3.contains(it);
    
    return ret;
}


bool check_intersect(const tmap& m1, const tmap& m2, const tmap& m3) {

    vector<key_type> e1, e2, e3;
    m1.keys(back_inserter(e1));
    m2.keys(back_inserter(e2));

    set_intersection(e1.begin(), e1.end(), e2.begin(), e2.end(), back_inserter(e3));

    bool ret = (m3.size() == e3.size());

    for (auto it : e3)
        ret &= m3.contains(it);

    return ret;
    
}

bool check_difference(const tmap& m1, const tmap& m2, const tmap& m3) {
    vector<key_type> e1, e2, e3;
    m1.keys(back_inserter(e1));
    m2.keys(back_inserter(e2));

    set_difference(e1.begin(), e1.end(), e2.begin(), e2.end(), back_inserter(e3));

    bool ret = (m3.size() == e3.size());

    for (auto it : e3)
        ret &= m3.contains(it);

    return ret;
}


bool check_split(key_type key, const par* v, const pair<tmap, tmap>& m) {
    size_t n = m.first.size() + m.second.size() + 1;

    bool ret = 1;
    for (size_t i = 0; i < n; ++i) {
        if (v[i].first < key) 
            ret &= m.first.contains(v[i].first);
        if (v[i].first > key) 
            ret &= m.second.contains(v[i].first);
    }

    return ret;
}


template<class T>
bool check_filter(const tmap& m, const tmap& f, T cond) {

    vector <par> e1, e2;
    m.content(back_inserter(e1));
    f.content(back_inserter(e2));

    vector<par> v;
    for (size_t i = 0; i < e1.size(); ++i) 
        if (cond(e1[i])) v.push_back(e1[i]);

    bool ret = (e2.size() == v.size());
    for (size_t i = 0; i < v.size(); ++i) 
        ret &= e2[i] == v[i];

    return ret;
}

bool contains(const tmap& m, const par* v) {
    bool ret = 1;
    for (size_t i = 0; i < m.size(); ++i) {
        ret &= m.contains(v[i].first);
    }

    return ret;
}


double test_union(size_t n, size_t m) {

    par* v1 = uniform_input(n, 20); 
    tmap m1(v1, v1 + n);
    
    par* v2 = uniform_input(m, (n/m) * 20); 
    tmap m2(v2, v2 + m);

    timer t;
    t.start();
    tmap m3 = map_union(m1, m2);
    double tm = t.stop();

    assert(m1.size() == n && "map size is wrong.");
    assert(m2.size() == m && "map size is wrong.");
    assert(check_union(m1, m2, m3) && "union is wrong");

    delete[] v1;
    delete[] v2;

    return tm;
}

double test_multi_insert_old(size_t n, size_t m) {

    tmap m1;
    par* v2 = uniform_input_unsorted(m, 1000000000); 

    timer t;
    t.start();
    m1.multi_insert(v2,v2+m,false);
    double tm = t.stop();

    delete[] v2;
    return tm;
}


double test_intersect(size_t n, size_t m) {    
    par* v1 = uniform_input(n, 20);
    tmap m1(v1, v1 + n);

    par* v2 = uniform_input(m, (n/m) * 20);
    tmap m2(v2, v2 + m);

    timer t;
    t.start();
    tmap m3 = map_intersect(m1, m2);
    double tm = t.stop();

    assert(m1.size() == n && "map size is wrong");
    assert(m2.size() == m && "map size is wrong");
    assert(check_intersect(m1, m2, m3) && "intersect is wrong");

    delete[] v1;
    delete[] v2;

    return tm;
}

double test_deletion(size_t n, size_t m) {
    par *v = uniform_input(n, 20);
    tmap m1(v, v + n);

    par *u = uniform_input(m, (n/m)*20, true);

    timer t;
    t.start();
    for (size_t i = 0; i < m; ++i) 
      m1.remove(u[i].first);
    double tm = t.stop();

    delete[] v;
    delete[] u;
    return tm;
}

double test_deletion_destroy(size_t n) {
    par *v = uniform_input(n, 20, true);
    tmap m1;
	for (size_t i = 0; i < n; ++i) 
		m1.insert(v[i]);
    pbbs::random_shuffle(sequence<par>(v, n));

    timer t;
    t.start();
    for (size_t i = 0; i < n; ++i) 
      m1.remove(v[i].first);
    double tm = t.stop();

    delete[] v;
    return tm;
}

double test_insertion_build(size_t n) {
  par *v = uniform_input(n, 20, true);
  tmap m1;

  timer t;
  t.start();
  for (size_t i = 0; i < n; ++i) 
    m1.insert(v[i]);
  double tm = t.stop();
  
  delete[] v;
  return tm;
}

double test_insertion(size_t n, size_t m) {
    par *v = uniform_input(n, 20);
    tmap m1(v, v + n);

    par *u = uniform_input(m, (n/m)*20, true);

    timer t;
    t.start();
    for (size_t i = 0; i < m; ++i) 
      m1.insert(u[i]);

    double tm = t.stop();
    delete[] v;
    delete[] u;
    return tm;
}

double test_multi_insert(size_t n, size_t m) {
    par *v = uniform_input(n, 20);
    par *u = uniform_input(m, (n/m)*20, true);
    tmap m1(v, v + n);
    tmap m2 = m1;

    timer t;
    t.start();
    m1.multi_insert(u,u+m,false);
    double tm = t.stop();

    delete[] v;
    delete[] u;
    return tm;
}


double test_dest_multi_insert(size_t n, size_t m) {
    par *v = uniform_input(n, 20);
    par *u = uniform_input(m, (n/m)*20, true);
    tmap m1(v, v + n);

    timer t;
    t.start();
    m1.multi_insert(u,u+m,false);
    double tm = t.stop();

    delete[] v;
    delete[] u;
    return tm;
}


double stl_insertion(size_t n, size_t m) {
    par *v = uniform_input(n, 20);
    par *u = uniform_input(m, (n/m)*20, true);

    map<key_type, key_type> m1;
    for (size_t i = 0; i < n; ++i) 
        m1.insert(v[i]);

    timer t;
    t.start();
    for (size_t i = 0; i < m; ++i) 
      m1.insert(u[i]);	
    double tm = t.stop();

    delete[] v;
    delete[] u;
    return tm;
}

double stl_insertion_build(size_t n) {
    par *v = uniform_input(n, 20, true);
    map<key_type, key_type> m1;

    timer t;
    t.start();
    for (size_t i = 0; i < n; ++i) 
        m1.insert(v[i]);
    double tm = t.stop();

    delete[] v;
    return tm;
}

double stl_delete_destroy(size_t n) {
    par *v = uniform_input(n, 20, true);
    map<key_type, key_type> m1;
    pbbs::random_shuffle(sequence<par>(v, n));

    timer t;
    t.start();
    for (size_t i = 0; i < n; ++i) 
      m1.erase(v[i].first);
    double tm = t.stop();

    delete[] v;
    return tm;
}


double test_build(size_t n) {
    par *v = uniform_input_unsorted(n, 1000000000);

    timer t;
    t.start();
    tmap m1(v, v + n, false, false, true);
    double tm = t.stop();

    assert(m1.size() == n && "map size is wrong");
    assert(contains(m1, v) && "build is wrong");

    delete[] v;
    return tm;
}

double test_filter(size_t n) {
    par* v = uniform_input(n, 20);
    tmap m1(v, v + n);
    tmap res = m1;

    timer t;
    t.start();
    auto cond = [] (par t) { return t.first < t.second; };
    res.filter(cond);
    double tm = t.stop();

    assert(m1.size() == n && "map size is wrong");
    assert(check_filter(m1, res, cond) && "filter is wrong");

    delete[] v;
    return tm;
}

double test_dest_union(size_t n, size_t m) {
    par* v1 = uniform_input(n, 20);
    tmap m1(v1, v1 + n);
    tmap m1_copy(v1, v1 + n);

    par* v2 = uniform_input(m, (n/m) * 20);
    tmap m2(v2, v2 + m);
    tmap m2_copy(v2, v2 + m);

    timer t;
    t.start();
    tmap m3 = map_union(move(m1), move(m2));
    double tm = t.stop();
    
    assert(m1.size() == 0 && "map size is wrong");
    assert(m2.size() == 0 && "map size is wrong");
    assert(check_union(m1_copy, m2_copy, m3) && "union is wrong");

    delete[] v1;
    delete[] v2;

    return tm;
}


double test_dest_intersect(size_t n, size_t m) {
    par* v1 = uniform_input(n, 20);
    tmap m1(v1, v1 + n);
    tmap m1_copy(v1, v1 + n);

    par* v2 = uniform_input(m, (n/m) * 20);
    tmap m2(v2, v2 + m);
    tmap m2_copy(v2, v2 + m);

    timer t;
    t.start();
    tmap m3 = map_intersect(move(m1), move(m2));
    double tm = t.stop();
    
    assert(m1.size() == 0 && "map size is wrong");
    assert(m2.size() == 0 && "map size is wrong");
    assert(check_intersect(m1_copy, m2_copy, m3) && "intersect is wrong");

    delete[] v1;
    delete[] v2;

    return tm;
}


double test_split(size_t n) {    
    par *v = uniform_input(n, 20);

    tmap m1(v, v+n);

    key_type key = v[n / 2].first;

    timer t;
    t.start();
    pair<tmap, tmap> res = m1.split(key);
    double tm = t.stop();

    assert(m1.size() == n && "map size is wrong");
    assert(res.first.size() + res.second.size() + 1 == n 
        && "splitted map size is wrong");
    assert(check_split(key, v, res) && "split is wrong");

    delete[] v;
    return tm;
}

double test_difference(size_t n, size_t m) {
    par* v1 = uniform_input(n, 20);
    tmap m1(v1, v1 + n);

    par *v2 = uniform_input(m, (n/m) * 20);
    tmap m2(v2, v2 + m);

    timer t;
    t.start();
    tmap m3 = map_difference(m1, m2);
    double tm = t.stop();

    assert(m1.size() == n && "map size is wrong");
    assert(m2.size() == m && "map size is wrong");
    assert(check_difference(m1, m2, m3));

    delete[] v1;
    delete[] v2;

    return tm;
}

double test_find(size_t n, size_t m) {
    par* v1 = uniform_input(n, 20);
    key_type max_key = v1[n-1].first;
    tmap m1(v1, v1 + n, true);

    par *v2 = uniform_input_unsorted(m, max_key);

    bool *v3 = new bool[m];

    timer t;
    t.start();
    parallel_for(size_t i=0; i < m; i++)
      v3[i] = m1.contains(v2[i].first);

    double tm = t.stop();

    delete[] v1;
    delete[] v2;
    delete[] v3;

    return tm;
}

double test_range(size_t n, size_t m) {
    par* v1 = uniform_input(n, 20);
    key_type max_key = v1[n-1].first;
    tmap m1(v1, v1 + n, true);

    par *v2 = uniform_input_unsorted(m, max_key);
    // window size is 1/1000 of total width
    size_t win = v1[n-1].first/1000;

    tmap *v3 = new tmap[m];

    timer t;
    t.start();
    parallel_for(size_t i=0; i < m; i++) {
      v3[i] = m1.range_fast(v2[i].first,v2[i].first+win);
    }
    
    double tm = t.stop();

    delete[] v1;
    delete[] v2;
    delete[] v3;

    return tm;
}

double test_aug_range(size_t n, size_t m) {
	par *v1 = uniform_input_unsorted(n, 1000000000);
    key_type max_key = 1000000000;
	tmap m1(v1, v1 + n, false, false, true);

    par *v2 = uniform_input_unsorted(m, max_key);
    // window size is 1/1000 of total width
    size_t win = max_key/1000;

    key_type *v3 = new key_type[m];

    timer t;
    t.start();
    parallel_for(size_t i=0; i < m; i++) {
      v3[i] = m1.aug_range(v2[i].first,v2[i].first+win);
    }
    
    double tm = t.stop();

    delete[] v1;
    delete[] v2;
    delete[] v3;

    return tm;
}

double test_aug_left(size_t n, size_t m) {
    par* v1 = uniform_input(n, 20);
    key_type max_key = v1[n-1].first;
    tmap m1(v1, v1 + n, true);

    par *v2 = uniform_input_unsorted(m, max_key);
    // window size is 1/1000 of total width

    key_type *v3 = new key_type[m];

    timer t;
    t.start();
    parallel_for(size_t i=0; i < m; i++) {
      v3[i] = m1.aug_left(v2[i].first);
    }
    
    double tm = t.stop();

    delete[] v1;
    delete[] v2;
    delete[] v3;

    return tm;
}


double stl_set_union(size_t n, size_t m) {
    par *v1 = uniform_input(n, 20);
    par *v2 = uniform_input(m, 20 * (n / m));

    set <key_type> s1, s2, sret;
    for (size_t i = 0; i < n; ++i) {
      s1.insert(v1[i].first);
    } 
    for (size_t i = 0; i < m; ++i) {
      s2.insert(v2[i].first);
    }
    
    timer t;
    t.start(); 
    set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(sret, sret.begin()));
    double tm = t.stop();

    delete[] v1;
    delete[] v2;

    return tm;
}

double stl_vector_union(size_t n, size_t m) {
    par *v1 = uniform_input(n, 20);
    par *v2 = uniform_input(m, 20 * (n / m));

    vector <mapped> s1, s2, sret;

    for (size_t i = 0; i < n; ++i) {
      s1.push_back(mapped(v1[i].first, v1[i].second));
    } 
    for (size_t i = 0; i < m; ++i) {
      s2.push_back(mapped(v2[i].first, v2[i].second));
    }
    
    timer t;
    t.start();
    set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), back_inserter(sret));
    double tm = t.stop();
    
    delete[] v1;
    delete[] v2;

    return tm;
}


string test_name[] = { 
    "persistent-union",      // 0
    "persistent-intersect",  // 1
    "insert",                // 2
    "std::map-insert",       // 3
    "build",                 // 4                
    "filter",                // 5
    "destructive-union",     // 6
    "destructive-intersect", // 7
    "split",                 // 8
    "difference",            // 9
    "std::set_union",        // 10
    "std::vector_union",     // 11
    "find",                  // 12
    "delete",                // 13
    "multi_insert",           // 14
    "destructive multi_insert",  //15
    "stl_insertion_build",   //16
    "test_deletion_destroy", //17
    "test_range", //18
    "test_aug_range", //19
	"test_insertion_build", //20
	"stl_delete_destroy", //21
	"test_aug_left", //22
    "nothing" 
};


double execute(size_t id, size_t n, size_t m) {

    switch (id) {
        case 0:  
            return test_union(n, m);
        case 1:  
            return test_intersect(n, m);
        case 2:  
            return test_insertion(n, m);
        case 3:  
            return stl_insertion(n, m);
        case 4:  
            return test_build(n);
        case 5:  
            return test_filter(n);
        case 6:  
            return test_dest_union(n, m);
        case 7:  
            return test_dest_intersect(n, m);
        case 8:  
            return test_split(n);
        case 9:  
            return test_difference(n, m);
        case 10: 
            return stl_set_union(n, m);
        case 11: 
            return stl_vector_union(n, m);
        case 12: 
            return test_find(n, m);
        case 13: 
            return test_deletion(n, m);
        case 14: 
			return test_multi_insert(n,m);
        case 15: 
			return test_dest_multi_insert(n,m);
        case 16:
			return stl_insertion_build(n);
        case 17:
			return test_deletion_destroy(n);
        case 18:
			return test_range(n,m);
        case 19:
			return test_aug_range(n,m);
	    case 20:
			return test_insertion_build(n);
	    case 21:
			return stl_delete_destroy(n);
	    case 22:
			return test_aug_left(n,m);
        default: 
            assert(false);
	    return 0.0;
    }
}

/*
 * argv[1] - test
 * argv[2] - n
 * argv[3] - m
 * argv[4] - repeat
 * argv[5] - isLarge (0 or 1) -- optional
 */
int main (int argc, char *argv[]) {
    if (argc < 5 || argc > 6) {
        fprintf(stderr, "Wrong number of arguments.\n");
		cout << argv[0] << " <test_id> <n> <m> <rounds>" << endl;
		cout << "test ids: " << endl;
		for (int i = 0; i < 23; i++) {
			cout << i << " - " << test_name[i] << endl;
		}
        exit(1);
    }

    size_t test_id = str_to_int(argv[1]);
    size_t n       = str_to_int(argv[2]);
    size_t m       = str_to_int(argv[3]);
    size_t repeat  = str_to_int(argv[4]);
	bool isLarge = 0;
    bool randomize = 0;
    if (argc > 5) isLarge = str_to_int(argv[5]);
    size_t threads = __cilkrts_get_nworkers();


    size_t reserve_size = (isLarge) ? n : 4 * n;
    for (size_t i = 0; i < repeat; ++i) {
        tmap::reserve(reserve_size, randomize);    
        double tm = execute(test_id, n, m);
        cout << "RESULT"  << fixed << setprecision(6)
             << "\ttest=" << test_name[test_id]
             << "\ttime=" << tm
             << "\tn=" << n
             << "\tm=" << m
             << "\titeration=" << i 
             << "\tp=" << threads << endl;

        tmap::finish();
     }

    return 0;
}
