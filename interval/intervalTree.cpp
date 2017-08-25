#include "augmented_map.h"
#include <iostream>
#include <algorithm>
#include "pbbs-include/random.h"


using namespace std;

using point = int;
typedef pair<point,point> par;

struct interval_map {

  using interval = pair<point, point>;
  struct aug {
  using aug_t = interval;
    static aug_t get_empty() {
      return interval(0, 0);}
    static aug_t from_entry(point k, point v) {
      return interval(k, v);}
    static aug_t combine(aug_t a, aug_t b) {
      return (a.second > b.second) ? a : b;}
  };

  using amap = augmented_map<point,point,aug>;
  amap m;

  interval_map(size_t n) {
    amap::reserve(n); 
  }

  void finish() {
    amap::finish();
  }

  interval_map(interval* A, size_t n) {
    m = amap(A,A+n); }

  bool stab(point p) {
    return (m.aug_left(p).second > p);}

  void insert(interval i) {m.insert(i);}

  vector<interval> report_all(point p) {
    vector<interval> vec;
    amap a = m;
    interval I = a.aug_left(p);
    while (I.second > p) {
      vec.push_back(I);
      a.remove(I.first); 
      I = a.aug_left(p); }
    return vec; }

  void remove_small(point l) {
    auto f = [&] (interval I) {
      return (I.second - I.first >= l);};
    m.filter(f); }
};

long str_to_long(char* str) {
    return strtol(str, NULL, 10);
}

int main(int argc, char** argv) {
  if (argc == 1) {
	  cout << argv[0] << " <n> <query_num> <rounds>" << endl;
	  cout << "Will run on n=query_num=100000000, rounds=5" << endl;
  }
  size_t n = 100000000;
  if (argc > 1) n = str_to_long(argv[1]);
  size_t q_num = n;
  if (argc > 2) q_num = str_to_long(argv[2]);
  size_t rounds = 5;
  if (argc > 3) rounds = str_to_long(argv[3]);

  par *v = pbbs::new_array<par>(n);
  par *vv = pbbs::new_array<par>(n);
  size_t max_size = (((size_t) 1) << 31)-1;

  pbbs::random r = pbbs::default_random;
  parallel_for (size_t i = 0; i < n; i++) {
    point start = r.ith_rand(2*i)%(max_size/2);
    point end = start + r.ith_rand(2*i+1)%(max_size-start);
    v[i] = make_pair(start, end);
  }

  bool* result = new bool[q_num];
  int* queries = new point[q_num];
  parallel_for (size_t i = 0; i < q_num; i++)
    queries[i] = r.ith_rand(6*i)%max_size;
  
  for (size_t i=0; i < rounds; i++) {
    parallel_for (size_t i = 0; i < n; i++) {
      vv[i] = v[i];
    }
    const size_t threads = __cilkrts_get_nworkers();
    interval_map xtree(n);
    timer t;
    t.start();
    interval_map itree(vv,n);
    double tm = t.stop();

    timer tq;
    tq.start();
    parallel_for (size_t i = 0; i < q_num; i++) 
      result[i] = itree.stab(queries[i]);
    double tm2 = tq.stop();

    cout << "n=" << n << "  threads=" << threads 
	 << "  build=" << tm << "  query=" << tm2 << endl;

    itree.finish();
  }
}
