//#define GLIBCXX_PARALLEL

#include <algorithm>
#include <iostream>
#include <vector>
#include <random>
#include <climits>
#include <cstdlib>
#include <iomanip>

#include <cilk/cilk_api.h>
#include "range_tree.h"
#include "pbbs-include/get_time.h"

using namespace std;

using data_type = int;
using point_type = Point<data_type, data_type>;
using tuple_type = tuple<data_type, data_type, data_type, data_type>;

struct Query_data {
	data_type x1, x2, y1, y2;
} query_data;

int str_to_int(char* str) {
    return strtol(str, NULL, 10);
}

int win;
int dist;

int random_hash(int seed, int a, int rangeL, int rangeR) {
	if (rangeL == rangeR) return rangeL;
	a = (a+seed) + (seed<<7);
	a = (a+0x7ed55d16) + (a<<12);
	a = (a^0xc761c23c) ^ (a>>19);
	a = (a+0x165667b1) + (a<<5);
	a = ((a+seed)>>5) ^ (a<<9);
	a = (a+0xfd7046c5) + (a<<3);
	a = (a^0xb55a4f09) ^ (a>>16);
	a = a % (rangeR-rangeL);
	if (a<0) a+= (rangeR-rangeL);
	a+=rangeL;
	return a;
}

vector<point_type> generate_points(size_t n, data_type a, data_type b) {
    vector<point_type> ret(n);

    cilk_for (size_t i = 0; i < n; ++i) {
        ret[i].x = random_hash('x', i, a, b);
        ret[i].y = random_hash('y', i, a, b);
        ret[i].w = random_hash('w', i, a, b);
    }

   return ret;
}

vector<Query_data> generate_queries(size_t q, data_type a, data_type b) {
    vector<Query_data> ret(q);

    cilk_for (size_t i = 0; i < q; ++i) {
        data_type x1 = random_hash('q'*'x', i*2, a, b);
        data_type y1 = random_hash('q'*'y', i*2, a, b);
		int dx = random_hash('d'*'x', i, 0, win);
		int dy = random_hash('d'*'y', i, 0, win);
		int x2 = x1+dx;
		int y2 = y1+dy;
		if (dist==0) {
			x2 = random_hash('q'*'x', i*2+1, a, b);
			y2 = random_hash('q'*'y', i*2+1, a, b);
		}
		if (x1 > x2) {
			data_type t = x1; x1 = x2; x2 = t;
		}
		if (y1 > y2) {
			data_type t = y1; y1 = y2; y2 = t;
		}
        ret[i].x1 = x1; ret[i].x2 = x2;
		ret[i].y1 = y1; ret[i].y2 = y2;
    }

    return ret;
}

void reset_timers() {
	reserve_tm.reset();
	init_tm.reset(); sort_tm.reset(); build_tm.reset(); total_tm.reset(); 	
}

void run_all(vector<point_type>& points, size_t iteration, data_type min_val, data_type max_val, size_t query_num) {
	string benchmark_name = "Query-All";
	RangeQuery<data_type, data_type> *r = new RangeQuery<data_type, data_type>(points);

	vector<Query_data> queries = generate_queries(query_num, min_val, max_val);
	
	size_t counts[query_num];

	timer t_query_total;
	t_query_total.start();
	cilk_for (size_t i = 0; i < query_num; i++) {
	  vector<pair<int,int>> out = r->query_range(queries[i].x1, queries[i].y1, queries[i].x2, queries[i].y2);
	  counts[i] = out.size();
	}

	t_query_total.stop();
	size_t total = pbbs::reduce_add(sequence<size_t>(counts,query_num));

  cout << "RESULT" << fixed << setprecision(3)
       << "\talgo=" << "RageTree"
       << "\tname=" << benchmark_name
       << "\tn=" << points.size()
       << "\tq=" << query_num
       << "\tp=" << __cilkrts_get_nworkers()
       << "\tmin-val=" << min_val
       << "\tmax-val=" << max_val
       << "\twin-mean=" << win
       << "\titeration=" << iteration
       << "\tbuild-time=" << total_tm.get_total()
       << "\treserve-time=" << reserve_tm.get_total()
       << "\tquery-time=" << t_query_total.get_total()
       << "\ttotal=" << total
       << endl;

     reset_timers();

     delete r;
}


void run_sum(vector<point_type>& points, size_t iteration, data_type min_val, data_type max_val, size_t query_num) {
	string benchmark_name = "Query-Sum";
	RangeQuery<data_type, data_type> *r = new RangeQuery<data_type, data_type>(points);

	vector<Query_data> queries = generate_queries(query_num, min_val, max_val);
	
  size_t counts[query_num];

	timer t_query_total;
	t_query_total.start();
	cilk_for (size_t i = 0; i < query_num; i++) {
	  counts[i] = r->query_sum(queries[i].x1, queries[i].y1,
	  			   queries[i].x2, queries[i].y2);
	}

	t_query_total.stop();
	size_t total = pbbs::reduce_add(sequence<size_t>(counts,query_num));

  cout << "RESULT" << fixed << setprecision(3)
       << "\talgo=" << "RageTree"
       << "\tname=" << benchmark_name
       << "\tn=" << points.size()
       << "\tq=" << query_num
       << "\tp=" << __cilkrts_get_nworkers()
       << "\tmin-val=" << min_val
       << "\tmax-val=" << max_val
       << "\twin-mean=" << win
       << "\titeration=" << iteration
       << "\tbuild-time=" << total_tm.get_total()
       << "\treserve-time=" << reserve_tm.get_total()
       << "\tquery-time=" << t_query_total.get_total()
       << "\ttotal=" << total
       << endl;

     reset_timers();
     delete r;
}


int main(int argc, char** argv) {

    if (argc != 9) {
      cout << argv[0] << " <n> <min_val> <max_val> <rounds> <queries> <dist> <window> <query_type>" << std::endl;
	  cout << "dist=1 to specify query window size" << endl << "query_type=1 for query-all, query_type=0 for query-sum" << endl;
      exit(1);
    }

	srand(2017);

    size_t n = str_to_int(argv[1]);
    data_type min_val  = str_to_int(argv[2]);
    data_type max_val  = str_to_int(argv[3]);

    size_t iterations = str_to_int(argv[4]);
  	dist = str_to_int(argv[6]);
  	win = str_to_int(argv[7]);
  	int type = str_to_int(argv[8]);
  	size_t query_num  = str_to_int(argv[5]);
	
	for (size_t i = 0; i < iterations; ++i) {
	    vector<point_type> points = generate_points(n, min_val, max_val);
	    if (type == 0) run_all(points, i, min_val, max_val, query_num);
	    else run_sum(points, i, min_val, max_val, query_num);
	}

    return 0;
}

