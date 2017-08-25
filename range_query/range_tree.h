#include <algorithm>
#include <iostream>
#include <cstring>
#include <vector>
#include <parallel/algorithm>
#include "tree_operations.h"
#include "augmented_map.h"
#include "pbbs-include/get_time.h"

using namespace std;

timer init_tm, build_tm, total_tm, aug_tm, sort_tm, reserve_tm, outmost_tm, globle_tm, g_tm;

template<typename c_type, typename w_type>
struct Point {
    c_type x, y;
	w_type w;
    Point(c_type _x, c_type _y, w_type _w) : x(_x), y(_y), w(_w){}
    Point(){}
};

template<typename value_type>
inline bool inRange(value_type x, value_type l, value_type r) {
	return ((x>=l) && (x<=r));
}

template<typename c_type, typename w_type>
struct RangeQuery {
	using x_type = c_type;
	using y_type = c_type;
	using point_type = Point<c_type, w_type>;
	using point_x = pair<x_type, y_type>;
	using point_y = pair<y_type, x_type>;
	
	struct aug_add {
		typedef w_type aug_t;
		inline static aug_t from_entry(point_y k, w_type v) { return v; }
		inline static aug_t combine(aug_t a, aug_t b) { return a+b; }
		static aug_t get_empty() { return 0;}
	};
	using sec_aug = augmented_map<point_y, w_type, aug_add>;
	
	struct aug_union {
		typedef sec_aug aug_t;
		inline static aug_t from_entry(point_x k, w_type v) {
			return aug_t(make_pair(make_pair(k.second, k.first), v));
		}
		inline static aug_t combine(aug_t a, aug_t b) {
			return map_union(a, b, [](w_type x, w_type y) {return x+y;});
		}
		static aug_t get_empty() { return aug_t();}
	};

	using main_aug = augmented_map<point_x, w_type, aug_union>;
	using aug_node = typename main_aug::node_type;
	using sec_node = typename sec_aug::node_type;
	using main_entry = pair<point_x, w_type>;
	typedef tree_operations<sec_node>   tree_ops;
	
    RangeQuery(vector<point_type>& points) {
		construct(points);
    }

    ~RangeQuery() {
        main_aug::finish();
		sec_aug::finish();
    }
	
    void construct(vector<point_type>& points) {
        const size_t n = points.size();
		
		reserve_tm.start();
		main_aug::reserve(n);
		sec_aug::reserve(24*n);
		reserve_tm.stop();
		
		pair<point_x, w_type> *pointsEle = new pair<point_x, w_type>[n];
		
		total_tm.start();
		
        cilk_for (size_t i = 0; i < n; ++i) 
			pointsEle[i] = make_pair(make_pair(points[i].x, points[i].y), points[i].w);
		
        range_tree = main_aug(pointsEle, pointsEle + n);

		delete[] pointsEle;
        total_tm.stop();
    }

  
	struct count_t {
		point_y y1, y2;
		int r;
		count_t(point_y y1, point_y y2) : y1(y1), y2(y2), r(0) {}
		void add_entry(point_x p, w_type wp) {
			if (p.second >= y1.first && p.second <= y2.first) r += 1;
		}
		void add_aug_val(sec_aug a) { 
		  r += (a.rank(y2) - a.rank(y1));}
	};
	  
	w_type query_count(x_type x1, y_type y1, x_type x2, y_type y2) {
		count_t qrc(make_pair(y1,x1), make_pair(y2,x2));
		range_tree.range_sum(make_pair(x1,y1), make_pair(x2,y2), qrc);
		return qrc.r;
	}

	struct sum_t {
	  point_y y1, y2;
	  int r;
	  sum_t(point_y y1, point_y y2) : y1(y1), y2(y2), r(0) {}
	  void add_entry(point_x p, w_type wp) {
		  if (p.second >= y1.first && p.second <= y2.first) r += wp;
	  }
	  void add_aug_val(sec_aug a) { r += a.aug_range(y1, y2); }
	};
	  
	w_type query_sum(x_type x1, y_type y1, x_type x2, y_type y2) {
		sum_t qrs(make_pair(y1,x1), make_pair(y2,x2));
		range_tree.range_sum(make_pair(x1,y1), make_pair(x2,y2), qrs);
		return qrs.r;
	}

	template<typename OutIter>
	struct range_t {
		point_y y1, y2;
		OutIter out;
		range_t(point_y y1, point_y y2, OutIter out) : y1(y1), y2(y2), out(out) {}
		void add_entry(point_x p, w_type wp) {
			if (p.second >= y1.first && p.second <= y2.first) {out = p; ++out;}
		}
		void add_aug_val(sec_aug a) { a.keys_in_range(y1, y2, out);}
	};

	template<typename OutIter>
	void query_range_iterate(x_type x1, y_type y1, x_type x2, y_type y2, OutIter out) {
		range_t<OutIter> qr(make_pair(y1,x1), make_pair(y2,x2), out);
		range_tree.range_sum(make_pair(x1,y1), make_pair(x2,y2), qr);
	}

	vector<point_y> query_range(x_type x1, y_type y1, x_type x2, y_type y2) {
		vector<point_y> out;
		size_t n = query_count(x1,y1,x2,y2);
		out.reserve(n);
		query_range_iterate(x1, y1, x2, y2, std::back_inserter(out));
		return out;
	}

	main_aug range_tree;
};
