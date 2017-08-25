#include <algorithm>
#include <iostream>
#include <cstring>
#include <vector>
#include <cctype>
#include <fstream>
#include <set>
#include "augmented_map.h"
#include "pbbs-include/get_time.h"
#include "pbbs-include/par_string.h"
#include "weighted_index.h"

using namespace std;
using namespace pbbs;

long str_to_long(char* str) {
    return strtol(str, NULL, 10);
}

using word = char*;
using post_elt = inv_index::post_elt;
using post_list = inv_index::post_list;
using index_elt = inv_index::index_elt;
using index_pair = pair<word,post_list>;

pair<pair<index_elt*,size_t>, char*> parse(string filename, size_t max_size) {
  startTime();
  sequence<char> Str = read_string_from_file(filename, 0, max_size);
  size_t n = Str.size();
  nextTime("read file");

  parallel_for (size_t i=0; i <n; i++) {
    Str[i] = tolower(Str[i]);
    if (!islower(Str[i])) Str[i] = 0;
  }

  sequence<char*> W = string_to_words(Str);
  size_t num_words = W.size();
  nextTime("find words");

  sequence<bool> start_flag(num_words);
  start_flag[num_words-1] = 0;
  char doc[] = "doc";
  char ids[] = "id";
  
  parallel_for(size_t i = 0; i < num_words-1; i++) {
    if ((strcmp(W[i],doc) == 0) && strcmp(W[i+1],ids) == 0)
      start_flag[i] = 1;
    else start_flag[i] = 0;
  }

  sequence<size_t> I = pack_index<size_t>(start_flag);
  size_t total_docs = I.size();
  cout << "Words = " << num_words << endl;
  cout << "Documents = " << total_docs << endl;
  nextTime("find documents");

  size_t header_size = 2;
  size_t total_pairs = num_words - total_docs * header_size;
  cout << "Pairs = " << total_pairs << endl;

  index_elt* KV = pbbs::new_array_no_init<index_elt>(total_pairs);
  parallel_for (size_t i = 0; i < total_docs; ++i) {
    size_t start = I[i];
    size_t end = (i == (total_docs-1)) ? num_words : I[i+1];
    size_t len = end - start - header_size;
    size_t start_out = start - (header_size * i);
    for (size_t j = 0; j < len; j++) {
      KV[start_out+j] = index_elt(W[start+header_size+j], 
				  post_elt(i,1.0));
    }
  }
  nextTime("create pairs");

  return make_pair(make_pair(KV,total_pairs),Str.as_array());
}

int main(int argc, char** argv) {

  pbbs::random r(0);
  bool write_output = 0;
  size_t max_chars = 1000000;
  size_t num_queries = 1000;
  if (argc > 1) max_chars = str_to_long(argv[1]);
  if (argc > 2) num_queries = str_to_long(argv[2]);
  if (argc > 3) write_output = 1;

  string fname = "/usr3/data/wikipedia/wikipedia.txt";
  //string fname = "/usr0/home/data/wikipedia/wikipedia.txt";
  //string fname = "/usr0/home/danielf/wikipedia.txt";
  pair<pair<index_elt*,size_t>,char*> X = parse(fname, max_chars);
  index_elt* KV = X.first.first;
  size_t n = X.first.second;

  pair<char*,char*>* test_words = new pair<char*,char*>[num_queries];
  size_t* size_in = new size_t[num_queries];
  size_t* size_out = new size_t[num_queries];

  inv_index test_idx(KV, KV +n);
  nextTime("build");

  // filters out any words which appear fewer than 100 times
  inv_index::index common_idx = test_idx.idx;
  common_idx.filter([n] (index_pair e) {
                       return (e.second.size() > 100); });
  cout << "filter size = " << common_idx.size() << endl;

  //size_t num_words = test_idx.idx.size();
  size_t num_common_words = common_idx.size();

  parallel_for(size_t i =0; i < num_queries; i++) {
    test_words[i].first = KV[r.ith_rand(i)%(n-1)].first;
    test_words[i].second = common_idx.select(r.ith_rand(num_queries+i)%num_common_words).first;
  }
  nextTime("waste");
  
  parallel_for(size_t i =0; i < num_queries; i++) {
    post_list l1 = test_idx.get_list(test_words[i].first);
    post_list l2 = test_idx.get_list(test_words[i].second);
    size_in[i] = l1.size() + l2.size();
    post_list l3 = test_idx.And(l1,l2);
    vector<post_elt> r = test_idx.top_k(l3,10);
    size_out[i] = l3.size() + r.size();
  }
  
  nextTime("query");
  size_t total_in = 0;
  size_t total_out = 0;
  
  for (size_t i =0; i < num_queries; i++) {
    if (i < 0) 
      cout << test_words[i].first << ", " << test_words[i].second << endl;
    total_in += size_in[i];
    total_out += size_out[i];
  }
  
  cout << "total in = " << total_in << endl;
  cout << "total out = " << total_out << endl;

  cout << "Unique words = " << test_idx.idx.size() << endl;
  cout << "Map: ";  inv_index::index::print_allocation_stats();
  cout << "Set: ";  post_list::print_allocation_stats();

  delete[] KV;

  if (write_output) {
    vector<pair<char*,post_list> >out;
    test_idx.idx.content(std::back_inserter(out));

    FILE* x = freopen("sol.out", "w", stdout);
    if (x == NULL) abort();

    for (size_t i = 0; i < out.size(); i++) { 
      std::cout << (out[i].first) << " " << (out[i].second).size() << std::endl;
    }
  }

  // free the original character array
  free(X.second);
  return 0;
}
