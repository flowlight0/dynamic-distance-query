#include <iostream>
#include "jlog.hpp"
#include "common.hpp"
#include "algorithm/distance_without_bit_parallel.hpp"
#include "algorithm/distance_with_bit_parallel.hpp"
#include "gflags/gflags.h"
#include <unordered_map>
#include <numeric>
#include <algorithm>
using namespace std;

DEFINE_string(graph_file, "-", "Input graph file");
DEFINE_int32(num_queries, 10000, "The number of vertex pairs.");
DEFINE_int32(num_updates,  100, "The number of dynamic updates.");
DEFINE_int32(num_trials ,  3, "The number of trials.");
DEFINE_int32(max_i, 5, "");


int main(int argc, char *argv[]) {
  JLOG_INIT(&argc, argv);
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  vector<pair<int, int> > es;
  vector<int> sources, targets, answers;
  int n;
  JLOG_PUT_BENCHMARK("setting.read_time") {
    n = ReadGraph(FLAGS_graph_file, es);
  }
  
  JLOG_PUT("setting.dataset", FLAGS_graph_file);
  JLOG_PUT("setting.num_vs", n);
  JLOG_PUT("setting.num_es", es.size());
  JLOG_PUT("setting.num_queries", FLAGS_num_queries);

  srand(0);
  for (int q = 0; q < FLAGS_num_queries; q++) {
    sources.push_back(rand() % n);
    targets.push_back(rand() % n);
  }

  // prepare distances for validation
  JLOG_PUT_BENCHMARK("time.prepare_answer") {
    GraphDistance<20> bps(n, es);
    for (int q = 0; q < FLAGS_num_queries; q++) {
      answers.push_back(bps.QueryDistance(sources[q], targets[q]));
    }
  }


  vector<int> dists(FLAGS_num_queries);
  vector<vector<pair<int, int > > > del_ess, ins_ess;
  for (int trial = 0; trial < FLAGS_num_trials; trial++) {
    random_shuffle(es.begin(), es.end());
    vector<pair<int, int> > upd_es(es.end() - FLAGS_num_updates, es.end());
    del_ess.push_back(upd_es);
    reverse(upd_es.begin(), upd_es.end());
    ins_ess.push_back(upd_es);
  }


  for (int i = 0; i <= FLAGS_max_i; i++) {
    string is = "reduction_" + to_string(i);
    
    JLOG_OPEN(is.c_str()) {
      GraphDistance<20, BitParallelLandmarks<20, DegreeLS<64>> > bpd;
      JLOG_PUT_BENCHMARK("time.construct") {
        bpd.Build(n, es, i);
      }
      JLOG_PUT("index_size", bpd.IndexSize());
      JLOG_PUT("graph_size", bpd.GraphSize());
      JLOG_PUT("memory_usage", jlog_internal::get_memory_usage());
      JLOG_PUT("initial_degree", bpd.DegreeSum());

      for (int trial = 0; trial < FLAGS_num_trials; trial++) {
        double q_start = jlog_internal::get_current_time_sec();
        JLOG_ADD_BENCHMARK("time.query.total") {
          for (int q = 0; q < FLAGS_num_queries; q++) {
            dists[q] = bpd.QueryDistance(sources[q], targets[q]);
          }
        }
        double q_time = jlog_internal::get_current_time_sec() - q_start;
        JLOG_ADD("time.query.average", q_time / FLAGS_num_queries);
        
        {
          // delete edges one by one. 
          double start = jlog_internal::get_current_time_sec();
          JLOG_ADD_BENCHMARK("time.delete.total") {
            for (const auto &e : del_ess[trial]) {
              bpd.DeleteEdge(e.fst, e.snd);
              bpd.DeleteEdge(e.snd, e.fst);
            }
          }
          double etime = jlog_internal::get_current_time_sec() - start;
          JLOG_ADD("time.delete.average", etime / FLAGS_num_updates);
        }
        {
          // insert edges in the reverse order. 
          double start = jlog_internal::get_current_time_sec();
          JLOG_ADD_BENCHMARK("time.insert.total") {
            for (const auto &e : ins_ess[trial]) {
              bpd.InsertEdge(e.fst, e.snd);
              bpd.InsertEdge(e.snd, e.fst);
            }
          }
          double etime = jlog_internal::get_current_time_sec() - start;
          JLOG_ADD("time.insert.average", etime / FLAGS_num_updates);
        }
        JLOG_ADD("result.degree_sum", bpd.DegreeSum());
        
        // verification.
        for (int q = 0; q < FLAGS_num_queries; q++) {
          CHECK_EQ(answers[q], dists[q]);
        }
      }
    }
  }
  return 0;
}
