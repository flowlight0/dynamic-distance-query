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
DEFINE_int32(num_bfs_queries, 100, "The number of vertex pairs.");
DEFINE_int32(num_bibfs_queries, 10000, "The number of vertex pairs.");
DEFINE_int32(num_trials, 3, "The number of vertex pairs.");

void OutputStats(const vector<int> dists, const vector<double> times, const vector<QueryStats> &stats){
  CHECK_EQ(dists.size(), times.size());
  CHECK_EQ(dists.size(), stats.size());
  for (size_t q = 0; q < dists.size(); q++) {
    JLOG_ADD("dists", dists[q]);
    JLOG_ADD("times", times[q]);
    JLOG_ADD("visit_vs", stats[q].visit_vs);
    JLOG_ADD("visit_es", stats[q].visit_es);
  }
  JLOG_PUT("time.average", accumulate(times.begin(), times.end(), 0.0) / times.size());
}

int main(int argc, char *argv[])
{
  JLOG_INIT(&argc, argv);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  int n = 0;
  vector<pair<int ,int> > es;

  JLOG_PUT_BENCHMARK("setting.readtime"){
    n = ReadGraph(FLAGS_graph_file, es);
  }
  JLOG_PUT("setting.num_vs", n);
  JLOG_PUT("setting.num_es", es.size());
  JLOG_PUT("setting.num_ls", 20);
  JLOG_PUT("setting.iset_param", 0);

  for (int trial = 0; trial < FLAGS_num_trials; trial++) {
    vector<int> dist1(FLAGS_num_bfs_queries);
    vector<int> dist2(FLAGS_num_bibfs_queries);
    vector<int> dist3(FLAGS_num_bibfs_queries);

    srand(trial);
    vector<int> sources;
    vector<int> targets;
    for (int i = 0; i < std::max(FLAGS_num_bfs_queries, FLAGS_num_bibfs_queries); i++) {
      sources.push_back(rand() % n);
      targets.push_back(rand() % n);
    }

    {
      Distance<0> naive(n, es);
      JLOG_ADD_OPEN("bfs") {
        JLOG_PUT("num_queries", FLAGS_num_bfs_queries);
        vector<double> times;
        vector<QueryStats> stats;
        
        JLOG_PUT_BENCHMARK("time.total") {
          for (int q = 0; q < FLAGS_num_bfs_queries; q++) {
            double start = jlog_internal::get_current_time_sec();
            dist1[q] = naive.QueryDistanceBFS(sources[q], targets[q]);
            double etime = jlog_internal::get_current_time_sec() - start;
            times.push_back(etime);
            stats.push_back(naive.GetQueryStats());
          }
        }
        OutputStats(dist1, times, stats);
      }

      JLOG_ADD_OPEN("bibfs.alternate") {
        
        JLOG_PUT("num_queries", FLAGS_num_bibfs_queries);
        vector<double> times;
        vector<QueryStats> stats;

        JLOG_PUT_BENCHMARK("time.total") {
          for (int q = 0; q < FLAGS_num_bibfs_queries; q++) {
            double start = jlog_internal::get_current_time_sec();
            dist2[q] = naive.QueryDistanceAlternate(sources[q], targets[q]);
            double etime = jlog_internal::get_current_time_sec() - start;
            times.push_back(etime);
            stats.push_back(naive.GetQueryStats());
          }
        }
        OutputStats(dist2, times, stats);

        for (int q = 0; q < min(FLAGS_num_bibfs_queries, FLAGS_num_bfs_queries); q++) {
          CHECK_EQ(dist1[q], dist2[q]);
        }
      }

      JLOG_ADD_OPEN("bibfs.immediate_return") {
        JLOG_PUT("num_queries", FLAGS_num_bibfs_queries);
        vector<double> times;
        vector<QueryStats> stats;

        JLOG_PUT_BENCHMARK("time.total") {
          for (int q = 0; q < FLAGS_num_bibfs_queries; q++) {
            double start = jlog_internal::get_current_time_sec();
            dist3[q] = naive.QueryDistance(sources[q], targets[q]);
            double etime = jlog_internal::get_current_time_sec() - start;
            times.push_back(etime);
            stats.push_back(naive.GetQueryStats());
          }
        }
        OutputStats(dist3, times, stats);

        for (int q = 0; q < FLAGS_num_bibfs_queries; q++) {
          CHECK_EQ(dist2[q], dist3[q]);
        }
      }
    }

    JLOG_ADD_OPEN("bibfs.bpbfs") {
      GraphDistance<20, BitParallelLandmarks<20, DegreeLS<64>>> bp;
      JLOG_PUT_BENCHMARK("time.construct") {
        bp.Build(n, es, 2);
      }

      JLOG_PUT("num_queries", FLAGS_num_bibfs_queries);
      vector<double> times;
      vector<pair<long long, long long> > visit;
      vector<QueryStats> stats;
      vector<bool> goal;
      vector<long long> prune;

      JLOG_PUT_BENCHMARK("time.total") {
        for (int q = 0; q < FLAGS_num_bibfs_queries; q++) {
          double start = jlog_internal::get_current_time_sec();
          dist3[q] = bp.QueryDistance(sources[q], targets[q]);
          double stop = jlog_internal::get_current_time_sec();
          times.push_back(stop - start);
          stats.push_back(bp.GetQueryStats());
          goal.push_back(bp.GetQueryStats().goal);
          prune.push_back(bp.GetQueryStats().num_prune);
        }
      }

      OutputStats(dist3, times, stats);
      JLOG_PUT("prune.average", accumulate(prune.begin(), prune.end(), 0.0) / prune.size());
      for (int q = 0; q < FLAGS_num_bibfs_queries; q++) {
        JLOG_ADD("goals", stats[q].goal);
        JLOG_ADD("sources", sources[q]);
        JLOG_ADD("targets", targets[q]);
        JLOG_ADD("dist_ub", bp.QueryDistanceUB(sources[q], targets[q]) );
      }

      for (int q = 0; q < FLAGS_num_bibfs_queries; q++) {
        CHECK_EQ(dist2[q], dist3[q]);
      }
    }
  }
  return 0;
}
