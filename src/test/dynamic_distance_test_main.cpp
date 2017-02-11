#include <iostream>
#include "common.hpp"
#include "algorithm/distance_without_bit_parallel.hpp"
#include "algorithm/distance_with_bit_parallel.hpp"
#include "algorithm/bit_parallel_landmarks.hpp"
#include <algorithm>
#include "gtest/gtest.h"
#include "../common.hpp"
#include <algorithm>
using namespace std;

namespace {
  const int NUM_VERTICES_TINY = 4;
  const int NUM_VERTICES_SMALL = 10;
  const int NUM_VERTICES_MEDIUM = 50;
  const int DEFAULT_NUM_TRIALS = 10;
  const int DEFAULT_NUM_QUERIES = 1000;
  const int DEFAULT_NUM_UPDATES = 10;
  
  string Log(const vector<string> &logs){
    ostringstream oss;
    for (const string &s : logs) oss << s << endl;
    return oss.str();
  }

  void DeleteEdge(vector<vector<int> > G[2], const pair<int, int> &e){
    G[0][e.fst].erase(find(ALL(G[0][e.fst]), e.snd));
    G[1][e.snd].erase(find(ALL(G[1][e.snd]), e.fst));
  }
}

void TestInsertBPLS(int n, vector<pair<int, int> > es, int iset_param){
  random_shuffle(es.begin(), es.end());

  vector<vector<int> > G[2];
  G[0].resize(n);
  G[1].resize(n);
  
  const auto mid = es.end() - min(size_t(DEFAULT_NUM_UPDATES), es.size());
  vector<pair<int, int> > tmp_es(es.begin(), mid);
  for (const auto &e : tmp_es){
    G[0][e.fst].push_back(e.snd);
    G[1][e.snd].push_back(e.fst);
  }
  vector<pair<int, int> > upd_es(mid, es.end());

  BitParallelLandmarks<1> bpls;
  bpls.Build(G, iset_param);

  vector<string> log;
  for (auto e : upd_es){
    const auto masks_before = bpls.GetNeighborSets();
    log.push_back("last: (" + to_string(e.fst) + ", " + to_string(e.snd) + ")");

    G[0][e.fst].push_back(e.snd);
    G[1][e.snd].push_back(e.fst);
    bpls.InsertEdge(e.fst, e.snd);
    ASSERT_EQ(bpls.GetDistancesFromScratch(), bpls.GetDistances())
      << G[0] << endl << G[1] << endl << Log(log) << endl;
    
    ASSERT_EQ(bpls.GetNeighborSetsFromScratch(), bpls.GetNeighborSets())
      << G[0] << endl << G[1] << endl << Log(log) << endl
      << bpls.GetRoots() << " " << bpls.GetRootNeighbors() << " " << bpls.GetDistances() << endl
      << "before: " << masks_before << endl;

    ASSERT_EQ(bpls.GetDistances(), bpls.GetDistancesFromScratch());
  }
}


void TestDeleteBPLS(int n, vector<pair<int, int> > es, int iset_param){
  random_shuffle(es.begin(), es.end());

  vector<vector<int> > G[2];
  G[0].resize(n);
  G[1].resize(n);
  vector<pair<int, int> > del_es(es.begin(), es.begin() + min(es.size(), size_t(DEFAULT_NUM_UPDATES)));

  for (const auto &e : es){
    G[0][e.fst].push_back(e.snd);
    G[1][e.snd].push_back(e.fst);
  }
  BitParallelLandmarks<1> bpls;
  bpls.Build(G, iset_param);

  vector<string> log;
  for (auto e : del_es){
    const auto masks_before = bpls.GetNeighborSets();
    log.push_back("last: (" + to_string(e.fst) + ", " + to_string(e.snd) + ")");
    DeleteEdge(G, e);
    bpls.DeleteEdge(e.fst, e.snd);
    ASSERT_EQ(bpls.GetDistancesFromScratch(), bpls.GetDistances())
      << G[0] << endl << G[1] << endl << Log(log) << endl;
    
    ASSERT_EQ(bpls.GetNeighborSetsFromScratch(), bpls.GetNeighborSets())
      << G[0] << endl << G[1] << endl << Log(log) << endl
      << bpls.GetRoots() << " " << bpls.GetRootNeighbors() << " " << bpls.GetDistances() << endl
      << "before: " << masks_before << endl;
  }
}


void TestInsertBPD(int n, vector<pair<int, int> > es, int iset_param){
  random_shuffle(es.begin(), es.end());

  vector<vector<int> > G(n);
  GraphDistance<1> dist_d;
  const auto mid = es.end() - min(size_t(DEFAULT_NUM_UPDATES), es.size());
  vector<pair<int, int> > tmp_es(es.begin(), mid);
  for (const auto &e : tmp_es){
    G[e.fst].push_back(e.snd);
  }
  for (int i = 0; i < n; i++) {
    tmp_es.emplace_back(i, i);
  }
  
  vector<pair<int, int> > upd_es(mid, es.end()); 
  dist_d.Build(n, tmp_es, iset_param);
  vector<string> log;

  for (const auto &e : upd_es){
    G[e.fst].push_back(e.snd);
    GraphDistance<1> dist_s;

    tmp_es.push_back(e);
    dist_s.Build(n, tmp_es, iset_param);
    dist_d.InsertEdge(e.fst, e.snd);

    log.push_back("last: (" + to_string(e.fst) + ", " + to_string(e.snd) + ")");
    log.push_back("static: " + dist_s.GetLandmarks());
    log.push_back("dynamic: " + dist_d.GetLandmarks());

    for (int u = 0; u < n; u++) {
      for (int v = 0; v < n; v++) {
        ASSERT_EQ(dist_s.QueryDistance(u, v), dist_d.QueryDistance(u, v))
          << u << " " << v << " " << G << endl
          << "upperbound: " << dist_s.QueryDistanceUB(u, v) << " " << dist_d.QueryDistanceUB(u, v) << endl << Log(log);
      }
    }
  }
}

void TestDeleteBPD(int n, vector<pair<int, int> > es, int iset_param){
  random_shuffle(es.begin(), es.end());
  const size_t num_es = es.size();
  for (int i = 0; i < n; i++) {
    es.emplace_back(i, i);
  }
  
  vector<pair<int, int> > del_es(es.begin(), es.begin() + min(num_es, size_t(DEFAULT_NUM_UPDATES)));
  GraphDistance<1> dist_d;
  dist_d.Build(n, es, iset_param);

  vector<string> log;
  for (size_t i = 0; i < del_es.size(); i++) {
    auto e = del_es[i];

    vector<pair<int, int> > tmp_es(es.begin() + i + 1, es.end());
    GraphDistance<1> dist_s;
    dist_s.Build(n, tmp_es, iset_param);
    dist_d.DeleteEdge(e.fst, e.snd);

    vector<vector<int> > G(n);
    for (const auto e_ : tmp_es){
      if (e_.fst != e_.snd) G[e_.fst].push_back(e_.snd);
    }

    log.push_back("last: (" + to_string(e.fst) + ", " + to_string(e.snd) + ")");
    log.push_back("static: " + dist_s.GetLandmarks());
    log.push_back("dynamic: " + dist_d.GetLandmarks());
    for (int u = 0; u < n; u++){
      for (int v = 0; v < n; v++){
        ASSERT_EQ(dist_s.QueryDistance(u, v), dist_d.QueryDistance(u, v))
          << u << " " << v << " " << G << endl
          << "upperbound: " << dist_s.QueryDistanceUB(u, v) << " " << dist_d.QueryDistanceUB(u, v) << endl << Log(log) << endl
          << "D: " << dist_d.GetNeighborSetsWithRecovery() << endl
          << "S: " << dist_s.GetNeighborSetsWithRecovery() << endl;
      }
    }
  }
}

void TestDeleteAndInsertBPD(int n, vector<pair<int, int> > es, int iset_param){
  random_shuffle(es.begin(), es.end());
  const size_t num_es = es.size();
  for (int i = 0; i < n; i++) {
    es.emplace_back(i, i);
  }
  
  vector<pair<int, int> > del_es(es.begin(), es.begin() + min(num_es, size_t(DEFAULT_NUM_UPDATES)));
  GraphDistance<1> dist_d;
  dist_d.Build(n, es, iset_param);

  vector<string> log;
  for (size_t i = 0; i < del_es.size(); i++) {
    const auto e = del_es[i];
    
    vector<pair<int, int> > tmp_es(es.begin() + i + 1, es.end());
    GraphDistance<1> dist_s;
    dist_s.Build(n, tmp_es, iset_param);
    dist_d.DeleteEdge(e.fst, e.snd);

    vector<vector<int> > G(n);
    for (const auto e_ : tmp_es){
      if (e_.fst != e_.snd) G[e_.fst].push_back(e_.snd);
    }

    log.push_back("last: (" + to_string(e.fst) + ", " + to_string(e.snd) + ")");
    log.push_back("static: " + dist_s.GetLandmarks());
    log.push_back("dynamic: " + dist_d.GetLandmarks());
    for (int u = 0; u < n; u++){
      for (int v = 0; v < n; v++){
        ASSERT_EQ(dist_s.QueryDistance(u, v), dist_d.QueryDistance(u, v))
          << u << " " << v << " " << G << endl
          << "upperbound: " << dist_s.QueryDistanceUB(u, v) << " " << dist_d.QueryDistanceUB(u, v) << endl
          << Log(log) << endl << dist_d.GetNeighborSetsWithRecovery() << endl;

      }
    }
  }
  
  for (int i = del_es.size() - 1; i >= 0; i--){
    auto e = del_es[i];
    vector<pair<int, int> > tmp_es(es.begin() + i, es.end());
    GraphDistance<1> dist_s;
    dist_s.Build(n, tmp_es, iset_param);
    dist_d.InsertEdge(e.fst, e.snd);

    vector<vector<int> > G(n);
    for (const auto e_ : tmp_es){
      if (e_.fst != e_.snd) G[e_.fst].push_back(e_.snd);
    }

    log.push_back("last: (" + to_string(e.fst) + ", " + to_string(e.snd) + ")");
    log.push_back("static: " + dist_s.GetLandmarks());
    log.push_back("dynamic: " + dist_d.GetLandmarks());
    for (int u = 0; u < n; u++){
      for (int v = 0; v < n; v++){
        ASSERT_EQ(dist_s.QueryDistance(u, v), dist_d.QueryDistance(u, v))
          << u << " " << v << " " << G << endl
          << "upperbound: " << dist_s.QueryDistanceUB(u, v) << " " << dist_d.QueryDistanceUB(u, v) << endl
          << Log(log) << endl << dist_d.GetNeighborSetsWithRecovery() << endl;
      }
    }
  }
}


template <int num_ts, int num_vs, int iset_param>
void ExamineInsertBPLSRandomGraph(){
  for (int trial = 0; trial < num_ts; trial++) {
    srand(trial);
    TestInsertBPLS(num_vs, GenerateErdosRenyi(num_vs, 0.3), iset_param);
  }
};

template <int num_ts, int num_vs, int iset_param>
void ExamineInsertBPDRandomGraph(){
  for (int trial = 0; trial < num_ts; trial++) {
    srand(trial);
    TestInsertBPD(num_vs, GenerateErdosRenyi(num_vs, 0.3), iset_param);
  }
};

template <int num_ts, int num_vs, int iset_param>
void ExamineDeleteBPLSRandomGraph(){
  for (int trial = 0; trial < num_ts; trial++) {
    srand(trial);
    TestDeleteBPLS(num_vs, GenerateErdosRenyi(num_vs, 0.4), iset_param);
  }
};


template <int num_ts, int num_vs, int iset_param>
void ExamineDeleteBPDRandomGraph(){
  for (int trial = 0; trial < num_ts; trial++) {
    srand(trial);
    TestDeleteBPD(num_vs, GenerateErdosRenyi(num_vs, 0.4), iset_param);
  }
};

template <int num_ts, int num_vs, int iset_param>
void ExamineDeleteAndInsertBPDRandomGraph(){
  for (int trial = 0; trial < num_ts; trial++) {
    srand(trial);
    TestDeleteAndInsertBPD(num_vs, GenerateErdosRenyi(num_vs, 0.4), iset_param);
  }
};


TEST(BPLS_EDGE_INSERT_WITHOUT_DEPENDENT_SET, TINY_RANDOM){
  ExamineInsertBPLSRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_TINY, 0>();
}

TEST(BPLS_EDGE_INSERT_WITHOUT_DEPENDENT_SET, SMALL_RANDOM) {
  ExamineInsertBPLSRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_SMALL, 0>();
}

TEST(BPLS_EDGE_INSERT_WITHOUT_DEPENDENT_SET, MIDDLE_RANDOM) {
  ExamineInsertBPLSRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_MEDIUM, 0>();
}


TEST(BPLS_EDGE_INSERT, TINY_RANDOM) {
  ExamineInsertBPLSRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_TINY, 3>();
}

TEST(BPLS_EDGE_INSERT, SMALL_RANDOM) {
  ExamineInsertBPLSRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_SMALL, 3>();
}

TEST(BPLS_EDGE_INSERT, MIDDLE_RANDOM) {
  ExamineInsertBPLSRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_MEDIUM, 3>();
}

TEST(BPD_EDGE_INSERT_WITHOUT_DEPENDENT_SET, TINY_RANDOM){
  ExamineInsertBPDRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_TINY, 0>();
}

TEST(BPD_EDGE_INSERT_WITHOUT_DEPENDENT_SET, SMALL_RANDOM){
  ExamineInsertBPDRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_SMALL, 0>();
}

TEST(BPD_EDGE_INSERT_WITHOUT_DEPENDENT_SET, MIDDLE_RANDOM){
  ExamineInsertBPDRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_MEDIUM, 0>();
}


TEST(BPD_EDGE_INSERT, TINY_RANDOM){
  ExamineInsertBPDRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_TINY, 3>();
}

TEST(BPD_EDGE_INSERT, SMALL_RANDOM){
  ExamineInsertBPDRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_SMALL, 3>();
}

TEST(BPD_EDGE_INSERT, MIDDLE_RANDOM){
  ExamineInsertBPDRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_MEDIUM, 3>();
}


TEST(BPLS_EDGE_DELETE_WITHOUT_DEPENDENT_SET, TINY_RANDOM){
  ExamineDeleteBPLSRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_TINY, 0>();
}

TEST(BPLS_EDGE_DELETE_WITHOUT_DEPENDENT_SET, SMALL_RANDOM) {
  ExamineDeleteBPLSRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_SMALL, 0>();
}

TEST(BPLS_EDGE_DELETE_WITHOUT_DEPENDENT_SET, MIDDLE_RANDOM) {
  ExamineDeleteBPLSRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_MEDIUM, 0>();
}

TEST(BPLS_EDGE_DELETE, TINY_RANDOM){
  ExamineDeleteBPLSRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_TINY, 3>();
}

TEST(BPLS_EDGE_DELETE, SMALL_RANDOM) {
  ExamineDeleteBPLSRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_SMALL, 3>();
}

TEST(BPLS_EDGE_DELETE, MIDDLE_RANDOM) {
  ExamineDeleteBPLSRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_MEDIUM, 3>();
}

TEST(BPD_EDGE_DELETE_WITHOUT_DEPENDENT_SET, TINY_RANDOM){
  ExamineDeleteBPLSRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_TINY, 0>();
}

TEST(BPD_EDGE_DELETE_WITHOUT_DEPENDENT_SET, SMALL_RANDOM) {
  ExamineDeleteBPDRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_SMALL, 0>();
}

TEST(BPD_EDGE_DELETE_WITHOUT_DEPENDENT_SET, MIDDLE_RANDOM) {
  ExamineDeleteBPDRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_MEDIUM, 0>();
}

TEST(BPD_EDGE_DELETE, TINY_RANDOM){
  ExamineDeleteBPDRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_TINY, 3>();
}

TEST(BPD_EDGE_DELETE, SMALL_RANDOM) {
  ExamineDeleteBPDRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_SMALL, 3>();
}

TEST(BPD_EDGE_DELETE, MIDDLE_RANDOM) {
  ExamineDeleteBPDRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_MEDIUM, 3>();
}


TEST(BPD_EDGE_DELETE_INSERT, TINY_RANDOM){
  ExamineDeleteAndInsertBPDRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_TINY, 3>();
}

TEST(BPD_EDGE_DELETE_INSERT, SMALL_RANDOM) {
  ExamineDeleteAndInsertBPDRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_SMALL, 3>();
}

TEST(BPD_EDGE_DELETE_INSERT, MIDDLE_RANDOM) {
  ExamineDeleteAndInsertBPDRandomGraph<DEFAULT_NUM_TRIALS, NUM_VERTICES_MEDIUM, 3>();
}


int main(int argc, char **argv){
  JLOG_INIT(&argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
