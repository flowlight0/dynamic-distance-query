#include <iostream>
#include "common.hpp"
#include "algorithm/distance_without_bit_parallel.hpp"
#include "algorithm/distance_with_bit_parallel.hpp"
#include "algorithm/bit_parallel_landmarks.hpp"
#include <algorithm>
#include "gtest/gtest.h"
using namespace std;

namespace {
  const int NUM_VERTICES_TINY = 5;
  const int NUM_VERTICES_SMALL = 10;
  const int NUM_VERTICES_MEDIUM = 50;
  const int DEFAULT_NUM_TRIALS = 10;
  const int DEFAULT_NUM_QUERIES = 1000;
  
  void GetUndirectedGraph(int n, const vector<pair<int, int> > &es, vector<vector<int> > G[]) {
    for (int i = 0; i < 2; i++) {
      G[i].resize(n);
      for (int v = 0; v < n; v++) {
        G[i][v].clear();
      }
    }
    
    for (const auto &e : es){
      for (int i = 0; i < 2; i++) {
        G[i][e.fst].push_back(e.snd);
        G[i][e.snd].push_back(e.fst);
      }
    }
  }

  int BFS(int s, int t, const vector<vector<int> > &G){
    int n = G.size();
    vector<uint8_t> dist(n, INF8);
    queue<int> que;
    que.push(s);
    dist[s] = 0;
    while (!que.empty()){
      int v = que.front(); que.pop();
      for (int w : G.at(v)){
        if (dist[w] == INF8){
          dist[w] = dist[v] + 1;
          que.push(w);
        }
      }
    }
    return dist[t];
  }
}

template <class DistanceOracle> 
void TestDistance(int n, const vector<pair<int, int> > &es){
  DistanceOracle search(n, es);
  vector<vector<int> > G(n);
  for (const auto &e : es){
    G[e.fst].push_back(e.snd);
  }
  
  for (int q = 0; q < DEFAULT_NUM_QUERIES; q++) {
    int s = rand() % n;
    int t = rand() % n;
    int td = BFS(s, t, G);
    ASSERT_EQ(td, search.QueryDistance(s, t))
      << s << " " << t << " "
      << G << " " << G[s] << " " << G[t] << endl;
  }
}

TEST(DISTANCE, RANDOM){
  for (int q = 0; q < DEFAULT_NUM_TRIALS; q++){
    vector<pair<int, int> > es = GenerateErdosRenyi(NUM_VERTICES_MEDIUM, 0.2);
    TestDistance<Distance<30> >(NUM_VERTICES_MEDIUM, es);
  }
}

// Verify that various bidirectional search methods that do not use landmarks return the same distances. 
TEST(DISTANCE_WITHOUT_LANDMARKS, RANDOM){
  vector<int> sources(DEFAULT_NUM_QUERIES), targets(DEFAULT_NUM_QUERIES);

  for (int i = 0; i < DEFAULT_NUM_TRIALS; i++){
    vector<pair<int, int> > es = GenerateErdosRenyi(NUM_VERTICES_MEDIUM, 0.15);
    Distance<0> naive(NUM_VERTICES_MEDIUM, es);

    for (int q = 0; q < DEFAULT_NUM_QUERIES; q++){
      sources[q] = rand() % NUM_VERTICES_MEDIUM;
      targets[q] = rand() % NUM_VERTICES_MEDIUM;
    }

    vector<vector<int> > G(NUM_VERTICES_MEDIUM);
    for (const auto &e : es){
      G[e.fst].push_back(e.snd);
    }

    for (int q = 0; q < DEFAULT_NUM_QUERIES; q++){
      int s = sources[q];
      int t = targets[q];
      int d1 = BFS(s, t, G);
      int d2 = naive.QueryDistanceBFS(s, t);
      int d3 = naive.QueryDistanceAlternate(s, t);
      int d4 = naive.QueryDistance(s, t);
      ASSERT_EQ(d1, d2);
      ASSERT_EQ(d2, d3);
      ASSERT_EQ(d3, d4);
    }
  }
}

void ExamineBitParallelLandmarksWithSizeReduction(int num_vs, int num_trials, int k){
  BitParallelLandmarks<3> bpls;
  BitParallelLandmarks<3> bpls_iset;
  
  for (int t = 0; t < num_trials; t++) {
    vector<pair<int, int> > es = GenerateErdosRenyi(num_vs, 0.2);
    vector<vector<int> > G[2];
    GetUndirectedGraph(num_vs, es, G);
    ASSERT_EQ(G[0].size(), num_vs);
    bpls     .Build(G, 0);
    bpls_iset.Build(G, k);
    for (int s = 0; s < num_vs; s++) {
      for (int t = 0; t < num_vs; t++) {
        int bpls_upperbound = bpls.GetDistanceUpperbound(s, t);
        int bpls_iset_upperbound = bpls_iset.GetDistanceUpperbound(s, t);
        ASSERT_EQ(bpls_upperbound, bpls_iset_upperbound)
          << G[0] << " " << s << " "<< t << " "
          << "bpls: " << bpls.GetRoots() << " " << bpls.GetDistances() << " "
          << "bpls_iset: " << bpls_iset.GetRoots() << " " << bpls_iset.GetDistances() << " " << endl;
      }
    }
  }
}

TEST(BIT_PARALLEL_LANDMARKS_WITH_SIZE_REDUCTION, RANDOM_TINY){
  for (int k = 0; k < 5; k++) {
    ExamineBitParallelLandmarksWithSizeReduction(NUM_VERTICES_TINY, DEFAULT_NUM_TRIALS, k);
  }
}


TEST(BIT_PARALLEL_LANDMARKS_WITH_SIZE_REDUCTION, RANDOM){
  for (int k = 0; k < 5; k++) {
    ExamineBitParallelLandmarksWithSizeReduction(NUM_VERTICES_MEDIUM, DEFAULT_NUM_TRIALS, k);
  }
}

TEST(BIT_PARALLEL_LANDMARKS_WITHOUT_LANDMARKS, RANDOM){
  for (int t = 0; t < DEFAULT_NUM_TRIALS; t++){
    vector<pair<int, int> > es = GenerateErdosRenyi(NUM_VERTICES_MEDIUM, 0.2);
    TestDistance<GraphDistance<0> > (NUM_VERTICES_MEDIUM, es);
  }
}

TEST(BP_DISTANCE, RANDOM){
  for (int t = 0; t < DEFAULT_NUM_TRIALS; t++) {
    vector<pair<int, int> > es = GenerateErdosRenyi(NUM_VERTICES_MEDIUM, 0.2);
    TestDistance<GraphDistance<30> >(NUM_VERTICES_MEDIUM, es);
  }
}

