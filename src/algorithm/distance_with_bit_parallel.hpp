#ifndef DISTANCE_WITH_BIT_PARALLEL_H
#define DISTANCE_WITH_BIT_PARALLEL_H

#include <vector>
#include <cstdlib>
#include <cstdint>
#include <limits>
#include <numeric>
#include <iostream>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <string>

#include "common.hpp"
#include "two_layer_queue.hpp"
#include "bit_parallel_landmarks.hpp"
#include "bit_parallel_landmarks_selector.hpp"
#include "gflags/gflags.h"
#include "queue.hpp"
#include "jlog.hpp"

using std::vector;
using std::pair;
using std::string;

struct QueryStats {
  long long visit_vs;
  long long visit_es;
  long long num_prune;
  long long num_prune_ub;
  int goal;
  int examined_trees;

  void Reset() { visit_vs = visit_es = num_prune = goal = examined_trees = num_prune_ub = 0; }
  friend std::ostream &operator<<(std::ostream &out, const QueryStats &stats) {
    out << "("
    << stats.visit_vs << ", "
    << stats.visit_es << ", "
    << stats.goal
    << ")";
    return out;
  }
};

template<int kNumBitParallelRoots, class upperbound_ls = BitParallelLandmarks<kNumBitParallelRoots>>
class GraphDistance {
  int V;
  vector<vector<int> > adj[2];
  vector<TwoLayerQueue<int> > que;
  upperbound_ls ubls;
  vector<uint8_t> dist[2];
  QueryStats query_stats;
private:

  void Clear() {
    V = 0;
    que.clear();
    ubls.Clear();
    for (int dir = 0; dir < 2; dir++) {
      adj[dir].clear();
      dist[dir].clear();
    }
  }

  void InitGraph(int n, const vector<pair<int, int> > &es);

  bool InsertEdgeIntoGraph(vector<vector<int> > &adj, int s, int t){
    auto iter = lower_bound(adj[s].begin(), adj[s].end(), t);
    if (iter != adj[s].end() && *iter == t){
      return false;
    } else {
      adj[s].insert(iter, t);
    }
    return true;
  }

  bool DeleteEdgeFromGraph(vector<vector<int> > &adj, int s, int t) {
    auto iter = lower_bound(adj[s].begin(), adj[s].end(), t);
    if (iter != adj[s].end() && *iter == t){
      adj[s].erase(iter);
      return true;
    } else {
      return false;
    }
  }

public:
  GraphDistance() {}

  GraphDistance(int n, const vector<pair<int, int> > &es);

  ~GraphDistance() {
    Clear();
  }

  void Build(int n, const vector<pair<int, int> > &es, int init_param);

  inline int QueryDistanceUB(int s, int t);

  inline int QueryDistance(int s, int t);

  void InsertEdge(int s, int t);

  void DeleteEdge(int s, int t);

  inline QueryStats GetQueryStats() const { return query_stats; }

  size_t IndexSize() const {
    size_t res = sizeof(adj);
    res += 2 * V * sizeof(int);
    res += 2 * V * sizeof(uint8_t);
    res += sizeof(query_stats);
    res += ubls.IndexSize();
    return res;
  }

  size_t GraphSize() const {
    return Vector2DSize(adj[0]) + Vector2DSize(adj[1]);
  }

  template<typename T>
  string to_string(const vector<T> &vec) const {
    std::ostringstream oss;
    oss << vec;
    return oss.str();
  }

  string GetLandmarks() const {
    string res;
    vector<int> roots = ubls.GetRoots();
    vector<vector<int> > dists_forward = ubls.GetDistances(0);
    vector<vector<int> > dists_backward = ubls.GetDistances(1);
    vector<vector<pair<uint64_t, uint64_t> > > masks = ubls.GetNeighborSets();
    vector<vector<int> > neighbors = ubls.GetRootNeighbors();
    res += "{";
    res += "roots: " + to_string(roots) + ", ";
    res += "neighbors: " + to_string(neighbors) + ", ";
    res += "dists: (" + to_string(dists_forward) + ", " + to_string(dists_backward) + "), ";
    res += "masks: " + to_string(masks);
    res += "}";
    return res;
  }
  vector<vector<pair<uint64_t, uint64_t> > > GetNeighborSetsWithRecovery() const {
    vector<vector<pair<uint64_t, uint64_t >>> res(kNumBitParallelRoots);
    return ubls.GetNeighborSetsWithRecovery();
  }

  vector<int> GetNeighbors(int v, int backward = 0){
    vector<int> res = adj[backward][v];
    sort(res.begin(), res.end());
    return res;
  }

  size_t DegreeSum() const {
    return ubls.DegreeSum();
  }
};

template<int kNumBitParallelRoots, class upperbound_ls>
void GraphDistance<kNumBitParallelRoots, upperbound_ls>::
InitGraph(int n, const vector<pair<int, int> > &es) {
  for (int dir = 0; dir < 2; dir++) {
    adj[dir].resize(n);
    dist[dir].resize(n, INF8);
  }
  for (const auto &e : es) {
    if (e.fst == e.snd) continue;
    adj[0][e.fst].push_back(e.snd);
    adj[1][e.snd].push_back(e.fst);
  }

  for (int dir = 0; dir < 2; dir++) {
    for (int v = 0; v < n; v++) {
      sort(adj[dir][v].begin(), adj[dir][v].end());
      MakeUnique(adj[dir][v]);
    }
  }
}

template<int kNumBitParallelRoots, class upperbound_ls>
GraphDistance<kNumBitParallelRoots, upperbound_ls>::
GraphDistance(int n, const vector<pair<int, int> > &es) {
  Build(n, es, 3);
}

template<int kNumBitParallelRoots, class upperbound_ls>
void GraphDistance<kNumBitParallelRoots, upperbound_ls>::
Build(int n, const vector<pair<int, int> > &es, int param) {
  Clear();
  V = n;
  que.push_back(TwoLayerQueue<int>(n));
  que.push_back(TwoLayerQueue<int>(n));
  InitGraph(n, es);
  ubls.Build(adj, param);
}



template<int kNumBitParallelRoots, class upperbound_ls>
int GraphDistance<kNumBitParallelRoots, upperbound_ls>::
QueryDistanceUB(int s, int t) {
  if (s == t) return 0;
  CHECK_LE(0, s); CHECK_LT(s, V);
  CHECK_LE(0, t); CHECK_LT(t, V);
  return ubls.GetDistanceUpperbound(s, t);
}

template <int kNumBitParallelRoots, class upperbound_ls>
int GraphDistance<kNumBitParallelRoots, upperbound_ls>::
QueryDistance(int s, int t) {
  if (s == t) {   return 0;  }
  query_stats.Reset();
  CHECK_LE(0, s); CHECK_LT(s, V);
  CHECK_LE(0, t); CHECK_LT(t, V);
  uint8_t dist_upper = QueryDistanceUB(s, t);
  if (ubls.SkipVertex(s) || ubls.SkipVertex(t)) {
    return dist_upper;
  }

  uint8_t res = dist_upper, dis[2] = {0, 0};
  for (int dir = 0; dir < 2; dir++){
    int v= dir == 0 ? s : t;
    que[dir].clear();
    que[dir].push(v);
    que[dir].next();
    dist[dir][v] = 0;
  }

  while (!que[0].empty() && !que[1].empty()) {
    int use = 0;
    use = (que[0].size() <= que[1].size()) ? 0 : 1;
    dis[use]++;

    if (dis[0] + dis[1] == dist_upper) {
      res = dis[0] + dis[1];
      goto LOOP_END;
    }

    while (!que[use].empty()) {
      query_stats.visit_vs++;

      int v = que[use].front();
      que[use].pop();

      if (ubls.SkipVertex(v)) continue;

      for (int w : adj[use][v]) {
        query_stats.visit_es++;

        uint8_t &src_d = dist[use][w];
        uint8_t &dst_d = dist[1 - use][w];
        if (src_d != INF8) continue;
        if (dst_d != INF8) {
          res = dist[use][v] + 1 + dst_d;
          query_stats.goal = 1;
          goto LOOP_END;
        } else {
          que[use].push(w);
          dist[use][w] = dist[use][v] + 1;
        }
      }
    }
    que[use].next();
  }
  LOOP_END:

  for (int dir = 0; dir < 2; dir++) {
    for (int v : que[dir]) {
      dist[dir][v] = INF8;
    }
    que[dir].clear();
  }
  return std::min(res, dist_upper);
}

template<int kNumBitParallelRoots, class upperbound_ls>
void GraphDistance<kNumBitParallelRoots, upperbound_ls>::
InsertEdge(int s, int t) {
  if (s == t) return;
  CHECK_LE(0, s); CHECK_LT(s, V);
  CHECK_LE(0, t); CHECK_LT(t, V);
  InsertEdgeIntoGraph(adj[0], s, t);
  InsertEdgeIntoGraph(adj[1], t, s);
  ubls.InsertEdge(s, t);
}

template<int kNumBitParallelRoots, class upperbound_ls>
void GraphDistance<kNumBitParallelRoots, upperbound_ls>::
DeleteEdge(int s, int t) {
  if (s == t) return;
  CHECK_LE(0, s); CHECK_LT(s, V);
  CHECK_LE(0, t); CHECK_LT(t, V);
  DeleteEdgeFromGraph(adj[0], s, t);
  DeleteEdgeFromGraph(adj[1], t, s);
  ubls.DeleteEdge(s, t);
}

#endif /* DISTANCE_WITH_BIT_PARALLEL_H */
