#ifndef DISTANCE_WITHOUT_BIT_PARALLEL_H
#define DISTANCE_WITHOUT_BIT_PARALLEL_H

#include "union_find.hpp"
#include "two_layer_queue.hpp"
#include "algorithm/distance_with_bit_parallel.hpp"
#include "common.hpp"
#include "queue.hpp"
#include <cstdint>
#include <algorithm>
#include <string>

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;


template <int num_landmarks = 0> class Distance {
  typedef std::pair<int, int> PI;
  int V;
  UnionFind uf;
  std::vector<std::vector<int> > adj[2];
  std::vector<uint8_t> dist[2];
  std::vector<int> stop;
  std::vector<int> ls;
  std::vector<std::vector<uint8_t> > tree_dist[2];
  vector<TwoLayerQueue<int> > que;
  vector<Queue<int> > alter_que;
  QueryStats query_stats;

public:
    
  Distance (int n, const std::vector<PI> &es) :
    V(n), uf(n), stop(n), que(2, TwoLayerQueue<int>(V)), alter_que(2, Queue<int>(V)){
    for (int dir = 0; dir < 2; dir++) {
      adj[dir].resize(n);
      dist[dir].resize(n, INF8);
    }
      
    for (const auto &e : es){
      CHECK_LT(e.fst, n);
      CHECK_LT(e.snd, n);
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
    
    std::vector<PI> ds;
    for (int v = 0; v < n; v++) {
      ds.emplace_back(-adj[0][v].size(), v);
    }
    sort(ds.begin(), ds.end());
      
    for (int i = 0; i < num_landmarks; i++) {
      stop[ds[i].second] = true;
      ls.push_back(ds[i].second);
      for (int dir = 0; dir < 2; dir++) {
        tree_dist[dir].push_back(std::vector<uint8_t>(n, INF8));
        tree_dist[dir][i][ls.back()] = 0;
        std::queue<int> que; 
        que.push(ls.back());
        while (!que.empty()){
          int v = que.front(); que.pop();
          for (int w : adj[dir][v]){
            if (tree_dist[dir][i][w] > tree_dist[dir][i][v] + 1){
              tree_dist[dir][i][w] = tree_dist[dir][i][v] + 1;
              que.push(w);
            }
          }
        }
      }
    }
    
    for (const auto &e : es){
      uf.unite(e.fst, e.snd);
    }

    for (int dir = 0; dir < 2; dir++) {
      for (int v = 0; v < n; v++) {
        if (stop[v]){
          adj[dir][v].clear();
        } else {
          size_t m = 0;
          for (size_t j = 0; j < adj[dir][v].size(); j++) {
            if (!stop[adj[dir][v][j]])  adj[dir][v][m++] = adj[dir][v][j];
          }
          adj[dir][v].resize(m);
        }
      }
    }
  }

  int QueryDistanceBFS(int s, int t){
    query_stats.Reset();
    
    vector<int> visit;
    std::queue<int> que;
    que.push(s);
    dist[0][s] = 0;
    visit.push_back(s);
    
    while (!que.empty()){
      int v = que.front(); que.pop();
      int d = dist[0][v];
      query_stats.visit_vs++;
      
      if (v == t){
        for (int w : visit){
          dist[0][w] = INF8;
        }
        return d;
      }

      for (int w : adj[0][v]){
        query_stats.visit_es++;
        if (dist[0][w] == INF8){
          dist[0][w] = d + 1;
          que.push(w);
          visit.push_back(w);
        }
      }
    }

    for (int w : visit){
      dist[0][w] = INF8;
    }
    return INF8;
  }
  
  int QueryDistance(int s, int t){
    query_stats.Reset();
    assert(0 <= s && s < V);
    assert(0 <= t && t < V);
    if (s == t) return 0;
    if (!uf.same(s, t)) return INF8;
    
    uint8_t dist_limit = INF8;
    for (size_t i = 0; i < ls.size(); i++) {
      dist_limit = std::min<uint8_t>(tree_dist[1][i][s]+tree_dist[0][i][t], dist_limit);
    }
    
    uint8_t res = INF8, dis[2] = {0, 0};
    que[0].clear(); que[1].clear();
    que[0].push(s); dist[0][s] = 0; 
    que[1].push(t); dist[1][t] = 0; 
    que[0].next(); que[1].next();
    
    
    while (!que[0].empty() && !que[1].empty()){
      int use = (que[0].size() <= que[1].size()) ? 0 : 1;
      dis[use]++;
      
      if (dis[0] + dis[1] == dist_limit) {
        res = dis[0] + dis[1];
        goto LOOP_END;
      }
      
      while (!que[use].empty()){
        query_stats.visit_vs++;
        int v = que[use].front(); que[use].pop();
        
        for (int w : adj[use][v]){
          query_stats.visit_es++;
          uint8_t &src_d = dist[    use][w];
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
      if (res != INF8) goto LOOP_END;
    }
  LOOP_END:
    for (int dir = 0; dir < 2; dir++) {
      for (int v : que[dir]) dist[dir][v] =  INF8;
      que[dir].clear();
    }
    return res == INF8 ? dist_limit : std::min(res, dist_limit);
  }

  int QueryDistanceAlternate(int s, int t) {
    query_stats.Reset();
    CHECK_LE(0, s); CHECK_LT(s, V);
    CHECK_LE(0, t); CHECK_LT(t, V);
    
    if (s == t) {
      return 0;
    } else if (!uf.same(s, t)){
      return INF8;
    }

    uint8_t res = INF8, dis[2] = {0, 0};
    for (size_t i = 0; i < ls.size(); i++) {
      res = std::min<uint8_t>(tree_dist[1][i][s]+tree_dist[0][i][t], res);
    }

    for (int dir = 0; dir < 2; dir++) {
      int v = dir == 0 ? s : t;
      alter_que[dir].clear();
      alter_que[dir].push(v);
      dist[dir][v] = 0;
    }

    int dir = 0;
    while (!alter_que[0].empty() && !alter_que[1].empty()) {
      query_stats.visit_vs++;
      int v = alter_que[dir].front(); alter_que[dir].pop();
      dis[dir] = std::max(dist[dir][v], dis[dir]);

      if (dis[0] + dis[1] + 1 >= res) {
        goto LOOP_END;
      }

      if (stop[v]) continue;

      for (int w : adj[dir][v]) {
        query_stats.visit_es++;
        uint8_t &src_d = dist[    dir][w];
        uint8_t &dst_d = dist[1 - dir][w];
        if (src_d != INF8) continue;
        if (dst_d != INF8) {
          res = std::min<uint8_t>(res, dist[dir][v] + 1 + dst_d);
        } else {
          alter_que[dir].push(w);
          dist[dir][w] = dist[dir][v] + 1;
        }
      }
      dir = 1 - dir;
    }
    LOOP_END:

    for (int dir = 0; dir < 2; dir++) {
      for (int v : alter_que[dir]) {
        dist[dir][v] = INF8;
      }
      alter_que[dir].clear();
    }
    return res;
  }

  inline QueryStats GetQueryStats() const { return query_stats; }
};

#endif /* DISTANCE_WITHOUT_BIT_PARALLEL_H */
