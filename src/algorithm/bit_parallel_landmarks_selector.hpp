#ifndef BIT_PARALLEL_LANDMARKS_SELECTOR_H
#define BIT_PARALLEL_LANDMARKS_SELECTOR_H

#include <vector>
#include <cstdlib>
#include <string>
#include <cassert>
#include <algorithm>
#include <queue>
#include <cmath>
#include <iostream>
#include <numeric>
#include "common.hpp"
#include "jlog.hpp"

using std::vector;
using std::pair;
using std::make_pair;

// Given a mutual graph, return the landmarks and their neighbors.
class LandmarkSelector {
protected:

  // Given adjacency lists consisting of forward and backward edges,
  // build an adjacency list consisting of vertex pairs where both of forward and backward edges exist between them. 
  virtual void ComputeMutualGraph(const vector<vector<int> > &forward_adj, const vector<vector<int> > &backward_adj, vector<vector<int> > &mutual_adj) {
    int V = forward_adj.size();
    vector<bool> has_edge(V, false);
    mutual_adj.resize(V);

    size_t mutual_es = 0;
    size_t es = 0;
    for (int v = 0; v < V; v++) {
      for (auto w : backward_adj[v]) has_edge[w] = true;
      for (auto w : forward_adj[v]) if (has_edge[w]) mutual_adj[v].push_back(w);
      for (auto w : backward_adj[v]) has_edge[w] = false;
      mutual_es += mutual_adj[v].size();
      es += forward_adj[v].size();
    }
    JLOG_PUT("setting.num_mutual_es", mutual_es / 2);
    JLOG_PUT("setting.frac_mutual_es", (double) mutual_es / es);
  }


public:
  virtual vector<pair<int, vector<int> > > SelectLandmarks(const vector<vector<int> > &forward_adj,
                                                           const vector<vector<int> > &backward_adj, int k) = 0;
};


template<int NumNeighbors>
class DegreeLS : public LandmarkSelector {
protected:
  vector<int> i_degree;
  vector<int> o_degree;
  vector<int> m_degree;
  vector<int> used;

  void InitializeDegree(const vector<vector<int> > &forward_adj, const vector<vector<int> > &backward_adj, const vector<vector<int> > &mutual_adj) {
    assert(i_degree.empty());
    for (size_t i = 0; i < forward_adj.size(); i++) {
      i_degree.push_back(backward_adj[i].size());
      o_degree.push_back(forward_adj[i].size());
      m_degree.push_back(mutual_adj[i].size());
    }
    used.resize(forward_adj.size(), false);
  }

  virtual void UseVertex(int v, const vector<vector<int> > &forward_adj, const vector<vector<int> > &backward_adj,
                         const vector<vector<int> > &mutual_adj) {
    assert(!used[v]);
    for (int w : forward_adj[v]) {
      i_degree[w]--;
      CHECK_GE(i_degree[w], 0);
    }

    for (int w : backward_adj[v]) {
      o_degree[w]--;
      CHECK_GE(o_degree[w], 0);
    }

    for (int w : mutual_adj[v]) {
      m_degree[w]--;
      CHECK_GE(m_degree[w], 0);
    }
    used[v] = true;
  }

  virtual double GetScore(int v) const {
    return this->i_degree.at(v) + this->o_degree.at(v);
  }

public:
  virtual vector<pair<int, vector<int> > > SelectLandmarks(const vector<vector<int> > &forward_adj,
                                                           const vector<vector<int> > &backward_adj, int k) {
    int V = forward_adj.size();
    vector<vector<int> > mutual_adj;
    this->ComputeMutualGraph(forward_adj, backward_adj, mutual_adj);

    std::for_each(mutual_adj.begin(), mutual_adj.end(), [&](vector<int> &vs) {
        std::sort(vs.begin(), vs.end());
        vs.erase(std::unique(vs.begin(), vs.end()), vs.end());
        std::sort(vs.begin(), vs.end(), [&mutual_adj](int a, int b) { return mutual_adj[a].size() > mutual_adj[b].size(); });
      });

    this->InitializeDegree(forward_adj, backward_adj, mutual_adj);


    // Use lazy-greedy algorithm to speed up landmark selection. 
    typedef pair<double, int> P;
    std::priority_queue<P> que;
    
    auto ComputeScoreSum = [&](int v) -> double {
      assert(!used[v]);
      double sum = this->GetScore(v);
      for (size_t i = 0, ns = 0; i < mutual_adj[v].size() && ns < NumNeighbors; i++) {
        int w = mutual_adj[v][i];
        if (!used[w]) {
          sum += this->GetScore(mutual_adj[v][i]);
          ns++;
        }
      }
      return sum;
    };

    for (int v = 0; v < V; v++) {
      que.push(P(ComputeScoreSum(v), v));
    }
    vector<pair<int, vector<int> > > res;
    double total = 0;

    while (int(res.size()) < k && !que.empty()) {
      P p = que.top();
      que.pop();
      if (used[p.second]) continue;
      
      int r = p.second;
      if (ComputeScoreSum(r) == p.first) { // add a new landmark.
        vector<int> neighbors;

        total += p.first;
        UseVertex(r, forward_adj, backward_adj, mutual_adj);

        for (size_t i = 0; i < mutual_adj[r].size() && neighbors.size() < NumNeighbors; i++) {
          int v = mutual_adj[r][i];
          if (!used[v]) {
            neighbors.push_back(v);
            UseVertex(v, forward_adj, backward_adj, mutual_adj);
          }
        }
        res.push_back(make_pair(r, neighbors));
      } else {
        mutual_adj[r].erase(std::remove_if(mutual_adj[r].begin(), mutual_adj[r].end(), [&](int v) { return used[v]; }), mutual_adj[r].end());
        std::sort(mutual_adj[r].begin(), mutual_adj[r].end(), [&](int a, int b) { return mutual_adj[a].size() > mutual_adj[b].size(); });
        que.push(make_pair(ComputeScoreSum(r), p.second));
      }
    }
    res.resize(k);

    JLOG_PUT("landmark_selection.total_score", total);
    JLOG_PUT("landmark_selection.size", res.size());
    while (int(res.size()) < k) {
      res.emplace_back(std::make_pair(0, vector<int>()));
    }
    return res;
  }
};

#endif /* BIT_PARALLEL_LANDMARKS_SELECTOR_H */
