#ifndef BIT_PARALLEL_LANDMARKS_H
#define BIT_PARALLEL_LANDMARKS_H

#include "common.hpp"
#include "distance_with_bit_parallel.hpp"
#include "bit_parallel_landmarks_selector.hpp"
#include <vector>
#include <queue>
#include <array>
#include <cstdlib>
#include <tuple>
#include <bitset>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include<memory>
#include "jlog.hpp"
#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>

using std::vector;
using std::array;
using std::pair;
using std::cout;
using std::endl;
using std::cerr;
using std::queue;
DEFINE_bool(output_iparam_log, false, "");
DEFINE_double(remove_threshold, 0.05, "");

typedef pair<uint64_t, uint64_t> neighbor_set_t;
template <typename T, typename E> using hash_map = google::dense_hash_map<T,E>;
template <typename T            > using hash_set = google::dense_hash_set<T>;

template <typename T, typename E> hash_map<T, E> CreateEmptyHashMap(){
  hash_map<T, E> h;
  h.set_empty_key(-1);
  h.set_deleted_key(-2);
  return h;
};

template <typename T> hash_set<T> CreateEmptyHashSet(){
  hash_set<T> h;
  h.set_empty_key(-1);
  h.set_deleted_key(-2);
  return h;
};


template<typename S, typename T>
std::ostream &operator<<(std::ostream &out, const std::pair<S, T> &p) {
  return out << "(" << p.first << ", " << p.second<< ")";
}

std::ostream &operator<<(std::ostream &out, const std::vector<uint8_t> &vec){
  out << "[";
  for (size_t i = 0; i < vec.size(); i++){
    out << (i > 0 ? ", ":  "") << int(vec[i]);
  }
  out << "]";
  return out;
}

neighbor_set_t operator&(const neighbor_set_t &m, uint64_t filter){
  return neighbor_set_t(m.fst & filter, m.snd & filter);
};

class BitParallelSPT {
private:
  const vector<vector<int> > *forward_adj, *backward_adj;
  vector<int> *covered_vs, *ns_offset;

  int V, root;
  vector<uint8_t> D;
  vector<neighbor_set_t> M;
  uint64_t neighbor_set_mask_filter;
  vector<int> id2neighbor;
  hash_map<int, int> neighbor2id;

public:
  //  variables for logging
  void RemoveLowDegreeNeighbors(double threshold){
    auto GetDegree = [&](int v) { return forward_adj->at(v).size() + backward_adj->at(v).size(); };
    
    vector<int> vs = {root};
    for (int v : id2neighbor) {
      if (InNeighborSet(v)) vs.push_back(v);
    }
    sort(vs.begin(), vs.end(), [&](int a, int b) {
        return GetDegree(a) > GetDegree(b);
      });
    
    JLOG_ADD("remove_low_degrees.maximum", GetDegree(*vs.begin()));
    JLOG_ADD("remove_low_degrees.minimum", GetDegree(*vs.rbegin()));
    
    for (int v : vs) {
      if (v != root && GetDegree(v) < GetDegree(*vs.begin()) * threshold){
        neighbor_set_mask_filter &= ~(uint64_t(1) << neighbor2id[v]);
        neighbor2id.erase(v);
      }
    }
    JLOG_ADD("remove_low_degrees.num_nvs", neighbor2id.size());
    assert(neighbor2id.size() == __builtin_popcountll(neighbor_set_mask_filter));
  }

private:
  inline bool HasNeighborSet(int v) const {
    return ns_offset->at(v) != -1;
  }

  inline neighbor_set_t GetStoredNeighborSet(int v) const {
    return M[ns_offset->at(v)] & neighbor_set_mask_filter;
  }

  inline uint64_t GetNeighborMask(int v) const {
    CHECK(InNeighborSet(v));
    return (uint64_t(1) << neighbor2id.find(v)->second) & neighbor_set_mask_filter;
  }

  inline uint64_t GetNeighborSetFilter() const {
    return neighbor_set_mask_filter;
  }

  inline void NeighborSetUnion(neighbor_set_t &mask_dst, const neighbor_set_t &mask_src) {
    mask_dst.first |= mask_src.first;
    mask_dst.second |= mask_src.second;
    mask_dst.second &= ~(mask_dst.first);
  }

  inline void ApplyNeighborSetFilter(neighbor_set_t &mask) {
    mask.first &= this->GetNeighborSetFilter();
    mask.second &= this->GetNeighborSetFilter();
  }

  inline uint64_t RecoverR(int v, hash_map<int, neighbor_set_t> &X, hash_map<int, uint64_t> &Y) const {
    if (v == root) return 0;
    else if (HasNeighborSet(v)) return GetStoredNeighborSet(v).fst;
    else if (InNeighborSet(v)) return GetNeighborMask(v);
    else if (X.count(v)) return X[v].fst;
    else if (Y.count(v)) return Y[v];

    uint64_t mask = 0;
    for (int w : backward_adj->at(v)) {
      if (D[v] == D[w] + 1) {
        mask |= RecoverR(w, X, Y);
      }
    }
    return Y[v] = mask;
  }

  inline neighbor_set_t Recover_(int v, hash_map<int, neighbor_set_t> &X, hash_map<int, uint64_t> &Y) const {
    if (v == root) return neighbor_set_t(uint64_t(0), uint64_t(0));
    else if (HasNeighborSet(v)) return GetStoredNeighborSet(v);
    else if (X.count(v)) return X[v];

    uint64_t mask0 = 0, mask1 = 0;
    if (InNeighborSet(v)) mask0 = GetNeighborMask(v);

    for (int w : backward_adj->at(v)) {
      CHECK(v != w)
        if (D[v] == D[w] + 1) {
          auto mask_ = Recover_(w, X, Y);
          mask0 |= mask_.fst;
          mask1 |= mask_.snd;
        }
    }

    for (int w : backward_adj->at(v)) {
      CHECK(v != w);
      if (D[v] == D[w]) {
        mask1 |= RecoverR(w, X, Y);
      }
    }
    mask1 &= ~mask0;
    return X[v] = std::make_pair(mask0, mask1);
  }

  inline neighbor_set_t Recover(int v) const {
    hash_map<int, neighbor_set_t> X = CreateEmptyHashMap<int, neighbor_set_t>();
    hash_map<int, uint64_t> Y = CreateEmptyHashMap<int, uint64_t>();
    return Recover_(v, X, Y);
  }
  
  bool UpdateI(int u, int v) {
    hash_map<int, neighbor_set_t> X = CreateEmptyHashMap<int, neighbor_set_t>();
    hash_map<int, uint64_t> Y = CreateEmptyHashMap<int, uint64_t>();

    bool has_update = false;
    neighbor_set_t mv, mu;

    if (D[u] + 1 < D[v]) {
      has_update = true;
      D[v] = D[u] + 1;
      // Calcluate v's neighbor sets information from v's neighbors. 
      mv = neighbor_set_t(uint64_t(0), uint64_t(0));
      for (int w : backward_adj->at(v)) {
        assert(v != w);
        if (D[w] + 1 == D[v]) {
          NeighborSetUnion(mv, Recover(w));
        } else if (D[w] == D[v]) {
          mv.snd |= Recover_(w, X, Y).first;
        }
      }
      X.clear();
      Y.clear();
      mv.snd &= ~mv.fst;
    } else if (D[u] + 1 == D[v]) {
      // We need to update v's neighbor sets if u's neighbor sets change it. 
      mu = Recover(u);
      mv = Recover(v);
      has_update = (mu.fst & ~mv.fst) > 0 || (mu.snd & ~(mv.fst | mv.snd)) > 0 || !HasNeighborSet(v);
      NeighborSetUnion(mv, mu);
    } else if (D[u] == D[v]) {
      // We need to update v's neighbor sets if u's neighbor sets change it. 
      mu = Recover(u);
      mv = Recover(v);
      has_update = (mu.fst & ~(mv.fst | mv.snd)) > 0 || !HasNeighborSet(v);
      mv.snd |= mu.fst;
      mv.snd &= ~mv.fst;
    }

    if (has_update && HasNeighborSet(v)) {
      M[ns_offset->at(v)] = mv;
    }
    return has_update;
  }

public:
  void Build(int root, const vector<int> root_neighbors,
             const vector<vector<int> > *forward_adj,
             const vector<vector<int> > *backward_adj,
             vector<int> *covered_vs, vector<int> *ns_offset, int dset_size) {
    
    this->forward_adj = forward_adj;
    this->backward_adj = backward_adj;
    this->root = root;
    this->V = forward_adj->size();
    this->covered_vs = covered_vs;
    this->ns_offset = ns_offset;
    
    int I = 0;
    id2neighbor.clear();
    neighbor2id.clear();
    neighbor2id = CreateEmptyHashMap<int, int>();
    for (int nv : root_neighbors) {
      id2neighbor.push_back(nv);
      neighbor2id[nv] = I++;
    }

    D.resize(V);
    M.resize(dset_size);
    fill(D.begin(), D.end(), INF8);
    fill(M.begin(), M.end(), std::make_pair(0, 0));
    this->neighbor_set_mask_filter = ~(unsigned long long) 0;
    // Set 0 for bits that do not correspond to any vertex
    for (int i = root_neighbors.size(); i < 64; i++){
      neighbor_set_mask_filter &= ~(uint64_t(1) << i);
    }
    
    vector<neighbor_set_t> tmp_s(V);
    BitParallelBFS(D, tmp_s);
    for (int v = 0; v < V; v++) {
      if (HasNeighborSet(v)) {
        tmp_s[v].snd &= ~tmp_s[v].fst;
        M[ns_offset->at(v)] = tmp_s[v];
      }
    }
  }

  inline uint8_t GetDistance(int v) const {
    return D.at(v);
  }

  inline neighbor_set_t GetNeighborSet(int v) const {
    return Recover(v) & neighbor_set_mask_filter;
  };

  inline bool InNeighborSet(int v) const {
    return neighbor2id.count(v);
  }

  inline size_t IndexSize() const {
    size_t total = 0;
    total += sizeof(forward_adj) + sizeof(backward_adj);
    total += sizeof(covered_vs) + sizeof(ns_offset) + sizeof(V) + sizeof(root);
    total += VectorSize(D);
    total += VectorSize(M);
    total += sizeof(neighbor_set_mask_filter);
    total += VectorSize(id2neighbor);
    total += sizeof(std::pair<int, int>) * neighbor2id.bucket_count();
    return total;
  }

  inline uint8_t GetLatestDistance(const int v, const hash_map<int, pair<uint8_t, neighbor_set_t> > &next_info) const {
    const auto iter = next_info.find(v);
    return iter != next_info.end() ? iter->second.first : this->D[v];
  }

  inline neighbor_set_t GetLatestNeighborSet(const int v, hash_map<int, neighbor_set_t> &X, hash_map<int, uint64_t> &Y,
                                             const hash_map<int, pair<uint8_t, neighbor_set_t> > &next_info){
    const auto iter = next_info.find(v);
    return iter != next_info.end() ? iter->second.second : Recover_(v, X, Y);
  }

  bool PropagateI(const int s, const int t, hash_map<int, neighbor_set_t> &X, hash_map<int, uint64_t> &Y,
                  hash_map<int, pair<uint8_t, neighbor_set_t> > &next_info, const int u, const int v)
  {
    auto new_u = next_info.count(u) ? next_info[u] : make_pair(this->D[u], this->Recover(u));
    uint8_t cur_vd = next_info.count(v) ? next_info[v].first : this->D[v];
    uint8_t new_ud = new_u.first;

    if (new_ud + 1 < cur_vd) {
      // obviously, the information of v should be updated. 
      uint8_t new_vd = new_ud + 1;
      neighbor_set_t new_vm(0, 0);
      
      // compute new neighbor sets with the updated distance. 
      for (int c : backward_adj->at(v)) {
        assert(v != c);
        auto new_c = next_info.count(c) ? next_info[c] : make_pair(this->D[c], this->Recover(c));
        if (new_c.first + 1 == new_vd) {
          NeighborSetUnion(new_vm, new_c.second);
        }
      }
      next_info[v] = make_pair(new_vd, new_vm);
      return true;
    } else if (new_ud  > cur_vd){
      return false;
    } else {
      // compute the current information of a vertex v. 
      pair<uint8_t, neighbor_set_t> cur_v;

      {
        if (next_info.count(v)){
          cur_v = next_info[v];
        } else if (this->HasNeighborSet(v)){
          cur_v.first  = this->D[v];
          cur_v.second = this->M[ns_offset->at(v)];
        } else {
          cur_v.first = this->D[v];
          if (InNeighborSet(v)) cur_v.second.first |= GetNeighborMask(v);

          for (int w : backward_adj->at(v)) {
            CHECK(v != w)
              if (w == s && v == t) continue;
            if (D[v] == D[w] + 1) {
              NeighborSetUnion(cur_v.second, Recover_(w, X, Y));
            } else if (D[v] == D[w] + 1){
              cur_v.second.second |= RecoverR(w, X, Y);
              cur_v.second.second &= ~cur_v.second.first;
            }
          }
        }
      }
      auto cur_vm = cur_v.second, new_vm = cur_v.second, new_um = new_u.second;

      bool has_update = false;
      if (new_u.first + 1 == cur_v.first) {
        NeighborSetUnion(new_vm, new_um);
        has_update = new_vm != cur_vm;
      } else if (new_u.first == cur_v.first) {
        new_vm.second |= new_um.first;
        new_vm.second &= ~new_vm.first;
        has_update = new_vm != cur_vm;
      }

      if (has_update){
        next_info[v] = make_pair(cur_vd, new_vm);
      }
      return has_update;
    }
    return false;
  }
  
  void InsertEdge(int s, int t) {
    hash_map<int, pair<uint8_t, neighbor_set_t> > next_info = CreateEmptyHashMap<int, pair<uint8_t, neighbor_set_t> > ();
    hash_map<int, neighbor_set_t> X = CreateEmptyHashMap<int, neighbor_set_t>();
    hash_map<int, uint64_t> Y = CreateEmptyHashMap<int, uint64_t>();
    hash_set<int> curr = CreateEmptyHashSet<int>();
    
    if (PropagateI(s, t, X, Y, next_info, s, t)) {
      curr.insert(t);
      CHECK(next_info.count(t));
    }

    while (!curr.empty()) {
      hash_set<int> next = CreateEmptyHashSet<int>();
      hash_set<int> temp = curr;

      for (int u : curr) {
        for (int v : forward_adj->at(u)) {
          CHECK(u != v);
          if (!PropagateI(s, t, X, Y, next_info, u, v)) continue;
          
          uint8_t du = GetLatestDistance(u, next_info), dv = GetLatestDistance(v, next_info);
          if (dv == du + 1) {
            next.insert(v);
          } else if (dv == du) {
            temp.insert(v);
          }
        }
        for (int v : backward_adj->at(u)) {
          CHECK(u != v);
          PropagateI(s, t, X, Y, next_info, v, u);
        }
      }
      
      curr = temp;
      for (int u : curr) {
        for (int v : forward_adj->at(u)) {
          CHECK(u != v);
          if (PropagateI(s, t, X, Y, next_info, u, v)) {
            if (GetLatestDistance(u, next_info) + 1 == GetLatestDistance(v, next_info)) {
              next.insert(v);
            }
          }
        }
      }
      curr = next;
    }

    for (const std::pair<int, pair<uint8_t, neighbor_set_t> > &p : next_info) {
      int v = p.first;
      this->D[v] = p.second.first;
      if (HasNeighborSet(v)) {
        this->M[ns_offset->at(v)] = p.second.second;
      }
    }
  }

  uint8_t FindMin(int v) const {
    uint8_t res = INF8;
    for (int u : backward_adj->at(v)) {
      res = std::min(res, uint8_t(D[u] + 1));
    }
    return res;
  }

  bool PropagateD(const int s, const int t, hash_map<int, neighbor_set_t> &X, hash_map<int, uint64_t> &Y,
                  hash_map<int, pair<uint8_t, neighbor_set_t> > &next_info, const int u, const int v) {
    bool has_update = false;
    if (D[u] > D[v]) {
      has_update = false;
    } else if (D[v] == D[u] || D[v] == D[u] + 1) {
        
      auto FindLatestMin = [&](int v) -> uint8_t {
        uint8_t res = INF8;
        for (int u : backward_adj->at(v)) {
          uint8_t du = GetLatestDistance(u, next_info);
          if (du == this->D[u] && du + 1 < res) {
            res = du + 1;
          }
        }
        return res;
      };
      
      uint8_t new_vd = FindLatestMin(v);
      neighbor_set_t  new_vm = neighbor_set_t(0, 0);
      if (InNeighborSet(v)) {
        new_vm.first |= GetNeighborMask(v);
      }

      if (new_vd > D[v]) {
        // If the distance from the root to v changes, we set the distance to infinity. 
        next_info[v] = make_pair(INF8, new_vm);
        has_update = true;
      } else {
        // Compute current mask new_vm
        // Some of v's neighbors' neighbor sets are not updated yet.
        // However, when these neighbor sets are updated, we try to propagate these changes to v again. 

        for (int w : backward_adj->at(v)) {
          if (GetLatestDistance(w, next_info) != D[w]) continue;

          if (D[w] + 1 == D[v]) {
            NeighborSetUnion(new_vm, GetLatestNeighborSet(w, X, Y, next_info));
          } else if (D[w] == D[v]) {
            new_vm.second |= GetLatestNeighborSet(w, X, Y, next_info).first & GetNeighborSetFilter();
            new_vm.second &= ~new_vm.first;
          }
        }
        // Compute previous mask old_vm
        neighbor_set_t old_vm = Recover_(v, X, Y);
        if (v == t){
          if (D[s] + 1 == D[t]) {
            NeighborSetUnion(old_vm, Recover_(s, X, Y));
          } else {
            old_vm.second |= RecoverR(s, X, Y);
            old_vm.second &= ~old_vm.first;
          }
        }
        
        has_update = new_vm != old_vm;
        if (has_update){
          next_info[v] = make_pair(D[v], new_vm);
        }
      }
    } else {
      CHECK(false);
    }
    return has_update;
  }
  
  hash_set<int> UpdatedVertexSet(int s, int t, hash_map<int, pair<uint8_t, neighbor_set_t> > &next_info) {
    hash_set<int> C = CreateEmptyHashSet<int>();
    hash_map<int, neighbor_set_t> X = CreateEmptyHashMap<int, neighbor_set_t>();
    hash_map<int, uint64_t> Y = CreateEmptyHashMap<int, uint64_t>();
    
    uint8_t d = D[t];
    hash_set<int> curr = CreateEmptyHashSet<int>();
    if (s != t && PropagateD(s, t, X, Y, next_info, s, t)) {
      curr.insert(t);
    }

    while (!curr.empty()) {
      hash_set<int> next = CreateEmptyHashSet<int>();
      hash_set<int> tmp_curr = curr;
      for (int u : tmp_curr) {
        C.insert(u);
      }
      
      // We use a queue to implement this propagation procedure in the paper, 
      // and the following procedure is doing the same thing. 
      for (int i = 0; i < 2; i++) {
        for (int u : tmp_curr) {
          for (int v : forward_adj->at(u)) {
            CHECK(u != v);
            if (PropagateD(s, t, X, Y, next_info, u, v)){
              CHECK_EQ(D[u], d);
              if (D[v] == d + 1) {
                next.insert(v);
              } else if (D[v] == d && i == 0) {
                curr.insert(v);
              }
            }
          }
        }
        tmp_curr = curr;
      }
      swap(curr, next);
      d++;
    }
    
    for (const pair<int, pair<uint8_t, neighbor_set_t> > &p : next_info){
      this->D[p.first] = p.second.first;
      if (HasNeighborSet(p.first)) {
        this->M[ns_offset->at(p.first)] = p.second.second;
      }
    }
    return C;
  }

  void DeleteEdge(int s, int t) {
    
    // First of all, we update a set of neighbor vertices. We need to delete a vertex s or t from a neighbor set
    // if one of them is the root and another one is in its neighbor set. 
    if (D[s] == 0 && InNeighborSet(t)){
      neighbor_set_mask_filter &= ~GetNeighborMask(t);
      covered_vs->at(t) = false;
      neighbor2id.erase(t);
    }
    if (D[t] == 0 && InNeighborSet(s)){
      neighbor_set_mask_filter &= ~GetNeighborMask(s);
      covered_vs->at(s) = false;
      neighbor2id.erase(s);
    }

    
    if (D[t] > 0) {
      // Compute a vertex set C whose distances from the root or neighbor sets changed. 
      hash_map<int, pair<uint8_t, neighbor_set_t> > next_info = CreateEmptyHashMap<int, pair<uint8_t, neighbor_set_t> > ();
      hash_set<int> C = UpdatedVertexSet(s, t, next_info);
      
      // 更新中の頂点集合を保持するための集合Uを容易
      static vector<hash_set<int> > U;
      if (U.empty()) {
        U.resize(INF8 + 1, CreateEmptyHashSet<int>());
      }
      
      uint8_t max_d = 0;
      for (int v : C) {
        uint8_t d = FindMin(v);
        CHECK(next_info.count(v));
        
        if (d < INF8) {
          // If a vertex v has neighbor sets in memory, we collect the imformation of neighbor sets of un-updated neighbors of v in advance.
          // The information of updated neighbors of v is propagated later. 
          D[v] = d;
          if (HasNeighborSet(v)){
            if (InNeighborSet(v)) {
              M[ns_offset->at(v)].fst = GetNeighborMask(v);
            }
            for (int u : backward_adj->at(v)){
              auto mask = Recover(u);
              if (D[u] + 1 == d){
                NeighborSetUnion(M[ns_offset->at(v)], mask);
              }
              if (D[u] == d){
                M[ns_offset->at(v)].snd |= mask.fst;
                M[ns_offset->at(v)].snd &= ~M[ns_offset->at(v)].fst;
              }
            }
          }
          assert(d < INF8);
          max_d = std::max(max_d, d);
          U[d].insert(v);
        }
      }

      for (uint8_t d = 0; d <= max_d; d++) {
        for (int i = 0; i < 2; i++) {
          hash_set<int> tmp_u = CreateEmptyHashSet<int>();
          
          for (int u : U[d]) {
            for (int v : forward_adj->at(u)) {
              CHECK(u != v);
              // We only have to update the information of vertices in C. 
              if (next_info.count(v) == 0 || !UpdateI(u, v)) {
                continue;
              } else if (D[v] == D[u] + 1) {
                max_d = std::max<uint8_t>(max_d, d + 1);
                if (d + 1 >= int(U.size())) {
                  U.resize(d + 2);
                }
                U[d + 1].insert(v);
              } else if (D[v] == D[u] && i == 1) {
                tmp_u.insert(v);
              }
            }
          }
          U[d] = tmp_u;
        }
        U[d].clear();
      }
      for (int d = 0; d < (int)U.size(); d++) {
        CHECK(U[d].empty());
      }
    }
  }

  void BitParallelBFS(vector<uint8_t> &tmp_d, vector<neighbor_set_t > &tmp_s) const {
    if (root == -1) return;

    std::queue<int> que;
    que.push(root);
    tmp_d[root] = 0;

    for (const auto &p : neighbor2id) {
      int v = p.fst;
      tmp_s[v].fst = 1ULL << p.snd;
    }

    while (!que.empty()) {
      int v = que.front();
      que.pop();
      uint8_t d = tmp_d[v];

      for (const int &w : backward_adj->at(v)) {
        if (tmp_d[v] == tmp_d[w] + 1) {
          tmp_s[v].fst |= tmp_s[w].fst;
          tmp_s[v].snd |= tmp_s[w].snd;
        }
      }
      
      for (const int &w : forward_adj->at(v)) {
        if (d > tmp_d[w]) { ;
        } else if (d == tmp_d[w]) {
          tmp_s[w].snd |= tmp_s[v].fst;
        } else if (tmp_d[w] == INF8) {
          que.push(w);
          tmp_d[w] = d + 1;
        }
      }
    }

    for (int v = 0; v < V; v++) {
      tmp_s[v].fst &= this->GetNeighborSetFilter();
      tmp_s[v].snd &= this->GetNeighborSetFilter();
    }
  }

  void Clear() {
    D.clear();
    M.clear();
    id2neighbor.clear();
    neighbor2id.clear();
  }
};

template<int kNumBitParallelRoots, class selector = DegreeLS<64>>
  class BitParallelLandmarks {
    const vector<vector<int> > *forward_adj;
    int V, I;
    int iset_param;
    array<BitParallelSPT, kNumBitParallelRoots> bpspts[2];
    vector<int> covered_vs;
    vector<int> ns_offset;
    array<int, kNumBitParallelRoots> roots;
    array<vector<int>, kNumBitParallelRoots> root_neighbors;

  public:
  
    size_t DegreeSum() const {
      size_t res = 0;
      for (int i = 0; i < kNumBitParallelRoots; i++) {
        int r = roots[i];
        vector<int> nrs = root_neighbors[i];
        res += forward_adj[0][r].size() + forward_adj[1][r].size();
        size_t max_deg = forward_adj[0][r].size() + forward_adj[1][r].size();
        size_t min_deg = forward_adj[0][r].size() + forward_adj[1][r].size();
        for (int nr : nrs) {
          if (bpspts[0][i].InNeighborSet(nr)) {
            res += forward_adj[0][nr].size() + forward_adj[1][nr].size();
            max_deg = std::max(max_deg, forward_adj[0][nr].size() + forward_adj[1][nr].size());
            min_deg = std::min(min_deg, forward_adj[0][nr].size() + forward_adj[1][nr].size());
          }
        }
        JLOG_ADD("result.max_degree", max_deg);
        JLOG_ADD("result.min_degree", min_deg);
      }
      return res;
    }

  private:

    void SelectBitParallelLandmarks() {
      std::unique_ptr<selector> ptr(new selector());
      auto landmarks = ptr->SelectLandmarks(forward_adj[0], forward_adj[1], kNumBitParallelRoots);
      CHECK_EQ(landmarks.size(), kNumBitParallelRoots)

        covered_vs.resize(V);
      fill(covered_vs.begin(), covered_vs.end(), false);
    
      for (int k = 0; k < kNumBitParallelRoots; k++) {
        roots[k] = landmarks[k].first;
        root_neighbors[k] = vector<int>(landmarks[k].second.begin(), landmarks[k].second.end());
        covered_vs[roots[k]] = true;
        for (int v : root_neighbors[k]){
          covered_vs[v] = true;
        }
      }
    }
  
    int SelectVerticesWithNS(int iset_param) {
      this->iset_param = iset_param;
      ns_offset.resize(V);
      fill(ns_offset.begin(), ns_offset.end(), -1);
    
      vector<bool> remove(V, false);
    
      vector<int> vs;
      {
        int n = forward_adj[0].size();
        vector<int> iset_in_degree_count(n);
        vector<int> iset_out_degree_count(n);
      
        vector<pair<int, int> > ds;
        for (int v = 0; v < n; v++){
          ds.emplace_back(forward_adj[0][v].size(), v);
        }
        sort(ds.begin(), ds.end());
        for (int vi = 0; vi < n; vi++){
          int v = ds[vi].second;
          if (std::max(iset_in_degree_count[v], iset_out_degree_count[v]) < iset_param){
            vs.push_back(v);
            for (int w : forward_adj[1][v]){
              iset_out_degree_count[w]++;
            }
            for (int w : forward_adj[0][v]){
              iset_in_degree_count[w]++;
            }
          }
        }
      
        for (int v : vs) {
          remove[v] = true;
        }
      }
    
      int dset_size = 0;
      size_t max_degree = 0;
      for (int v = 0; v < V; v++) {
        if (!remove[v]) {
          ns_offset[v] = dset_size++;
        } else {
          max_degree = std::max(max_degree, forward_adj[0][v].size());
        }
      }

      if (FLAGS_output_iparam_log) {
        JLOG_OPEN(("reduction.param_" + std::to_string(iset_param)).c_str()) {
          JLOG_PUT("num_vs", dset_size);
          JLOG_PUT("num_vs_ratio", (double) dset_size / V);
          JLOG_PUT("max_degree", max_degree);
        }
      }
      return dset_size;
    }

    void BuildLandmarks() {
      for (int k = 0; k < kNumBitParallelRoots; k++) {
        cerr << "Constructing bit-parallel spanning trees. ["
             << k + 1 << " / " << kNumBitParallelRoots << "]" << endl;
        bpspts[0][k].Build(roots[k], root_neighbors[k], &forward_adj[0], &forward_adj[1], &covered_vs, &ns_offset, I);
        bpspts[1][k].Build(roots[k], root_neighbors[k], &forward_adj[1], &forward_adj[0], &covered_vs, &ns_offset, I);
      }
    }
  
  public:
    BitParallelLandmarks() { }

    BitParallelLandmarks(const vector<vector<int> > forward_adj_[2], int iset_param) {
      Build(forward_adj_, iset_param);
    }

    void Build(const vector<vector<int> > forward_adj_[2], int iset_param) {
      Clear();
      CHECK_GT(forward_adj_[0].size(), 0);
      CHECK_GT(forward_adj_[1].size(), 0);
      this->forward_adj = forward_adj_;
      this->V = forward_adj[0].size();
      CHECK_GT(V, 0);
      SelectBitParallelLandmarks();
      I = SelectVerticesWithNS(iset_param);
      BuildLandmarks();

      if (FLAGS_remove_threshold > 0) {
        JLOG_ADD_OPEN("before") {
          JLOG_ADD("total_degree", DegreeSum());
        }
        for (int k = 0; k < kNumBitParallelRoots; k++) {
          bpspts[0][k].RemoveLowDegreeNeighbors(FLAGS_remove_threshold);
          bpspts[1][k].RemoveLowDegreeNeighbors(FLAGS_remove_threshold);
        }
        JLOG_ADD_OPEN("after") {
          JLOG_ADD("total_degree", DegreeSum());
        }
      }
    }

    void Clear() {
      V = I = iset_param = 0;
      for (int k = 0; k < kNumBitParallelRoots; k++) {
        bpspts[0][k].Clear();
        bpspts[1][k].Clear();
        covered_vs.clear();
        ns_offset.clear();
        root_neighbors[k].clear();
      }
    }

    inline int GetDistanceUpperbound(int s, int t) {
      CHECK_LE(0, s);
      CHECK_LT(s, V);
      CHECK_LE(0, t);
      CHECK_LT(t, V);
      int dist_upper = INF8;
      for (int k = 0; k < kNumBitParallelRoots; k++) {
        dist_upper = std::min<int>(dist_upper, bpspts[0][k].GetDistance(t) + bpspts[1][k].GetDistance(s));
      }

      for (int k = 0; k < kNumBitParallelRoots; k++) {
        int d = bpspts[0][k].GetDistance(t) + bpspts[1][k].GetDistance(s);
        if (d >= dist_upper + 2) {
          continue;
        }

        auto mask_s = bpspts[1][k].GetNeighborSet(s);
        auto neighbor_set_t = bpspts[0][k].GetNeighborSet(t);
        int td = d +
          ((neighbor_set_t.fst & mask_s.fst) ? -2 :
           (neighbor_set_t.snd & mask_s.fst) ? -1 :
           (neighbor_set_t.fst & mask_s.snd) ? -1 : 0);
        dist_upper = std::min(dist_upper, td);
      }
      return dist_upper;
    }

    // Return if a vertex v is covered by any of bit-parallel landmarks. 
    inline int SkipVertex(int v) const {
      CHECK_LE(0, v);
      CHECK_LT(v, V);
      return covered_vs[v];
    }
  
    // Insert an edge (s, t)
    void InsertEdge(int s, int t) {
      if (s == t) {
        return;
      }
      for (int k = 0; k < kNumBitParallelRoots; k++) {
        bpspts[0][k].InsertEdge(s, t);
        bpspts[1][k].InsertEdge(t, s);
      }
    }

    // Delete an edge (s, t)
    void DeleteEdge(int s, int t) {
      if (s == t) {
        return;
      }
      for (int k = 0; k < kNumBitParallelRoots; k++) {
        bpspts[0][k].DeleteEdge(s, t);
        bpspts[1][k].DeleteEdge(t, s);
      }
    }

    // Return the total size of data structures.
    size_t IndexSize() const {
      size_t total = 0;
      // since we don't include the size of the original graph, here we only add the size of pointers of the graph. 
      total += sizeof(forward_adj) + sizeof(I) + sizeof(V) + sizeof(iset_param);;

      // add each bit-prallel landmark's size. 
      for (int k = 0; k < kNumBitParallelRoots; k++) {
        total += bpspts[0][k].IndexSize();
        total += bpspts[1][k].IndexSize();
      }

      // add the sizes of each member variable. 
      total += VectorSize(covered_vs);
      total += VectorSize(ns_offset);
      total += sizeof(int) * roots.size();
      total += root_neighbors.size() * (sizeof(vector<int>) + 64 * sizeof(int));
      return total;
    }

    // Return the landmarks. 
    vector<int> GetRoots() const {
      return vector<int>(roots.begin(), roots.end());
    }

    // Return the selected neighbors for all landmarks. 
    vector<vector<int> > GetRootNeighbors() const {
      vector<vector<int> > res(kNumBitParallelRoots);
      for (int k = 0; k < kNumBitParallelRoots; k++) {
        for (size_t i = 0; i < root_neighbors[k].size(); i++) {
          res[k].push_back(root_neighbors[k][i]);
        }
      }
      return res;
    }
  
    // A function for debugging distance update procedures, which returns stored distances from landmarks. 
    vector<vector<int> > GetDistances(int dir = 0) const {
      vector<vector<int> > ds(kNumBitParallelRoots);
      for (int k = 0; k < kNumBitParallelRoots; k++) {
        for (int v = 0; v < V; v++) {
          ds[k].push_back(bpspts[dir][k].GetDistance(v));
        }
      }
      return ds;
    }
  
    // A function for debugging distance update procedures, which computes the distances from landmarks from scratch. 
    vector<vector<int> > GetDistancesFromScratch() const {
      vector<vector<int> > ds(kNumBitParallelRoots);
      for (int k = 0; k < kNumBitParallelRoots; k++) {
        vector<uint8_t> tmp_d(V, INF8);
        vector<neighbor_set_t> tmp_s(V), tmp_m(V);
        bpspts[0][k].BitParallelBFS(tmp_d, tmp_s);
        for(int d : tmp_d){
          ds[k].push_back(d);
        }
      }
      return ds;
    }

    //  A function for debugging mask recovery, which computes neighor sets without recovering them. 
    vector<vector<neighbor_set_t> >  GetNeighborSets() const {
      vector<vector<neighbor_set_t> >  res(kNumBitParallelRoots);
      for (int k = 0; k < kNumBitParallelRoots; k++){
        for (int v = 0; v < V; v++) {
          if (ns_offset[v] != -1){
            res[k].push_back(bpspts[0][k].GetNeighborSet(v));
          } else {
            res[k].push_back(std::make_pair(-1, 0));
          }
        }
      }

      return res;
    };

    //  A function for debugging mask recovery, which computes recovered neighbor sets for all the vertices. 
    vector<vector<neighbor_set_t> >  GetNeighborSetsWithRecovery() const {
      vector<vector<neighbor_set_t > >  res(kNumBitParallelRoots);
      for (int k = 0; k < kNumBitParallelRoots; k++) {
        for (int v = 0; v < V; v++) {
          res[k].push_back(bpspts[0][k].GetNeighborSet(v));
        }
      }
      return res;
    };

    //  A function for debugging neighbor set updates, which computes neighbor sets for all the vertices from scratch.
    vector<vector<neighbor_set_t> > GetNeighborSetsFromScratch() const {
      vector<vector<neighbor_set_t> >  res;
      for (int k = 0; k < kNumBitParallelRoots; k++){
        vector<uint8_t> tmp_d(V, INF8);
        vector<neighbor_set_t> tmp_s(V), tmp_m(V);
      
        bpspts[0][k].BitParallelBFS(tmp_d, tmp_s);
        for (int v = 0; v < V; v++) {
          if (ns_offset[v] != -1) {
            tmp_m[v] = tmp_s[v];
            tmp_m[v].snd &= ~tmp_m[v].fst;
          } else {
            tmp_m[v] = std::make_pair(-1, 0);
          }
        }
        res.push_back(tmp_m);
      }
      return res;
    };
  };

#endif /* BIT_PARALLEL_LANDMARKS_H */
