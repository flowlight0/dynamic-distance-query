#ifndef UNION_FIND_H
#define UNION_FIND_H

#include <vector>

class UnionFind{
  int num_components;
  std::vector<int> parent;
  std::vector<int> weight;
  std::vector<int> rank;
public:
  UnionFind(int N) : num_components(N),
                     parent(std::vector<int>(N)),
                     weight(std::vector<int>(N, 1)),
                     rank(std::vector<int>(N, 0)){
    for(int i = 0; i < N; i++) parent[i] = i;
  }
  
  int find(int x){
    if(x == parent[x]) return x;
    else return parent[x] = find(parent[x]);
  }
  
  int size(int x){
    return weight[find(x)];
  }

  
  bool same(int x, int y){
    return find(x) == find(y);
  }
    
  void unite(int x, int y){
    x = find(x);
    y = find(y);
    if(x == y) return;

    num_components--;
    if(rank[x] < rank[y]){
      weight[y] += weight[x];
      parent[x] = y;
    }else{
      weight[x] += weight[y];
      parent[y] = x;
      if(rank[x] == rank[y]) rank[y]++;
    }
  }
  
  int count(){
    return num_components;
  }
};


#endif /* UNION_FIND_H */
