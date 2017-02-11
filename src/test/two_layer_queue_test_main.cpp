#include <iostream>
#include "common.hpp"
#include "algorithm/two_layer_queue.hpp"
#include "gtest/gtest.h"
#include <algorithm>
using namespace std;

// Verify that BFSs using std::queue and TwoLayerQueue show the same results. 
void CheckTwoLayerQueue(int n, const vector<pair<int, int>> &es){
  vector<vector<int> > G(n);
  for (const auto &e : es){
    G[e.fst].push_back(e.snd);
    G[e.snd].push_back(e.fst);
  }

  TwoLayerQueue<int> two_layer_queue(n);

  for (int s = 0; s < n; s++){
    vector<int> dist1(n, -1), dist2(n, -1);
    
    queue<int> que;
    que.push(s);
    dist1[s] = 0;
    while (!que.empty()){
      int v = que.front(); que.pop();
      for (int w : G[v]){
        if (dist1[w] == -1){
          dist1[w] = dist1[v] + 1;
          que.push(w);
        }
      }
    }

    dist2[s] = 0;
    two_layer_queue.push(s);
    two_layer_queue.next();
    while (!two_layer_queue.empty()){
      while (!two_layer_queue.empty()){
        int v = two_layer_queue.front();
        two_layer_queue.pop();
        for (int w : G[v]){
          if (dist2[w] == -1){
            dist2[w] = dist2[v] + 1;
            two_layer_queue.push(w);
          }
        }
      }
      two_layer_queue.next();
    }
    two_layer_queue.clear();
    ASSERT_EQ(dist1, dist2) << dist1 << "  " << dist2 << endl;
  }
}



TEST(TWO_LAYER_QUEUE, RANDOM){
  const int num_vs = 50;
  const int num_ts = 10;
  for (int q = 0; q < num_ts; q++){
    vector<pair<int, int> > es = GenerateErdosRenyi(num_vs, 0.2);
    CheckTwoLayerQueue(num_vs, es);
  }
}
