#include "common.hpp"

using namespace std;

DEFINE_bool(undirected, false, "");
DEFINE_bool(read_binary_file, false, "");

vector<pair<int, int> > GenerateErdosRenyi(int n, double p){
  vector<pair<int, int> > es;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      if (i != j && (double)rand() / RAND_MAX < p){
        es.emplace_back(i, j);
      }
    }
  }
  return es;
}

int ReadGraphBinary(const string &graph_file, vector<pair<int, int> > &es){
  es.clear();
  FILE *fp = fopen(graph_file.c_str(), "rb");
  
  if (fp == NULL){
    cerr << "Binary graph_file doesn't exist." << endl;
    return 0;
  }
  
  uint64_t V, E;
  
  CHECK(fread(&V, sizeof(V), 1, fp) == 1);
  CHECK(fread(&E, sizeof(E), 1, fp) == 1);
  for (size_t i = 0; i < E; i++){
    int u, v;
    CHECK(fread(&u, sizeof(u), 1, fp) == 1);
    CHECK(fread(&v, sizeof(v), 1, fp) == 1);
    CHECK_LT(u, int(V));
    CHECK_LT(v, int(V));
    es.emplace_back(u, v);
  }
  
  fclose(fp);
  return int(V);
}

int ReadGraph(const string &graph_file, vector<pair<int, int> > &es){
  
  int V = 0;

  if (FLAGS_read_binary_file && (V = ReadGraphBinary(graph_file + ".bin", es)) > 0) {
    return V;
  }
  
  ifstream ifs(graph_file.c_str());
  if (!ifs.good()){
    cerr << "Error: open graph_file." << endl;
    exit(EXIT_FAILURE);
  }
  
  unordered_map<int, int> vertex2id;
  
  for (int u, v; ifs >> u >> v;){
    if (vertex2id.count(u) == 0) vertex2id[u] = V++;
    if (vertex2id.count(v) == 0) vertex2id[v] = V++;
    u = vertex2id[u];
    v = vertex2id[v];
    if (u != v){
      es.emplace_back(u, v);
    }
  }

  if (FLAGS_undirected){
    size_t m = es.size();
    for (size_t i = 0;i < m; i++){
      es.emplace_back(es[i].snd, es[i].fst);
    }
  }
  JLOG_PUT("setting.undirected", FLAGS_undirected);


  sort(es.begin(), es.end());
  es.erase(unique(es.begin(), es.end()), es.end());
  ifs.close();
  return V;
}
