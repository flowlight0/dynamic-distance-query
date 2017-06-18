#include "common.hpp"
#include "gflags/gflags.h"
#include "algorithm/distance_with_bit_parallel.hpp"
using namespace std;

DEFINE_string(graph_file, "graph.txt", "Input graph file");
DEFINE_string(query_file, "query.txt", "Input query file");
DEFINE_int32(index_reduction_level, 1, "");

int main(int argc, char *argv[])
{
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  if (!FLAGS_undirected) {
    cerr << "The option --undirected must be specified since this CUI can handle only undirected graphs." << endl;
    return 1;
  }

  vector<pair<int ,int> > es;
  const int V = ReadGraph(FLAGS_graph_file, es);
  
  GraphDistance<10> bpd;
  bpd.Build(V, es, FLAGS_index_reduction_level);

  ifstream ifs;
  ifs.open(FLAGS_query_file);
  if (ifs.fail()) {
    cerr << "Failed to open query input file: " << FLAGS_query_file << endl;
    exit(-1);
  }

  string query, op;
  int u, v;
  while (getline(ifs, query)) {
    istringstream iss(query);
    iss >> op >> u >> v;
    if (op == "Q") {
      cout << bpd.QueryDistance(u, v) << endl;
    } else if (op == "EI") {
      bpd.InsertEdge(u, v);
      bpd.InsertEdge(v, u);
    } else if (op == "ED") {
      bpd.DeleteEdge(u, v);
      bpd.DeleteEdge(v, u);
    } else {
      cerr << "Invalid op " << op << ", skipping it" << endl;
    }
  }
  ifs.close();
  return 0;
}
