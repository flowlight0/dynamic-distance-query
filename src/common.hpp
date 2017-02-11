#ifndef COMMON_H
#define COMMON_H

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <queue>
#include <stack>
#include <math.h>
#include <string>
#include <vector>
#include <limits>
#include <sstream>
#include <set>
#include <unordered_map>
#include "jlog.hpp"
#include "gflags/gflags.h"


#define CHECK(expr)                                 \
  if (expr) {                                       \
  } else {                                          \
    fprintf(stderr, "CHECK Failed (%s:%d): %s\n",   \
            __FILE__, __LINE__, #expr);             \
    exit(EXIT_FAILURE);                             \
  }

#define CHECK_PRED(expr1, expr2, pred)                        \
  if ((expr1) pred (expr2)) {                                 \
  } else {                                                    \
    fprintf(stderr, "CHECK Failed (%s:%d): %s %s %s\n", \
            __FILE__, __LINE__, #expr1, #pred, #expr2);       \
    std::cerr << #expr1 << " = " << expr1 << ".\n"            \
              << #expr2 << " = " << expr2 << ".\n";           \
      exit(EXIT_FAILURE);                                     \
  }

#define CHECK_EQ(expr1, expr2)  CHECK_PRED(expr1, expr2, ==)
#define CHECK_LE(expr1, expr2)  CHECK_PRED(expr1, expr2, <=)
#define CHECK_LT(expr1, expr2)  CHECK_PRED(expr1, expr2, < )
#define CHECK_GE(expr1, expr2)  CHECK_PRED(expr1, expr2, >=)
#define CHECK_GT(expr1, expr2)  CHECK_PRED(expr1, expr2, > )


#define ALL(S) (S).begin(), (S).end()
#define fst first
#define snd second


const uint8_t INF8 = 100;

DECLARE_bool(undirected);

template <typename T> size_t VectorSize(const std::vector<T> &vec){
  return sizeof(vec) + sizeof(T) * vec.size();
}

template <typename T> size_t Vector2DSize(const std::vector<std::vector<T> > &vecs){
  size_t res = sizeof(vecs);
  for (const auto &vec : vecs){
    res += VectorSize(vec);
  }
  return res;
}

template <typename T> void MakeUnique(std::vector<T> &vec){
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}


template <typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &vec){
  out << "[";
  for (size_t i = 0; i < vec.size(); i++){
    if (i > 0) out << ", ";
    out << vec[i];
  }
  out << "]";
  return out;
}

template <typename T, size_t nmemb>
std::ostream &operator<<(std::ostream &out, const std::array<T, nmemb> &vec){
  out << "[";
  for (size_t i = 0; i < nmemb; i++){
    if (i > 0) out << ", ";
    out << vec[i];
  }
  out << "]";
  return out;
}

template <typename S, typename  T>
std::ostream &operator<<(std::ostream &out, const std::unordered_map<S, T> &m){
  int c = 0;
  out << "[";
  for (const auto &e : m){
    if (c++ > 0) out << ", ";
    out << e;
  }
  out << "]";
  return out;
}


template<typename T>
void ParseSpaceSeparatedString(const std::string &str, std::vector<T> *res) {
  std::istringstream ss(str);
  res->clear();
  for (T t; ss >> t; ) res->push_back(t);
}

template<typename T>
void ParseCommaSeparaedString(std::string str, std::vector<T> *res) {
  replace(str.begin(), str.end(), ',', ' ');
  ParseSpaceSeparatedString(str, res);
}

int ReadGraph(const std::string &graph_file, std::vector<std::pair<int, int> > &es);
std::vector<std::pair<int, int> > GenerateErdosRenyi(int n, double p);


#endif /* COMMON_H */







