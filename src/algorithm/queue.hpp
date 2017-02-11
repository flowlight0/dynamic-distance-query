#ifndef QUEUE_H
#define QUEUE_H

#include "common.hpp"

template <typename T> class Queue {
  std::vector<T> data;
  size_t V;
  size_t curr_;
  size_t end_;
public:
  explicit Queue(size_t V) : data(V), V(V), curr_(0), end_(0){

  };

  inline bool empty() const { return curr_ == end_; }
  inline bool full() const { return end_ == V; }
  inline T &front() { return data[curr_];}
  inline size_t size() const { return end_; }
  inline void pop() { ++curr_; assert(curr_ <= end_);}
  inline void push(const T &val){ data[end_++] = val; assert(end_ <= V);}
  inline void clear() { curr_ = end_ = 0; }

  inline typename std::vector<T>::iterator begin() { return data.begin();}
  inline typename std::vector<T>::iterator end() { return data.begin() + end_;}
};

#endif /* QUEUE_H */
