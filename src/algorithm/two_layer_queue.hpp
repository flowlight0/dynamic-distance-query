#ifndef TWO_LAYER_QUEUE_H
#define TWO_LAYER_QUEUE_H

#include <cstdlib>
#include <iostream>
#include <vector>

template <typename T> class TwoLayerQueue {
  std::vector<T> data;
  size_t V;
  size_t curr_;
  size_t next_;
  size_t end_;
    
public:
  explicit TwoLayerQueue(size_t V) : data(V), V(V), curr_(0), next_(0), end_(0){
  }
  inline bool empty() const { return curr_ == next_; }
  inline bool full() const { return end_ == V; }
  inline T &front() { return data[curr_];}
  inline size_t size() const { return end_; }
  inline void pop() { ++curr_; assert(curr_ <= end_);}
  inline void push(const T &val){ data[end_++] = val; assert(end_ <= V);}
  inline void next() { assert(curr_ == next_); next_ = end_; }
  inline void clear() { curr_ = next_ = end_ = 0; }
  
  inline typename std::vector<T>::iterator begin() { return data.begin();}
  inline typename std::vector<T>::iterator end() { return data.begin() + end_;}
};

#endif /* TWO_LAYER_QUEUE_H */
