// This code is part of the project "Fast Parallel Algorithms for Euclidean
// Minimum Spanning Tree and Hierarchical Spatial Clustering"
// Copyright (c) 2021 Yiqiu Wang, Shangdi Yu, Yan Gu, Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma once

namespace pargeo {

  template<typename T> struct parBuf {
    using intT = size_t;

    // current buffer
    T *m_arr;
    intT m_currentSize;
    intT m_count;

    // array of buffers
    T **m_parent;
    intT *m_parentSizes;
    intT m_parentUsed;
    intT m_parentTotal;
    const intT m_defaultParent = 10;

    parBuf(intT t_initSize) {
      m_arr = (T *) malloc(sizeof(T) *t_initSize);
      m_parent = (T **) malloc(sizeof(T*) *m_defaultParent);
      m_parentSizes = (intT *) malloc(sizeof(intT) *(m_defaultParent+1));

      m_currentSize = t_initSize;
      m_count = 0;

      m_parent[0] = m_arr;
      m_parentSizes[0] = t_initSize;

      m_parentUsed = 1;
      m_parentTotal = m_defaultParent;
    }

    intT size() {
      intT tmp = 0;
      for(intT i = 0; i < m_parentUsed-1; ++ i) {
	//cout << "size: parent " << i << ": " << m_parentSizes[i] << endl;
	tmp += m_parentSizes[i];
      }
      //cout << "size: last parent " << ": " << m_count << endl;
      tmp += m_count;
      return tmp;
    }

    intT incrementParent() {
      T **parent1 = (T **) malloc(sizeof(T) * m_parentTotal * 2);
      intT *parentSizes1 = (intT *) malloc(sizeof(intT) * (m_parentTotal * 2+1));

      for(intT i = 0; i < m_parentUsed; ++ i) {
	parent1[i] = m_parent[i];
	parentSizes1[i] = m_parentSizes[i];
      }

      free(m_parent);
      free(m_parentSizes);

      m_parent = parent1;
      m_parentSizes = parentSizes1;
      m_parentTotal *= 2;
      return m_parentUsed ++;
    }

    inline void finalize() {
      m_parentSizes[m_parentUsed-1] = m_count;
    }

    T* increment() {
      if (m_count + 1 > m_currentSize) {
	// book keeping for current arr
	finalize();

	// allocate new arr
	T *arr1 = (T *) malloc(sizeof(T) * m_currentSize * 2);

	if (m_parentUsed < m_parentTotal) {
	  m_parent[m_parentUsed++] = arr1;
	} else {
	  m_parent[incrementParent()] = arr1;
	}

	m_arr = arr1;
	m_count = 0;
	m_currentSize *= 2;
      }
      return m_arr + (m_count++);
    }

    ~parBuf() {
      for (intT i = 0; i < m_parentUsed; ++ i) {
	free(m_parent[i]);
      }

      free(m_parent);
    }
  };

  template <typename T>
  T prefixSumSerial(T* data, size_t s, size_t e) {
    T res = 0;
    for (size_t i = s; i < e; ++i) {
      res += data[i];
      data[i] = res - data[i];
    }
    return res;
  }

  template<typename T>
  parlay::sequence<T> parBufCollect(parBuf<T> **t_threadVecs, size_t P) {
    using namespace parlay;
    using intT = size_t;

    parlay::sequence<intT> vecSizes(P);
    for(int p = 0; p < P; ++ p) {
      t_threadVecs[p]->finalize();
      vecSizes[p] = t_threadVecs[p]->size();
    }

    intT total = scan_inplace(make_slice(vecSizes));

    T* all= (T*) malloc(sizeof(T) * total);

    parallel_for(0, P,
		 [&](intT p) {
		   auto buf = t_threadVecs[p];
		   intT threadTotal = prefixSumSerial<intT>(buf->m_parentSizes, 0, buf->m_parentUsed);
		   buf->m_parentSizes[buf->m_parentUsed] = threadTotal; // there's 1 extra space at the end of array

		   parallel_for (0, buf->m_parentUsed, [&](intT parent) {
		       for (intT elem = 0; elem < buf->m_parentSizes[parent+1]-buf->m_parentSizes[parent]; ++ elem) {
			 all[vecSizes[p] + buf->m_parentSizes[parent] + elem] = buf->m_parent[parent][elem];
		       }}, 1);
		 }, 1);

    return sequence<T>(all, all + total);
  }

} // End namespace
