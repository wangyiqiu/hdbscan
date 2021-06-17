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

#include "limits.h"
#include <tuple>

namespace pargeo {

  using namespace std;

  namespace bccpInternal {

    template <typename objT>
    struct bcp {
      using floatT = typename objT::floatT;

      objT* u;
      objT* v;
      floatT dist;

      bcp(): u(NULL), v(NULL), dist(std::numeric_limits<floatT>::max()) {}

      void update(objT* _u, objT* _v, floatT _dist) {
	if (_dist < dist) {
	  u = _u; v = _v; dist = _dist;
	}
      }
    };

    template <typename nodeT>
    bcp<typename nodeT::objT> bcpBruteforce(nodeT* n1, nodeT* n2) {
      bcp<typename nodeT::objT> r;
      for (size_t i = 0; i < n1->size(); ++ i) {
	for (size_t j = 0; j < n2->size(); ++ j) {
	  double tmp = n1->at(i)->dist( *(n2->at(j)) );
	  r.update(n1->at(i), n2->at(j), tmp);
	}
      }
      return r;
    }

    template <typename nodeT>
    inline void bcpHelper(nodeT* n1, nodeT* n2, bcp<typename nodeT::objT>* r) {
      if (nodeDistance(n1, n2) > r->dist) return;

      if (n1->isLeaf() && n2->isLeaf()) {

	for (size_t i=0; i<n1->size(); ++i) {
	  for (size_t j=0; j<n2->size(); ++j) {
	    r->update(n1->at(i), n2->at(j),
		      n1->at(i)->dist(*n2->at(j)));
	  }
	}

      } else {

	if (n1->isLeaf()) {

	  if (nodeDistance(n1, n2->L()) < nodeDistance(n1, n2->R())) {
	    bcpHelper(n1, n2->L(), r); bcpHelper(n1, n2->R(), r);
	  } else {
	    bcpHelper(n1, n2->R(), r); bcpHelper(n1, n2->L(), r);
	  }

	} else if (n2->isLeaf()) {

	  if (nodeDistance(n2, n1->L()) < nodeDistance(n2, n1->R())) {
	    bcpHelper(n1->L(), n2, r); bcpHelper(n1->R(), n2, r);
	  } else {
	    bcpHelper(n1->R(), n2, r); bcpHelper(n1->L(), n2, r);
	  }

	} else {

	  pair<nodeT*, nodeT*> ordering[4]; //todo change to tuple
	  ordering[0] = make_pair(n1->L(), n2->L());
	  ordering[1] = make_pair(n1->L(), n2->R());
	  ordering[2] = make_pair(n1->R(), n2->L());
	  ordering[3] = make_pair(n1->R(), n2->R());

	  auto cmp = [&](pair<nodeT*,nodeT*> p1, pair<nodeT*,nodeT*> p2) {
									  return nodeDistance(p1.first, p1.second) < nodeDistance(p2.first, p2.second);};
	  sort(ordering, ordering + 4, cmp);

	  for (int o=0; o<4; ++o) {
	    bcpHelper(ordering[o].first, ordering[o].second, r);}

	}

      }
    }

  } // End namespace bcpInternal

  template <typename nodeT>
  tuple<typename nodeT::objT*,
	typename nodeT::objT*,
	typename nodeT::objT::floatT> bccp(nodeT* n1, nodeT* n2) {
    using namespace bccpInternal;
    using floatT = double;

    auto r = bcp<typename nodeT::objT>();

    bcpHelper(n1, n2, &r);

    // auto verify = bcpBruteforce(n1, n2);
    // if (r.u != verify.u || r.v != verify.v) {
    //   throw std::runtime_error("bcp wrong");
    // }
    return tuple(r.u, r.v, r.dist);
  }

} // End namespace
