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

#include <atomic>
#include <tuple>
#include <limits>
#include "wspd.h"
#include "parBuf.h"
#include "getTime.h"
#include "parlay/utilities.h"
#include "hdbscanBccp.h"

namespace pargeo {
  namespace hdbscanInternal {

    template<class nodeT, class objT, typename floatT>
    inline void nodeCD(nodeT *nd,
		       sequence<floatT> &coreDist,
		       sequence<floatT> &cdMin,
		       sequence<floatT> &cdMax,
		       nodeT* root,
		       objT* P
		       ) {
      if(nd->isLeaf()){
	for (intT i = 0; i < nd->size(); ++i) {
	  if (coreDist[nd->getItem(i)-P] > cdMax[nd-root]) {
	    cdMax[nd-root] = coreDist[nd->getItem(i)-P];
	  }
	  if (coreDist[nd->getItem(i)-P] < cdMin[nd-root]) {
	    cdMin[nd-root] = coreDist[nd->getItem(i)-P];
	  }
	}
      } else {
	if (nd->size() > 2000) {
	  parlay::par_do([&](){nodeCD(nd->L(), coreDist, cdMin, cdMax, root, P);},
			 [&](){nodeCD(nd->R(), coreDist, cdMin, cdMax, root, P);});
	} else {
	  nodeCD(nd->L(), coreDist, cdMin, cdMax, root, P);
	  nodeCD(nd->R(), coreDist, cdMin, cdMax, root, P);
	}
	cdMax[nd-root] = max(cdMax[nd->L()-root], cdMax[nd->R()-root]);
	cdMin[nd-root] = min(cdMin[nd->L()-root], cdMin[nd->R()-root]);
      }
    }

    //geo-separated || unreachable
    template<class nodeT, typename floatT>
    inline bool unreachable(nodeT *u,
			    nodeT *v,
			    sequence<floatT> &cdMin,
			    sequence<floatT> &cdMax,
			    nodeT *root
			    ) {
      floatT circleDiam_u = 0;
      floatT circleDiam_v = 0;
      floatT circleDistance = 0;
      for (int d = 0; d < u->dim; ++ d) {
	floatT uTmpDiff = u->getMax(d) - u->getMin(d);
	floatT vTmpDiff = v->getMax(d) - v->getMin(d);
	floatT uTmpAvg = (u->getMax(d) + u->getMin(d))/2;
	floatT vTmpAvg = (v->getMax(d) + v->getMin(d))/2;
	circleDistance += (uTmpAvg - vTmpAvg) * (uTmpAvg - vTmpAvg);
	circleDiam_u += uTmpDiff * uTmpDiff;
	circleDiam_v += vTmpDiff * vTmpDiff;
      }
      circleDiam_u = sqrt(circleDiam_u);
      circleDiam_v = sqrt(circleDiam_v);

      floatT myRadius = max(circleDiam_u, circleDiam_v)/2;
      floatT myDiam = max(2*myRadius, cdMax[u-root]);
      myDiam = max(myDiam, cdMax[v-root]);

      circleDistance = sqrt(circleDistance) - circleDiam_u/2 - circleDiam_v/2;
      bool geoSep = circleDistance >= 2 * myRadius;
      circleDistance = max(circleDistance, cdMin[u-root]);
      circleDistance = max(circleDistance, cdMin[v-root]);

      if (circleDistance >= myDiam) {
	return true || geoSep;
      } else {
	return false || geoSep;
      }
    }

    template<class nodeT, class UF>
    struct rhoUpdateParallel {
      using floatT = double;

      std::atomic<floatT> rho;
      floatT beta;
      UF *uf;
      nodeT* tree;
      sequence<floatT>& coreDist;
      sequence<floatT>& cdMin;
      sequence<floatT>& cdMax;

      rhoUpdateParallel(floatT _beta,
			UF *_uf,
			nodeT* _tree,
			sequence<floatT>& _coreDist,
			sequence<floatT>& _cdMin,
			sequence<floatT>& _cdMax
			) :
	beta(_beta), uf(_uf), coreDist(_coreDist), cdMin(_cdMin), cdMax(_cdMax), tree(_tree) {
	rho = std::numeric_limits<floatT>::max();}

      void run(nodeT *u, nodeT *v) {
	floatT myDist = max(nodeDistance(u, v), cdMin[u - tree]);
	myDist = max(myDist, cdMin[v - tree]);
	parlay::write_min(&rho, myDist, std::less<floatT>()); //check
      }

      bool moveon(nodeT *u, nodeT *v) {
	if (u->hasId() && u->getId() == v->getId()) return false; // filtering
	if (u->size()+v->size() <= beta) return false; // not in E_u, not considered for rho
	floatT myDist = max(nodeDistance(u, v), cdMin[u - tree]);
	myDist = max(myDist, cdMin[v - tree]);
	if (myDist >= rho) return false; // no subsequent finds can update rho
	return true;
      }

      bool start(nodeT *u) {
	if (u->size() > beta) {
	  return true;
	} else {
	  return false;// if node size < beta, so would children
	}
      }

      floatT getRho() { return rho;}

      bool wellSeparated(nodeT* u, nodeT* v, floatT s) {
	return unreachable(u, v, cdMin, cdMax, tree);}
    };

    template<class nodeT, class UF>
    struct wspGetParallel {
      using floatT = double;
      using bcpT = std::tuple<typename nodeT::objT*,
			      typename nodeT::objT*,
			      typename nodeT::objT::floatT>;
      using bufT = parBuf<bcpT>;

      floatT rhoLo;
      floatT rhoHi;
      floatT beta;
      nodeT *tree;
      UF *uf;
      bufT **out;
      sequence<floatT>& coreDist;
      sequence<floatT>& cdMin;
      sequence<floatT>& cdMax;
      typename nodeT::objT* P;

      wspGetParallel(floatT betaa,
		     floatT rhoLoo,
		     floatT rhoHii,
		     nodeT *treee,
		     UF *uff,
		     sequence<floatT>& _coreDist,
		     sequence<floatT>& _cdMin,
		     sequence<floatT>& _cdMax,
		     typename nodeT::objT* _P
		     ) :
	beta(betaa), rhoLo(rhoLoo), rhoHi(rhoHii),
	tree(treee), uf(uff),
	coreDist(_coreDist), cdMin(_cdMin), cdMax(_cdMax), P(_P) {
	size_t procs = num_workers();
	out = (bufT**) malloc(sizeof(bufT*)*procs);
	parallel_for(0, procs, [&](size_t p) {
				 out[p] = new bufT(tree->size()/procs);
			       });
      }

      ~wspGetParallel() {
	size_t procs = num_workers();
	parallel_for(0, procs, [&](size_t p) {
				 delete out[p];});
	free(out);
      }

      sequence<bcpT> collect() {
	int procs = num_workers();
	return parBufCollect<bcpT>(out, procs);
      }

      void run(nodeT *u, nodeT *v) {
	auto bcp = hdbscanBccp(u, v, coreDist, P);
	if (u->size() + v->size() <= beta &&
	    std::get<2>(bcp) >= rhoLo &&
	    std::get<2>(bcp) < rhoHi) {
	  auto tmp = out[worker_id()]->increment();
	  get<0>(*tmp) = get<0>(bcp);
	  get<1>(*tmp) = get<1>(bcp);
	  get<2>(*tmp) = get<2>(bcp);
	}
      }

      bool moveon(nodeT *u, nodeT *v) {
	if (u->hasId() && u->getId() == v->getId()) {return false;}
	floatT dist = max(nodeDistance(u, v), cdMin[u-tree]);
	dist = max(dist, cdMin[v-tree]);
	if (dist >= rhoHi) return false; // too separated to be considered
	dist = max(nodeFarDistance(u, v), cdMax[u-tree]);
	dist = max(dist, cdMax[v-tree]);
	if (dist < rhoLo) return false; // too close to be considered, bug!!
	return true;
      }

      bool start(nodeT *u) {
	if (max(u->diag(), cdMax[u-tree]) >= rhoLo) {
	  return true;
	} else {
	  return false;
	}
      }
      bool wellSeparated(nodeT* u, nodeT* v, floatT s) {
	return unreachable(u, v, cdMin, cdMax, tree);}
    };

    template <class nodeT, class UF>
    sequence<std::tuple<typename nodeT::objT*,
			typename nodeT::objT*,
			typename nodeT::objT::floatT>>
    filterWspdParallel(double t_beta,
		       double t_rho_lo,
		       double& t_rho_hi,
		       nodeT *t_kdTree,
		       UF *t_mst,
		       sequence<floatT> &coreDist,
		       sequence<floatT> &cdMin,
		       sequence<floatT> &cdMax,
		       typename nodeT::objT* P
		       ) {
      using floatT = double;
      using objT = typename nodeT::objT;
      using bcpT = std::tuple<objT*, objT*, floatT>;

      auto myRho = rhoUpdateParallel<nodeT, UF>(t_beta, t_mst, t_kdTree, coreDist, cdMin, cdMax);

      pargeo::computeWspdParallel<nodeT, rhoUpdateParallel<nodeT, UF>>(t_kdTree, &myRho);

      auto mySplitter = wspGetParallel<nodeT, UF>(t_beta, t_rho_lo, myRho.getRho(), t_kdTree, t_mst, coreDist, cdMin, cdMax, P);

      pargeo::computeWspdParallel<nodeT, wspGetParallel<nodeT, UF>>(t_kdTree, &mySplitter);

      t_rho_hi = myRho.getRho();
      return mySplitter.collect();
    }

  } // End namespace hdbscanInternal
} // End namespace pargeo
