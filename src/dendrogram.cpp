#include <tuple>
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "unionFind.h"
#include "getTime.h"

#include "hdbscan/point.h"
#include "hdbscan/edge.h"
#include "hdbscan/hdbscan.h"

using namespace std;
using namespace parlay;
using namespace pargeo;

parlay::sequence<pargeo::dendroNode> pargeo::dendrogram(parlay::sequence<pargeo::wghEdge> &edges, size_t n) {
  timer t; t.start();

  sequence<pargeo::wghEdge> edgesSorted =
    parlay::sort(make_slice(edges), [&](wghEdge e1, wghEdge e2) {
				      return e1.weight < e2.weight;
				    });

  unionFind uf = unionFind<size_t>(n);

  size_t idx = n;

  sequence<size_t> idxMap(n);

  sequence<size_t> sizes(n);

  parallel_for(0, n,[&](int i) {
		      idxMap[i] = i;
		      sizes[i] = 1;
		    });

  auto dendro = parlay::sequence<dendroNode>(edgesSorted.size());

  for(size_t i=0; i<n-1; ++i){
    size_t u = uf.find(edgesSorted[i].u);
    size_t v = uf.find(edgesSorted[i].v);
    dendro[i] = tuple(idxMap[u], idxMap[v], edgesSorted[i].weight, sizes[u]+sizes[v]);
    uf.link(u, v);
    size_t newIdx = uf.find(u);
    idxMap[newIdx] = idx;
    sizes[newIdx] = sizes[u]+sizes[v];
    idx++;
  }

  cout << "dendrogram-time = " << t.stop() << endl;

  // for (auto d: dendro) {
  //   cout << get<0>(d) << " " << get<1>(d) << " ";
  //   cout << get<2>(d) << " " << get<3>(d) << endl;
  // }
  return dendro;
}
