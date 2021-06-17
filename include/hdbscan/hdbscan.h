#pragma once

#include "parlay/sequence.h"
#include "edge.h"
#include "point.h"

namespace pargeo {

  template<int dim>
  parlay::sequence<pargeo::wghEdge> hdbscan(parlay::sequence<pargeo::point<dim>> &, size_t);

  typedef tuple<size_t, size_t, double, size_t> dendroNode;

  parlay::sequence<dendroNode> dendrogram(parlay::sequence<pargeo::wghEdge> &, size_t);

}
