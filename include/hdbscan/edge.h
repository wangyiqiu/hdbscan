#pragma once

namespace pargeo {

  template<bool _directed, typename idxT>
  class _edge {
  public:
    static constexpr bool directed = _directed;
    idxT u,v;

    _edge(idxT _u, idxT _v): u(_u), v(_v) {
      if (!directed && u > v) std::swap(u,v);
    }

    _edge(): u(-1), v(-1) { }

    bool isEmpty() { return u == -1; }

    bool operator==(_edge e2) {
      return e2.u == u && e2.v == v;
    }

    bool operator!=(_edge e2) {
      return e2.u != u || e2.v != v;
    }
  };

  using edge = pargeo::_edge<false, size_t>;
  using dirEdge = pargeo::_edge<true, size_t>;

  template<bool _directed, typename idxT, typename weightT>
  class _weightedEdge {
  public:
    static constexpr bool directed = _directed;
    idxT u,v;
    weightT weight;

    _weightedEdge(idxT _u, idxT _v): u(_u), v(_v), weight(-1) {
      if (!directed && u > v) std::swap(u,v);
    }

    _weightedEdge(idxT _u, idxT _v, weightT _w): u(_u), v(_v), weight(_w) {
      if (!directed && u > v) std::swap(u,v);
    }

    _weightedEdge(): u(-1), v(-1), weight(-1) { }

    bool isEmpty() { return u == -1; }

    bool operator==(_weightedEdge e2) {
      return e2.u == u && e2.v == v;
    }

    bool operator!=(_weightedEdge e2) {
      return e2.u != u || e2.v != v;
    }
  };

  using wghEdge = pargeo::_weightedEdge<false, size_t, double>;
  using wghDirEdge = pargeo::_weightedEdge<true, size_t, double>;

} // END namespace pargeo
