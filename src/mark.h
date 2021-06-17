#pragma once

#include "parlay/parallel.h"

namespace pargeo {
namespace hdbscanInternal {

template<class nodeT>
void markAll(nodeT *nd, size_t id){
  if(!nd->isLeaf() && nd->getId() != id){
    if (nd->size() > 2000) {
      parlay::par_do([&](){markAll(nd->L(), id);},
	     [&](){markAll(nd->R(), id);});
    } else {
      markAll(nd->L(), id);
      markAll(nd->R(), id);}
  }
  nd->setId(id);
}

template<class nodeT, class pointT, class UF>
void mark(nodeT *nd, UF *uf, pointT* s){

  if(nd->hasId()) {
    markAll(nd, uf->find(nd->getItem(0)-s));
    return;
  }

  nd->setId(uf->find(nd->getItem(0)-s));

  if(nd->isLeaf()){
    for (size_t i=1; i < nd->size(); ++i) {
      if (nd->getId() != uf->find(nd->getItem(i)-s)) {
        nd->resetId();
	return;}
    }
  } else {
    if (nd->size() > 2000) {
      parlay::par_do([&](){mark(nd->L(), uf, s);},
	     [&](){mark(nd->R(), uf, s);});
    } else {
      mark(nd->L(), uf, s);
      mark(nd->R(), uf, s);
    }
    if (nd->getId() != nd->L()->getId()) {
      nd->resetId();
      return;
    }
    if (nd->getId() != nd->R()->getId()) {
      nd->resetId();
      return;
    }
  }
}

} // End namespace emstInternal
} // End namespace pargeo
