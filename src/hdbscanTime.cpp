#include <iostream>
#include <algorithm>
#include "parlay/utilities.h"
#include "parlay/parallel.h"
#include "pointIO.h"
#include "getTime.h"
#include "parseCommandLine.h"

#include "hdbscan/point.h"
#include "hdbscan/hdbscan.h"

using namespace std;
using namespace parlay;
using namespace pargeo;
using namespace pargeo::pointIO;

// *************************************************************
//  TIMING
// *************************************************************

template<int dim>
void timeHdbscan(sequence<point<dim>>& P, size_t minPts, int rounds, char* outFile, int perturb) {
  if (perturb) {
    parallel_for(0, P.size(), [&](size_t i) {
	for (int j = 0; j < dim; ++ j) {
	  double myRand = P[i][j] / 10000;
	  P[i][j] += -myRand + 2*myRand*hash64(i)/numeric_limits<size_t>::max();
	}});
  }

  for (int i=0; i < rounds; i++) {
    timer t0; t0.start();
    sequence<wghEdge> E = hdbscan<dim>(P, minPts);
    dendrogram(E, P.size());
    cout << "timing = " << t0.stop() << endl;
  }
}

int main(int argc, char* argv[]) {
  //commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] [-p <0/1 perturb points>] <inFile>");
  commandLine P(argc,argv,"[-m <minPts value>] <inFile>");

  char* iFile = P.getArgument(0);
  size_t minPts = P.getOptionIntValue("-m",1);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);
  int perturb = P.getOptionIntValue("-p",0);

  int dim = readHeader(iFile);

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> Points = readPointsFromFile<pargeo::point<2>>(iFile);
    timeHdbscan<2>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> Points = readPointsFromFile<pargeo::point<3>>(iFile);
    timeHdbscan<3>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> Points = readPointsFromFile<pargeo::point<4>>(iFile);
    timeHdbscan<4>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> Points = readPointsFromFile<pargeo::point<5>>(iFile);
    timeHdbscan<5>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> Points = readPointsFromFile<pargeo::point<6>>(iFile);
    timeHdbscan<6>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> Points = readPointsFromFile<pargeo::point<7>>(iFile);
    timeHdbscan<7>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 8) {
    parlay::sequence<pargeo::point<8>> Points = readPointsFromFile<pargeo::point<8>>(iFile);
    timeHdbscan<8>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 9) {
    parlay::sequence<pargeo::point<9>> Points = readPointsFromFile<pargeo::point<9>>(iFile);
    timeHdbscan<9>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 10) {
    parlay::sequence<pargeo::point<10>> Points = readPointsFromFile<pargeo::point<10>>(iFile);
    timeHdbscan<10>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 11) {
    parlay::sequence<pargeo::point<11>> Points = readPointsFromFile<pargeo::point<11>>(iFile);
    timeHdbscan<11>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 12) {
    parlay::sequence<pargeo::point<12>> Points = readPointsFromFile<pargeo::point<12>>(iFile);
    timeHdbscan<12>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 13) {
    parlay::sequence<pargeo::point<13>> Points = readPointsFromFile<pargeo::point<13>>(iFile);
    timeHdbscan<13>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 14) {
    parlay::sequence<pargeo::point<14>> Points = readPointsFromFile<pargeo::point<14>>(iFile);
    timeHdbscan<14>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 15) {
    parlay::sequence<pargeo::point<15>> Points = readPointsFromFile<pargeo::point<15>>(iFile);
    timeHdbscan<15>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 16) {
    parlay::sequence<pargeo::point<16>> Points = readPointsFromFile<pargeo::point<16>>(iFile);
    timeHdbscan<16>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 17) {
    parlay::sequence<pargeo::point<17>> Points = readPointsFromFile<pargeo::point<17>>(iFile);
    timeHdbscan<17>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 18) {
    parlay::sequence<pargeo::point<18>> Points = readPointsFromFile<pargeo::point<18>>(iFile);
    timeHdbscan<18>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 19) {
    parlay::sequence<pargeo::point<19>> Points = readPointsFromFile<pargeo::point<19>>(iFile);
    timeHdbscan<19>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 20) {
    parlay::sequence<pargeo::point<20>> Points = readPointsFromFile<pargeo::point<20>>(iFile);
    timeHdbscan<20>(Points, minPts, rounds, oFile, perturb);}
  else {
    std::cout << "dim = " << dim << "\n";
    throw std::runtime_error("dimension not supported yet (2-20 is supported at the moment)");
  }
}
