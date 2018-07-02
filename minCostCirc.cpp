//----------------------------------------------------------------------------
// Minimum cost circulation problem, cast as min-cost flow. Uses lemon package
// to solve.
//
// Compilation:  g++ -std=c++11 minCostCirc.cpp -lemon
// Test: a.out < circtest_1.txt
//
// 2018-07-01 22:33:26 leo 
//----------------------------------------------------------------------------
#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <lemon/list_graph.h>
#include <lemon/capacity_scaling.h>
#include <lemon/lgf_writer.h>

using namespace std;
using namespace lemon;

#define INF_CAP 1000

// input format (Ref circtest_1.txt):
//    numVerts
//    u supply nArcs (v cost)+
//    u supply nArcs (v cost)+
//    ...
//
// Reads problem definition and builds Digraph. 
// Adds master src and sink nodes and edges sNode -> {U}, {U} -> tNode. 
//    Capacity of sNode -> {U} edges are the supply at each vertex
//    Capacity of {U} -> tNode edges are the demands are each vertex (hardcoded to 1)
// Capacity of other edges are set to INF_CAP. Should be set to max flow possible.
pair <ListDigraph::Node, ListDigraph::Node> 
readGraph (ListDigraph& g, 
           ListDigraph::ArcMap <int>& arcCap, 
           ListDigraph::ArcMap <int>& arcCosts)
{
  int numVerts; cin >> numVerts;
  for (int i=0; i < numVerts; ++i) {
    g. addNode ();
  }
  ListDigraph::Node sNode = g.addNode (); // src node
  ListDigraph::Node tNode = g.addNode (); // snk node
  if (0) cout << "Graph added " << numVerts << " verts + s " << g.id (sNode) 
              << " t " << g.id (tNode) << endl;

  for (int i=0; i<numVerts; ++i) {
    int u, s, nArcs; 
    cin >> u >> s >> nArcs;
    
    auto uNode =  g. nodeFromId (u);
    auto s_u = g. addArc (sNode, uNode); 
    if (0) cout << "addArc " << g.id (sNode) << " -> " << g.id (uNode) << endl;

    auto u_t = g. addArc (uNode, tNode); 
    if (0) cout << "addArc " << g.id (uNode) << " -> " << g.id (tNode) << endl;

    arcCap [s_u] = s;
    arcCap [u_t] = 1;
    arcCosts [s_u] = 0;
    arcCosts [u_t] = 0;

    for ( ; nArcs; --nArcs) {
      int v, cost; cin >> v >> cost;
      auto vNode =  g. nodeFromId (v);
      auto u_v = g. addArc (uNode, vNode); 
      if (0) cout << "addArc " << g.id (uNode) << " -> " << g.id (vNode) << endl;

      arcCosts [u_v] = cost;
      arcCap [u_v] = INF_CAP;
    }
  }
  return make_pair (sNode, tNode);
}

// read flow assignments from solver and print only for real edges
void readFlow (const CapacityScaling <ListDigraph>& solver, 
               const ListDigraph& g,
               const ListDigraph::Node& sNode,
               const ListDigraph::Node& tNode)
{
  for (ListDigraph::ArcIt a (g); a != INVALID; ++a) {
    int flow = solver. flow (a);
    if (flow && 
        g. source (a) != sNode && 
        g. target (a) != tNode) 
      cout << g.id (g. source (a)) << " -> " << 
              g.id (g. target (a)) << " flow " << flow << endl;
  }
}

main ()
{
  int t; cin >> t;
  for ( ; t; --t) {
    cout << "Problem " << t << endl;

    ListDigraph g;
    ListDigraph::ArcMap <int> arcCap (g);
    ListDigraph::ArcMap <int> arcCosts (g);

    ListDigraph::Node sNode, tNode;

    tie (sNode, tNode) = readGraph (g, arcCap, arcCosts);

    int sNodeSupply = 0;
    for (ListDigraph::OutArcIt aId (g, sNode); aId != INVALID; ++aId) {
      sNodeSupply += arcCap [aId];
    }
    cout << "Supply at sNode " << sNodeSupply << endl;
    // using this as max flow that needs to be pushed

    if (0) {
      DigraphWriter <ListDigraph> gWriter (g, std::cout);
      gWriter. run ();
    }

    CapacityScaling <ListDigraph> solver (g);
    solver. costMap (arcCosts);
    solver. upperMap (arcCap);
    solver. stSupply (sNode, tNode, sNodeSupply);
    auto ret = solver. run ();
    cout << "Solver status: ";

    switch (ret) {
      case CapacityScaling <ListDigraph>::INFEASIBLE: 
        cout << "Infeasible" << endl; break;
      case CapacityScaling <ListDigraph>::OPTIMAL: 
        cout << "Optimal. Total cost: " << solver. totalCost () << endl; 
        readFlow (solver, g, sNode, tNode);
        break;
      case CapacityScaling <ListDigraph>::UNBOUNDED:
        cout << "Unbounded" << endl; break;
      default: cout << "Unknown" << endl;
    }
  }
}
