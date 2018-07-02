//----------------------------------------------------------------------------
// Generic solver for Minimum Cost Flow problem. Uses successive shortest path
// algorithm. Uses Dijkstra to find the shortest paths and a cool trick using
// node-potentials to avoid negative costs.
//
// 2018-06-27 21:15:32 leo
//----------------------------------------------------------------------------
#include <iostream>
#include <vector>
#include <limits>
#include <queue>
#include <algorithm>
#include <tuple>

using namespace std;

int g_dbg = 0;

template <typename T>
class Graph {
public:
  struct EdgeData {
    int v_;     // "to" vertex
    int c_;     // capacity
    int r_;     // residual capacity
    T w_;       // weight
    bool be_;   // back edge
    EdgeData (int v, int c, T w, bool be)
    : v_ (v), c_ (c), w_ (w), be_ (be), r_(c) {}
  };

public:
  Graph (int nVerts) 
    : m_adj (nVerts), m_badEdge (-1, 0, 0, false) {}

  ~Graph () {}

  void addEdge (int u, int v, int c, T w, bool be) {
    m_adj [u]. push_back (EdgeData (v, c, w, be));
  }

  int numVerts () const {
    return m_adj. size ();
  }

  const vector <EdgeData>& adj (int u) const {
    return m_adj [u];
  }

  vector <EdgeData>& adj (int u) {
    return m_adj [u];
  }

  EdgeData& getEdge (int u, int offset) {
    return m_adj [u] [offset];
  }

  EdgeData& getEdge (const pair <int, int>& nodeAndOffset) {
    return m_adj [nodeAndOffset.first] [nodeAndOffset.second];
  }

  EdgeData& findEdge (int u, int v) {
    for (auto& e : m_adj [u]) {
      if (e.v_ == v) {
        return e;
      }
    }
    return m_badEdge;
  }

  void resetResiduals () {
    for (auto& a : m_adj) {
      for (auto& e : a) {
        if (e.be_) 
          e.r_ = 0;
        else
          e.r_ = e.c_;
      }
    }
  }

  void printGraph () const {
    for (int u=0; u < m_adj. size (); ++u) {
      cout << u << " -> ";
      for (const auto& e : adj (u)) {
        cout << e.v_ << "(" << e.c_ << "," << e.w_ << "," << e.r_ << ") ";
      }
      cout << endl;
    }
  }

  void printGraph (const vector <T>& nodePot) const {
    for (int u=0; u < m_adj. size (); ++u) {
      cout << u << " (" << nodePot [u] << ") -> ";
      for (const auto& e : adj (u)) {
        cout << e.v_ << "(" << e.c_ << "," << e.w_ << ") ";
      }
      cout << endl;
    }
  }

private:
  vector < vector <EdgeData> > m_adj;
  EdgeData m_badEdge;
};

//----------------------------------------------------------------------------

#define NO_PARENT -1
#define FAKE_PARENT -2
#define BAD_OFF -1

template <typename T>
class MinCostFlow {

public:
  MinCostFlow (Graph <T>& g, int s, int t) 
    : m_g (g), m_potential (m_g. numVerts (), 0), m_s (s), m_t (t) {}
  ~MinCostFlow () {}

  pair <int, T> findMcf ();

private:

  struct ExpNode {
    int id_;  // current vertex
    int pid_; // parent vertex
    int off_; // offset of edge on parent
    T cost_;

    ExpNode (int id, int pid, int off, T cost) 
      : id_(id), pid_(pid), off_ (off), cost_(cost) {}

    // for mincost priority que. 
    bool operator< (const ExpNode& rhs) const {
      if (cost_ != rhs.cost_) return cost_ > rhs.cost_;
      else if (id_ != rhs.id_) return id_ > rhs.id_;
      else return pid_ > rhs.pid_;
    }
  };

  pair <int, T> getSpFlow (const vector < pair<int,int> >& pa);
  void updateRes (int flow, const vector < pair<int,int> >& pa);
  int st_dijkstra ();
  T computeCost ();

private:
  Graph <T>& m_g;
  vector <T> m_potential;
  int m_s, m_t;
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template <typename T>
pair <int, T> MinCostFlow<T>::findMcf () 
{
  int flow = 0;
  while (int f = st_dijkstra ()) {
    flow += f;
    if (g_dbg >1) cout << "MCF: flow " << f << endl;
    if (g_dbg) m_g. printGraph (m_potential);
  }
  T cost = computeCost ();
  return make_pair (flow, cost);
}

template <typename T>
T MinCostFlow<T>::computeCost ()
{
  T cost = 0;
  for (int u=0; u<m_g. numVerts (); ++u) {
    for (const auto& e : m_g. adj (u)) {
      if (e.be_) continue; 
      cost += (e.c_ - e.r_) * e.w_;
    }
  }
  return cost;
}

template <typename T>
int MinCostFlow<T>::st_dijkstra ()
{
  vector < pair<int,int> > par (m_g. numVerts (), make_pair (NO_PARENT, BAD_OFF));
  vector <T> cost (m_g. numVerts (), numeric_limits <T>::max ());
  std::priority_queue <ExpNode> q;
  cost [m_s] = 0;

  q. push (ExpNode (m_s, FAKE_PARENT, BAD_OFF, 0));
  if (g_dbg >1) cout << "Push: " << m_s << " " << 0 << endl;

  while (! q.empty ()) {
    int curNode = q.top ().id_;
    int parNode = q.top ().pid_;
    int parOff  = q.top ().off_;
    T curCost = q.top ().cost_;
    q. pop ();
    if (par [curNode].first != NO_PARENT) {
      continue; // already popped
    }
    par [curNode] = make_pair (parNode, parOff);
    cost [curNode] = curCost;
    if (g_dbg > 1) cout << "Pop: " << curNode << " " << curCost << "(" <<
      parNode << "," << parOff << ")" << endl;

    // what the...
    typename vector <typename Graph<T>::EdgeData>::const_iterator 
           eItr = m_g.  adj (curNode). begin (),
        eItrEnd = m_g.  adj (curNode). end ();
    // much better
    //auto eItr = m_g. adj (curNode). begin ();
    //auto eItrEnd = m_g. adj (curNode). end ();

    for (int off=0; eItr != eItrEnd; ++eItr, ++off) {
      const auto& e = *eItr;
      if (e.r_ == 0) {
        continue;  // no capacity
      }
      int childNode = e.v_;
      if (par [childNode].first != NO_PARENT) {
        continue; // child already popped
      }
      T edgeCost = m_potential [curNode] + e.w_ - m_potential [childNode];
      T newCost = curCost + edgeCost;
      if (newCost >= cost [childNode]) {
        continue; // prune
      }
      if (g_dbg > 1) cout << "Push: " << childNode << "(" << off << ") " << newCost << endl;
      q. push (ExpNode (childNode, curNode, off, newCost));
    }
  }

  if (par [m_t].first == NO_PARENT) {
    if (g_dbg > 1) cout << "No path" << endl;
    return false; // no path to target
  }

  int f; T c;
  tie (f, c) = getSpFlow (par);
  if (g_dbg > 1) cout << "getSpFlow " << f << " cost " << c << endl;

  updateRes (f, par);

  if (g_dbg > 1) {
    cout << "dists: "; for (auto t : cost) cout << t << " "; cout << endl;
  }
  copy (cost. begin (), cost.end (), m_potential. begin ());

  return f;
}

template <typename T>
pair <int, T> MinCostFlow<T>::getSpFlow (const vector < pair<int,int> >& par)
{
  int f = numeric_limits <int>::max ();
  T pcost = 0;
  int idx = m_t;
  while (par [idx].first != FAKE_PARENT) {
    if (idx == m_s) {
      cout << "ERROR: getSpFlow hit src " << idx << " par " << par [idx].first << endl;
      break;
    }
    typename Graph<T>::EdgeData& e1 = m_g. getEdge (par [idx]);
    f = min (f, e1.r_);

    if (g_dbg > 1) cout << idx << " <- ";

    if (e1.be_)
      pcost -= e1.w_;
    else
      pcost += e1.w_;
    idx = par [idx].first;
  }
  if (g_dbg > 1) cout << endl;

  return make_pair (f, pcost);
}

template <typename T>
void MinCostFlow<T>::updateRes (int flow, const vector < pair<int,int> >& par)
{
  int idx = m_t;
  while (par [idx].first != FAKE_PARENT) {
    typename Graph<T>::EdgeData& e1 = m_g. getEdge (par [idx]);
    e1.r_ -= flow;
    typename Graph<T>::EdgeData& e2 = m_g. findEdge (idx, par [idx].first);
    e2.r_ += flow;

    idx = par [idx].first;
  }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

main ()
{
  int t; cin >> t;
  for ( ; t; --t) {
    int n; cin >> n;  // n is nSrcs/nSnks in bpg
    int N = 2*n + 2;  // nSrcs + nSnks + s + t

    Graph <int> g (N);
    for (int u=1; u <= n; ++ u) {
      for (int v=n+1; v <= 2*n; ++ v) {
        int c=1;
        int w; cin >> w;
        g. addEdge (u, v, c, w, false); // regular edge
        g. addEdge (v, u, 0, -w, true); // back edge, zero cap initially
      }
    }

    int s = 0, t = 2*n+1;
    for (int u=1; u <= n; ++ u) {
      g. addEdge (s, u, 1, 0, false);
      g. addEdge (u, s, 0, 0, true);
    }
    for (int v=n+1; v <= 2*n; ++ v) {
      g. addEdge (v, t, 1, 0, false);
      g. addEdge (t, v, 0, 0, true);
    }

    if (g_dbg) { 
      g. printGraph ();
      cout << "Graph " << g. numVerts () << " vertices" << endl;
    }

    MinCostFlow <int> mcf (g, s, t);
    int flow, cost;
    tie (flow, cost) = mcf. findMcf ();
    cout << "MinCost: " << cost << " flow: " << flow << endl;

    if (0) {
      // testing residual reset
      g. resetResiduals ();
      MinCostFlow <int> mcf2 (g, s, t);
      int flow, cost;
      tie (flow, cost) = mcf2. findMcf ();
      cout << "MinCost2: " << cost << " flow: " << flow << endl;
    }

  }

}
