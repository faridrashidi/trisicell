#ifndef UNWEIGHTED_GRAPH_H
#define UNWEIGHTED_GRAPH_H

#include <iostream>
#include <list>
#include <set>
#include <string>
#include <vector>
//#include <limits>
using namespace std;

#include "Utils.h"

void OutputQuotedString(ofstream &outFile, const char *buf);

// ***************************************************************************
// Buneman graph utilities
// ***************************************************************************
class BGVertex {
  friend class BunemanGraph;
  friend class BGEdge;
  friend class UnWeightedGraph;

public:
  BGVertex() : id(0), speciesID("") {}
  BGVertex(const string &nm) : name(nm), id(0), speciesID("") {}
  void AddBlock(bool blk) { blocks.push_back(blk); }
  void SetSpeciesID(const string &id) { speciesID = id; }
  friend ostream &operator<<(ostream &out, const BGVertex &v);

private:
  string name; // name of vertex, such as "v1"
  int id;      // unique name for node
  string speciesID;
  vector<bool> blocks; // false: first block, true: second block
};
typedef vector<BGVertex> LIST_VERTEX;

class BGEdge {
  friend class BunemanGraph;
  friend class UnWeightedGraph;

public:
  BGEdge() : v1Pos(-1), v2Pos(-1) { pv1 = pv2 = NULL; }
  BGEdge(string nm) : name(nm), v1Pos(-1), v2Pos(-1) { pv1 = pv2 = NULL; }
  BGEdge(string nm, int v1p, int v2p, LIST_VERTEX &listVerts)
      : name(nm), v1Pos(v1p), v2Pos(v2p) {
    pv1 = &listVerts[v1p];
    pv2 = &listVerts[v2p];
  }

private:
  string name; // name of edge, such as "s1"
  int v1Pos;
  int v2Pos;
  BGVertex *pv1; // end vertex1 of edge
  BGVertex *pv2; // end vertex2 of edge
};

// ***************************************************************************
typedef vector<BGEdge> LIST_EDGE;
typedef set<int> INCOMPATIBLE_SET; // each integer is for indexing
                                   // into matrix columns
typedef set<int> SPLIT_BLOCK_SET;
typedef list<INCOMPATIBLE_SET> INC_CLIQUE_LIST;
// typedef vector<int> SEQUENCE;			// represent sequence of
// bits

#define NIL_VERTEX 0x7FFFFFFF

// ***************************************************************************
// UnWightedGraph class
// ***************************************************************************
class UnWeightedGraph {
public:
  UnWeightedGraph() {}
  UnWeightedGraph(LIST_VERTEX &listVerts, LIST_EDGE &listEs)
      : listVertices(listVerts), listEdges(listEs) {}
  ~UnWeightedGraph() {}
  void SetVertices(LIST_VERTEX &verts) { listVertices = verts; }
  void SetEdges(LIST_EDGE &edges) { listEdges = edges; }
  int GetNumVertices() const { return listVertices.size(); }
  int GetAdjVert(int src, int lastAdj);
  bool IsConnected();
  void OutputGML(const char *fileName);
  bool IsNeighour(int i, int j);
  LIST_VERTEX &GetListVerts() { return listVertices; }

private:
  LIST_VERTEX listVertices;
  LIST_EDGE listEdges;
};

// ***************************************************************************
// UndirectedGraph class
// ***************************************************************************
class GraphVertex {
public:
  GraphVertex() {
    value = id = 0;
    visited = false;
  }
  GraphVertex(int id1) : visited(false) { id = id1; }
  GraphVertex(int id1, int val) : visited(false) {
    id = id1;
    value = val;
  }
  GraphVertex(const GraphVertex &rhs) {
    value = rhs.value;
    id = rhs.id;
    visited = rhs.visited;
  }
  //  GraphVertex&  operator=(const GraphVertex &rhs) {value = rhs.value; id =
  //  rhs.id; visited = rhs.visited; return this;}

  void SetVisited(bool f) { visited = f; }
  bool IsVisited() { return visited; }
  int GetID() { return id; }
  void SetValue(int v) { value = v; }
  int GetValue() { return value; }
  string GetLabel() const { return label; }
  void SetLabel(string lbl) { label = lbl; }

private:
  int value; // a vertex can have a value
  int id;    // id is unique, when removal no reuse
  bool visited;
  string label;
};

class GraphEdge {
public:
  GraphEdge() {
    vid1 = -1;
    vid2 = -1;
    value = -1;
  }
  GraphEdge(int id1, int id2) {
    vid1 = id1;
    vid2 = id2;
  }
  GraphEdge(int id1, int id2, int v) {
    vid1 = id1;
    vid2 = id2;
    value = v;
  }
  GraphEdge(const GraphEdge &rhs) {
    vid1 = rhs.vid1;
    vid2 = rhs.vid2;
    value = rhs.value;
    label = rhs.label;
  }

  void GetVertexIDs(int &v1, int &v2) {
    v1 = vid1;
    v2 = vid2;
  }
  int GetValue() { return value; }
  void SetValue(int v) { value = v; }
  string GetLabel() const { return label; }
  void SetLabel(string lbl) { label = lbl; }

private:
  int vid1;     // id of the vertex #1
  int vid2;     // id of vertex #2
  int value;    // an edge can have a value
  string label; // label of the edge
};

#if 0
class GraphEdgeIterator
{
public:
    GraphEdgeIterator( );

private:
    int curPos;
};
#endif

typedef int GRAPH_TRAV_POSITION;

// ***************************************************************************
// Graph classes
// ***************************************************************************
// This is a hopefully generic class for directed graph
// we store, at each vertex, the list of
class GenericGraph {

public:
  typedef vector<GraphEdge> EDGE_LIST;

  GenericGraph();
  virtual ~GenericGraph() {}

  // Basic graph functions
  virtual int AddVertex(int val);
  virtual bool RemoveVertex(int id);
  virtual bool AddEdge(int vid1, int vid2, int val) = 0;
  int GetNumVertices() const { return vertices.size(); }
  virtual int GetNumEdges() const;
  virtual int GetEdgeNum(int vid);
  virtual GraphEdge *GetEdgeByIndex(int vid, int index);
  virtual GraphEdge *GetEdge(int vid, int uid);
  virtual bool IsEdge(int vid, int uid);
  virtual bool FindVertexByID(int id, GraphVertex &v);
  virtual GraphVertex *FindVertex(int id);
  virtual void SetVertexVisited(int vid, bool flag);
  virtual bool IsVertexVisited(int vid);
  virtual void SetVertexLabel(int vid, string lbl);
  virtual GraphVertex *GetVertexByLabel(string lbl);
  virtual void SetEdgeLabel(int vid1, int vid2, string lbl);

protected:
  map<int, GraphVertex> vertices; // for faster access, use a map, indexed by id
  map<int, EDGE_LIST> adjacencyList; // Indexed by id of the vertex
  int nextId;
};

class UndirectedGraph : public GenericGraph {
  //    typedef vector<GraphEdge> EDGE_LIST;
public:
  UndirectedGraph();
  virtual ~UndirectedGraph();

  // Basic graph functions
  bool AddEdge(int vid1, int vid2, int val);
  int GetNumEdges() const;

  // bool GetFirstEdge (int vid, GraphEdge &e);
  // bool GetNextEdge( int vid, GraphEdge &e );

  GraphEdge *GetEdge(int vid, int uid);
  int GetFirstNode(GraphVertex &v); // return -1 when done
  int GetNextNode(GraphVertex &v);  // return id of the node

  // Some traversal utility here, need improvements later
  void InitTraversal(); // Set visited flag to false
  int TraversalFrom(int id, set<int> &listOfCCVertices); // Start traversal from
                                                         // node's id = id
  int FindUnvisitedNode(); // return the id of an unvisited node, if none return
                           // -1

  // Some basic functions here
  bool IsBipartitie();
  void FindComponents(set<set<int> > &comps);

  // DFS functions
  void InitPrevConfig();
  void DFSSetPrevNode(int u, int uprev);
  int DFSGetPrevNode(int u);
  void DFSSetNextNode(int u, int unext);
  int DFSGetNextNode(int u);

protected:
  // map<int, GraphVertex> vertices;      // for faster access, use a map,
  // indexed by id map<int, EDGE_LIST> adjacencyList;  // Indexed by id of the
  // vertex int nextId;

  // private member here
private:
  // Mainly used temporiliy for traversal
  map<int, int> prevMap;
  map<int, int> nextMap; // not really neccessary, but to make it simple....
  map<int, GraphVertex>::iterator itCurrent;
  //    GRAPH_TRAV_POSITION                     curpos;
  int *prevArray;
};

// This is a hopefully generic class for directed graph
// we store, at each vertex, the list of
class DirectedGraph : public GenericGraph {
public:
  // Basic graph functions
  bool AddEdge(int vid1, int vid2, int val);
  bool IsNodeSink(int vid);
  bool IsNodeSource(int vid);
  bool IsAcyclic();
  void TopologicalSort(vector<int> &listNodesIds);

  // Output
  void OutputGML(const char *fileName);
  void TrimTreeArcs(); // recursivly remove all nodes as sinks

private:
  void DFSVisitAcyclic(int nid, int &time, map<int, int> &nodesColor,
                       map<int, int> &nodesdval, map<int, int> &nodesfval,
                       vector<int> *plistFinishedNodes = NULL);
};

#endif // UNWEIGHTED_GRAPH_H
