#ifndef RBT_H
#define RBT_H

//
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

#include "BinaryMatrix.h"
#include "Utils.h"
#include "Utils2.h"

// define a leaf-labeled rooted binary tree
// note that we do not store the leaf label explicitly
// rather, we just assume leafs are not the same
// and leafs are simply labeled as 0, 1, 2, 3, 4 etc
// thus, it is enough to ORDER leaves
// IMPORTANT: we assume left node is always smaller (in
// the sense of its minimum left leave is always smaller
// than minimum right leaves

// a class for tree node
class RBTNode {
public:
  // create a leaf node
  RBTNode(int lvid)
      : pLeft(NULL), pRight(NULL), pParent(NULL), lvid(lvid), height(0.0) {}
  // create an internal node with two
  RBTNode(RBTNode *pLeft, RBTNode *pRight);
  ~RBTNode() { Clear(); }

  // operation
  void SetHeight(double ht) { height = ht; }
  double GetHeight() const { return height; }
  RBTNode *CopySubTree();
  void AddToLeftEdge(int lvid);
  void AddToRightEdge(int lvid);
  RBTNode *AddSibling(int lvid);
  void AddSiblingToLeaf(int lvid);
  RBTNode *FindLeaf(int lvid,
                    int &ponid); // IMPORTANT, in traversal,
                                 // assume post-order search, and return the how
                                 // many nodes visited so far
  bool RemoveLeafSelf();         // only remove self if it is a leaf
  void DetachSubtree();          // detach this node from the rest of the tree
  RBTNode *AttachSubtree(RBTNode *pSib);

  // access
  RBTNode *GetLeftChild() const { return pLeft; }
  RBTNode *GetRightChild() const { return pRight; }
  RBTNode *GetParent() { return pParent; }
  RBTNode *GetSibling();
  void SetLeftChild(RBTNode *pLeft) { this->pLeft = pLeft; }
  void SetRightChild(RBTNode *pRight) { this->pRight = pRight; }
  void SetParent(RBTNode *pParent) { this->pParent = pParent; }
  int GetLeafId() { return lvid; }
  void SetLeafId(int idNew) { this->lvid = idNew; }
  RBTNode *GetLeftMostChild();
  int GetMinLeaveId();
  void GetLeaves(set<int> &lvs);
  bool IsLeaf() const;
  int GetNumLeavesUnder();
  bool IsLeftChild();
  bool IsRoot() { return this->pParent == NULL; }
  void Dump() const;
  string GetNewick() const;
  void OutputNodeGML(ofstream &ofs);
  void OutputEdgeGML(ofstream &ofs);

  // memory. free recursively
  void Clear();

private:
  void AdjustLRChildUpwards();

  // two children
  RBTNode *pLeft;
  RBTNode *pRight;
  RBTNode *pParent;
  int lvid;
  double height; // useful in some situations, normalized to between 0-1

  // utility
  static int idNodeNextToUse;
};

// define triplets
// Triplets are important for rooted tree, since the set of triplets
// uniquely define a RBT
typedef struct {
  // note by convention, a < b. But c is on the other side of partition (a,b), c
  int a;
  int b;
  int c;
} TripleLeaves;

// define for traversal
typedef struct {
  RBTNode *pCurNode;
} TraversRecord;

// sometimes, we want to an ID for the tree
// for the use say HMM states. We just simply
// use integer. Note that we need a way of interpreting this state
// and convert back and force between the real tree and its ID
typedef int RBT_ID;

// main class
class RBT {
public:
  // different ways of initializing a tree
  // it can be by a supplied id
  RBT(int numLeaves, RBT_ID tid);
  RBT(const RBT &rhs);
  // interop with simple representation
  RBT(int numLeaves, const vector<int> &listNodeLabels,
      const vector<int> &listParentNodePos, const vector<double> &listEdgeDist);
  RBT &operator=(const RBT &rhs);
  // bool operator == (const RBT &rhs) { return IsSame(rhs); }
  ~RBT();

  // ID functions
  RBT_ID GetId();
  RBT_ID MapToId();
  bool IsSame(const RBT &tr) const;

  // splits functions
  bool IsSplitContained(const set<int> &split); // test whether a split is in
                                                // the tree
  void GetAllSplits(vector<set<int> > &listSplits);

  // SPR function
  void FindSPRDistOneNgbrs(set<int> &ngbrIds);
  void FindSPRDistOneNgbrs(vector<RBT *> &ngbrTrees);
  void FindSPRDistOneNgbrsRestricted(vector<RBT *> &ngbrTrees,
                                     const vector<RBT *> &ConstraintTrees);
  bool IsOneSPRAway(const RBT &rbt) const; // testing whether it is one or two
                                           // SPR away
  bool IsTwoSPRAway(const RBT &rbt) const;
  static void Consolidate(RBT &treeOpt, RBT &treeCmp);

  // editing
  bool RemoveLeaf(int lvid);
  void ReconstructById(RBT_ID tid);
  bool ReconstructNewick(const string &strNewick);
  void PruneLargeIdNodes(int idThres);
  void DeleteLeaves(set<int> &lvids); // delete leaves designated
  void RealignLeaves(); // sometimes, say after leave is deleted, leaves are no
                        // longer contiguous, this op sets it back to contiguous
  void AugamentDupRows(
      const vector<REMOVED_ROWS_INFO> &rmLvsStage); // restore the leaves
                                                    // removed during matrix
                                                    // preprocessing
  void SetRoot(RBTNode *pRootNew);
  RBTNode *GetRoot() { return pRoot; }

  // dynamic functions: allow adding new nodes
  bool AddLeaf(int pos);

  // access
  void GetLeaves(set<int> &lvs);
  void Dump() const;
  void OutputGML(const char *fileName);

  // Int-op with another format
  void RetrievePlainDesc(int &numLeaves, vector<int> &listNodeLabels,
                         vector<int> &listParentNodePos,
                         vector<double> &listEdgeDist);
  int GetNodesNum() { return 2 * numLeaves - 1; }
  string GetNewick() const;
  int GetLeafNum() { return numLeaves; }
  bool IsEmpty() const { return pRoot == NULL && numLeaves == 0; }

  // compare
  int Compare(RBT &rhs);
  bool IsSameUnrootedTree(RBT &rhs);
  void CollectTips();
  RBTNode *GetTip(int id);
  void GetAllTips(vector<RBTNode *> &tips);

private:
  RBT() {}     // do not allow default construction
  void Init(); // common initialization
  // void ConsTripleMap();   // save all the triples
  // support traversal
  bool InitPostorderTranvers(TraversRecord &tr);
  bool NextPostorderTranvers(TraversRecord &tr);
  void RetrieveSplits();
  RBTNode *FindLeaf(int lvid, int &ponid);
  RBTNode *ReconstructNewickInternal(const string &strNewick);
  bool InternalAddleaf(int lvid, int pos);
  bool ReconstructByPlainDesc(const vector<int> &listNodeLabels,
                              const vector<int> &listParentNodePos,
                              const vector<double> &listEdgeDist);
  void SetLvids(const vector<int> &mapLvids); // configure the name for the
                                              // leaves

  // save a dynamic root node
  RBTNode *pRoot;

  // we also save the splits
  map<set<int>, bool> mapSplitsInTree;
  map<int, RBTNode *> mapTipPtrs;

  // note we do not normally allow morhping the tree
  // EXCEPT during initialtion. Since convert to id
  // can be slow, we cache it
  int numLeaves;
  RBT_ID tid;

  // collect of triples
  // map< TripleLeaves, bool > mapTriples;
};

///////////////////////////////////////////////////////////////////////////////////////
// useful functions

int GetNumRBT(int nlv);

#endif // RBT_H
