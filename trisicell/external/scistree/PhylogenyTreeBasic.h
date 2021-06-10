#ifndef PHYLOGENY_TREE_BASIC_H
#define PHYLOGENY_TREE_BASIC_H

#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
#include <stack>
#include <string>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "Utils.h"

using namespace std;

//*****************************************************************************
// Utility classes
//*****************************************************************************

// utilities for Newick format
class NewickUtils {
public:
  NewickUtils() {}

  static void RetrieveLabelSet(const string &strNW,
                               multiset<string> &setLabels);
  static bool FindSplitIn(const string &strNW, string &strPart1,
                          string &strPart2);
  static void UpdateLabells(string &strNW,
                            const map<string, string> &mapOldLabelToNew);
  static string RemoveBrLenFromTree(string &strNW);
  static void ConsolidateSinglChildChain(string &strNW);
  static double GetLenAt(const string &strNW, int posLen);
};

// map between string-based taxa to integer based id (used internally by the
// code)
class TaxaMapper {
public:
  //
  TaxaMapper();

  // utility
  bool IsInitialized() { return fInit; }
  void SetInitialized(bool f) { fInit = f; }
  void InitToDec1Mode(int numTaxa);
  bool IsEmpty();
  bool IsIdIn(int id);
  int AddTaxaString(const string &str);
  void AddTaxaStringWithId(int tid, const string &str);
  int GetId(const string &str);
  string GetString(const int id);
  string ConvIdStringWithOrigTaxa(const string &strId);
  int GetNumTaxaInMapper() const { return mapIdToStr.size(); }
  void GetAllTaxaIds(set<int> &taxaIndices) const;
  void GetAllTaxaStrs(set<string> &setStrs) const;
  void Dump() const;
  static string ExtractIdPartFromStr(const string &strIdNW);
  static int GetIdFromStr(const string &strPart, TaxaMapper *pTMapper);

private:
  map<string, int> mapStrToId;
  map<int, string> mapIdToStr;
  int curId;
  bool fInit;
};

//*****************************************************************************
// Defintions and utilties class, not for external use.
// Myabe I should create a separate file for these implementation-only stuff.
// Later
// ****************************************************************************
typedef enum { PHY_TN_DEFAULT_SHAPE = 0, PHY_TN_RECTANGLE = 1 } TREE_NODE_SHAPE;

class TreeNode {
  friend class PhylogenyTreeBasic;
  friend class PhylogenyTree;

public:
  TreeNode();
  TreeNode(int iid);
  ~TreeNode();

  TreeNode *Copy();
  void AddChild(TreeNode *pChild, const vector<int> &labels);
  void AddEdgeLabelToChild(int cIndex, int lbl);
  void RemoveChild(TreeNode *pChild);
  void RemoveAllChildren();
  void DetachAllChildren();
  void DetachSelf();
  void SetLength(double len) { lenBranchAbove = len; }
  double GetLength() const { return lenBranchAbove; }
  void SetLabel(const string str) { label = str; }
  bool IsLeaf() const { return listChildren.size() == 0; }
  void AddNodeValue(int val) { nodeValues.push_back(val); }
  int GetChildrenNum() const { return listChildren.size(); }
  int GetNumNodesUnder(bool fInternalOnly,
                       bool fAddNonBinary) const; // include itself if this is
                                                  // an internal node
  int GetLevel() const; // level: leaf at 0, internal: longest path to some leaf
                        // under
  TreeNode *GetChild(int i) { return listChildren[i]; }
  void GetDescendentLabelSet(set<int> &labelSet);
  bool IsAncesterOf(TreeNode *pAssumedDescend, int &branchIndex);
  int GetNumEdgesToAncestor(TreeNode *pAssumedAncestor);
  int GetID() const { return id; }
  void SetID(int i) { id = i; }
  string GetLabel() const { return label; }
  void SetUserLabel(const string &str) { labelUserProvided = str; }
  string GetUserLabel() const { return labelUserProvided; }
  void RemoveLabels();
  void RemoveLabelsPar();
  void IncEdgeLabelsBy(int offset, bool fSub);
  int GetIntLabel() const;
  void SetIntLabel(int lbl);
  TREE_NODE_SHAPE GetShape() { return shape; }
  void SetShape(TREE_NODE_SHAPE param) { shape = param; }
  void GetEdgeLabelsAtBranch(int i, vector<int> &labels) {
    labels = listEdgeLabels[i];
  }
  void GetEdgeLabelsToChild(TreeNode *pChild, vector<int> &lbls);
  TreeNode *GetParent() { return parent; }
  void SetParent(TreeNode *ppar) { parent = ppar; }
  TreeNode *GetRoot() const;
  void GetSiblings(vector<TreeNode *> &listSibs);
  void GetAllChildren(set<TreeNode *> &setChildren) const;
  void GetAllDescendents(set<TreeNode *> &setDescendents);
  void GetAllLeavesUnder(set<TreeNode *> &setDescendents);
  void GetAllLeavesIdUnder(set<int> &setDescendents);
  void GetAllLeafLabeles(vector<string> &listLeafLabels);
  void GetAllLeafIntLabeles(vector<int> &listLeafLabels);
  void GetAllDistinctLeafLabeles(set<string> &setLeafLabels);
  void GetAllDescendIntLbls(set<int> &setIntLbs);
  void GetAllAncestors(set<TreeNode *> &listAncestors);
  string GetShapeLabel(const set<int> &idTerms,
                       map<int, int> &mapNodeLabel) const;
  string GetShapeLabel(const set<int> &idTerms, bool fSort = true) const;
  // string GetShapeLabelDistinct(const set<int> &idTerms) const;
  string
  GetShapeLabelNodeBrNum(map<TreeNode *, pair<int, int> > &mapNodeNumBrannches,
                         vector<int> &listORderedLeaves);
  TreeNode *GetMRCA(TreeNode *pOther);
  void Order();
  bool IsMulfurcate();
  bool IsCheryNode() {
    return (GetChildrenNum() == 2 && GetChild(0)->IsLeaf() == true &&
            GetChild(1)->IsLeaf());
  }
  bool IsRoot() const { return parent == NULL; }
  int GetChildIndex(TreeNode *pchild) const;
  void Binarize(int &idToUseNext);
  int GetMaxIdWithinSubtree() const;
  void Dump() const;

private:
  vector<TreeNode *> listChildren;
  vector<vector<int> > listEdgeLabels; // What labels is used in the edge
  TreeNode *parent;
  int id;                 // id of this node, should be UNIQUE
  vector<int> nodeValues; // A node can have several values, for example, nodes
                          // labeling CAUTION: we assume node value is >=0 !!!!!
  string label;
  string labelUserProvided; // this ist he label before any conversion
  TREE_NODE_SHAPE shape;
  double lenBranchAbove;
};

// ***************************************************************************
// Utilities
// ***************************************************************************
class PhylogenyTreeBasic;

class PhylogenyTreeIteratorBacktrack {
public:
  PhylogenyTreeIteratorBacktrack(PhylogenyTreeBasic &pt) : phyTree(pt) {}
  void Init();
  void Next();
  void Back(); // do not continue going downwards (i.e. do not explore its
               // descendent)
  bool IsDone();
  TreeNode *GetCurrNode();

private:
  PhylogenyTreeBasic &phyTree;
  stack<TreeNode *> stackNodesToExplore;
  // TreeNode *pCurr;
};

class PhylogenyTreeIterator {
public:
  PhylogenyTreeIterator(PhylogenyTreeBasic &pt) : phyTree(pt) {}
  void Init();
  void Next();
  bool IsDone();
  TreeNode *GetCurrNode();

private:
  PhylogenyTreeBasic &phyTree;
  stack<TreeNode *> stackPostorder;
  // TreeNode *pCurr;
};

// ***************************************************************************
// Define phylogeny tree class
// ***************************************************************************

class PhylogenyTreeBasic {
  friend class PhylogenyTreeIterator;
  friend class PhylogenyTreeIteratorBacktrack;

public:
  PhylogenyTreeBasic(); // Empty tree
  virtual ~PhylogenyTreeBasic();
  PhylogenyTreeBasic *Copy();
  void InitPostorderWalk(); // when walk, return the value of the node if any
  TreeNode *NextPostorderWalk();
  void OutputGML(const char *inFileName);
  void OutputGMLNoLabel(const char *inFileName);
  void ConsNewick(string &strNewick, bool wGridLen = false,
                  double gridWidth = 1.0, bool fUseCurLbl = false);
  void ConsNewickSorted(string &strNewick, bool wGridLen = false,
                        double gridWidth = 1.0, bool fUseCurLbl = false);
  void ConsNewickEdgeLabel(string &strNewick);
  TreeNode *AddTreeNode(TreeNode *parNode, int id);
  void ConsOnNewick(const string &nwString, int numLeaves = -1,
                    bool fBottomUp = false, TaxaMapper *pTMapper = NULL);
  void ConsOnNewickDupLabels(const string &nwString,
                             TaxaMapper *pTMapper = NULL);
  void ConsOnNewickEdgeLabelTree(const string &nwString);
  int GetNumVertices() const;
  int GetNumLeaves();
  int GetNumInternalNodes();
  void GetNodeParInfo(vector<int> &nodeIds, vector<int> &parPos);
  void GetNodeParInfoNew(vector<int> &nodeIds, vector<int> &parPos);
  bool ConsOnParPosList(const vector<int> &parPos, int numLeaves = -1,
                        bool fBottupUpLabel = false);
  void GetLeaveIds(set<int> &lvids);
  void GetLeafIntLabels(set<int> &setIntLabels);
  void GetLeavesIdsWithLabel(const string &label, set<int> &lvids);
  void GetLeavesWithLabels(const set<string> &setLabels,
                           set<TreeNode *> &setLvNodes);
  void UpdateIntLabel(const vector<int> &listLabels);
  TreeNode *GetRoot() const { return rootNode; }
  void SetRoot(TreeNode *rn) {
    if (rootNode != NULL)
      delete rootNode;
    rootNode = rn;
  }
  void SetRootPlain(TreeNode *rn) { rootNode = rn; }
  void GetAllLeafLabeles(vector<string> &listLeafLabels) {
    rootNode->GetAllLeafLabeles(listLeafLabels);
  }
  void GetAllLeafIntLabeles(vector<int> &listLeafLabels) {
    rootNode->GetAllLeafIntLabeles(listLeafLabels);
  }
  string GetShapeLabel(const set<int> &idTerms,
                       map<int, int> &mapNodeLabel) const {
    return rootNode->GetShapeLabel(idTerms, mapNodeLabel);
  }
  string GetShapeLabel(const set<int> &idTerms, bool fSort = true) const {
    return rootNode->GetShapeLabel(idTerms, fSort);
  }
  // string GetShapeLabelDistinct(const set<int> &idTerms ) const { return
  // rootNode->GetShapeLabelDistinct(idTerms); }
  string
  GetShapeLabelNodeBrNum(map<TreeNode *, pair<int, int> > &mapNodeNumBrannches,
                         vector<int> &listORderedLeaves);
  bool TestIsomorphic(PhylogenyTreeBasic &treeOther,
                      map<TreeNode *, TreeNode *> &mapOldNodeToNew) const;
  void Reroot(TreeNode *pRootDesc); // pRootDesc: the node in the current tree
                                    // (must be, but we will not check) which
                                    // will be root's descendent
  void GetAllLeafNodes(vector<TreeNode *> &listLeafNodes) const;
  void GetAllNodes(vector<TreeNode *> &listLeafNodes) const;
  void Order() { rootNode->Order(); }
  bool IsMulfurcate() { return GetRoot()->IsMulfurcate(); }
  void CleanNonLabeledLeaves();
  void RemoveNode(TreeNode *pn);
  void RemoveNodeKeepChildren(TreeNode *pn);
  void RemoveDegreeOneNodeAt(TreeNode *pn);
  void RemoveDegreeOneNodes();
  void RemoveEdgeLabels();
  void RemoveEdgeLabelsToLeaves();
  void IncEdgeLabelsBy(int offset);
  void ConsPhyTreeFromClusters(const set<set<int> > &setClusters);
  static void RemoveDescendentsFrom(set<TreeNode *> &setTreeNodes);
  void FindCladeOfSubsetLeaves(const set<TreeNode *> &setLeaves,
                               set<set<TreeNode *> > &setSubtreeClades);
  void FindCladeOfSubsetLeavesExact(const set<TreeNode *> &setLeaves,
                                    set<set<TreeNode *> > &setSubtreeClades);
  static void
  GroupLeavesToSubtrees(const set<TreeNode *> &setLeaves,
                        const set<set<TreeNode *> > &cladeNodesToProc,
                        set<set<TreeNode *> > &setSubtreeClades);
  static void
  GroupLeavesToSubtreesSamePar(const set<TreeNode *> &setLeaves,
                               const set<set<TreeNode *> > &cladeNodesToProc,
                               set<set<TreeNode *> > &setSubtreeClades);
  static void GroupNodesWithCommonPars(
      const set<TreeNode *> &setNodes,
      map<TreeNode *, set<TreeNode *> > &mapNodesWithSamePar);
  void GetAllClades(set<set<int> > &setClades);
  void GetAllCladesList(vector<set<int> > &listClades);
  void GetAllCladesById(set<set<int> > &setClades);
  void GetAllCladeNodess(set<set<TreeNode *> > &setClades);
  void GetAllCladeGroupsIntLabel(
      multiset<multiset<multiset<int> > > &setCladeGroupsDupLabels,
      multiset<int> &rootClade);
  TreeNode *GetSubtreeRootForLeaves(const set<TreeNode *> &setLvNodes);
  void GetSubtreesWithMaxSize(set<TreeNode *> &setSTRoots,
                              int maxSzSubtree) const;
  void GetMaxSubtrees(set<TreeNode *> &setSTRootsIdents);
  void MakeSubtreeUnrefined(TreeNode *pSubtree);
  void Binarize();
  void CreatePhyTreeFromLeavesWithLabels(const set<string> &setLeafLabels,
                                         PhylogenyTreeBasic &treeToProc,
                                         bool fUseOldTaxonName);
  void AssignLeafLabels(const map<int, string> &mapLeafLbls);
  void ReassignLeafLabels(const map<string, string> &mapLeafLbls);
  void SetUserLabelToCurrLabels();
  void SetLabelsToCurrUserLabels();
  int GetMaxDegree() const;
  static bool GetSiblingsPairFrom(const set<TreeNode *> &setNodesToChoose,
                                  pair<TreeNode *, TreeNode *> &pairSibs);
  static bool GetSiblingsNodesFrom(const set<TreeNode *> &setNodesToChoose,
                                   set<TreeNode *> &setSibs);
  static void FindAllLabelsInSubtrees(const set<TreeNode *> &setSTRoots,
                                      set<string> &setLabels);
  static void
  FindDescendentsOfNodeWithin(TreeNode *pAnc,
                              const set<TreeNode *> &setNodesToChoose,
                              set<TreeNode *> &setDescendents);
  void Dump() const;

protected:
  void PostOrderPushStack(TreeNode *treeNode,
                          stack<TreeNode *> &stackPostorder);
  string ConsNewickTreeNode(TreeNode *pNode, bool wGridLen = false,
                            double gridWidth = 1.0, bool fUseCurLbl = false,
                            bool fSort = false, bool fOutEdgeLabel = false);
  TreeNode *ConsOnNewickSubtree(const string &nwStringPart, int &leafId,
                                int &invId, int numLeaves = -1,
                                bool fBottomUp = false,
                                TaxaMapper *pTMapper = NULL);
  bool ConvParPosToNewick(const vector<int> &parPos, string &strNewick);
  void ConvParPosToNewickSubtree(int nodeInd, const vector<int> &parPos,
                                 string &strNewick);
  TreeNode *ConsOnNewickSubtreeDupLabels(const string &nwStringPart, int &invId,
                                         int &leafId,
                                         TaxaMapper *pTMapper = NULL);
  // void GetSubtreesWithMaxSizeExcludeTaxa(set<TreeNode *> &setSTRoots, int
  // maxSzSubtree, const set<string> &setTaxaAllowed) const; int GetIdFromStr(
  // const string &strPart, TaxaMapper *pTMapper );

  // Privaet data members
  TreeNode *rootNode;

  // Postoder traversal
  stack<TreeNode *> stackPostorder;
  int numLeaves;
};

//*****************************************************************************
string GetStringFromId(int id);
int GetNewickNumLeaves(const string &strNewick, char chSepLeft = '(',
                       char chSepRight = ')', char midSep = ',');
bool GetTripleType(TreeNode *pn1, TreeNode *pn2, TreeNode *pn3,
                   pair<pair<TreeNode *, TreeNode *>, TreeNode *> &triple);
bool ReadinPhyloTreesNewick(ifstream &inFile, int numLeaves,
                            vector<PhylogenyTreeBasic *> &treePtrList,
                            TaxaMapper *pTMapper = NULL);
void InitRandomTree(PhylogenyTreeBasic &treeToInit, int numTaxa,
                    int rndSeed = -1);
void CreatePhyTreeWithRootedSplits(PhylogenyTreeBasic &treeToProc, int numTaxa,
                                   const set<set<int> > &setGivenSplits);
void DumpAllSubtreesWithTaxaSize(
    const vector<PhylogenyTreeBasic *> &listPtrGTrees, int numTaxonSubtree,
    const char *fileNameOut);
void DumpAllSubtreesWithBoundedSize(
    const vector<PhylogenyTreeBasic *> &listPtrGTrees, int maxSzSubtree,
    int maxIdentSubtreeSz, const char *fileNameOut);
PhylogenyTreeBasic *ConsPhyTreeShrinkIdentSubtrees(PhylogenyTreeBasic *ptreeIn,
                                                   int maxIdentSubtreeSz,
                                                   bool fIdConsecutive = false);
void ChangebackLeafLabelForTreeWithZeroBaseId(PhylogenyTreeBasic *ptree,
                                              TaxaMapper *pTMapper);
void ChangeLeafIntLabelOfTree(PhylogenyTreeBasic &treeToChange,
                              const map<int, int> &mapOldIntLblToNewIntLbl,
                              bool fSetUserLblToo = false);
void AssignConsecutiveIdsForTree(PhylogenyTreeBasic &treeToChange);
bool ConvPhyloTreesToZeroBasedId(vector<PhylogenyTreeBasic *> &treePtrList,
                                 TaxaMapper *pTMapper);
void RandTrimLeavesFromTree(PhylogenyTreeBasic *ptreeToTrim,
                            int numLeavesRemain);
PhylogenyTreeBasic *ConsPhyTreeSubsetTaxa(PhylogenyTreeBasic *ptreeIn,
                                          const set<int> &setTaxaKept);
string ConsEdgeLabeTree(const string &strNWWithLabels);

#endif // PHYLOGENY_TREE_H
