#ifndef MARGINAL_TREE_H
#define MARGINAL_TREE_H

#include <stack>
#include <vector>
using namespace std;

#include "Utils.h"
#include "Utils2.h"
#include "Utils3.h"

//////////////////////////////////////////////////////////////////////////////
// Define a simple coalescent tree. My experience shows that such a data
// structure can be quite useful

// yet another structure to represent marginal tree

class MarginalTree {
public:
  MarginalTree();
  void Clear();
  void Binarize();
  void Consolidate();
  void BuildDescendantInfo();
  void InitDefaultEdgeLen();
  void InitUnitEdgelen();
  double GetDefaultEdgeLen(int child);
  void SetParent(int child, int par, bool fAdjLen = true);
  int GetParent(int child) const;
  int GetLeftDescendant(int node) const;
  int GetRightDescendant(int node) const;
  double GetEdgeLen(int childNodeIndex) const;
  double GetTotEdgeLen() const;
  int GetTotNodesNum() const { return listNodeLabels.size(); }
  int GetNumLeaves() const { return numLeaves; }
  void SetNumLeaves(int nl) { numLeaves = nl; }
  void ConsDecedentInfo(vector<vector<int> > &descNodes) const;
  void ConsAllDecedentInfo(vector<set<int> > &descNodes,
                           bool fIncSelf = true) const;
  void ConsDecedentLeavesInfo(vector<set<int> > &descNodes) const;
  void ConsDecedentLeavesInfoLabels(vector<set<int> > &leafNodeLabels) const;
  void ConsHeightsInfo(vector<int> &nodesHt) const;
  void Dump() const;
  int GetLabel(int r) const {
    YW_ASSERT_INFO(r >= 0 && r < (int)listNodeLabels.size(), "wrong3");
    return listNodeLabels[r];
  }
  void SetLabel(int node, int lbl) {
    YW_ASSERT_INFO(node >= 0 && node < (int)listNodeLabels.size(), "wrong4");
    listNodeLabels[node] = lbl;
  }
  int GetPosForLabel(int lbl) const;
  void GetlabelsFor(const set<int> &setPos, set<int> &setLbls) const;
  bool IsLeaf(int node) const { return node >= 0 && node < numLeaves; }
  bool IsToplogicSame(const MarginalTree &tree) const;
  int GetMRCA(int v1, int v2) const;
  int GetFirstNonselfAnces(int v, const set<int> &setAnces) const;
  void GetChildren(int node, set<int> &listChildren) const;
  int GetMaxHt() const;
  void RemoveLeafNodeFromBinaryTree(int lfn);
  bool AreTwoPathsDisjoint(int sn1, int en1, int sn2, int en2) const;
  int GetPath(int sn, int en, set<int> &edgesOnPath) const;
  double GetPathLen(int sn, int en);
  void OutputGML(const char *fileName) const;
  string GetNewick() const;
  string GetNewickSorted(bool fLen) const;
  string GetNewickNoBrLen() const;
  string GetNewickNoBrLen2() const;
  void GetLeavesUnder(int nn, set<int> &leavesUnder) const;
  void GetLeafSetsForCuts(const vector<int> &listCuts,
                          vector<set<int> > &listLeafSets) const;
  int GetMRCAForNodes(const set<int> &listNodes) const;
  bool IsNodeUnder(int nn, int ancesNode) const;
  void RandPermuateLeaves();
  int GetTriple(int a, int b, int c) const;
  int GetSibling(int a) const;
  bool AreNodesSibling(int a, int b) const;
  void SetBranchLen(int b, double len) {
    YW_ASSERT_INFO(b < (int)listEdgeDist.size(), "Branch wrong");
    listEdgeDist[b] = len;
  }
  void SetLabelList(const vector<int> &listLbls) { listNodeLabels = listLbls; }
  void GetLabelList(vector<int> &listLbls) const { listLbls = listNodeLabels; }
  void GetLabelListForLeaf(vector<int> &listLbls) const;
  void SetParList(const vector<int> &listPars) { listParentNodePos = listPars; }
  void SetBranchLenList(const vector<double> &listLens) {
    listEdgeDist = listLens;
  }
  void SortByLeafId();
  void FixDupIds();
  double GetHeight() const;
  int GetRoot() const { return GetTotNodesNum() - 1; }
  void SwapBranches(int nodeBranch1, int nodeBranch2);
  void RearrangeParIncOrder();
  void ResetIncLabel();
  void IncLabels();
  void GetTreeEdgeLen(vector<double> &listEdgeDistOut) const {
    listEdgeDistOut = this->listEdgeDist;
  }
  void MapLeafLblConsecutiveOrder(vector<int> &listLeafLblsOld);
  void RemapLeafLabels(const map<int, int> &mapLeafLblsToNew);
  void FindAllSplits(vector<set<int> > &listSplits) const;
  void FindSibLeafPairs(vector<pair<int, int> > &listSibPairs) const;
  void MakeLeafSubtreeOfTwo(int posLeaf, int lblChild1, int lblChild2,
                            double len1, double len2);
  void FindDiffSubtreesFrom(const MarginalTree &mtreeRef, set<int> &setDiffBrs,
                            set<int> &setDiffBrsOrigOnly) const;
  bool IsOutgroup(int lvid) const;

public:
  int CalcNormHeight(int node);
  void GetParPosInfo(vector<int> &parPosList) {
    parPosList = listParentNodePos;
  }
  double GetHeightOfNode(int node) const;

  // Use an array to store  leaves
  int numLeaves;
  // assume the first numLeaves nodes are leaves
  vector<int> listNodeLabels;
  vector<int> listParentNodePos;
  vector<double> listEdgeDist;
  vector<int> listLeftDescs;
  vector<int> listRightDescs;

private:
  string GetNewickAt(int node, bool fSort = false, bool fLen = true) const;
  void MapLeafLblConsecutiveOrderAt(int rootST, int &idNext,
                                    vector<int> &listLeafLblsOld);
};

////////////////////////////////////////////////////////////////////////////////////
// Global Utilities
class TaxaMapper;

bool ReadinMarginalTrees(ifstream &inFile, vector<MarginalTree> &treeList);
bool ReadinMarginalTreesNewick(ifstream &inFile, int numLeaves,
                               vector<MarginalTree> &treeList,
                               TaxaMapper *pTMapper = NULL, bool fDup = false);
bool ReadinMarginalTreesNewickWLen(ifstream &inFile, int numLeaves,
                                   vector<MarginalTree> &treeList,
                                   TaxaMapper *pTMapper = NULL);
void AddRootAsLeafToTree(MarginalTree &tree1, bool fIdNonNeg = false);
void GenRandBinaryTree(int numLeaves, MarginalTree &tree1);
void GenRandBinaryTreeClock(int numLeaves, double totHt, MarginalTree &tree1);
// vector<int>: list of leaves in the order from top down, int = top node of
// chain
void FindChainsInTree(const MarginalTree &tree1,
                      map<vector<int>, int> &foundChains);
void InitMarginalTree(MarginalTree &mTree, int numLeaves,
                      const vector<int> &listLabels,
                      const vector<int> &listParentNodePos);
bool ReadinMarginalTreesNewickWLenString(const string &strNewick, int numLeaves,
                                         MarginalTree &treeOut,
                                         bool fStartFromZero = true,
                                         TaxaMapper *pTMapper = NULL);
void CollapseEquivTrees(const vector<MarginalTree> &listOrigTrees,
                        vector<MarginalTree> &listUniqTrees,
                        vector<int> &listMultiplicity);
void FindOneNNIMTreesFrom(
    MarginalTree &mTreeSrc, vector<MarginalTree> &listNNITrees,
    vector<pair<int, int> > *pListPairEdgesSwapped = NULL);
void CreateSubtreeFromLeaves(MarginalTree &mTreeOrig,
                             const set<int> &setLeafLabels,
                             MarginalTree &mTreeSub,
                             map<int, int> &mapNewNodeToOldNode);
void UpdateBranchLenInSubtree(MarginalTree &mTreeOrig,
                              map<int, int> &mapNewNodeToOldNode,
                              MarginalTree &mTreeSub);
void RemapLeafIntLabelsTaxaMap(MarginalTree &mtree,
                               map<string, string> &mapper);
void RemapMargTree(MarginalTree &mtree, TaxaMapper &refTMapper);
void FindMatchedSubtrees(MarginalTree &mtreeNew, MarginalTree &mtreeRef,
                         map<int, int> &mapSTNewToRef);

#endif // MARGINAL_TREE_H
