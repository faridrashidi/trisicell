//
//  ScistPerfPhyUtils.hpp
//
//
//  Created by Yufeng Wu on 5/25/18.
//
//

#ifndef ScistPerfPhyUtils_hpp
#define ScistPerfPhyUtils_hpp

#include <map>
#include <set>
#include <vector>

class ScistGenGenotypeMat;
class PhylogenyTree;

// *************************************************************************************
// Cluster

class ScistPerfPhyCluster;

class ScistPerfPhyClusterItor {
public:
  ScistPerfPhyClusterItor(const ScistPerfPhyCluster &clusIn) : clus(clusIn) {
    First();
  }
  void First();
  void Next();
  bool IsDone();
  int GetCurrentSC() const;

private:
  const ScistPerfPhyCluster &clus;
  std::set<int>::const_iterator it;
};

class ScistPerfPhyCluster {
  friend class ScistPerfPhyClusterItor;

public:
  ScistPerfPhyCluster();
  ScistPerfPhyCluster(const std::set<int> &clus);
  ScistPerfPhyCluster(const ScistPerfPhyCluster &rhs);
  ScistPerfPhyCluster &operator=(const ScistPerfPhyCluster &rhs);

  bool operator<(const ScistPerfPhyCluster &rhs) const;
  int GetSize() const { return setMutSCs.size(); }
  void IntersectWith(const ScistPerfPhyCluster &rhs,
                     ScistPerfPhyCluster &clusInt,
                     ScistPerfPhyCluster &clusThisOnly,
                     ScistPerfPhyCluster &clusRHSOnly) const;
  void SubtractFrom(const ScistPerfPhyCluster &rhs);
  void UnionWith(const ScistPerfPhyCluster &rhs);
  void Clear() { setMutSCs.clear(); }
  void GetGenoBinVec(int numSCs, std::vector<int> &vecGeno) const;
  void GetClus(std::set<int> &clus) const { clus = setMutSCs; }
  bool IsCompatibleWith(const ScistPerfPhyCluster &rhs) const;
  bool IsCompatibleWith(const std::set<ScistPerfPhyCluster> &setClus) const;
  void GetSplitPartsWith(const ScistPerfPhyCluster &rhs,
                         std::vector<std::set<int> > &listParts) const;
  void AddMutSC(int r) { setMutSCs.insert(r); }
  void FlipAlleleAt(int r);
  int GetAlleleAt(int r) const;
  void Dump() const;

private:
  std::set<int> setMutSCs;
};

// *************************************************************************************
// Cluster partial order tree node

class ScistPerfPhyClusTreeNode {
public:
  ScistPerfPhyClusTreeNode(const ScistPerfPhyCluster *pClusIn)
      : pClus(pClusIn), pParent(NULL) {}
  ~ScistPerfPhyClusTreeNode();
  static ScistPerfPhyClusTreeNode *
  ConsClusterTree(const std::map<int, ScistPerfPhyCluster> &setSeedSites,
                  bool fNoDup = false);
  static ScistPerfPhyClusTreeNode *
  ConsClusterTree(const std::set<ScistPerfPhyCluster> &setSeedSites);
  void SetParent(ScistPerfPhyClusTreeNode *pParentIn) { pParent = pParentIn; }
  ScistPerfPhyClusTreeNode *GetParent() { return pParent; }
  int GetNumChildren() const { return listChildren.size(); }
  ScistPerfPhyClusTreeNode *GetChild(int i) const { return listChildren[i]; }
  const ScistPerfPhyCluster *GetClus() const { return pClus; }
  void AddChild(ScistPerfPhyClusTreeNode *pChild);
  void RemoveChild(ScistPerfPhyClusTreeNode *pChild);
  void InsertNode(ScistPerfPhyClusTreeNode *pNode);
  bool IsRoot() const { return pParent == NULL; }
  bool IsLeaf() const { return GetNumChildren() == 0; }
  void Dump() const;

private:
  const ScistPerfPhyCluster *pClus;
  ScistPerfPhyClusTreeNode *pParent;
  std::vector<ScistPerfPhyClusTreeNode *> listChildren;
};

// *************************************************************************************
// Guide tree

class ScistPerfPhyGuideTree {
public:
  ScistPerfPhyGuideTree();
  void Init(const std::string &strGuideTree);
  void InitDecAll(const std::string &strGuideTree1Base);
  double EvalClus(const ScistPerfPhyCluster &clus) const;
  void GetAllClusters(std::set<ScistPerfPhyCluster> &clusAll) const {
    clusAll = this->setGuideTreeClus;
  }

private:
  static int EvalClusWith(const ScistPerfPhyCluster &clus,
                          const ScistPerfPhyCluster &clusInTree);

  std::set<ScistPerfPhyCluster> setGuideTreeClus;
};

// *************************************************************************************
// Inf perfect phylogeny from genotypes

class ScistInfPerfPhyUtils {
public:
  ScistInfPerfPhyUtils();
  ~ScistInfPerfPhyUtils();
  static void FillClusterFromMat(const ScistGenGenotypeMat &genos, int site,
                                 ScistPerfPhyCluster &clus);
  std::string
  ConsTreeWCombDistClus(const ScistGenGenotypeMat &genos,
                        const std::map<int, ScistPerfPhyCluster> &setClus,
                        bool fUseGenoName = true) const;

private:
  void ClearClusTree();

  ScistPerfPhyClusTreeNode *pClusTreeRoot;
};

// *************************************************************************************

void ScistInfPerfPhyTest();

#endif /* ScistPerfPhy_hpp */
