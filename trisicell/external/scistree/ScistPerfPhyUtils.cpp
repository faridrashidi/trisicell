//
//  ScistPerfPhyUtils.cpp
//
//
//  Created by Yufeng Wu on 5/25/18.
//
//

#include "ScistPerfPhyUtils.hpp"
#include "PhylogenyTree.h"
#include "ScistGenotype.hpp"
#include "TreeBuilder.h"
#include "Utils3.h"
#include "Utils4.h"
#include "UtilsNumerical.h"
#include <iomanip>

// *************************************************************************************
// Cluster

void ScistPerfPhyClusterItor ::First() { it = clus.setMutSCs.begin(); }
void ScistPerfPhyClusterItor ::Next() { ++it; }
bool ScistPerfPhyClusterItor ::IsDone() { return it == clus.setMutSCs.end(); }
int ScistPerfPhyClusterItor ::GetCurrentSC() const { return *it; }

// *************************************************************************************
// Cluster

ScistPerfPhyCluster ::ScistPerfPhyCluster() {}

ScistPerfPhyCluster ::ScistPerfPhyCluster(const std::set<int> &clus)
    : setMutSCs(clus) {}

ScistPerfPhyCluster ::ScistPerfPhyCluster(const ScistPerfPhyCluster &rhs)
    : setMutSCs(rhs.setMutSCs) {}

ScistPerfPhyCluster &
ScistPerfPhyCluster ::operator=(const ScistPerfPhyCluster &rhs) {
  setMutSCs = rhs.setMutSCs;
  return *this;
}

bool ScistPerfPhyCluster ::operator<(const ScistPerfPhyCluster &rhs) const {
  return this->setMutSCs < rhs.setMutSCs;
}

void ScistPerfPhyCluster ::IntersectWith(
    const ScistPerfPhyCluster &rhs, ScistPerfPhyCluster &clusInt,
    ScistPerfPhyCluster &clusThisOnly, ScistPerfPhyCluster &clusRHSOnly) const {
  //
  JoinSets(this->setMutSCs, rhs.setMutSCs, clusInt.setMutSCs);
  clusThisOnly = *this;
  clusThisOnly.SubtractFrom(rhs);
  clusRHSOnly = rhs;
  clusRHSOnly.SubtractFrom(*this);
}

void ScistPerfPhyCluster ::SubtractFrom(const ScistPerfPhyCluster &rhs) {
  //
  SubtractSets(this->setMutSCs, rhs.setMutSCs);
}

void ScistPerfPhyCluster ::UnionWith(const ScistPerfPhyCluster &rhs) {
  UnionSets(this->setMutSCs, rhs.setMutSCs);
}

void ScistPerfPhyCluster ::GetGenoBinVec(int numHaps,
                                         vector<int> &vecGeno) const {
  vecGeno.clear();
  for (int i = 0; i < numHaps; ++i) {
    int g = 0;
    if (this->setMutSCs.find(i) != this->setMutSCs.end()) {
      g = 1;
    }
    vecGeno.push_back(g);
  }
}

bool ScistPerfPhyCluster ::IsCompatibleWith(
    const ScistPerfPhyCluster &rhs) const {
  // YW: assume rooted compatibility (i.e. three gamates test)
  ScistPerfPhyCluster clusInt, clusThisOnly, clusRHSOnly;
  IntersectWith(rhs, clusInt, clusThisOnly, clusRHSOnly);
  return clusInt.GetSize() == 0 || clusThisOnly.GetSize() == 0 ||
         clusRHSOnly.GetSize() == 0;
}

bool ScistPerfPhyCluster ::IsCompatibleWith(
    const std::set<ScistPerfPhyCluster> &setClus) const {
  for (set<ScistPerfPhyCluster>::const_iterator it = setClus.begin();
       it != setClus.end(); ++it) {
    if (IsCompatibleWith(*it) == false) {
      return false;
    }
  }
  return true;
}

void ScistPerfPhyCluster ::GetSplitPartsWith(
    const ScistPerfPhyCluster &rhs,
    std::vector<std::set<int> > &listParts) const {
  // get 10, 01, and 11
  ScistPerfPhyCluster clusInt, clusThisOnly, clusRHSOnly;
  IntersectWith(rhs, clusInt, clusThisOnly, clusRHSOnly);
  //
  listParts.push_back(clusThisOnly.setMutSCs);
  listParts.push_back(clusRHSOnly.setMutSCs);
  listParts.push_back(clusInt.setMutSCs);
}

void ScistPerfPhyCluster ::FlipAlleleAt(int r) {
  // if row r is in the cluster, remove it; otherwise add it
  if (setMutSCs.find(r) != setMutSCs.end()) {
    setMutSCs.erase(r);
  } else {
    setMutSCs.insert(r);
  }
}

int ScistPerfPhyCluster ::GetAlleleAt(int r) const {
  if (setMutSCs.find(r) != setMutSCs.end()) {
    return 1;
  } else {
    return 0;
  }
}

void ScistPerfPhyCluster ::Dump() const { DumpIntSet(setMutSCs); }

// *************************************************************************************
// Cluster partial order tree node

ScistPerfPhyClusTreeNode ::~ScistPerfPhyClusTreeNode() {}

ScistPerfPhyClusTreeNode *ScistPerfPhyClusTreeNode ::ConsClusterTree(
    const std::map<int, ScistPerfPhyCluster> &setSeedSites, bool fNoDup) {
  // the root has no clus attached (i.e. contains everything)
  ScistPerfPhyClusTreeNode *pTreeRoot = new ScistPerfPhyClusTreeNode(NULL);
  set<ScistPerfPhyCluster> setClusDone;

  for (map<int, ScistPerfPhyCluster>::const_iterator it = setSeedSites.begin();
       it != setSeedSites.end(); ++it) {
    if (fNoDup) {
      if (setClusDone.find(it->second) != setClusDone.end()) {
        continue;
      }
    }

    // cout << "Init cluster tree: node: ";
    // it->second.Dump();
    ScistPerfPhyClusTreeNode *pNode =
        new ScistPerfPhyClusTreeNode(&(it->second));
    pTreeRoot->InsertNode(pNode);

    setClusDone.insert(it->second);
  }
  return pTreeRoot;
}
ScistPerfPhyClusTreeNode *ScistPerfPhyClusTreeNode ::ConsClusterTree(
    const std::set<ScistPerfPhyCluster> &setSeedSites) {
  // the root has no clus attached (i.e. contains everything)
  ScistPerfPhyClusTreeNode *pTreeRoot = new ScistPerfPhyClusTreeNode(NULL);

  for (set<ScistPerfPhyCluster>::const_iterator it = setSeedSites.begin();
       it != setSeedSites.end(); ++it) {
    // cout << "Init cluster tree: node: ";
    // it->second.Dump();
    ScistPerfPhyClusTreeNode *pNode = new ScistPerfPhyClusTreeNode(&(*it));
    pTreeRoot->InsertNode(pNode);
  }
  return pTreeRoot;
}

void ScistPerfPhyClusTreeNode ::AddChild(ScistPerfPhyClusTreeNode *pChild) {
  //
  listChildren.push_back(pChild);
  pChild->SetParent(this);
}

void ScistPerfPhyClusTreeNode ::RemoveChild(ScistPerfPhyClusTreeNode *pChild) {
  //
  pChild->SetParent(NULL);
  listChildren.erase(
      std::remove(listChildren.begin(), listChildren.end(), pChild),
      listChildren.end());
}

void ScistPerfPhyClusTreeNode ::InsertNode(ScistPerfPhyClusTreeNode *pNode) {
  // cout << "Insert node: ";
  // pNode->Dump();
  // cout << " under parent node: ";
  // Dump();
  // insert this node below it; may need to split children if there are multiple
  // ones we assume there is no incompatibility occuring here
  vector<ScistPerfPhyClusTreeNode *> listChildrenContained;
  for (int i = 0; i < GetNumChildren(); ++i) {
    //
    ScistPerfPhyCluster clusInt, clusThisOnly, clusRHSOnly;
    pNode->GetClus()->IntersectWith(*GetChild(i)->GetClus(), clusInt,
                                    clusThisOnly, clusRHSOnly);

    // test if contained by one subtree; if so, add it to it inteaad
    bool fContained = (clusThisOnly.GetSize() == 0);
    if (fContained) {
      GetChild(i)->InsertNode(pNode);
      return;
    }
    bool fContaining = (clusRHSOnly.GetSize() == 0);
    bool fDisjoint = (clusInt.GetSize() == 0);
    if (fContaining) {
      listChildrenContained.push_back(GetChild(i));
    } else {
      YW_ASSERT_INFO(fDisjoint == true,
                     "Wrong: the site is not compatible with the tree");
    }
  }
  // if( listChildrenContained.size() == 0 )
  //{
  // just add it below
  //    AddChild(pNode);
  //}
  // else
  //{
  //
  for (int i = 0; i < (int)listChildrenContained.size(); ++i) {
    RemoveChild(listChildrenContained[i]);
    pNode->AddChild(listChildrenContained[i]);
  }
  AddChild(pNode);
  //}
}

void ScistPerfPhyClusTreeNode ::Dump() const {
  cout << "Node: "
       << ", num of children: " << GetNumChildren() << ": ";
  if (GetClus() == NULL) {
    cout << "root. \n";
  } else {
    GetClus()->Dump();
  }
}

// *************************************************************************************
// Guide tree

ScistPerfPhyGuideTree ::ScistPerfPhyGuideTree() {}

void ScistPerfPhyGuideTree ::Init(const std::string &strGuideTree) {
  this->setGuideTreeClus.clear();
  // cout << "INIT guide tree...\n";
  // extract clusters in the tree
  PhylogenyTreeBasic treeGuide;
  treeGuide.ConsOnNewick(strGuideTree);

  // get all clusters that don't have zero length (i.e. not non-informative)
  PhylogenyTreeIterator itorTree(treeGuide);
  itorTree.Init();
  while (itorTree.IsDone() == false) {
    TreeNode *pn = itorTree.GetCurrNode();
    if (pn->IsLeaf() == false && pn->IsRoot() == false)
    // if( pn->GetLength() >= MIN_POS_VAL && pn->IsLeaf() == false &&
    // pn->IsRoot() == false)
    {
      set<int> ss;
      pn->GetAllDescendIntLbls(ss);
      DecAllNumInSet(ss);
      ScistPerfPhyCluster clus(ss);
      this->setGuideTreeClus.insert(clus);
      // cout << "guide tree cluster: ";
      // clus.Dump();
    }

    itorTree.Next();
  }

  // set<set<int> > setClades;
  // treeGuide.GetAllClades(setClades);
  // for(set<set<int> > :: iterator it = setClades.begin(); it !=
  // setClades.end(); ++it)
  //{
  //    set<int> ss=*it;
  //    DecAllNumInSet(ss);
  //    ScistPerfPhyCluster clus(ss);
  //    this->setGuideTreeClus.insert(clus);
  // cout << "guide tree cluster: ";
  // clus.Dump();
  //}
}

void ScistPerfPhyGuideTree ::InitDecAll(const std::string &strGuideTree1Base) {
  //
  this->setGuideTreeClus.clear();
  // cout << "INIT guide tree...\n";
  // extract clusters in the tree
  PhylogenyTreeBasic treeGuide;
  treeGuide.ConsOnNewick(strGuideTree1Base);

  // dec by one
  // map<int,int> mapOldToNew;
  // for(int i=0; i<treeGuide.GetNumLeaves(); ++i)
  //{
  //    mapOldToNew[i+1] = i;
  //}
  // ChangeLeafIntLabelOfTree(treeGuide, mapOldToNew, true);

  // get all clusters that don't have zero length (i.e. not non-informative)
  PhylogenyTreeIterator itorTree(treeGuide);
  itorTree.Init();
  while (itorTree.IsDone() == false) {
    TreeNode *pn = itorTree.GetCurrNode();
    if (pn->IsLeaf() == false && pn->IsRoot() == false) {
      set<int> ss;
      pn->GetAllDescendIntLbls(ss);
      DecAllNumInSet(ss);
      ScistPerfPhyCluster clus(ss);
      this->setGuideTreeClus.insert(clus);
      // cout << "guide tree cluster: ";
      // clus.Dump();
    }

    itorTree.Next();
  }
}

double ScistPerfPhyGuideTree ::EvalClus(const ScistPerfPhyCluster &clus) const {
  // find the best match
  double res = 0.0;
  if (setGuideTreeClus.size() == 0) {
    return res;
  }
  for (set<ScistPerfPhyCluster>::const_iterator it = setGuideTreeClus.begin();
       it != setGuideTreeClus.end(); ++it) {
    //
    int score = EvalClusWith(clus, *it);
    res += score;
    // if( res < score )
    //{
    //    res = score;
    //}
  }
  // return res;
  return res / setGuideTreeClus.size();
}

int ScistPerfPhyGuideTree ::EvalClusWith(
    const ScistPerfPhyCluster &clus, const ScistPerfPhyCluster &clusInTree) {
  // score (dissimlarity): high means the better fit. Score: percentage of
  // smallest diff; use Jaccard distance that is, size of intersection over size
  // of union YW: try compat (1.0) and incompat (0.0) and take average
  int res = 1;
  if (clus.IsCompatibleWith(clusInTree) == true) {
    res = 0;
  }
  return res;

  // ScistPerfPhyCluster clusUnion = clus;
  // clusUnion.UnionWith(clusInTree);
  // ScistPerfPhyCluster clusInt, clus1, clus2;
  // clus.IntersectWith(clusInTree, clusInt, clus1, clus2);
  // return ((double) clusInt.GetSize() )/ clusUnion.GetSize();
}

// *************************************************************************************
// Inf perfect phylogeny from genotypes

ScistInfPerfPhyUtils ::ScistInfPerfPhyUtils() {}

ScistInfPerfPhyUtils ::~ScistInfPerfPhyUtils() {}

std::string ScistInfPerfPhyUtils ::ConsTreeWCombDistClus(
    const ScistGenGenotypeMat &genos,
    const std::map<int, ScistPerfPhyCluster> &setClus,
    bool fUseGenoName) const {
  // not only construct a tree, but also consider the distance
  set<set<int> > setClustersMustHave;
  for (map<int, ScistPerfPhyCluster>::const_iterator it = setClus.begin();
       it != setClus.end(); ++it) {
    set<int> setOnes;
    it->second.GetClus(setOnes);
    setClustersMustHave.insert(setOnes);
  }

#if 0
    // add all NJ clusters when compatible
    std::set< ScistPerfPhyCluster > clusAll;
    this->treeGuide.GetAllClusters( clusAll );
    //string strGuideTree = genos.ConsNJTree();
//cout << "Neighbor joining tree from corrected genotypes: " << strGuideTree << endl;
    //PhylogenyTreeBasic treeGuide;
    //treeGuide.ConsOnNewick(strGuideTree);
    //set<set<int> > setClades;
    //treeGuide.GetAllClades(setClades);
    //for(set<set<int> > :: iterator it = setClades.begin(); it != setClades.end(); ++it)
    for(set<ScistPerfPhyCluster> :: iterator it = clusAll.begin(); it != clusAll.end(); ++it)
    {
        //set<int> ss=*it;
        set<int> ss;
        it->GetClus(ss);
        //DecAllNumInSet(ss);
        //ScistPerfPhyCluster clus(ss);
        ScistPerfPhyCluster clus = *it;
        if( clus.GetSize() <=1 )
        {
            continue;
        }
        // make sure compatible with exisitng clusters
        bool fCompat = true;
        for( map<int, ScistPerfPhyCluster> :: const_iterator it2 = setClus.begin(); it2 != setClus.end(); ++it2 )
        {
            if( clus.IsCompatibleWith(it2->second) == false )
            {
                fCompat = false;
                break;
            }
        }
        if( fCompat)
        {
//cout << "ADDING guide tree cluster: ";
//clus.Dump();
            setClustersMustHave.insert( ss );
        }

    }
#endif

  //
  PhyloDistance phyDist;
#if 0
    BinaryMatrix matClus;
    matClus.SetSize( genosInput.GetNumHaps(), genosInput.GetNumSites() );
    //
    for( map<int, ScistPerfPhyCluster> :: const_iterator it = setClus.begin(); it != setClus.end(); ++it )
    {
        for(int i=0; i<genosInput.GetNumHaps(); ++i)
        {
            matClus.SetValAt( i, it->first, 0 );
        }
        set<int> setOnes;
        it->second.GetClus(setOnes);
        for( set<int> :: iterator it2 = setOnes.begin(); it2 != setOnes.end(); ++it2 )
        {
            matClus.SetValAt(*it2, it->first, 1);
        }
    }
//cout << "matClus: ";
//matClus.Dump();
#endif
  for (int r1 = 0; r1 < genos.GetNumHaps(); ++r1) {
    phyDist.SetDistance(r1, r1, 0.0);
    for (int r2 = r1 + 1; r2 < genos.GetNumHaps(); ++r2) {
      // set<int> setDiffs;
      // matClus.GetSequencesDiffSites(r1,r2, setDiffs);
      // double dist = ((double)setDiffs.size())/genosInput.GetNumSites();
      double dist = genos.CalcHammingDistBetwHaps(r1, r2);
      phyDist.SetDistance(r1, r2, dist);
      phyDist.SetDistance(r2, r1, dist);
    }
  }

  // build tree
  set<set<int> > setClustersForbiddenEmpty;
  ConstrainedUPGMATreeBuilder treeBuilder(phyDist, setClustersMustHave,
                                          setClustersForbiddenEmpty);
  while (treeBuilder.IsDone() == false) {
    set<int> st1, st2;
    double minDist = treeBuilder.GetMinCoalSubtrees(st1, st2);
    // cout << "Merging subtrees ht " << minDist << ": ";
    // DumpIntSet(st1);
    // DumpIntSet(st2);
    treeBuilder.MergeSubtrees(st1, st2, minDist);
  }
  string strTreeRaw = treeBuilder.GetTree();
  // cout << "Constructed tree with clusters and distance: " << strTreeRaw <<
  // endl;

  // convert labels
  PhylogenyTreeBasic phTree;
  phTree.ConsOnNewick(strTreeRaw);

  map<string, string> mapIdToLabels;
  for (int i = 0; i < genos.GetNumHaps(); ++i) {
    // cout << "i: " << i << ", name: " << this->genosInput.GetGenotypeName(i)
    // << endl; string str = "(" + std::to_string(i) + ")";
    string str = std::to_string(i);
    if (fUseGenoName) {
      mapIdToLabels[str] = genos.GetGenotypeName(i);
    } else {
      mapIdToLabels[str] = std::to_string(i + 1);
    }
  }
  phTree.ReassignLeafLabels(mapIdToLabels);
  // use base-1 label
  phTree.IncEdgeLabelsBy(1);

  string res;
  phTree.ConsNewickSorted(res);
  // phTree.ConsNewick(res, false, 0.0, true);

  // output a tree in GML format
  // if( this->fOutput )
  //{
  //    string fileNameOut =  genosInput.GetFileName() + ".tree.gml";
  //    phTree.OutputGML(fileNameOut.c_str());
  //}

  return res;
}

void ScistInfPerfPhyUtils ::FillClusterFromMat(const ScistGenGenotypeMat &genos,
                                               int site,
                                               ScistPerfPhyCluster &clus) {
  //
  for (int r = 0; r < genos.GetNumHaps(); ++r) {
    int v = genos.GetGenotypeAt(r, site);
    if (v != 0) {
      clus.AddMutSC(r);
    }
  }
}

// *************************************************************************************

void ScistInfPerfPhyTest() {
  ScistHaplotypeMat genoMat;
  const int numSCs = 4, numSites = 3;
  genoMat.SetSize(numSCs, numSites);
  genoMat.SetGenotypeAt(0, 0, 1);
  genoMat.SetGenotypeAt(0, 1, 0);
  genoMat.SetGenotypeAt(0, 2, 1);
  genoMat.SetGenotypeAt(1, 0, 1);
  genoMat.SetGenotypeAt(1, 1, 1);
  genoMat.SetGenotypeAt(1, 2, 0);
  genoMat.SetGenotypeAt(2, 0, 0);
  genoMat.SetGenotypeAt(2, 1, 1);
  genoMat.SetGenotypeAt(2, 2, 1);
  genoMat.SetGenotypeAt(3, 0, 0);
  genoMat.SetGenotypeAt(3, 1, 1);
  genoMat.SetGenotypeAt(3, 2, 0);

  // ScistInfPerfPhy ppInf( genoMat );
  // ppInf.InferGreedy();
}
