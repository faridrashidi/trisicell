//
//  TreeBuilder.h
//
//
//  Created by Yufeng Wu on 9/11/15.
//
//

#ifndef ____TreeBuilder__
#define ____TreeBuilder__

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
using namespace std;

//***********************************************************************
void TestNJ();

//***********************************************************************
// implement various methods to build a phylogenetic tree

// define distances between taxa
class PhyloDistance {
public:
  void SetDistance(int node1, int node2, double dist);
  double GetDistance(int node1, int node2) const;
  void GetAllNodes(set<int> &nodesAll) const;
  double GetDistanceNonNeg(int node1, int node2) const;
  double CalcAveDistBtwClusters(const set<set<int> > &setClusters) const;
  void Dump() const;

private:
  map<pair<int, int>, double> mapDists;
};

// distance based tree builder
class DistanceTreeBuilder {
public:
  DistanceTreeBuilder(PhyloDistance &distPairwiseTaxaIn);
  string NJ();
  string ConstrainedUPGMA(const set<set<int> > &setClustersMustHave,
                          const set<set<int> > &setClustersForbidden,
                          map<set<int>, double> &mapSTHts, int numTotElem = -1);
  string ConstrainedUPGMA(const set<set<int> > &setClustersMustHave,
                          const set<set<int> > &setClustersDesired,
                          int numTopCandidates,
                          const set<set<int> > &setClustersForbidden,
                          map<set<int>, double> &mapSTHts, int numTotElem);
  void SetTaxonName(int id, const string &tname) { mapIndexToName[id] = tname; }
  void SetOutgroup(int og) { taxonOutgroup = og; }

private:
  void NJFindNgbrs(int nodeIdNew, set<int> &nodesToSearch, int &ngbr1,
                   int &ngbr2);
  double NJCalcAveDist(int nodecur, const set<int> &nodesToSearch);
  bool IsClusterIncompatible(const set<int> &clus1, const set<int> &clus2,
                             int numTotElem = -1) const;
  bool IsClusterIncompatibleWithSetofClus(const set<int> &clus1,
                                          const set<set<int> > &setClus,
                                          int numTotElem = -1) const;
  void UpdateDistUPGMA(const pair<set<int>, set<int> > &pairClus,
                       const map<set<int>, pair<string, double> > &mapSubtree,
                       map<pair<set<int>, set<int> >, double> &distMapCur);
  string GetTaxonNameFor(int index) const;
  int GetNumCompatCladesIn(const set<int> &clus1,
                           const set<set<int> > &setCladesTest,
                           int numTotElem) const;

  PhyloDistance &distPairwiseTaxa;
  map<int, string> mapIndexToName;
  int taxonOutgroup;
};

//***********************************************************************
// tool for building UPGMA tree

class ConstrainedUPGMATreeBuilder {
public:
  ConstrainedUPGMATreeBuilder(PhyloDistance &distPairwiseTaxaIn,
                              const set<set<int> > &setClustersMustHave,
                              const set<set<int> > &setClustersForbidden,
                              int numTotElemIn = -1);
  ConstrainedUPGMATreeBuilder(const ConstrainedUPGMATreeBuilder &rhs);
  string GetTree() const;
  string GetPartialConsTree() const;
  double GetMinCoalSubtrees(set<int> &st1, set<int> &st2) const;
  void GetCoalSubtreesHtBound(
      double htBound,
      set<pair<pair<set<int>, set<int> >, double> > &setCandidates) const;
  void MergeSubtrees(const set<int> &st1, const set<int> &st2,
                     double htMergedST);
  void GetMergeCandidates(
      map<pair<set<int>, set<int> >, double> &setCandidates) const;
  double GetCurDistForTwoClusters(const set<int> &clus1,
                                  const set<int> &clus2) const;
  void SetDistForTwoClusters(const set<int> &clus1, const set<int> &clus2,
                             double dist);
  int GetNumSubtrees() const;
  void GetAllSubtrees(map<set<int>, string> &mapSTs) const;
  void GetActiveSubtrees(set<set<int> > &setActiveSTs) const;
  bool IsDone() const;
  void Dump() const;

private:
  void Init();
  bool IsClusterIncompatible(const set<int> &clus1,
                             const set<int> &clus2) const;
  bool IsClusterIncompatibleWithSetofClus(const set<int> &clus1,
                                          const set<set<int> > &setClus) const;
  void UpdateDistUPGMA(const set<int> &st1, const set<int> &st2);

  PhyloDistance &distPairwiseTaxa;
  const set<set<int> > &setClustersMustHave;
  const set<set<int> > &setClustersForbidden;
  int numTotElem;
  map<pair<set<int>, set<int> >, double> distMapActivePair;
  map<set<int>, pair<string, double> > mapClusSubtree;
  vector<pair<set<int>, set<int> > > histSTMerge;
};

//***********************************************************************
// tool for building near-optimal UPGMA tree

class ConstrainedNearUPGMATreesBuilder {
public:
  ConstrainedNearUPGMATreesBuilder(PhyloDistance &distPairwiseTaxaIn,
                                   const set<set<int> > &setClustersMustHave,
                                   const set<set<int> > &setClustersForbidden,
                                   int numTotElem);
  void Construct(int maxNumTrees, double thresMaxDistRatio);
  void GetTrees(set<string> &setConsTrees) const { setConsTrees = setTreeCons; }

private:
  PhyloDistance &distPairwiseTaxa;
  const set<set<int> > &setClustersMustHave;
  const set<set<int> > &setClustersForbidden;
  int numTotElem;
  set<string> setTreeCons;
};

#endif /* defined(____TreeBuilder__) */
