//
//  TreeBuilder.cpp
//
//
//  Created by Yufeng Wu on 9/11/15.
//
//

#include "TreeBuilder.h"
#include "Utils.h"
#include <string>

//***********************************************************************
void TestNJ() {
  //
  PhyloDistance distNJ;
  distNJ.SetDistance(1, 2, 5.0);
  distNJ.SetDistance(1, 3, 9.0);
  distNJ.SetDistance(1, 4, 9.0);
  distNJ.SetDistance(1, 5, 8.0);
  distNJ.SetDistance(2, 3, 10.0);
  distNJ.SetDistance(2, 4, 10.0);
  distNJ.SetDistance(2, 5, 9.0);
  distNJ.SetDistance(3, 4, 8.0);
  distNJ.SetDistance(3, 5, 7.0);
  distNJ.SetDistance(4, 5, 3.0);

  DistanceTreeBuilder builder(distNJ);
  string treeNW = builder.NJ();
  cout << "Constructed NJ tree: " << treeNW << endl;
  // distNJ.Dump();
}

//***********************************************************************
// define distances between taxa

void PhyloDistance ::SetDistance(int node1, int node2, double dist) {
  //
  pair<int, int> pp(node1, node2);
  mapDists.insert(map<pair<int, int>, double>::value_type(pp, dist));
}

double PhyloDistance ::GetDistance(int node1, int node2) const {
  //
  PhyloDistance *pthis = const_cast<PhyloDistance *>(this);
  pair<int, int> pp1(node1, node2), pp2(node2, node1);
  if (mapDists.find(pp1) != mapDists.end()) {
    //
    return pthis->mapDists[pp1];
  }
  if (mapDists.find(pp2) != mapDists.end()) {
    //
    return pthis->mapDists[pp2];
  }
  YW_ASSERT_INFO(false, "Fail to find");
  return 0.0;
}

double PhyloDistance ::GetDistanceNonNeg(int node1, int node2) const {
  //
  double dist = GetDistance(node1, node2);
  if (dist < 0.0) {
    dist = 0.0;
  }
  return dist;
}

void PhyloDistance ::GetAllNodes(set<int> &nodesAll) const {
  // cout << "PhyloDistance :: GetAllNodes: dump: ";
  // this->Dump();
  //
  nodesAll.clear();
  for (map<pair<int, int>, double>::const_iterator it = mapDists.begin();
       it != mapDists.end(); ++it) {
    nodesAll.insert(it->first.first);
    nodesAll.insert(it->first.second);
  }
}

double PhyloDistance ::CalcAveDistBtwClusters(
    const set<set<int> > &setClusters) const {
  //
  double res = 0.0;
  int numDist = 0;

  for (set<set<int> >::const_iterator it1 = setClusters.begin();
       it1 != setClusters.end(); ++it1) {
    set<set<int> >::const_iterator it2 = it1;
    ++it2;
    for (; it2 != setClusters.end(); ++it2) {
      // now sum over all dist
      for (set<int>::const_iterator it3 = it1->begin(); it3 != it1->end();
           ++it3) {
        for (set<int>::const_iterator it4 = it2->begin(); it4 != it2->end();
             ++it4) {
          res += GetDistance(*it3, *it4);
          ++numDist;
        }
      }
    }
  }

  return res / numDist;
}

void PhyloDistance ::Dump() const {
  //
  for (map<pair<int, int>, double>::const_iterator it = mapDists.begin();
       it != mapDists.end(); ++it) {
    cout << "[" << it->first.first << "," << it->first.second
         << "]: " << it->second << endl;
  }
}

// distance based tree builder

DistanceTreeBuilder ::DistanceTreeBuilder(PhyloDistance &distPairwiseTaxaIn)
    : distPairwiseTaxa(distPairwiseTaxaIn), taxonOutgroup(-1) {
  //
}

// build tree using neighbor joining
string DistanceTreeBuilder ::NJ() {
  // get all the things into the search set
  set<int> nodesToSearch;
  distPairwiseTaxa.GetAllNodes(nodesToSearch);

  // must have at least two nodes
  YW_ASSERT_INFO(nodesToSearch.size() >= 2, "Must have two nodes at least");

  // get the next largest one
  int nodeToUse = 1 + (*nodesToSearch.rbegin());

  // build a Newick string
  string strNW;
  map<int, string> mapSubtreeStr;
  for (set<int>::iterator it = nodesToSearch.begin(); it != nodesToSearch.end();
       ++it) {
    //
    // char buf[100];
    // sprintf(buf, "%d", *it);
    // string strName = buf;
    string strName = GetTaxonNameFor(*it);
    mapSubtreeStr.insert(map<int, string>::value_type(*it, strName));
    // cout << "Init node: " << *it << ": string: " << strName << endl;
  }

  int ngbr1 = -1, ngbr2 = -1;
  while (nodesToSearch.size() >= 3) {
    NJFindNgbrs(nodeToUse, nodesToSearch, ngbr1, ngbr2);
    // cout << "Neighbors found: " << ngbr1 << ", " << ngbr2 << ", and merged
    // into node: " << nodeToUse << endl;

    char buf1[100];
    sprintf(buf1, "%f", distPairwiseTaxa.GetDistanceNonNeg(nodeToUse, ngbr1));
    string strDist1 = buf1;
    char buf2[100];
    sprintf(buf2, "%f", distPairwiseTaxa.GetDistanceNonNeg(nodeToUse, ngbr2));
    string strDist2 = buf2;

    string strSubtree = "(" + mapSubtreeStr[ngbr1] + ":" + strDist1 + "," +
                        mapSubtreeStr[ngbr2] + ":" + strDist2 + ")";
    mapSubtreeStr.insert(map<int, string>::value_type(nodeToUse, strSubtree));
    // cout << "For node: " << nodeToUse << ": string: " << strSubtree << endl;
    ++nodeToUse;
  }

  // create a root
  int rootNode = nodeToUse;
  int ngbr1Final = *(nodesToSearch.begin());
  int ngbr2Final = *(nodesToSearch.rbegin());
  double distRootBranch =
      distPairwiseTaxa.GetDistanceNonNeg(ngbr1Final, ngbr2Final);
  double distNew = 0.5 * distRootBranch;
  distPairwiseTaxa.SetDistance(rootNode, ngbr1Final, distNew);
  distPairwiseTaxa.SetDistance(rootNode, ngbr2Final, distNew);

  // final subtree
  char buf1[100];
  sprintf(buf1, "%f", distNew);
  string strDist = buf1;
  string strSubtree = "(" + mapSubtreeStr[ngbr1Final] + ":" + strDist + "," +
                      mapSubtreeStr[ngbr2Final] + ":" + strDist + ")";
  // cout << "Final neighbor joining tree: " << strSubtree << endl;
  // now dump out all distances
  // distPairwiseTaxa.Dump();

  return strSubtree;
}

void DistanceTreeBuilder ::NJFindNgbrs(int nodeIdNew, set<int> &nodesToSearch,
                                       int &ngbr1, int &ngbr2) {
  // cout << "set of nodes to search: ";
  // DumpIntSet( nodesToSearch);
  // find two best ngbrs from the nodes to search (which are
  // ngbr1 and ngbr2, and remove these two from the nodes to search)
  // create a new node with the given id, add into nodestosearch, update dists

  // first compute ave distances of all current nodes
  map<int, double> mapAveDists;
  for (set<int>::iterator it = nodesToSearch.begin(); it != nodesToSearch.end();
       ++it) {
    //
    double dist = NJCalcAveDist(*it, nodesToSearch);
    // cout << "single node distance for " << *it << ": " << dist << endl;
    mapAveDists.insert(map<int, double>::value_type(*it, dist));
  }

  // search all pair to find the best one to merge
  double distNJMin = HAP_MAX_INT * 1.0;
  int node1Min = -1;
  int node2Min = -1;
  double dist12Min = 0.0;
  for (set<int>::iterator it1 = nodesToSearch.begin();
       it1 != nodesToSearch.end(); ++it1) {
    int node1cur = *it1;
    // don't consider outgroup
    if (node1cur == taxonOutgroup) {
      continue;
    }

    double dist1 = mapAveDists[node1cur];

    set<int>::iterator it2 = it1;
    ++it2;

    for (; it2 != nodesToSearch.end(); ++it2) {
      int node2cur = *it2;
      if (node2cur == taxonOutgroup) {
        continue;
      }
      double dist12 = distPairwiseTaxa.GetDistance(node1cur, node2cur);
      double dist2 = mapAveDists[node2cur];
      double distNJ = dist12 - dist1 - dist2;
      // cout << "For nodes: " << node1cur << ", " << node2cur << ", dist1=" <<
      // dist1 << ", dist2: " << dist2 << ", dist12: " << dist12 << ", distNJ: "
      // << distNJ << endl;
      if (distNJ < distNJMin) {
        distNJMin = distNJ;
        node1Min = node1cur;
        node2Min = node2cur;
        dist12Min = dist12;
      }
    }
  }

  // add the new node with right dist
  YW_ASSERT_INFO(node1Min >= 0 && node2Min >= 0, "Wrong");
  double dist1toNew =
      0.5 * dist12Min + 0.5 * (mapAveDists[node1Min] - mapAveDists[node2Min]);
  double dist2toNew =
      0.5 * dist12Min + 0.5 * (mapAveDists[node2Min] - mapAveDists[node1Min]);
  distPairwiseTaxa.SetDistance(nodeIdNew, node1Min, dist1toNew);
  distPairwiseTaxa.SetDistance(nodeIdNew, node2Min, dist2toNew);

  // calc remaining distances
  for (set<int>::iterator it = nodesToSearch.begin(); it != nodesToSearch.end();
       ++it) {
    int nodecur = *it;
    if (nodecur == node1Min || nodecur == node2Min) {
      continue;
    }
    double distNew = 0.5 * (distPairwiseTaxa.GetDistance(nodecur, node1Min) +
                            distPairwiseTaxa.GetDistance(nodecur, node2Min) -
                            distPairwiseTaxa.GetDistance(node1Min, node2Min));
    // dist12Min);
    distPairwiseTaxa.SetDistance(nodeIdNew, nodecur, distNew);
  }

  // maintain the search set
  nodesToSearch.insert(nodeIdNew);
  nodesToSearch.erase(node1Min);
  nodesToSearch.erase(node2Min);

  ngbr1 = node1Min;
  ngbr2 = node2Min;
}

double DistanceTreeBuilder ::NJCalcAveDist(int nodecur,
                                           const set<int> &nodesToSearch) {
  // calc average distance from nodecur to all nodes in the search set
  // must have at least three nodes
  YW_ASSERT_INFO(nodesToSearch.size() >= 3, "Too few nodes");
  // YW_ASSERT_INFO( nodesToSearch.find(nodecur) != nodesToSearch.end(),
  // "current node must be in the set");
  double res = 0.0;
  for (set<int>::const_iterator it = nodesToSearch.begin();
       it != nodesToSearch.end(); ++it) {
    if (*it != nodecur) {
      res += distPairwiseTaxa.GetDistance(nodecur, *it);
    }
  }
  return res / (nodesToSearch.size() - 2);
}

string DistanceTreeBuilder ::GetTaxonNameFor(int index) const {
  //
  map<int, string>::const_iterator it = mapIndexToName.find(index);
  if (it == mapIndexToName.end()) {
    // juse use the index itself
    char buf[100];
    sprintf(buf, "%d", index);
    string res(buf);
    return res;
  } else {
    return it->second;
  }
}

//********************************************************************************************************
// UPGMA utilities

string DistanceTreeBuilder ::ConstrainedUPGMA(
    const set<set<int> > &setClustersMustHave,
    const set<set<int> > &setClustersForbidden, map<set<int>, double> &mapSTHts,
    int numTotElem) {
  // construct UPGMA trees with constraints that exclude some clusters and must
  // have some clusters
  map<pair<set<int>, set<int> >, double> mapClusDist;
  map<set<int>, pair<string, double> > mapClusSubtree; // subtree with height
  // init all singleton
  set<int> nodesAll;
  distPairwiseTaxa.GetAllNodes(nodesAll);
  // cout << "nodesAll: ";
  // DumpIntSet(nodesAll);
  for (set<int>::const_iterator it1 = nodesAll.begin(); it1 != nodesAll.end();
       ++it1) {
    set<int> ss1;
    ss1.insert(*it1);
    string strLeaf = std::to_string(*it1);
    pair<string, double> sp(strLeaf, 0.0);
    mapClusSubtree.insert(
        map<set<int>, pair<string, double> >::value_type(ss1, sp));
    // cout << "Process leaf: " << strLeaf << endl;

    set<int>::const_iterator it2 = it1;
    ++it2;
    for (; it2 != nodesAll.end(); ++it2) {
      set<int> ss2;
      ss2.insert(*it2);
      pair<set<int>, set<int> > ss(ss1, ss2);
      mapClusDist.insert(map<pair<set<int>, set<int> >, double>::value_type(
          ss, distPairwiseTaxa.GetDistance(*it1, *it2)));
      // cout << "init pairwise distance with leaf " << *it2 << " dist=" <<
      // distPairwiseTaxa.GetDistance(*it1, *it2) << endl;
    }
  }
  // now start UPGMA procedure
  while (mapClusDist.size() >= 1) {
    // cout << "size of mapClusDist: " << mapClusDist.size() << endl;
    // find the smallest dist
    map<pair<set<int>, set<int> >, double>::iterator itOpt = mapClusDist.end();
    for (map<pair<set<int>, set<int> >, double>::iterator it =
             mapClusDist.begin();
         it != mapClusDist.end(); ++it) {
      // cout << "Dist=" << it->second << ", subtree1: ";
      // DumpIntSet(it->first.first);
      // cout << "subtree2: ";
      // DumpIntSet(it->first.second);
      set<int> scoal = it->first.first;
      UnionSets(scoal, it->first.second);

      if (itOpt != mapClusDist.end() && itOpt->second <= it->second) {
        // cout << "not optimal\n";
        continue;
      }

      bool fForbid =
          setClustersForbidden.find(scoal) != setClustersForbidden.end();
      if (fForbid == true) {
        // cout << "Not allowed\n";
        continue;
      }
      bool fCompat = IsClusterIncompatibleWithSetofClus(
          scoal, setClustersMustHave, numTotElem);

      if (fCompat == true) {
        itOpt = it;
      } else {
        // cout << "Not compatible\n";
      }
    }
    // must find something
    if (itOpt == mapClusDist.end()) {
      YW_ASSERT_INFO(false, "Fail to construct the tree");
    }
    // cout << "Best pair to merge: ";
    // DumpIntSet(itOpt->first.first);
    // cout << " and ";
    // DumpIntSet(itOpt->first.second);
    // now merge the two
    set<int> ssNew = itOpt->first.first;
    UnionSets(ssNew, itOpt->first.second);
    YW_ASSERT_INFO(
        mapClusSubtree.find(itOpt->first.first) != mapClusSubtree.end() &&
            mapClusSubtree.find(itOpt->first.second) != mapClusSubtree.end(),
        "Clusters: not found");
    double htSt1 = mapClusSubtree[itOpt->first.first].second;
    double htSt2 = mapClusSubtree[itOpt->first.second].second;
    double distSt1 = itOpt->second / 2 - htSt1;
    double distSt2 = itOpt->second / 2 - htSt2;
    // YW_ASSERT_INFO( distSt1 >= 0.0 && distSt2 >= 0.0, "Distance: should be
    // positive" );
    string strDist1 = std::to_string(distSt1);
    string strDist2 = std::to_string(distSt2);
    string strST = "(";
    strST += mapClusSubtree[itOpt->first.first].first;
    strST += ":";
    strST += strDist1;
    strST += ",";
    strST += mapClusSubtree[itOpt->first.second].first;
    strST += ":";
    strST += strDist2;
    strST += ")";
    pair<string, double> sp2(strST, itOpt->second / 2);
    mapClusSubtree.insert(
        map<set<int>, pair<string, double> >::value_type(ssNew, sp2));
    // cout << "subtree: " << strST << ", height: " << itOpt->second/2 << ", for
    // subtree: "; DumpIntSet( ssNew );
    // update the distance map
    UpdateDistUPGMA(itOpt->first, mapClusSubtree, mapClusDist);
    // cout << "mapClusDist: size = " << mapClusDist.size() << endl;
  }
  //
  YW_ASSERT_INFO(mapClusSubtree.find(nodesAll) != mapClusSubtree.end(),
                 "Not fully constructed yet");
  string strNWHt = mapClusSubtree[nodesAll].first;

  // record subtree ht
  mapSTHts.clear();
  for (map<set<int>, pair<string, double> >::iterator it =
           mapClusSubtree.begin();
       it != mapClusSubtree.end(); ++it) {
    mapSTHts.insert(
        map<set<int>, double>::value_type(it->first, it->second.second));
  }

  // strNWHt += ":" + std::to_string( mapClusSubtree[nodesAll].second );
  return strNWHt;
}

bool DistanceTreeBuilder ::IsClusterIncompatible(const set<int> &clus1,
                                                 const set<int> &clus2,
                                                 int numTotElem) const {
  // cout << "Clus1: ";
  // DumpIntSet(clus1);
  // cout << "clus2: ";
  // DumpIntSet(clus2);
  // four gamate test
  set<int> sint;
  JoinSets(clus1, clus2, sint);
  if (sint.size() == 0) {
    return true;
  }
  // set<int> sdiff1 = clus1;
  // SubtractSets(sdiff1, clus2);
  if (sint == clus1 || sint == clus2) {
    return true;
  }
  if (numTotElem > 0) {
    set<int> sunion = clus1;
    UnionSets(sunion, clus2);
    if ((int)sunion.size() == numTotElem) {
      return true;
    }
  }
  // set<int> sdiff2 = clus2;
  // SubtractSets(sdiff2, clus1);
  // if( sdiff2.size() == clus1.size() )
  //{
  //    return true;
  //}
  return false;
}

bool DistanceTreeBuilder ::IsClusterIncompatibleWithSetofClus(
    const set<int> &clus1, const set<set<int> > &setClus,
    int numTotElem) const {
  //
  for (set<set<int> >::const_iterator it = setClus.begin(); it != setClus.end();
       ++it) {
    if (IsClusterIncompatible(clus1, *it, numTotElem) == false) {
      return false;
    }
  }
  return true;
}

void DistanceTreeBuilder ::UpdateDistUPGMA(
    const pair<set<int>, set<int> > &pairClus,
    const map<set<int>, pair<string, double> > &mapSubtree,
    map<pair<set<int>, set<int> >, double> &distMapCur) {
  // remove all entries with one components as the newly merged subtree
  map<pair<set<int>, set<int> >, double> distMapUpdated;
  // set<set<int> > setClusCurr;
  for (map<pair<set<int>, set<int> >, double>::iterator it = distMapCur.begin();
       it != distMapCur.end(); ++it) {
    if (it->first.first != pairClus.first &&
        it->first.first != pairClus.second &&
        it->first.second != pairClus.first &&
        it->first.second != pairClus.second) {
      distMapUpdated.insert(*it);
      // setClusCurr.insert( it->first.first );
      // setClusCurr.insert( it->first.second );
    }
  }

  set<int> snew = pairClus.first;
  UnionSets(snew, pairClus.second);
  YW_ASSERT_INFO(mapSubtree.find(snew) != mapSubtree.end(), "Fail to find223");

  // collect all subsets that are not done yet
  set<set<int> > setsToProc;
  for (map<pair<set<int>, set<int> >, double>::const_iterator it =
           distMapCur.begin();
       it != distMapCur.end(); ++it) {
    setsToProc.insert(it->first.first);
    setsToProc.insert(it->first.second);
  }

  // now update distance with the new one
  // set< set<int> > setsDone;
  for (set<set<int> >::const_iterator it = setsToProc.begin();
       it != setsToProc.end(); ++it) {
    // if( setsDone.find(*it) != setsDone.end() )
    //{
    //    continue;
    //}
    // setsDone.insert( *it );
    // cout << "process cluster: ";
    // DumpIntSet(*it);
    if (*it == pairClus.first || *it == pairClus.second || *it == snew) {
      // cout << "Skipped\n";
      continue;
    }

    set<int> s1 = snew;
    set<int> s2 = *it;
    if (s2 < s1) {
      s1 = *it;
      s2 = snew;
    }

    pair<set<int>, set<int> > pp1(pairClus.first, *it);
    if (*it < pairClus.first) {
      pp1.first = *it;
      pp1.second = pairClus.first;
    }
    // cout << "pp1: ";
    // DumpIntSet(pp1.first);
    // DumpIntSet(pp1.second);
    YW_ASSERT_INFO(distMapCur.find(pp1) != distMapCur.end(), "Fail to find111");
    double htSt1 = distMapCur[pp1];
    pair<set<int>, set<int> > pp2(pairClus.second, *it);
    if (*it < pairClus.second) {
      pp2.first = *it;
      pp2.second = pairClus.second;
    }
    YW_ASSERT_INFO(distMapCur.find(pp2) != distMapCur.end(), "Fail to find112");
    double htSt2 = distMapCur[pp2];
    double distNew =
        (pairClus.first.size() * htSt1 + pairClus.second.size() * htSt2) /
        (pairClus.first.size() + pairClus.second.size());
    // cout << "htSt1: " << htSt1 << ", htSt2: " << htSt2 << ", distNew: " <<
    // distNew << ", for clusters: " << endl; DumpIntSet(s1); DumpIntSet(s2);
    pair<set<int>, set<int> > sp(s1, s2);
    distMapUpdated.insert(
        map<pair<set<int>, set<int> >, double>::value_type(sp, distNew));
  }

  // update map
  distMapCur = distMapUpdated;
}

string DistanceTreeBuilder ::ConstrainedUPGMA(
    const set<set<int> > &setClustersMustHave,
    const set<set<int> > &setClustersDesired, int numTopCandidates,
    const set<set<int> > &setClustersForbidden, map<set<int>, double> &mapSTHts,
    int numTotElem) {
  // picking the top-k candidates that matches the best of the desired splits
  // construct UPGMA trees with constraints that exclude some clusters and must
  // have some clusters
  map<pair<set<int>, set<int> >, double> mapClusDist;
  map<set<int>, pair<string, double> > mapClusSubtree; // subtree with height
  // init all singleton
  set<int> nodesAll;
  distPairwiseTaxa.GetAllNodes(nodesAll);
  // cout << "nodesAll: ";
  // DumpIntSet(nodesAll);
  for (set<int>::const_iterator it1 = nodesAll.begin(); it1 != nodesAll.end();
       ++it1) {
    set<int> ss1;
    ss1.insert(*it1);
    string strLeaf = std::to_string(*it1);
    pair<string, double> sp(strLeaf, 0.0);
    mapClusSubtree.insert(
        map<set<int>, pair<string, double> >::value_type(ss1, sp));
    // cout << "Process leaf: " << strLeaf << endl;

    set<int>::const_iterator it2 = it1;
    ++it2;
    for (; it2 != nodesAll.end(); ++it2) {
      set<int> ss2;
      ss2.insert(*it2);
      pair<set<int>, set<int> > ss(ss1, ss2);
      mapClusDist.insert(map<pair<set<int>, set<int> >, double>::value_type(
          ss, distPairwiseTaxa.GetDistance(*it1, *it2)));
      // cout << "init pairwise distance with leaf " << *it2 << " dist=" <<
      // distPairwiseTaxa.GetDistance(*it1, *it2) << endl;
    }
  }
  // now start UPGMA procedure
  while (mapClusDist.size() >= 1) {
    // cout << "size of mapClusDist: " << mapClusDist.size() << endl;
    // find the smallest dist
    map<double, set<pair<set<int>, set<int> > > > mapScoredPairs;
    int index = 0;
    const double MIN_DIST_INC = 0.00000000000000000000000001;
    for (map<pair<set<int>, set<int> >, double>::iterator it =
             mapClusDist.begin();
         it != mapClusDist.end(); ++it) {
      ++index;
      // cout << "Dist=" << it->second << ", subtree1: ";
      // DumpIntSet(it->first.first);
      // cout << "subtree2: ";
      // DumpIntSet(it->first.second);
      set<int> scoal = it->first.first;
      UnionSets(scoal, it->first.second);
      double distUse = it->second + index * MIN_DIST_INC;

      // if( (int)mapScoredPairs.size() >= numTopCandidates &&
      // mapScoredPairs.rbegin()->first < distUse )
      //{
      //    //cout << "not optimal\n";
      //    continue;
      //}

      bool fForbid =
          setClustersForbidden.find(scoal) != setClustersForbidden.end();
      if (fForbid == true) {
        // cout << "Not allowed\n";
        continue;
      }
      bool fCompat = IsClusterIncompatibleWithSetofClus(
          scoal, setClustersMustHave, numTotElem);

      if (fCompat == true) {
        // add it to the list
        mapScoredPairs[distUse].insert(it->first);
      } else {
        // cout << "Not compatible\n";
      }
    }
    // find the best
    int maxDesired = -1;
    pair<set<int>, set<int> > ppBest;
    double htBest = 0.0;
    YW_ASSERT_INFO(mapScoredPairs.size() > 0, "Must have some candidates");
    const double THRES_DIST_RATIO = 1.05;
    double minHit = mapScoredPairs.begin()->first;
    int index2 = 0;
    for (map<double, set<pair<set<int>, set<int> > > >::iterator it =
             mapScoredPairs.begin();
         it != mapScoredPairs.end(); ++it, ++index2) {
      if (index2 >= numTopCandidates || it->first > minHit * THRES_DIST_RATIO) {
        break;
      }

      for (set<pair<set<int>, set<int> > >::const_iterator it2 =
               it->second.begin();
           it2 != it->second.end(); ++it2) {
        // cout << "Choice: ";
        // DumpIntSet(it2->first);
        // cout << "   with ";
        // DumpIntSet(it2->second);
        set<int> ssCombo = it2->first;
        UnionSets(ssCombo, it2->second);

        // YW: desired clade must have exact match
        int numDesired = 0;
        if (setClustersDesired.find(ssCombo) != setClustersDesired.end()) {
          numDesired = 1;
        }

        // int numDesired = GetNumCompatCladesIn( ssCombo, setClustersDesired,
        // numTotElem );
        if (numDesired > maxDesired) {
          maxDesired = numDesired;
          ppBest = *it2;
          htBest = it->first;
        }
        // cout << "Hitting number of desired one: " << numDesired << endl;
      }
    }
    // cout << "Chosen one: ";
    // DumpIntSet(ppBest.first);
    // cout << "    with ";
    // DumpIntSet(ppBest.second);

    // cout << "Best pair to merge: ";
    // DumpIntSet(itOpt->first.first);
    // cout << " and ";
    // DumpIntSet(itOpt->first.second);
    // now merge the two
    set<int> ssNew = ppBest.first;
    UnionSets(ssNew, ppBest.second);
    YW_ASSERT_INFO(mapClusSubtree.find(ppBest.first) != mapClusSubtree.end() &&
                       mapClusSubtree.find(ppBest.second) !=
                           mapClusSubtree.end(),
                   "Clusters: not found");
    double htCurr = htBest;
    double htSt1 = mapClusSubtree[ppBest.first].second;
    double htSt2 = mapClusSubtree[ppBest.second].second;
    double distSt1 = htBest / 2 - htSt1;
    double distSt2 = htBest / 2 - htSt2;
    // YW_ASSERT_INFO( distSt1 >= 0.0 && distSt2 >= 0.0, "Distance: should be
    // positive" );
    string strDist1 = std::to_string(distSt1);
    string strDist2 = std::to_string(distSt2);
    string strST = "(";
    strST += mapClusSubtree[ppBest.first].first;
    strST += ":";
    strST += strDist1;
    strST += ",";
    strST += mapClusSubtree[ppBest.second].first;
    strST += ":";
    strST += strDist2;
    strST += ")";
    pair<string, double> sp2(strST, htCurr / 2);
    mapClusSubtree.insert(
        map<set<int>, pair<string, double> >::value_type(ssNew, sp2));
    // cout << "subtree: " << strST << ", height: " << itOpt->second/2 << ", for
    // subtree: "; DumpIntSet( ssNew );
    // update the distance map
    UpdateDistUPGMA(ppBest, mapClusSubtree, mapClusDist);
    // cout << "mapClusDist: size = " << mapClusDist.size() << endl;
  }
  //
  YW_ASSERT_INFO(mapClusSubtree.find(nodesAll) != mapClusSubtree.end(),
                 "Not fully constructed yet");
  string strNWHt = mapClusSubtree[nodesAll].first;

  // record subtree ht
  mapSTHts.clear();
  for (map<set<int>, pair<string, double> >::iterator it =
           mapClusSubtree.begin();
       it != mapClusSubtree.end(); ++it) {
    mapSTHts.insert(
        map<set<int>, double>::value_type(it->first, it->second.second));
  }

  // strNWHt += ":" + std::to_string( mapClusSubtree[nodesAll].second );
  return strNWHt;
}

int DistanceTreeBuilder ::GetNumCompatCladesIn(
    const set<int> &clus1, const set<set<int> > &setCladesTest,
    int numTotElem) const {
  //
  int res = 0;
  for (set<set<int> >::const_iterator it = setCladesTest.begin();
       it != setCladesTest.end(); ++it) {
    if (IsClusterIncompatible(clus1, *it, numTotElem) == true) {
      ++res;
    }
  }
  return res;
}

//***********************************************************************
// tool for building UPGMA tree

ConstrainedUPGMATreeBuilder ::ConstrainedUPGMATreeBuilder(
    PhyloDistance &distPairwiseTaxaIn,
    const set<set<int> > &setClustersMustHaveIn,
    const set<set<int> > &setClustersForbiddenIn, int numTotElemIn)
    : distPairwiseTaxa(distPairwiseTaxaIn),
      setClustersMustHave(setClustersMustHaveIn),
      setClustersForbidden(setClustersForbiddenIn), numTotElem(numTotElemIn) {
  Init();
}

ConstrainedUPGMATreeBuilder ::ConstrainedUPGMATreeBuilder(
    const ConstrainedUPGMATreeBuilder &rhs)
    : distPairwiseTaxa(rhs.distPairwiseTaxa),
      setClustersMustHave(rhs.setClustersMustHave),
      setClustersForbidden(rhs.setClustersForbidden),
      numTotElem(rhs.numTotElem), distMapActivePair(rhs.distMapActivePair),
      mapClusSubtree(rhs.mapClusSubtree), histSTMerge(rhs.histSTMerge) {
  //
}

string ConstrainedUPGMATreeBuilder ::GetTree() const {
  set<int> nodesAll;
  distPairwiseTaxa.GetAllNodes(nodesAll);
  map<set<int>, pair<string, double> >::const_iterator it =
      mapClusSubtree.find(nodesAll);
  YW_ASSERT_INFO(it != mapClusSubtree.end(), "Not fully constructed yet");
  string strNWHt = it->second.first;
  // strNWHt += ":" + std::to_string( mapClusSubtree[nodesAll].second );
  return strNWHt;
}

string ConstrainedUPGMATreeBuilder ::GetPartialConsTree() const {
  // get partially constructed tree for now; only consider those merged; that
  // is, if nothing occur, empty
  map<set<int>, string> mapSTs;
  //
  for (int i = 0; i < (int)histSTMerge.size(); ++i) {
    //
    map<set<int>, string>::iterator it1 = mapSTs.find(histSTMerge[i].first);
    map<set<int>, string>::iterator it2 = mapSTs.find(histSTMerge[i].second);

    //
    string strLeft, strRight;
    if (it1 == mapSTs.end()) {
      YW_ASSERT_INFO(histSTMerge[i].first.size() == 1, "Singleton");
      char buf[10000];
      sprintf(buf, "%d", *histSTMerge[i].first.begin());
      strLeft = buf;
    } else {
      strLeft = it1->second;
    }
    if (it2 == mapSTs.end()) {
      YW_ASSERT_INFO(histSTMerge[i].second.size() == 1, "Singleton");
      char buf[10000];
      sprintf(buf, "%d", *histSTMerge[i].second.begin());
      strRight = buf;
    } else {
      strRight = it2->second;
    }
    string strLeftUse = strLeft;
    string strRightUse = strRight;
    if (strRight < strLeft) {
      strLeftUse = strRight;
      strRightUse = strLeft;
    }

    string strMerge = "(" + strLeftUse + "," + strRightUse + ")";
    set<int> ss = histSTMerge[i].first;
    UnionSets(ss, histSTMerge[i].second);
    mapSTs[ss] = strMerge;
    if (it1 != mapSTs.end()) {
      mapSTs.erase(it1);
    }
    if (it2 != mapSTs.end()) {
      mapSTs.erase(it2);
    }
  }
  // result is concatnation of all the remaining stuff
  string res = "(";
  for (map<set<int>, string>::iterator it = mapSTs.begin(); it != mapSTs.end();
       ++it) {
    if (it != mapSTs.begin()) {
      res += ",";
    }
    res += it->second;
  }
  res += ")";

  return res;
}

double ConstrainedUPGMATreeBuilder ::GetMinCoalSubtrees(set<int> &st1,
                                                        set<int> &st2) const {
  // cout << "*GetMinCoalSubtrees\n";
  map<pair<set<int>, set<int> >, double>::const_iterator itOpt =
      distMapActivePair.end();
  for (map<pair<set<int>, set<int> >, double>::const_iterator it =
           distMapActivePair.begin();
       it != distMapActivePair.end(); ++it) {
    // cout << "Dist=" << it->second << ", subtree1: ";
    // DumpIntSet(it->first.first);
    // cout << "subtree2: ";
    // DumpIntSet(it->first.second);
    set<int> scoal = it->first.first;
    UnionSets(scoal, it->first.second);

    if (itOpt != distMapActivePair.end() && itOpt->second <= it->second) {
      // cout << "not optimal\n";
      continue;
    }

    bool fForbid =
        setClustersForbidden.find(scoal) != setClustersForbidden.end();
    if (fForbid == true) {
      // cout << "Not allowed: forbidden\n";
      continue;
    }
    bool fCompat =
        IsClusterIncompatibleWithSetofClus(scoal, setClustersMustHave);

    if (fCompat == false) {
      // cout << "Not allowed: incomaptible\n";
      continue;
    }

    //
    itOpt = it;
  }
  // must find something
  if (itOpt == distMapActivePair.end()) {
    YW_ASSERT_INFO(false, "Fail to construct the tree");
  }
  // cout << "here..\n";
  st1 = itOpt->first.first;
  st2 = itOpt->first.second;
  // cout << "Min dist: " << itOpt->second << ", subtrees: ";
  // DumpIntSet(st1);
  // cout << "  and ";
  // DumpIntSet(st2);
  return itOpt->second;
}

void ConstrainedUPGMATreeBuilder ::GetCoalSubtreesHtBound(
    double htBound,
    set<pair<pair<set<int>, set<int> >, double> > &setCandidates) const {
  //
  for (map<pair<set<int>, set<int> >, double>::const_iterator it =
           distMapActivePair.begin();
       it != distMapActivePair.end(); ++it) {
    // cout << "Dist=" << it->second << ", subtree1: ";
    // DumpIntSet(it->first.first);
    // cout << "subtree2: ";
    // DumpIntSet(it->first.second);

    if (it->second > htBound) {
      // cout << "not optimal\n";
      continue;
    }

    set<int> scoal = it->first.first;
    UnionSets(scoal, it->first.second);

    bool fForbid =
        setClustersForbidden.find(scoal) != setClustersForbidden.end();
    if (fForbid == true) {
      // cout << "Not allowed\n";
      continue;
    }
    bool fCompat =
        IsClusterIncompatibleWithSetofClus(scoal, setClustersMustHave);

    if (fCompat == false) {
      continue;
    }

    //
    pair<pair<set<int>, set<int> >, double> pp;
    pp.first = it->first;
    pp.second = it->second;
    setCandidates.insert(pp);
  }
}

void ConstrainedUPGMATreeBuilder ::MergeSubtrees(const set<int> &st1,
                                                 const set<int> &st2,
                                                 double htMergedST) {
  // now merge the two
  set<int> ssNew = st1;
  UnionSets(ssNew, st2);
  YW_ASSERT_INFO(mapClusSubtree.find(st1) != mapClusSubtree.end() &&
                     mapClusSubtree.find(st2) != mapClusSubtree.end(),
                 "Clusters: not found");
  double htSt1 = mapClusSubtree[st1].second;
  double htSt2 = mapClusSubtree[st2].second;
  double distSt1 = htMergedST / 2 - htSt1;
  double distSt2 = htMergedST / 2 - htSt2;
  // YW_ASSERT_INFO( distSt1 >= 0.0 && distSt2 >= 0.0, "Distance: should be
  // positive" );
  string strDist1 = std::to_string(distSt1);
  string strDist2 = std::to_string(distSt2);
  string strST = "(";
  strST += mapClusSubtree[st1].first;
  strST += ":";
  strST += strDist1;
  strST += ",";
  strST += mapClusSubtree[st2].first;
  strST += ":";
  strST += strDist2;
  strST += ")";
  double distSet = htMergedST / 2;
  pair<string, double> sp2(strST, distSet);
  mapClusSubtree.insert(
      map<set<int>, pair<string, double> >::value_type(ssNew, sp2));
  // cout << "MergeSubtrees: subtree: " << strST << ", height: " << htMergedST
  // << ", for subtree: "; DumpIntSet( ssNew );
  // update the distance map
  UpdateDistUPGMA(st1, st2);

  pair<set<int>, set<int> > pp(st1, st2);
  histSTMerge.push_back(pp);
}

void ConstrainedUPGMATreeBuilder ::GetMergeCandidates(
    map<pair<set<int>, set<int> >, double> &setCandidates) const {
  setCandidates.clear();
  for (map<pair<set<int>, set<int> >, double>::const_iterator it =
           distMapActivePair.begin();
       it != distMapActivePair.end(); ++it) {
    // cout << "GetMergeCandidates: candidate clades: ";
    // DumpIntSet(it->first.first);
    // cout << "  ";
    // DumpIntSet(it->first.second);
    set<int> scoal = it->first.first;
    UnionSets(scoal, it->first.second);
    bool fForbid =
        setClustersForbidden.find(scoal) != setClustersForbidden.end();
    if (fForbid == true) {
      // cout << "Not allowed\n";
      continue;
    }
    // cout << "socal: not forbidden\n";
    // DumpIntSet(scoal);
    bool fCompat =
        IsClusterIncompatibleWithSetofClus(scoal, setClustersMustHave);

    if (fCompat == false) {
      continue;
    }
    // cout << "A good candidate: ";
    // DumpIntSet(it->first.first);
    // cout << "  ";
    // DumpIntSet(it->first.second);

    setCandidates.insert(map<pair<set<int>, set<int> >, double>::value_type(
        it->first, it->second));
  }
  // cout << "Done: GetMergeCandidates\n";
}

double ConstrainedUPGMATreeBuilder ::GetCurDistForTwoClusters(
    const set<int> &clus1, const set<int> &clus2) const {
  //
  pair<set<int>, set<int> > ss(clus1, clus2);
  map<pair<set<int>, set<int> >, double>::const_iterator it =
      distMapActivePair.find(ss);
  YW_ASSERT_INFO(it != distMapActivePair.end(), "Fail to find");
  return it->second;
}

void ConstrainedUPGMATreeBuilder ::SetDistForTwoClusters(const set<int> &clus1,
                                                         const set<int> &clus2,
                                                         double dist) {
  pair<set<int>, set<int> > ss(clus1, clus2);
  map<pair<set<int>, set<int> >, double>::const_iterator it =
      distMapActivePair.find(ss);
  YW_ASSERT_INFO(it != distMapActivePair.end(), "Fail to find");
  distMapActivePair[ss] = dist;
}

bool ConstrainedUPGMATreeBuilder ::IsDone() const {
  return distMapActivePair.size() == 0;
}

void ConstrainedUPGMATreeBuilder ::Init() {
  set<int> nodesAll;
  distPairwiseTaxa.GetAllNodes(nodesAll);
  // cout << "nodesAll: ";
  // DumpIntSet(nodesAll);
  for (set<int>::const_iterator it1 = nodesAll.begin(); it1 != nodesAll.end();
       ++it1) {
    set<int> ss1;
    ss1.insert(*it1);
    string strLeaf = std::to_string(*it1);
    pair<string, double> sp(strLeaf, 0.0);
    mapClusSubtree.insert(
        map<set<int>, pair<string, double> >::value_type(ss1, sp));
    // cout << "Process leaf: " << strLeaf << endl;

    set<int>::const_iterator it2 = it1;
    ++it2;
    for (; it2 != nodesAll.end(); ++it2) {
      set<int> ss2;
      ss2.insert(*it2);

      pair<set<int>, set<int> > ss(ss1, ss2);
      distMapActivePair.insert(
          map<pair<set<int>, set<int> >, double>::value_type(
              ss, distPairwiseTaxa.GetDistance(*it1, *it2)));
      // cout << "init pairwise distance with leaf " << *it2 << " dist=" <<
      // distPairwiseTaxa.GetDistance(*it1, *it2) << endl;
    }
  }
}

bool ConstrainedUPGMATreeBuilder ::IsClusterIncompatible(
    const set<int> &clus1, const set<int> &clus2) const {
  // cout << "Clus1: ";
  // DumpIntSet(clus1);
  // cout << "clus2: ";
  // DumpIntSet(clus2);
  // four gamate test
  set<int> sint;
  JoinSets(clus1, clus2, sint);
  if (sint.size() == 0) {
    return true;
  }
  // set<int> sdiff1 = clus1;
  // SubtractSets(sdiff1, clus2);
  if (sint == clus1 || sint == clus2) {
    return true;
  }
  // set<int> sdiff2 = clus2;
  // SubtractSets(sdiff2, clus1);
  // if( sdiff2.size() == clus1.size() )
  //{
  //    return true;
  //}
  if (this->numTotElem > 0) {
    set<int> sunion = clus1;
    UnionSets(sunion, clus2);
    if ((int)sunion.size() == numTotElem) {
      return true;
    }
  }
  return false;
}

bool ConstrainedUPGMATreeBuilder ::IsClusterIncompatibleWithSetofClus(
    const set<int> &clus1, const set<set<int> > &setClus) const {
  //
  for (set<set<int> >::const_iterator it = setClus.begin(); it != setClus.end();
       ++it) {
    if (IsClusterIncompatible(clus1, *it) == false) {
      return false;
    }
  }
  return true;
}

void ConstrainedUPGMATreeBuilder ::UpdateDistUPGMA(const set<int> &st1,
                                                   const set<int> &st2) {
  // remove all entries with one components as the newly merged subtree
  map<pair<set<int>, set<int> >, double> distMapUpdated;
  // set<set<int> > setClusCurr;
  for (map<pair<set<int>, set<int> >, double>::iterator it =
           distMapActivePair.begin();
       it != distMapActivePair.end(); ++it) {
    if (it->first.first != st1 && it->first.first != st2 &&
        it->first.second != st1 && it->first.second != st2) {
      distMapUpdated.insert(*it);
      // setClusCurr.insert( it->first.first );
      // setClusCurr.insert( it->first.second );
    }
  }

  set<int> snew = st1;
  UnionSets(snew, st2);
  YW_ASSERT_INFO(mapClusSubtree.find(snew) != mapClusSubtree.end(),
                 "Fail to find223");

  // collect all subsets that are not done yet
  set<set<int> > setsToProc;
  for (map<pair<set<int>, set<int> >, double>::const_iterator it =
           distMapActivePair.begin();
       it != distMapActivePair.end(); ++it) {
    setsToProc.insert(it->first.first);
    setsToProc.insert(it->first.second);
  }

  // now update distance with the new one
  // set< set<int> > setsDone;
  for (set<set<int> >::const_iterator it = setsToProc.begin();
       it != setsToProc.end(); ++it) {
    // if( setsDone.find(*it) != setsDone.end() )
    //{
    //    continue;
    //}
    // setsDone.insert( *it );
    // cout << "process cluster: ";
    // DumpIntSet(*it);
    if (*it == st1 || *it == st2 || *it == snew) {
      // cout << "Skipped\n";
      continue;
    }

    // make sure this is allowed
    // set<int> scoal = snew;
    // UnionSets( scoal, *it);
    // bool fForbid = setClustersForbidden.find(scoal) !=
    // setClustersForbidden.end(); if( fForbid == true )
    //{
    //    //cout << "Not allowed\n";
    //    continue;
    //}
    // bool fCompat = IsClusterIncompatibleWithSetofClus( scoal,
    // setClustersMustHave );
    //
    // if( fCompat == false)
    //{
    //    continue;
    //}

    set<int> s1 = snew;
    set<int> s2 = *it;
    if (s2 < s1) {
      s1 = *it;
      s2 = snew;
    }

    pair<set<int>, set<int> > pp1(st1, *it);
    if (*it < st1) {
      pp1.first = *it;
      pp1.second = st1;
    }
    // cout << "pp1: ";
    // DumpIntSet(pp1.first);
    // DumpIntSet(pp1.second);
    YW_ASSERT_INFO(distMapActivePair.find(pp1) != distMapActivePair.end(),
                   "Fail to find111");
    double htSt1 = distMapActivePair[pp1];
    pair<set<int>, set<int> > pp2(st2, *it);
    if (*it < st2) {
      pp2.first = *it;
      pp2.second = st2;
    }
    YW_ASSERT_INFO(distMapActivePair.find(pp2) != distMapActivePair.end(),
                   "Fail to find112");
    double htSt2 = distMapActivePair[pp2];
    double distNew =
        (st1.size() * htSt1 + st2.size() * htSt2) / (st1.size() + st2.size());
    // cout << "htSt1: " << htSt1 << ", htSt2: " << htSt2 << ", distNew: " <<
    // distNew << ", for clusters: " << endl; DumpIntSet(s1); DumpIntSet(s2);
    pair<set<int>, set<int> > sp(s1, s2);
    distMapUpdated.insert(
        map<pair<set<int>, set<int> >, double>::value_type(sp, distNew));
  }

  // update map
  distMapActivePair = distMapUpdated;

  // cout <<"After update, ";
  // Dump();
}

int ConstrainedUPGMATreeBuilder ::GetNumSubtrees() const {
  //
  return mapClusSubtree.size();
}

void ConstrainedUPGMATreeBuilder ::GetAllSubtrees(
    map<set<int>, string> &mapSTs) const {
  //
  mapSTs.clear();
  for (map<set<int>, pair<string, double> >::const_iterator it =
           mapClusSubtree.begin();
       it != mapClusSubtree.end(); ++it) {
    //
    mapSTs.insert(
        map<set<int>, string>::value_type(it->first, it->second.first));
  }
}

void ConstrainedUPGMATreeBuilder ::GetActiveSubtrees(
    set<set<int> > &setActiveSTs) const {
  //
  for (map<pair<set<int>, set<int> >, double>::const_iterator it =
           distMapActivePair.begin();
       it != distMapActivePair.end(); ++it) {
    setActiveSTs.insert(it->first.first);
    setActiveSTs.insert(it->first.second);
  }
}

void ConstrainedUPGMATreeBuilder::Dump() const {
  cout << "List of coalescent pairs: \n";
  for (map<pair<set<int>, set<int> >, double>::const_iterator it =
           distMapActivePair.begin();
       it != distMapActivePair.end(); ++it) {
    cout << "[" << it->second << "] ";
    DumpIntSet(it->first.first);
    DumpIntSet(it->first.second);
  }
}

//***********************************************************************
// tool for building near-optimal UPGMA tree

ConstrainedNearUPGMATreesBuilder ::ConstrainedNearUPGMATreesBuilder(
    PhyloDistance &distPairwiseTaxaIn,
    const set<set<int> > &setClustersMustHaveIn,
    const set<set<int> > &setClustersForbiddenIn, int numTotElemIn)
    : distPairwiseTaxa(distPairwiseTaxaIn),
      setClustersMustHave(setClustersMustHaveIn),
      setClustersForbidden(setClustersForbiddenIn), numTotElem(numTotElemIn) {}

void ConstrainedNearUPGMATreesBuilder ::Construct(int maxNumTrees,
                                                  double thresMaxDistRatio) {
  // thresMaxDistRatio: say 1.2, meaning consdiering 1.2*min distance to use as
  // candidate
  YW_ASSERT_INFO(thresMaxDistRatio >= 1.0,
                 "Threshold: cannot be less than 1.0");

  map<string, ConstrainedUPGMATreeBuilder *> listTreeBuilders;
  // start with a single tree
  ConstrainedUPGMATreeBuilder *pBuild0 = new ConstrainedUPGMATreeBuilder(
      this->distPairwiseTaxa, this->setClustersMustHave,
      this->setClustersForbidden, this->numTotElem);
  string strDummy;
  listTreeBuilders[strDummy] = pBuild0;

  // start to build near-upgma trees
  while (true) {
    // process each builder
    map<string, ConstrainedUPGMATreeBuilder *> listTreeBuildersNext;
    bool fDone = false;
    for (map<string, ConstrainedUPGMATreeBuilder *>::iterator it =
             listTreeBuilders.begin();
         it != listTreeBuilders.end(); ++it) {
      ConstrainedUPGMATreeBuilder *pCurr = it->second;
      // perform
      if (pCurr->IsDone() == true) {
        fDone = true;
        break;
      }
      set<int> st1, st2;
      double minDist = pCurr->GetMinCoalSubtrees(st1, st2);

      // get near-min dist
      double distUse = thresMaxDistRatio * minDist;
      // set< pair<pair<set<int>, set<int> >, double> > setCandidates;
      // listTreeBuilders[i]->GetCoalSubtreesHtBound(distUse, setCandidates );
      map<pair<set<int>, set<int> >, double> setCandidates;
      pCurr->GetMergeCandidates(setCandidates);

      YW_ASSERT_INFO(setCandidates.size() > 0, "Fail to find candidates");

      // process if room for more trees
      for (map<pair<set<int>, set<int> >, double>::iterator it =
               setCandidates.begin();
           it != setCandidates.end(); ++it) {
        if (it->second > distUse) {
          continue;
        }

        // make sure it is not the mimimum one found before
        if ((it->first.first != st1 || it->first.second != st2) &&
            (it->first.first != st2 || it->first.second != st1))

        {
          if ((int)listTreeBuildersNext.size() < maxNumTrees) {
            ConstrainedUPGMATreeBuilder *pBuildCopy =
                new ConstrainedUPGMATreeBuilder(*pCurr);
            pBuildCopy->MergeSubtrees(it->first.first, it->first.second,
                                      it->second);
            string strTreeCons = pBuildCopy->GetPartialConsTree();
            if (listTreeBuildersNext.find(strTreeCons) ==
                listTreeBuildersNext.end()) {
              listTreeBuildersNext[strTreeCons] = pBuildCopy;
              // cout << "Candidate merge: ";
              // DumpIntSet(it->first.first);
              // cout << "  ";
              // DumpIntSet(it->first.second);
            } else {
              delete pBuildCopy;
            }
          }
        }
      }

      // do the merge of the optimal one
      pCurr->MergeSubtrees(st1, st2, minDist);
      string strTreeCons2 = pCurr->GetPartialConsTree();
      if (listTreeBuildersNext.find(strTreeCons2) ==
          listTreeBuildersNext.end()) {
        listTreeBuildersNext[strTreeCons2] = pCurr;
        // cout << "Candidate (minimum) merge: ";
        // DumpIntSet(st1);
        // cout << "  ";
        // DumpIntSet(st2);
      } else {
        delete pCurr;
      }
    }
    if (fDone) {
      break;
    }

    // add if there is no duplicate
    listTreeBuilders = listTreeBuildersNext;
  }

  // collect trees
  setTreeCons.clear();
  for (map<string, ConstrainedUPGMATreeBuilder *>::iterator it =
           listTreeBuilders.begin();
       it != listTreeBuilders.end(); ++it) {
    string treres = it->second->GetTree();
    setTreeCons.insert(treres);
    // cout << "Tree constructed: " << treres << endl;
  }

  // clean
  for (map<string, ConstrainedUPGMATreeBuilder *>::iterator it =
           listTreeBuilders.begin();
       it != listTreeBuilders.end(); ++it) {
    delete it->second;
  }
  listTreeBuilders.clear();
}
