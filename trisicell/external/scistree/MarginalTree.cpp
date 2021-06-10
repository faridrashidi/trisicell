#include "MarginalTree.h"
#include "PhylogenyTreeBasic.h"
#include "UnWeightedGraph.h"
#include "Utils4.h"
#include <queue>
#include <stack>

//////////////////////////////////////////////////////////////////////////////
// Define a utility class
// GLobal utility function

// static void OutputQuotedString(ofstream &outFile, const char *buf)
//{
//	outFile << '"';
//	outFile << buf;
//	outFile << '"';
//}

void RemapLeafIntLabelsTaxaMap(MarginalTree &mtree,
                               map<string, string> &mapper) {
  // map the leaf labels to new integer labels
  for (int i = 0; i < mtree.GetNumLeaves(); ++i) {
    int lbl = mtree.GetLabel(i);
    char buf[100];
    sprintf(buf, "%d", lbl);
    string strbuf = buf;
    YW_ASSERT_INFO(mapper.find(strbuf) != mapper.end(), "Fail to find");
    string strLbl = mapper[strbuf];
    int lblNewInt;
    sscanf(strLbl.c_str(), "%d", &lblNewInt);
    mtree.SetLabel(i, lblNewInt);
  }
}

void RemapMargTree(MarginalTree &mtree, TaxaMapper &refTMapper) {
  //
  // map the leaf labels to new integer labels
  // cout << "RemapMargTree: mtree:" << mtree.GetNewick() << endl;
  // mtree.Dump();
  for (int i = 0; i < mtree.GetNumLeaves(); ++i) {
    int lbl = mtree.GetLabel(i);
    string strlbl = refTMapper.GetString(lbl);
    int lblNew = lbl;
    sscanf(strlbl.c_str(), "%d", &lblNew);
    mtree.SetLabel(i, lblNew);
  }
}

static bool ReadinOneMarginalTree(ifstream &inFile, int numNodes,
                                  MarginalTree &tree) {
  // first read in the node ids
  for (int i = 0; i < numNodes; ++i) {
    int tmp;
    inFile >> tmp;
    tree.listNodeLabels.push_back(tmp);
  }
  for (int i = 0; i < numNodes; ++i) {
    int tmp;
    inFile >> tmp;
    tree.listParentNodePos.push_back(tmp);
  }
  for (int i = 0; i < numNodes; ++i) {
    double tmp;
    inFile >> tmp;
    tree.listEdgeDist.push_back(tmp);
  }

  return true;
}

static void ReadNewickLen(const string &strNewick,
                          map<set<int>, double> &mapClusterLen,
                          TaxaMapper *pTMapper) {
  // cout << "ReadNewickLen: strNewick = " << strNewick << endl;
  // the first letter must be (
  // YW_ASSERT_INFO( strNewick.length() > 0 && strNewick[0] == '(', "Bad Newick
  // format" );

  const char *strNwBuf = strNewick.c_str();

  // reverse and find the last ) to get dist
  int posLastState = -1;
  double bLen = 1.0;
  for (int i = (int)strNewick.length() - 1; i >= 0; --i) {
    if (strNewick[i] == ':') {
      float fLen = 1.0;
      sscanf(strNwBuf + i + 1, "%f", &fLen);
      bLen = fLen;
      if (strNewick[i] != ')') {
        posLastState = i - 2;
        break;
      }
    } else if (strNewick[i] == ')') {
      // should also stop
      posLastState = i - 1;
      break;
    }
  }
  // accumlate all the labels
  PhylogenyTreeBasic phTree;
  phTree.ConsOnNewick(strNewick, -1, false, pTMapper);

  // see if zero is in, if not, must have 1 and decrease by 1
  set<int> lvids;
  phTree.GetLeaveIds(lvids);
  // cout << "ReadNewickLen: lvids = ";
  // DumpIntSet( lvids );

  // add a record
  mapClusterLen.insert(map<set<int>, double>::value_type(lvids, bLen));
  // cout << "Subtree len = " << bLen << ", for leaf set = ";
  // DumpIntSet( lvids );

  // given newick format, read in the edge length of the clusters
  // we will perform this recursively
  // first find the position where it is the first ,
  int posSplit = -1;
  int netParen = 0;
  for (int i = 0; i < (int)strNewick.length(); ++i) {
    if (strNewick[i] == '(') {
      netParen++;
    } else if (strNewick[i] == ')') {
      netParen--;
    }
    if (netParen == 1 && strNewick[i] == ',') {
      posSplit = i;
      break;
    }
  }
  // YW_ASSERT_INFO( netParen >= 0 && posSplit >= 1, "Bad Newick format" );

  // now recurisvely to two children (if needed)
  if (posSplit >= 0) {
    YW_ASSERT_INFO(posSplit - 1 >= 1, "Newick format wrong");
    string strLeft = strNewick.substr(1, posSplit - 1);
    ReadNewickLen(strLeft, mapClusterLen, pTMapper);
    YW_ASSERT_INFO(posSplit + 1 <= posLastState, "Newick format wrong");
    string strRight = strNewick.substr(posSplit + 1, posLastState - posSplit);
    ReadNewickLen(strRight, mapClusterLen, pTMapper);
  }
}

static int UpdateMTreeWithNWString(MarginalTree &treeToChange, int &leafNext,
                                   int &nodeIntNext, string &strNewick,
                                   TaxaMapper *pTMapper) {
  // cout << "UpdateMTreeWithNWString: strNewick = " << strNewick  << ",
  // leafNext = " << leafNext << ", nodeIntNext = " << nodeIntNext << endl;
  // a recursive call to change all the nodes from nodeToChnage to the correct
  // length specified by the string strNewick (and all the underlying nodes)
  // return the current node

  // conslidate the newick string first
  string strNewickUse = strNewick;
  NewickUtils ::ConsolidateSinglChildChain(strNewickUse);
  if (strNewickUse != strNewick) {
    // cout << "**Newick: " << strNewick << ", after consolidate: " <<
    // strNewickUse << endl;
  }
  // first find current length by finding the rightmost : (outside any ))
  // now find the separator in order to proceed recurisvely
  string strNW1, strNW2;
  bool fNonAtom = NewickUtils ::FindSplitIn(strNewickUse, strNW1, strNW2);
  int nodeCurrent;
  if (fNonAtom == true) {
    // now recursive
    if (nodeIntNext < treeToChange.GetNumLeaves()) {
      treeToChange.Dump();
      cout << "nodeIntNext: " << nodeIntNext << ", ";
      cout << "Tree to chagne: " << treeToChange.GetNewick() << endl;
    }
    YW_ASSERT_INFO(nodeIntNext >= treeToChange.GetNumLeaves(),
                   "UpdateBranchLenInfo: internal node out of range");
    nodeCurrent = nodeIntNext--;
  } else {
    YW_ASSERT_INFO(leafNext < treeToChange.GetNumLeaves(),
                   "UpdateBranchLenInfo: Leaf out of range");
    // this node is a leaf
    nodeCurrent = leafNext++;

    // int idNew;
    // sscanf(strNewick.c_str(), "%d", &idNew);
    int idNew = TaxaMapper ::GetIdFromStr(strNewickUse, pTMapper);
    treeToChange.SetLabel(nodeCurrent, idNew);
  }
  // cout << "nodeCurrent = " << nodeCurrent << endl;

  // by default, set branch length to be 1.0
  float lenCur = 1.0;
  size_t posSep1 = strNewickUse.rfind(':');
  size_t posSep2 = strNewickUse.rfind(')');
  if (posSep1 != string::npos &&
      (posSep1 > posSep2 || posSep2 == string::npos)) {
    // yes, there is a length specified for this
    sscanf(strNewickUse.c_str() + (int)posSep1 + 1, "%f", &lenCur);
  }
  double lenCurUse = lenCur;
  treeToChange.SetBranchLen(nodeCurrent, lenCurUse);
  // cout << "Found length: " << lenCurUse << " for node " << nodeCurrent <<
  // endl;
  // now recurse
  // cout << "In UpdateMTreeWithNWString: \n";
  // treeToChange.Dump();
  if (fNonAtom == true) {
    // now recursive
    int nodeChild1 = UpdateMTreeWithNWString(treeToChange, leafNext,
                                             nodeIntNext, strNW1, pTMapper);
    int nodeChild2 = UpdateMTreeWithNWString(treeToChange, leafNext,
                                             nodeIntNext, strNW2, pTMapper);

    // update the par pos
    treeToChange.SetParent(nodeChild1, nodeCurrent, false);
    treeToChange.SetParent(nodeChild2, nodeCurrent, false);
  }

  return nodeCurrent;
}

bool ReadinMarginalTrees(ifstream &inFile, vector<MarginalTree> &treeList) {
  // first read in the number of chrom
  int numLeaves;
  inFile >> numLeaves;
  int nTreeNodes;
  inFile >> nTreeNodes;
  int nTrees;
  inFile >> nTrees;
  treeList.clear();
  for (int i = 0; i < nTrees; ++i) {
    // cout << "Reading TREE " << i << endl;
    MarginalTree tree;
    ReadinOneMarginalTree(inFile, nTreeNodes, tree);
    tree.numLeaves = numLeaves;
    treeList.push_back(tree);

    // tree.Dump();
    // YW_ASSERT_INFO(false, "early abort");
  }
  return true;
}

void CollapseEquivTrees(const vector<MarginalTree> &listOrigTrees,
                        vector<MarginalTree> &listUniqTrees,
                        vector<int> &listMultiplicity) {
  // collect ordered-leaf tree newick strings
  vector<string> listRepOrderedLeafList;
  for (int tr = 0; tr < (int)listOrigTrees.size(); ++tr) {
    int numLeaves = listOrigTrees[tr].GetNumLeaves();
    // make the tree binary
    // listTreesMT[tr].Binarize();
    // cout << "Processing gene tree: ";
    // listTreesMT[tr].Dump();
    PhylogenyTreeBasic *pphtree = new PhylogenyTreeBasic;
    pphtree->ConsOnParPosList(listOrigTrees[tr].listParentNodePos, numLeaves,
                              true);
    pphtree->UpdateIntLabel(listOrigTrees[tr].listNodeLabels);
    pphtree->Order();
    string strNewick1;
    pphtree->ConsNewick(strNewick1);
    delete pphtree;
    //		cout << "Constructed one gene Tree = " << strNewick1 << endl;
    listRepOrderedLeafList.push_back(strNewick1);
    //		const int ROOT_LABLE = 7;
    // TestReroot(pphtree, ROOT_LABLE);
  }

  listUniqTrees.clear();
  listMultiplicity.clear();
  vector<string> listStoredLeafList;

  // check each tree in orig list
  for (int tr = 0; tr < (int)listRepOrderedLeafList.size(); ++tr) {
    bool fFound = false;
    for (int trNew = 0; trNew < (int)listStoredLeafList.size(); ++trNew) {
      if (listStoredLeafList[trNew] == listRepOrderedLeafList[tr]) {
        // add it
        listMultiplicity[trNew]++;
        fFound = true;
        break;
      }
    }
    if (fFound == false) {
      listUniqTrees.push_back(listOrigTrees[tr]);
      listMultiplicity.push_back(1);
      listStoredLeafList.push_back(listRepOrderedLeafList[tr]);
    }
  }
}

bool ReadinMarginalTreesNewick(ifstream &inFile, int numLeaves,
                               vector<MarginalTree> &treeList,
                               TaxaMapper *pTMapper, bool fDup) {
  // NOTE: RETURN TRUE IF NO LABEL ADJUSTMENT IS DONE
  // RETURN FALSE IF WE SWITCHED LABEL BY DECREASING BY ONE
  // figure out leave num
  bool fNoChange = true;
  int nLvs = numLeaves;

  // read marginal trees in newick format
  // here there is no preamble, one line per tree
  while (inFile.eof() == false) {
    string treeNewick;
    inFile >> treeNewick;
    if (treeNewick.size() == 0) {
      break;
    }
    // cout << "newick tree = " << treeNewick << endl;

    //#if 0
    // update numleaves
    multiset<string> setLabels;
    NewickUtils ::RetrieveLabelSet(treeNewick, setLabels);
    nLvs = setLabels.size();
    //#endif
    //
    PhylogenyTreeBasic phTree;
    // if( fDup == false )
    //{
    phTree.ConsOnNewick(treeNewick, -1, false, pTMapper);
    //}
    // else
    //{
    //	phTree.ConsOnNewickDupLabels(treeNewick, pTMapper);
    //}

    if (pTMapper != NULL) {
      pTMapper->SetInitialized(true);
    }
    // string strTr;
    // phTree.ConsNewick(strTr);
    // cout << "After reconstruction: strTr = " << strTr << endl;
    // see if zero is in, if not, must have 1 and decrease by 1
    set<int> lvids;
    phTree.GetLeaveIds(lvids);
    if (lvids.find(0) == lvids.end()) {
      YW_ASSERT_INFO(lvids.find(1) != lvids.end(), "Wrong");

      // decrease by one
      phTree.InitPostorderWalk();
      while (true) {
        TreeNode *pn = phTree.NextPostorderWalk();
        if (pn == NULL) {
          break; // done with all nodes
        }
        if (pn->IsLeaf() == true) {
          // cout << "Found leaf id: " << pn->GetID() << endl;
          pn->SetID(pn->GetID() - 1);
          // YW: 8/18/11, changed. NEED VERIFICATION
          // char buf[1000];
          // sprintf(buf, "%d", pn->GetID() );
          // string lblNew = buf;
          // pn->SetLabel( lblNew );
        }
      }

      // mark the change
      fNoChange = false;
    }

    vector<int> nidsList, nparsList;
    phTree.GetNodeParInfo(nidsList, nparsList);
    // cout << "nidsList: ";
    // DumpIntVec( nidsList );
    // cout << "nparsList" << ": ";
    // DumpIntVec(nparsList);
    // phTree.GetNodeParInfoNew(nidsList, nparsList);
    // phTree.GetNodeParInfo(nidsList, nparsList);
    // if( nLvs <= 0 )
    //{
    // YW: 09072010, ASSUME the tree is binary tree
    nLvs = (phTree.GetNumVertices() + 1) / 2;
    // cout << "nlvs = " << nLvs << endl;
    //}
    MarginalTree tree;
    InitMarginalTree(tree, nLvs, nidsList, nparsList);
    // cout << "After init, mtree = ";
    // tree.Dump();
    // YW: 01/30/12, sort the leaf first
    tree.SortByLeafId();
    // cout << "After sorting, tree = ";
    // tree.Dump();

    // cout << "Initialize a tree: ";
    // tree.Dump();
    treeList.push_back(tree);
    // cout << "Newick format of this marginal tree: ";
    // cout << tree.GetNewick() << endl;
  }
  return fNoChange;
}

bool ReadinMarginalTreesNewickWLenString(const string &strNewick, int numLeaves,
                                         MarginalTree &treeOut,
                                         bool fStartFromZero,
                                         TaxaMapper *pTMapper) {
  // YW_ASSERT_INFO(pTMapper != NULL, "Stop here2");
  // mark the change
  bool fNoChange = true;
  // NOTE: RETURN TRUE IF NO LABEL ADJUSTMENT IS DONE
  // RETURN FALSE IF WE SWITCHED LABEL BY DECREASING BY ONE
  // figure out leave num

  if (strNewick.size() == 0) {
    return fNoChange;
  }
  // cout << "newick tree = " << strNewick << endl;

  // make sure leave num is correct
  if (numLeaves < 0) {
    //
    multiset<string> setLabels;
    NewickUtils ::RetrieveLabelSet(strNewick, setLabels);
    numLeaves = setLabels.size();
    // cout << "Set number of leaves of marginal tree to: " << numLeaves <<
    // endl;
  }

  int nLvs = numLeaves;

  // assume binary tree for now
  int numTotNodes = 2 * nLvs - 1;
  // init Marginal tree for now
  vector<int> trLbls, trPos;
  vector<double> trDist;
  for (int i = 0; i < numTotNodes; ++i) {
    trLbls.push_back(i);
    trPos.push_back(-1);
    trDist.push_back(0.0);
  }
  treeOut.SetNumLeaves(nLvs);
  treeOut.SetLabelList(trLbls);
  treeOut.SetParList(trPos);
  treeOut.SetBranchLenList(trDist);
  // InitMarginalTree(treeOut, nLvs, trLbls, trPos);

  // now update tree
  int leafNext = 0;
  int nodeIntNext = numTotNodes - 1;
  string strNewickUse = strNewick;
  UpdateMTreeWithNWString(treeOut, leafNext, nodeIntNext, strNewickUse,
                          pTMapper);
  // cout << "Immediate after UpdateMTreeWithNWString: treeOut: \n";
  // treeOut.Dump();

  // finally prepare marginal tree for query
  treeOut.BuildDescendantInfo();
  // cout << "ReadinMarginalTreesNewickWLenString: newick string = \n" <<
  // treeOut.GetNewick() << endl;

  if (pTMapper != NULL) {
    pTMapper->SetInitialized(true);
  }

  return fNoChange;
}

bool ReadinMarginalTreesNewickWLen(ifstream &inFile, int numLeaves,
                                   vector<MarginalTree> &treeList,
                                   TaxaMapper *pTMapper) {
  // YW_ASSERT_INFO(pTMapper != NULL, "Stop here");
  // NOTE: RETURN TRUE IF NO LABEL ADJUSTMENT IS DONE
  // RETURN FALSE IF WE SWITCHED LABEL BY DECREASING BY ONE
  // figure out leave num
  bool fNoChange = true;
  // int nLvs = numLeaves;

  // read marginal trees in newick format
  // here there is no preamble, one line per tree
  while (inFile.eof() == false) {
    string treeNewick;
    inFile >> treeNewick;
    if (treeNewick.size() == 0) {
      break;
    }
    MarginalTree tree;
    bool fres = ReadinMarginalTreesNewickWLenString(treeNewick, numLeaves, tree,
                                                    true, pTMapper);
    if (fres == false) {
      fNoChange = false;
    }
    if (pTMapper != NULL) {
      pTMapper->SetInitialized(true);
    }

    // cout << "Initialize a tree: ";
    // tree.Dump();
    treeList.push_back(tree);
  }
  return fNoChange;
}

void AddRootAsLeafToTree(MarginalTree &tree1, bool fIdNonNeg) {
  // cout << "AddRootAsLeafToTree: tree1 = \n";
  // tree1.Dump();
  // we now add the root to the tree as a special leaf
  vector<int> nodesIdNew, nodesParsNew;
  for (int i = 0; i < tree1.GetNumLeaves(); ++i) {
    nodesIdNew.push_back(tree1.listNodeLabels[i]);
    nodesParsNew.push_back(tree1.GetParent(i) + 1);
  }
  // add the new special leaf
  int idLeafNew = -2; // -2 is the default unique id for this speical node
  int idNewStart = 3 * tree1.GetNumLeaves() + 1;
  if (fIdNonNeg == true) {
    // use  continuous id
    idLeafNew = tree1.GetNumLeaves();
  }
  nodesIdNew.push_back(idLeafNew);
  nodesParsNew.push_back(tree1.GetTotNodesNum() + 1);
  // add the rest
  for (int i = tree1.GetNumLeaves(); i < tree1.GetTotNodesNum(); ++i) {
    nodesIdNew.push_back(tree1.listNodeLabels[i]);
    int oldpar = tree1.GetParent(i);
    if (oldpar < 0) {
      // this is the old root
      nodesParsNew.push_back(tree1.GetTotNodesNum() + 1);
    } else {
      nodesParsNew.push_back(oldpar + 1);
    }
  }
  // finally the new root
  int idRootId = -3; // -3 is the unique id for this speical root
  if (fIdNonNeg == true) {
    // use  it
    idRootId = ++idNewStart;
  }
  nodesIdNew.push_back(idRootId);
  nodesParsNew.push_back(-1);

  // finally increment the number of leaves
  tree1.listNodeLabels = nodesIdNew;
  tree1.listParentNodePos = nodesParsNew;
  tree1.numLeaves++;
  // cout << "After adding the root, now tree1 = \n";
  // tree1.Dump();
}

void GenRandBinaryTree(int numLeaves, MarginalTree &tree1) {
  // generate a binary marginal tree with certain number of leaves
  // we do this by random pick two active nodes (a leave without assiging
  // parents)
  tree1.Clear();
  tree1.numLeaves = numLeaves;

  // first add a list of leaves
  set<int> activeNodes;
  for (int i = 0; i < numLeaves; ++i) {
    tree1.listNodeLabels.push_back(i);
    tree1.listParentNodePos.push_back(-1); // for now, set to -1
                                           // (un-initialized)
    tree1.listEdgeDist.push_back(0.0);
    activeNodes.insert(i);
  }

  // now start to setup new internal nodes (and assign parents)
  while (activeNodes.size() >= 2) {
    // cout << "activeNodes = ";
    // DumpIntSet( activeNodes );

    // uniformly pick two nodes
    int node1 = GetRandItemInSet(activeNodes);
    activeNodes.erase(node1);
    int node2 = GetRandItemInSet(activeNodes);
    activeNodes.erase(node2);
    // cout << "Select node1 = " << node1 << ", node2 = " << node2 << endl;
    // now create a new node
    int nodeNew = tree1.listNodeLabels.size();
    tree1.listNodeLabels.push_back(nodeNew);
    tree1.listParentNodePos.push_back(-1); // for now, set to -1
                                           // (un-initialized)
    tree1.listEdgeDist.push_back(0.0);
    activeNodes.insert(nodeNew);
    // cout << "nodeNew = " << nodeNew << endl;
    // setup parent of two children to it
    tree1.SetParent(node1, nodeNew);
    tree1.SetParent(node2, nodeNew);
  }
}

void GenRandBinaryTreeClock(int numLeaves, double totHt, MarginalTree &tree1) {
  // generate a binary marginal tree with certain number of leaves and have
  // clock property we do this by random pick two active nodes (a leave without
  // assiging parents)
  map<int, double> mapNodeHeights;

  tree1.Clear();
  tree1.numLeaves = numLeaves;

  // first add a list of leaves
  set<int> activeNodes;
  for (int i = 0; i < numLeaves; ++i) {
    tree1.listNodeLabels.push_back(i);
    tree1.listParentNodePos.push_back(-1); // for now, set to -1
                                           // (un-initialized)
    tree1.listEdgeDist.push_back(0.0);
    activeNodes.insert(i);
    mapNodeHeights.insert(map<int, double>::value_type(i, 0.0));
  }

  // now start to setup new internal nodes (and assign parents)
  while (activeNodes.size() >= 2) {
    // cout << "activeNodes = ";
    // DumpIntSet( activeNodes );

    // uniformly pick two nodes
    int node1 = GetRandItemInSet(activeNodes);
    activeNodes.erase(node1);
    int node2 = GetRandItemInSet(activeNodes);
    activeNodes.erase(node2);
    // cout << "Select node1 = " << node1 << ", node2 = " << node2 << endl;
    // now create a new node
    int nodeNew = tree1.listNodeLabels.size();
    tree1.listNodeLabels.push_back(nodeNew);
    tree1.listParentNodePos.push_back(-1); // for now, set to -1
                                           // (un-initialized)
    tree1.listEdgeDist.push_back(0.0);
    activeNodes.insert(nodeNew);
    double htNodeCur =
        totHt * (numLeaves - (double)activeNodes.size()) / (numLeaves - 1);
    // cout << "Node: " << nodeNew << ", ht = " << htNodeCur << endl;
    // set branches
    mapNodeHeights.insert(map<int, double>::value_type(nodeNew, htNodeCur));
    YW_ASSERT_INFO(mapNodeHeights.find(node1) != mapNodeHeights.end(),
                   "Not found");
    YW_ASSERT_INFO(node1 < (int)tree1.listEdgeDist.size(), "Wrong");
    tree1.listEdgeDist[node1] = htNodeCur - mapNodeHeights[node1];
    // cout << "Setting edge " << node1 << " to " <<  tree1.listEdgeDist[node1]
    // << endl;
    YW_ASSERT_INFO(tree1.listEdgeDist[node1] >= 0.0, "Negative");
    YW_ASSERT_INFO(mapNodeHeights.find(node2) != mapNodeHeights.end(),
                   "Not found");
    YW_ASSERT_INFO(node2 < (int)tree1.listEdgeDist.size(), "Wrong");
    tree1.listEdgeDist[node2] = htNodeCur - mapNodeHeights[node2];
    YW_ASSERT_INFO(tree1.listEdgeDist[node2] >= 0.0, "Negative");
    // cout << "Setting edge " << node2 << " to " <<  tree1.listEdgeDist[node2]
    // << endl; cout << "nodeNew = " << nodeNew << endl;
    // setup parent of two children to it
    tree1.SetParent(node1, nodeNew, false);
    tree1.SetParent(node2, nodeNew, false);
  }

  // cout << "Edge dist list: ";
  // DumpDoubleVec(tree1.listEdgeDist);
}

// find a chain with specified length
static bool FindChainAtNodeInTree(const MarginalTree &tree1, int nodeHead,
                                  int lenChain, vector<int> &leaves,
                                  vector<int> &leaves2) {
  leaves.clear();
  leaves2.clear();

  int curn = nodeHead;
  for (int i = 0; i < lenChain; ++i) {
    // cout << "curn = " << curn << endl;
    bool fLeftLeave = tree1.IsLeaf(tree1.GetLeftDescendant(curn));
    bool fRightLeave = tree1.IsLeaf(tree1.GetRightDescendant(curn));
    if ((fLeftLeave == true && fRightLeave == true && i < lenChain - 2) ||
        (fLeftLeave == false && fRightLeave == false)) {
      // cout << "Fail\n";
      return false;
    }
    // now move down
    if (fLeftLeave == false) {
      // put right to store
      int child = tree1.GetRightDescendant(curn);
      YW_ASSERT_INFO(tree1.IsLeaf(child), "Not a leaf");
      leaves.push_back(child);
      curn = tree1.GetLeftDescendant(curn);
    } else if (fRightLeave == false) {
      int child = tree1.GetLeftDescendant(curn);
      YW_ASSERT_INFO(tree1.IsLeaf(child), "Not a leaf");
      leaves.push_back(child);
      curn = tree1.GetRightDescendant(curn);
    } else {
      YW_ASSERT_INFO(i >= lenChain - 2, "wrong1");
      // YW_ASSERT_INFO( tree1.IsLeaf( curn ) == false, " a leaf" );
      // in this case, we just save itself
      if (i == lenChain - 2) {
        leaves2 = leaves;
        leaves.push_back(tree1.GetLeftDescendant(curn));
        leaves.push_back(tree1.GetRightDescendant(curn));
        leaves2.push_back(tree1.GetRightDescendant(curn));
        leaves2.push_back(tree1.GetLeftDescendant(curn));
        break;
      } else {
        // only save one, try both possibilities
        leaves2 = leaves;
        leaves.push_back(tree1.GetLeftDescendant(curn));
        // save leaves2
        leaves2.push_back(tree1.GetRightDescendant(curn));
        break;
      }
    }
  }
  return true;
}

// here, we simply found chain of fixed length (i.e. 4)
void FindChainsInTree(const MarginalTree &tree1,
                      map<vector<int>, int> &foundChains) {
  //
  foundChains.clear();

  // we simply enumerate all internal node and trace down from it
  for (int nn = tree1.GetNumLeaves(); nn < tree1.GetTotNodesNum(); ++nn) {
    // is this nn a chain-hnead of 4 leaves on one side?
    vector<int> listLeaves, listLeaves2;
    if (FindChainAtNodeInTree(tree1, nn, 4, listLeaves, listLeaves2) == true) {
      // cout << "Found one chain at node nn = " << nn << ", with leaves = ";
      // DumpIntVec( listLeaves);
      foundChains.insert(map<vector<int>, int>::value_type(listLeaves, nn));

      if (listLeaves2.size() > 0) {
        // cout << "Found one chain at node nn = " << nn << ", with leaves = ";
        // DumpIntVec( listLeaves2);
        foundChains.insert(map<vector<int>, int>::value_type(listLeaves2, nn));
      }
    }
  }
}

// construct a marginal tree from nodes and parent info
// NOTE: this function does not take distance. Therefore, we arbitarily assign
// nodes to their respective heights and thus also assign branch length ALSO
// NOTE: when we assign branch length, the branch length are set uniformly
// distributed within [0-1].
void InitMarginalTree(MarginalTree &mTree, int numLeaves,
                      const vector<int> &listLabels,
                      const vector<int> &listParentNodePos) {
  // cout << "numLeaves = " << numLeaves << endl;
  // cout << "InitMarginalTree: numLeaves = " << numLeaves << endl;
  // cout << "listLabels = ";
  // DumpIntVec(listLabels);
  // cout << "listParentNodePos = ";
  // DumpIntVec( listParentNodePos );
  //
  mTree.numLeaves = numLeaves;
  mTree.listNodeLabels = listLabels;
  mTree.listParentNodePos = listParentNodePos;

  // now init edge dist
  mTree.listEdgeDist.clear();
  int numNonLeafNodes = listLabels.size() - numLeaves;
  double unitLen = 1.0 / numNonLeafNodes;
  for (int i = 0; i < (int)listLabels.size() - 1; ++i) {
    int parPos = listParentNodePos[i] - numLeaves + 1;
    // cout << "par = " << listParentNodePos[i]  << " for node i = " << i <<
    // endl; cout << "normalized par pos = " << parPos << endl;
    YW_ASSERT_INFO(parPos > 0, "Fatal error in InitMarginalTree");
    if (i < numLeaves) {
      // leaf
      mTree.listEdgeDist.push_back(parPos * unitLen);
    } else {
      // need to subtract current pos
      int curpos = i - numLeaves + 1;
      // cout << "curpos = " << curpos << endl;
      YW_ASSERT_INFO(curpos < parPos, "Trouble in InitMarginalTree");
      mTree.listEdgeDist.push_back((parPos - curpos) * unitLen);
    }
  }
  // the root has length-0 by default
  mTree.listEdgeDist.push_back(0.0);
  // also build up descendents
  mTree.BuildDescendantInfo();
}

// find the neighborhood of marginal trees within one NNI operation away (incl.
// the current tree)
void FindOneNNIMTreesFrom(MarginalTree &mTreeSrc,
                          vector<MarginalTree> &listNNITrees,
                          vector<pair<int, int> > *pListPairEdgesSwapped) {
  //
  listNNITrees.clear();

  // process each internal node (w/ at least three leaves below) of the mtree,
  // and
  for (int node = mTreeSrc.GetNumLeaves(); node < mTreeSrc.GetTotNodesNum();
       ++node) {
    //
    int nodeLeft = mTreeSrc.GetLeftDescendant(node);
    int nodeRight = mTreeSrc.GetRightDescendant(node);
    if (mTreeSrc.IsLeaf(nodeLeft) == true &&
        mTreeSrc.IsLeaf(nodeRight) == true) {
      // skip if both children are leaves since in this case swapping has no
      // effect
      continue;
    }
    // now swap its two children's subtree in up to four ways
    int nodesProc1[2], nodesProc2[2];
    nodesProc1[0] = nodeLeft;
    nodesProc1[1] = nodeRight;
    nodesProc2[1] = nodeLeft;
    nodesProc2[0] = nodeRight;
    for (int ii = 0; ii < 2; ++ii) {
      int n1Proc = nodesProc1[ii];
      int n2Proc = nodesProc2[ii];
      if (mTreeSrc.IsLeaf(n1Proc) == false) {
        int node1Left = mTreeSrc.GetLeftDescendant(n1Proc);
        int node1Right = mTreeSrc.GetRightDescendant(n1Proc);
        YW_ASSERT_INFO(node1Left >= 0 && node1Right >= 0, "Can not miss");

        // two choices to swap: n2Proc with one of the descendents
        int nodesProc1Child[2];
        nodesProc1Child[0] = node1Left;
        nodesProc1Child[1] = node1Right;
        for (int jj = 0; jj < 2; ++jj) {
          MarginalTree mtreeNNI1 = mTreeSrc;
          mtreeNNI1.SwapBranches(nodesProc1Child[jj], n2Proc);
          mtreeNNI1.BuildDescendantInfo();
          // cout << "After swap: \n";
          // mtreeNNI1.Dump();
          mtreeNNI1.RearrangeParIncOrder();
          // cout << "Found a new mtreeNNI1: " << mtreeNNI1.GetNewick() << endl;
          // mtreeNNI1.Dump();
          mtreeNNI1.BuildDescendantInfo();
          // sort by leaf id: YW: Feb 19,2016
          mtreeNNI1.SortByLeafId();
          mtreeNNI1.BuildDescendantInfo();
          listNNITrees.push_back(mtreeNNI1);

          if (pListPairEdgesSwapped != NULL) {
            pair<int, int> pp(nodesProc1Child[jj], n2Proc);
            pListPairEdgesSwapped->push_back(pp);
          }
          // cout << "After descendent rebult, " << mtreeNNI1.GetNewick() <<
          // endl; mtreeNNI1.Dump();
        }
      }
    }
  }
  // finally add self
  listNNITrees.push_back(mTreeSrc);
  // exit(1);
}

void CreateSubtreeFromLeaves(MarginalTree &mTreeOrig,
                             const set<int> &setLeafLabels,
                             MarginalTree &mTreeSub,
                             map<int, int> &mapNewNodeToOldNode) {
  // cout << "Original tree: " << mTreeOrig.GetNewick() << ": set of leaves to
  // process: "; DumpIntSet( setLeafLabels );

  // find a subset of trees with the desired leaves (as matching the given
  // labels) mapNewNodeToOldNode: new node index ==> old node index
  map<pair<int, set<int> >, int> mapShrunkLeavesWithNum;

  // get all the clades
  for (int i = 0; i < mTreeOrig.GetTotNodesNum(); ++i) {
    //
    set<int> setGetDesc;
    mTreeOrig.GetLeavesUnder(i, setGetDesc);
    set<int> setGetDescLbls;
    for (set<int>::iterator it = setGetDesc.begin(); it != setGetDesc.end();
         ++it) {
      int lbl = mTreeOrig.GetLabel(*it);
      setGetDescLbls.insert(lbl);
    }
    set<int> sIntsect;
    JoinSets(setGetDescLbls, setLeafLabels, sIntsect);

    // ignore empty nodes
    if (sIntsect.size() <= 0) {
      //
      continue;
    }

    // save it
    pair<int, set<int> > ss(sIntsect.size(), sIntsect);
    if (mapShrunkLeavesWithNum.find(ss) == mapShrunkLeavesWithNum.end()) {
      mapShrunkLeavesWithNum.insert(
          map<pair<int, set<int> >, int>::value_type(ss, i));
    } else {
      // save the lower (smaller)
      if (mapShrunkLeavesWithNum[ss] > i) {
        mapShrunkLeavesWithNum[ss] = i;
      }
    }
  }
#if 0
cout << "mapShrunkLeavesWithNum: ";
for( map< pair<int,set<int> >, int > :: iterator it = mapShrunkLeavesWithNum.begin(); it != mapShrunkLeavesWithNum.end(); ++it )
{
cout << "Size: " << it->first.first << ", orig. node = " << it->second << ", set of leaves: ";
DumpIntSet( it->first.second);
}
#endif

  // set up the old and new node position map
  map<int, int> mapNewToOldPos, mapOldToNewPos;
  set<int> setNewParsPosOld;
  int index = 0;
  for (map<pair<int, set<int> >, int>::iterator it =
           mapShrunkLeavesWithNum.begin();
       it != mapShrunkLeavesWithNum.end(); ++it, ++index) {
    //
    mapNewToOldPos.insert(map<int, int>::value_type(index, it->second));
    mapOldToNewPos.insert(map<int, int>::value_type(it->second, index));

    setNewParsPosOld.insert(it->second);
  }

  // now init the tree: note edge labels are ignored!
  mTreeSub.Clear();
  mTreeSub.SetNumLeaves(setLeafLabels.size());
  vector<int> listLbls;
  PopulateVecBySet(listLbls, setLeafLabels);
  for (int i = (int)setLeafLabels.size();
       i < (int)mapShrunkLeavesWithNum.size(); ++i) {
    // these are internal nodes
    listLbls.push_back(-1);
  }
  mTreeSub.SetLabelList(listLbls);
  vector<int> listParPos;
  // now set up parent
  for (int i = 0; i < (int)listLbls.size(); ++i) {
    YW_ASSERT_INFO(mapNewToOldPos.find(i) != mapNewToOldPos.end(),
                   "Fail to find2");
    int posOrig = mapNewToOldPos[i];
    int anc = mTreeOrig.GetFirstNonselfAnces(posOrig, setNewParsPosOld);
    int posNewAnc = -1;
    if (anc >= 0) {
      YW_ASSERT_INFO(mapOldToNewPos.find(anc) != mapOldToNewPos.end(),
                     "Fail to find3");
      posNewAnc = mapOldToNewPos[anc];
    }
    listParPos.push_back(posNewAnc);
  }
  mTreeSub.SetParList(listParPos);

  // create nodes mapping
  mapNewNodeToOldNode = mapNewToOldPos;

  //
  mTreeSub.BuildDescendantInfo();

  // YW: how do we assign branch length
  UpdateBranchLenInSubtree(mTreeOrig, mapNewNodeToOldNode, mTreeSub);
#if 0
cout << "Constructed subtree: " << mTreeSub.GetNewick() << endl;
mTreeSub.Dump();
cout << "mapNewNodeToOldNode: ";
for(map<int,int> :: iterator it=mapNewNodeToOldNode.begin(); it != mapNewNodeToOldNode.end(); ++it)
{
cout << "[" << it->first << "," << it->second << "]  ";
}
cout << endl;
#endif
}

void UpdateBranchLenInSubtree(MarginalTree &mTreeOrig,
                              map<int, int> &mapNewNodeToOldNode,
                              MarginalTree &mTreeSub) {
  // inverse map
  // map<int,int> mapOldNodeToNewNode;
  // for( map<int,int> :: iterator it = mapNewNodeToOldNode.begin(); it !=
  // mapNewNodeToOldNode.end(); ++it  )
  //{
  //    //
  //    YW_ASSERT_INFO( mapOldNodeToNewNode.find(it->second) ==
  //    mapOldNodeToNewNode.end(), "Wrong" ); mapOldNodeToNewNode.insert(
  //    map<int,int> :: value_type(it->second, it->first) );
  //}

  //
  vector<double> listBrLens;
  for (map<int, int>::iterator it = mapNewNodeToOldNode.begin();
       it != mapNewNodeToOldNode.end(); ++it) {
    double distcur = 0.0;
    //
    int pnew = it->first;
    int pold = it->second;
    int pnewpar = mTreeSub.GetParent(pnew);
    if (pnewpar >= 0) {
      YW_ASSERT_INFO(mapNewNodeToOldNode.find(pnewpar) !=
                         mapNewNodeToOldNode.end(),
                     "Fail to find");
      int poldpar = mapNewNodeToOldNode[pnewpar];
      distcur = mTreeOrig.GetPathLen(pold, poldpar);
    }

    listBrLens.push_back(distcur);
  }
  mTreeSub.SetBranchLenList(listBrLens);
}

void FindMatchedSubtrees(MarginalTree &mtreeNew, MarginalTree &mtreeRef,
                         map<int, int> &mapSTNewToRef) {
  // find the shared subtrees that are in both trees, then create a map: map the
  // subtree index in mtreeNew to mtreeRef find all branches (subtrees below
  // them) that are not in the reference tree setDiffBrs: in this tree but not
  // in reference tree setDiffRefMissed: in reference tree but not in this tree
  vector<set<int> > listSubtreesNew, listSubtreesRef;
  mtreeNew.ConsDecedentLeavesInfoLabels(listSubtreesNew);
  mtreeRef.ConsDecedentLeavesInfoLabels(listSubtreesRef);

  // create fast searching
  map<set<int>, int> mapIndexSTRef;
  for (int i = 0; i < (int)listSubtreesRef.size(); ++i) {
    mapIndexSTRef.insert(map<set<int>, int>::value_type(listSubtreesRef[i], i));
  }

  //
  mapSTNewToRef.clear();
  for (int i = 0; i < (int)listSubtreesNew.size(); ++i) {
    if (mapIndexSTRef.find(listSubtreesNew[i]) == mapIndexSTRef.end()) {
      mapSTNewToRef.insert(
          map<int, int>::value_type(i, mapIndexSTRef[listSubtreesNew[i]]));
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// Define a utility class

MarginalTree ::MarginalTree() : numLeaves(0) {
  // here, we initiailize distance
  // TBD
}

void MarginalTree ::Clear() {
  numLeaves = 0;
  listNodeLabels.clear();
  listParentNodePos.clear();
  listEdgeDist.clear();
  listLeftDescs.clear();
  listRightDescs.clear();
}

void MarginalTree ::BuildDescendantInfo() {
  // Note, this only works for binary tree!!!!!
  listLeftDescs.clear();
  listRightDescs.clear();
  int numNodes = GetTotNodesNum();
  // cout << "BuildDescendantInfo: numNodes: " << numNodes << endl;
  listLeftDescs.resize(numNodes);
  listRightDescs.resize(numNodes);

  // for leaves, there is no children
  for (int i = 0; i < numNodes; ++i) {
    listLeftDescs[i] = -1;
    listRightDescs[i] = -1;
  }

  // handle other cases
  for (int i = 0; i < numNodes; ++i) {
    int p = GetParent(i);
    // cout << "paret of " << i << " is " << p << endl;
    if (p < 0) {
      continue;
    }
    // setup p's child to i
    if (listLeftDescs[p] < 0) {
      listLeftDescs[p] = i;
      // cout << "Set left descendent of " << p << " to be " << i << endl;
    } else {
      if (listRightDescs[p] >= 0) {
        cout << "Something wrong: the current tree:";
        Dump();
      }

      // make sure this is binary tree
      YW_ASSERT_INFO(listRightDescs[p] < 0, "Not a binary tree2");
      listRightDescs[p] = i;
      // cout << "Set right descendent of " << p << " to be " << i << endl;
    }
  }
}

bool MarginalTree ::IsToplogicSame(const MarginalTree &tree) const {
  // this function test whether two things are topologically the same
  if (GetTotNodesNum() != tree.GetTotNodesNum()) {
    // cout << "Tree node numbers are not equal\n";
    // nodes number are different, then different
    return false;
  }
  //
  if (GetNumLeaves() != tree.GetNumLeaves()) {
    // cout << "Tree leaf numbers are not equal\n";
    return false;
  }
  // make sure node id the same
  // if( listNodeLabels != tree.listNodeLabels )
  //{
// cout << "Tree node ids are not equal\n";
//        return false;
//    }
#if 0
    if( listParentNodePos != tree.listParentNodePos  )
    {
//cout << "Tree node parents are not equal\n";
        return false;
    }
#endif

  // sort the leaves
  // MarginalTree t1 = *this;
  // MarginalTree t2 = tree;
  // t1.SortByLeafId();
  // t2.SortByLeafId();
  vector<set<int> > t1splits, t2splits;
  ConsDecedentLeavesInfo(t1splits);
  tree.ConsDecedentLeavesInfo(t2splits);
  set<set<int> > st1splits, st2splits;
  for (int i = 0; i < (int)t1splits.size(); ++i) {
    st1splits.insert(t1splits[i]);
  }
  for (int i = 0; i < (int)t2splits.size(); ++i) {
    st2splits.insert(t2splits[i]);
  }

  if (st1splits != st2splits) {
    // cout << "Tree node parents are not equal\n";
    // cout << "** tree 1: \n";
    // for(int i=0; i<(int)t1splits.size(); ++i )
    //{
    // DumpIntSet(t1splits[i]);
    //}
    // cout << "** tree 2: \n";
    // for(int i=0; i<(int)t2splits.size(); ++i )
    //{
    // DumpIntSet(t2splits[i]);
    //}
    return false;
  }

  return true;
}

int MarginalTree ::GetLeftDescendant(int node) const {
  YW_ASSERT_INFO((int)listLeftDescs.size() == GetTotNodesNum() &&
                     (int)listRightDescs.size() == GetTotNodesNum(),
                 "descendant info not set");
  return listLeftDescs[node];
}
int MarginalTree ::GetRightDescendant(int node) const {
  YW_ASSERT_INFO((int)listLeftDescs.size() == GetTotNodesNum() &&
                     (int)listRightDescs.size() == GetTotNodesNum(),
                 "descendant info not set");
  return listRightDescs[node];
}

int MarginalTree ::GetFirstNonselfAnces(int v, const set<int> &setAnces) const {
  // find the first non-self ancestor from the list; if not found return -1
  int res = -1;

  int ncv = v;
  while (ncv >= 0) {
    // get parent
    ncv = GetParent(ncv);
    if (setAnces.find(ncv) != setAnces.end()) {
      res = ncv;
      break;
    }
  }

  return res;
}

void MarginalTree ::InitDefaultEdgeLen() {
  listEdgeDist.clear();

  // the default assume the following:
  // (a) all leaves are on the same level
  // (b) the rest of tree nodes are orgnized uniformly in distance
  for (int i = 0; i < GetTotNodesNum() - 1; ++i) {
    double distRel = GetDefaultEdgeLen(i);
    listEdgeDist.push_back(distRel);
  }
  // the root has no edge here
  listEdgeDist.push_back(0.0);
}

void MarginalTree ::InitUnitEdgelen() {
  //
  listEdgeDist.clear();

  // the default assume the following:
  // (a) all leaves are on the same level
  // (b) the rest of tree nodes are orgnized uniformly in distance
  for (int i = 0; i < GetTotNodesNum() - 1; ++i) {
    listEdgeDist.push_back(1.0);
  }
  // the root has no edge here
  listEdgeDist.push_back(0.0);
}

double MarginalTree ::GetDefaultEdgeLen(int child) {
  int curpos = child;
  int parpos = listParentNodePos[child];

  int punorm = CalcNormHeight(parpos);
  int plnorm = CalcNormHeight(curpos);
  int numLeaves = GetNumLeaves();

  if (punorm >= numLeaves) {
    punorm = numLeaves - 1;
  }
  if (plnorm >= numLeaves) {
    plnorm = numLeaves - 1;
  }
  // YW: changed back to old distance, 082306, to see if this matters
  double res =
      2.0 * (1.0 / (numLeaves - punorm) - 1.0 / (numLeaves - plnorm + 1));
  // cout << "numLeaves = " << numLeaves << ", punorm = " <<  punorm << ",
  // plnorm = " << plnorm << ", res = " << res << endl;
  // here we assume the distrbution of time is according to exponential
  // distibution of mean 2.0/k(k+1) waiting time
  return res;
}

void MarginalTree ::SetParent(int child, int par, bool fAdjLen) {
  YW_ASSERT_INFO(child < GetTotNodesNum() && par < GetTotNodesNum(),
                 "Wrong here");
  listParentNodePos[child] = par;
  // also setup height
  if (fAdjLen == true) {
    listEdgeDist[child] = GetDefaultEdgeLen(child);
  }
}

void MarginalTree ::SwapBranches(int nodeBranch1, int nodeBranch2) {
  // cout << "Swapping nodes: " << nodeBranch1 << ", " << nodeBranch2 << endl;
  // swap two branches ending at the two nodes passed in; here assume the branch
  // length will not change note: may need to reset some other descendents' info
  // after this
  int p1 = GetParent(nodeBranch1);
  int p2 = GetParent(nodeBranch2);
  SetParent(nodeBranch1, p2, false);
  SetParent(nodeBranch2, p1, false);
}

int MarginalTree ::CalcNormHeight(int node) {
  int normHt = node - (GetNumLeaves() - 1);
  if (normHt < 0) {
    normHt = 0;
  }
  return normHt;
}

void MarginalTree ::Binarize() {
  // first initialize distance if not yet
  if (listEdgeDist.size() == 0) {
    InitDefaultEdgeLen();
  }

  // assume distance has been set properly
  YW_ASSERT_INFO(listEdgeDist.size() > 0, "Tree edge length not set");

  // This function makes this marginal binary
  vector<int> updatedLabels, updatedPars;
  vector<double> updatedDist;

  // find out the current largest label, for the purpose of adding new labels
  int maxLabel = -1;
  for (int i = 0; i < (int)listNodeLabels.size(); ++i) {
    if (listNodeLabels[i] > maxLabel) {
      maxLabel = listNodeLabels[i];
    }
  }
  int labelNextToUse = maxLabel + 1;

  // before doing anything, get the descendent info for each tree node
  vector<vector<int> > listDescendentsVec;
  ConsDecedentInfo(listDescendentsVec);
  // vector< set<int> > listDescendents;
  // for( unsigned int i=0; i<listDescendents.size(); ++i )
  //{
  //    set<int> tmpSet;
  //    PopulateSetByVec( tmpSet,listDescendentsVec[i]  );
  //    listDescendents.push_back(tmpSet);
  //}

  // we need another auxilary data structure to map old position to new position
  // we need this because we are adding some new nodes between two old nodes
  vector<int> mapOldPosToNewPos(GetTotNodesNum());

  // first copy every thing up to the leaves
  for (int i = 0; i < numLeaves; ++i) {
    updatedLabels.push_back(listNodeLabels[i]);
    updatedPars.push_back(listParentNodePos[i]);
    updatedDist.push_back(listEdgeDist[i]);

    // leaf is never changed position
    mapOldPosToNewPos[i] = i;
  }
  // now we treat each internal node one by one, and split it when needed
  for (int i = numLeaves; i < GetTotNodesNum(); ++i) {
    // the first thing to do is: find children from the constructed portion of
    // tree
    vector<int> &listChildren = listDescendentsVec[i];
    // cout << "IN node = " << i << ", children num = " << listChildren.size()
    // << endl;

    // do nothing if there is no mor than 2 children
    // it is possible that an internal node does not have any children
    // Then what to do here? TBD
    if (listChildren.size() == 2 || listChildren.size() == 0) {
      // cout << "Simply go over the originals...\n";
      updatedLabels.push_back(listNodeLabels[i]);
      updatedPars.push_back(listParentNodePos[i]); // do it for now, will update
                                                   // later
      updatedDist.push_back(listEdgeDist[i]);

      // record current position
      mapOldPosToNewPos[i] = (int)updatedLabels.size() - 1;

      // now update its children's parent to this new location
      for (int jjj = 0; jjj < (int)listChildren.size(); ++jjj) {
        int oldpos = listChildren[jjj];
        int newpos = mapOldPosToNewPos[oldpos];
        updatedPars[newpos] = mapOldPosToNewPos[i];
      }

      continue;
    }
    if (listChildren.size() == 1) {
      // we should remove this node
      int childOldPos = listChildren[0];
      // skip this node, but update the node
      // let its (only) child points to its parent
      // cout << "childOldPos = " << childOldPos << endl;
      listParentNodePos[childOldPos] = listParentNodePos[i];
      // cout << "childOldPos's parent set to  = " << listParentNodePos[i] <<
      // endl;

      // also update listChildren
      if (listParentNodePos[i] >= 0) {
        int pppos = listParentNodePos[i];
        vector<int> listNewChildAtIParent;
        for (int ii = 0; ii < (int)listDescendentsVec[pppos].size(); ++ii) {
          if (i != listDescendentsVec[pppos][ii]) {
            // do not append i anymore
            listNewChildAtIParent.push_back(listDescendentsVec[pppos][ii]);
          }
        }
        //
        YW_ASSERT_INFO((int)listNewChildAtIParent.size() ==
                           (int)listDescendentsVec[pppos].size() - 1,
                       "Something wrong");
        // append a new thing
        listNewChildAtIParent.push_back(childOldPos);
        // update the orginal list
        listDescendentsVec[pppos] = listNewChildAtIParent;
      } else {
        int newpos = mapOldPosToNewPos[childOldPos];
        updatedPars[newpos] = -1;
        updatedDist[newpos] = 0.0;
      }
      continue;
    }

    // otherwise, we have to split the node
    for (int jjj = 0; jjj < (int)listChildren.size() - 2; ++jjj) {
      updatedLabels.push_back(labelNextToUse++); // new IN is assigned an
                                                 // arbitary label
      updatedPars.push_back(-1); // do it for now, will update later
      // for any new internal node, edge length (out of it) is 0
      updatedDist.push_back(0.0);

      // now update children
      int curINPos = (int)updatedLabels.size() - 1;
      if (jjj == 0) {
        // Then we use the first original child
        int oldpos = listChildren[0];
        int newpos = mapOldPosToNewPos[oldpos];
        updatedPars[newpos] = curINPos;
      } else {
        // otherwise, we use the previous IN
        updatedPars[curINPos - 1] = curINPos;
      }
      // the right branch is always an original branch
      int oldpos = listChildren[jjj + 1];
      int newpos = mapOldPosToNewPos[oldpos];
      updatedPars[newpos] = curINPos;
    }
    // now we append the original internal node in
    updatedLabels.push_back(listNodeLabels[i]);
    updatedPars.push_back(listParentNodePos[i]); // do it for now, will update
                                                 // later
    updatedDist.push_back(listEdgeDist[i]);

    // record current position
    mapOldPosToNewPos[i] = (int)updatedLabels.size() - 1;

    // update its two children, one of them is the last new node to add
    updatedPars[(int)updatedPars.size() - 2] = mapOldPosToNewPos[i];
    int oldpos = listChildren[(int)listChildren.size() - 1];
    int newpos = mapOldPosToNewPos[oldpos];
    updatedPars[newpos] = mapOldPosToNewPos[i];
  }
  // finally, we update the mtree
  this->listNodeLabels = updatedLabels;
  this->listParentNodePos = updatedPars;
  this->listEdgeDist = updatedDist;

  // check to make sure this is indeed binary
  YW_ASSERT_INFO(this->listNodeLabels.size() == this->listParentNodePos.size(),
                 "In binaralize: size wrong1");
  YW_ASSERT_INFO(this->listNodeLabels.size() == this->listEdgeDist.size(),
                 "In binaralize: size wrong1");
  // now iterator the degree
#if 0
    vector<int> nodeOutDegrees;
    for(int i=0; i<(int)this->listNodeLabels.size(); ++i)
    {
        nodeOutDegrees.push_back( 0 );
    }
    for(int i=0; i<(int)this->listNodeLabels.size(); ++i)
    {
        int ppos = listParentNodePos[i] ;
        YW_ASSERT_INFO( ppos < (int)listParentNodePos.size(), "pos wrong" );
        if( ppos >= 0 )
        {
            nodeOutDegrees[ ppos ]++;
            if( nodeOutDegrees[ ppos ] >= 3 )
            {
                YW_ASSERT_INFO( false, "Error in binarinize." );
            }
        }
    }
#endif
  // Dump();
}

void MarginalTree ::Consolidate() {
  // cout << "Before consolidate, tree = ";
  // this->Dump();
  // Remove degree-2 intermediate nodes
  // first find out which nodes are those to be removed
  set<int> nodesToDel;
  // this is very simple: scan parent list
  // if a node (non-leaf) only appears at most once of them, then remove it
  vector<int> occurTimes;
  vector<bool> nodeVisitedFlags;
  for (int i = 0; i < GetTotNodesNum(); ++i) {
    occurTimes.push_back(0);
    nodeVisitedFlags.push_back(false);
  }
  stack<int> nodesToExplore;
  for (int i = 0; i < GetNumLeaves(); ++i) {
    nodesToExplore.push(i);
  }
  while (nodesToExplore.empty() == false) {
    // find one node
    int node = nodesToExplore.top();
    nodesToExplore.pop();

    // if this is already visited, skip
    if (nodeVisitedFlags[node] == true) {
      continue;
    }
    // this is a new node, so explore it
    nodeVisitedFlags[node] = true;
    int pp = GetParent(node);
    if (pp >= 0) {
      nodesToExplore.push(pp);
      occurTimes[pp]++;
    }
  }
  // now figure out how many to remove up to a point
  vector<int> listNumDelItems;
  for (int i = 0; i < GetNumLeaves(); ++i) {
    listNumDelItems.push_back(0);
  }
  int numToDelete = 0;
  for (int i = GetNumLeaves(); i < GetTotNodesNum(); ++i) {
    if (occurTimes[i] <= 1 && i != GetTotNodesNum() - 1) {
      numToDelete++;
    }
    listNumDelItems.push_back(numToDelete);
  }

  // now store a new set of items
  vector<int> listNodeLabelsNew;
  vector<int> listParentNodePosNew;
  vector<double> listEdgeDistNew;
  // now mark those with at most once to be deleted
  for (int i = 0; i < GetTotNodesNum(); ++i) {
    // leaves and the root is always there
    if (occurTimes[i] > 1 || i < GetNumLeaves() || i == GetTotNodesNum() - 1) {
      listNodeLabelsNew.push_back(listNodeLabels[i]);

      // for parent, we trace upwards until either find a occur time > 1 or root
      double distNew = listEdgeDist[i];
      int parNew = GetParent(i);
      // now trace back to see if we need them
      while (occurTimes[parNew] <= 1 && parNew >= 0) {
        int parNext = GetParent(parNew);
        if (parNext < 0) {
          break;
        }
        distNew += listEdgeDist[parNew];
        parNew = parNext;
      }

      // save this (and make adjustment)
      int parToSet = parNew - listNumDelItems[parNew];
      if (parToSet < 0) {
        parToSet = -1;
      }
      listParentNodePosNew.push_back(parToSet);
      listEdgeDistNew.push_back(distNew);
    }
  }

  // finally store this
  listNodeLabels = listNodeLabelsNew;
  listParentNodePos = listParentNodePosNew;
  listEdgeDist = listEdgeDistNew;

  // cout << "After consolidate, tree = ";
  // this->Dump();
}

double MarginalTree ::GetEdgeLen(int childNodeIndex) const {
  YW_ASSERT_INFO(childNodeIndex < (int)listEdgeDist.size(), "List overflow");
  return listEdgeDist[childNodeIndex];
}

double MarginalTree ::GetTotEdgeLen() const {
  //
  double res = 0.0;
  for (int i = 0; i < GetTotNodesNum(); ++i) {
    if (i != GetRoot()) {
      res += GetEdgeLen(i);
    }
  }
  return res;
}

void MarginalTree ::ConsDecedentInfo(vector<vector<int> > &descNodes) const {
  descNodes.clear();
  int numNodes = GetTotNodesNum();
  // vector< vector<int> > listDescendents;
  for (int i = 0; i < numNodes; ++i) {
    vector<int> emptyVec;
    descNodes.push_back(emptyVec);
  }
  for (int i = 0; i < numNodes; ++i) {
    int parpos = listParentNodePos[i];
    if (parpos >= 0) {
      descNodes[parpos].push_back(i);
    }
  }
  // cout << "Descedents info:\n";
  // for( unsigned int i=0; i<descNodes.size(); ++i)
  //{
  // DumpIntVec( descNodes[i] );
  //}
}

void MarginalTree ::ConsAllDecedentInfo(vector<set<int> > &descNodes,
                                        bool fIncSelf) const {
  descNodes.clear();
  int numNodes = GetTotNodesNum();
  // vector< vector<int> > listDescendents;
  for (int i = 0; i < numNodes; ++i) {
    set<int> emptySet;
    descNodes.push_back(emptySet);
  }
  for (int i = 0; i < numNodes; ++i) {
    // Always contain itself if set
    if (fIncSelf == true) {
      descNodes[i].insert(i);
    }

    int parpos = listParentNodePos[i];
    if (parpos >= 0) {
      UnionSets(descNodes[parpos], descNodes[i]);
      if (fIncSelf == false) {
        // otherwise, we need to append this current node to
        descNodes[parpos].insert(i);
      }
    }
  }
}

void MarginalTree ::ConsDecedentLeavesInfo(vector<set<int> > &descLaves) const {
  descLaves.clear();
  // vector< vector<int> > listDescendents;
  int numNodes = GetTotNodesNum();
  for (int i = 0; i < numNodes; ++i) {
    set<int> emptyVec;
    descLaves.push_back(emptyVec);
  }
  for (int i = 0; i < numNodes; ++i) {
    // If this is a leave, push itself into
    if (i < numLeaves) {
      descLaves[i].insert(i);
    }

    int parpos = listParentNodePos[i];
    if (parpos >= 0) {
      UnionSets(descLaves[parpos], descLaves[i]);
    }
  }
  // cout << "Descedents info:\n";
  // for( unsigned int i=0; i<descLaves.size(); ++i)
  //{
  // DumpIntSet( descLaves[i] );
  //}
}

void MarginalTree ::ConsDecedentLeavesInfoLabels(
    vector<set<int> > &leafNodeLabels) const {
  //
  leafNodeLabels.clear();
  vector<set<int> > leafNodePos;
  ConsDecedentLeavesInfo(leafNodePos);
  for (int i = 0; i < (int)leafNodePos.size(); ++i) {
    set<int> ss;
    for (set<int>::const_iterator it = leafNodePos[i].begin();
         it != leafNodePos[i].end(); ++it) {
      ss.insert(GetLabel(*it));
    }
    leafNodeLabels.push_back(ss);
  }
}

void MarginalTree ::FindAllSplits(vector<set<int> > &listSplits) const {
  //
  listSplits.clear();
  // vector< vector<int> > listDescendents;
  int numNodes = GetTotNodesNum();
  for (int i = 0; i < numNodes; ++i) {
    set<int> emptyVec;
    listSplits.push_back(emptyVec);
  }
  for (int i = 0; i < numNodes; ++i) {
    // If this is a leave, push itself into
    if (i < numLeaves) {
      listSplits[i].insert(GetLabel(i));
    }

    int parpos = listParentNodePos[i];
    if (parpos >= 0) {
      UnionSets(listSplits[parpos], listSplits[i]);
    }
  }
}

int MarginalTree ::GetParent(int child) const {
  if (child >= GetTotNodesNum()) {
    cout << "child = " << child << ", tot num of nodes = " << GetTotNodesNum()
         << endl;
  }
  YW_ASSERT_INFO(child < GetTotNodesNum(), "Range bug");
  return listParentNodePos[child];
}

void MarginalTree ::ConsHeightsInfo(vector<int> &nodesHt) const {
  nodesHt.clear();
  int numNodes = GetTotNodesNum();
  for (int i = 0; i < numNodes; ++i) {
    nodesHt.push_back(0);
  }
  for (int i = 0; i < numNodes; ++i) {
    // test whether the parent node should be updated its height
    int parpos = listParentNodePos[i];
    if (parpos >= 0 && nodesHt[parpos] < nodesHt[i] + 1) {
      nodesHt[parpos] = nodesHt[i] + 1;
    }
  }
}

void MarginalTree ::Dump() const {
  // Output marginal tree states
  cout << "Tree: number of leaves: " << numLeaves << endl;
  cout << "Node list = ";
  DumpIntVec(this->listNodeLabels);
  cout << "Parent list = ";
  DumpIntVec(this->listParentNodePos);
  cout << "Tree dist = ";
  DumpDoubleVec(this->listEdgeDist);
}

int MarginalTree ::GetPosForLabel(int lbl) const {
  //
  int res = -1;
  for (int i = 0; i < (int)listNodeLabels.size(); ++i) {
    if (listNodeLabels[i] == lbl) {
      res = i;
      break;
    }
  }
  return res;
}

int MarginalTree ::GetMRCA(int v1, int v2) const {
  // retrieve MRCA from it
  // cout << "v1 = " << v1 << ", v2= " << v2 << endl;
  int n1 = v1, n2 = v2;
  while (n1 != n2) {
    // we alternatively move up, depend on which one is smaller
    if (n1 < n2) {
      // move n1
      n1 = GetParent(n1);
    } else {
      // move n2
      n2 = GetParent(n2);
    }
    // cout << "GetMRCA1: n1 = " << n1 << ", n2 = " << n2 << endl;
  }
  // n1 (or n2) is the result)
  return n1;
}

void MarginalTree ::GetChildren(int node, set<int> &listChildren) const {
  listChildren.clear();

  // we just search parent list to see who has entry equal to node
  for (int i = 0; i < (int)listParentNodePos.size(); ++i) {
    if (listParentNodePos[i] == node) {
      listChildren.insert(i);
    }
  }
}

int MarginalTree ::GetMaxHt() const {
  vector<int> heights;
  ConsHeightsInfo(heights);
  int maxHt = 0;
  for (int i = 0; i < (int)heights.size(); ++i) {
    if (maxHt < heights[i]) {
      maxHt = heights[i];
    }
  }
  return maxHt;
}

double MarginalTree ::GetHeight() const {
  int root = GetRoot();
  return GetHeightOfNode(root);
}
double MarginalTree ::GetHeightOfNode(int node) const {
  // get descendent
  int lchild = GetLeftDescendant(node);
  int rchild = GetRightDescendant(node);
  if (lchild < 0 || rchild < 0) {
    return 0.0;
  }
  return max(GetEdgeLen(lchild) + GetHeightOfNode(lchild),
             GetEdgeLen(rchild) + GetHeightOfNode(rchild));
}

void MarginalTree ::RemoveLeafNodeFromBinaryTree(int lfn) {
  YW_ASSERT_INFO(IsLeaf(lfn) == true, "Not a leaf");
  // rmeove a leaf node (and suppress the degree-2 node if so
  // first fill in leaves
  vector<int> listNodeLabelsNew;
  vector<int> listParentNodePosNew;
  int pp = GetParent(lfn);
  for (int i = 0; i < GetTotNodesNum(); ++i) {
    if (i != lfn && i != pp) {
      listNodeLabelsNew.push_back(this->listNodeLabels[i]);

      int parNew;
      int oldPar = GetParent(i);
      if (oldPar < pp) {
        // just minus 1
        parNew = oldPar - 1;
      } else if (oldPar > pp) {
        // otherwise, we lost two
        parNew = oldPar - 2;
      } else {
        // In this case, we are pointing to pp, since pp is removed, we need to
        // move up by one
        parNew = GetParent(pp) - 2;
      }
      if (parNew < 0) {
        parNew = -1;
      }
      listParentNodePosNew.push_back(parNew);
    }
  }
  //
  this->listNodeLabels = listNodeLabelsNew;
  this->listParentNodePos = listParentNodePosNew;

  this->numLeaves--;
}

bool MarginalTree ::AreTwoPathsDisjoint(int sn1, int en1, int sn2,
                                        int en2) const {
  // test whether two path (sn1, en1) and (sn2, en2) are (vertex) disjoint
  // note that for binary tree, this is also checking for edge disjoint
  // we use a dumb method here
  set<int> nodesVisitedTree1;

  int n1 = sn1, n2 = en1;
  nodesVisitedTree1.insert(n1);
  nodesVisitedTree1.insert(n2);
  while (n1 != n2) {
    // we alternatively move up, depend on which one is smaller
    int nodeNew;
    if (n1 < n2) {
      // move n1
      n1 = GetParent(n1);
      nodeNew = n1;
    } else {
      // move n2
      n2 = GetParent(n2);
      nodeNew = n2;
    }

    //
    nodesVisitedTree1.insert(nodeNew);
  }
  // cout << "Path 1=";
  // DumpIntSet( nodesVisitedTree1 );
  // now we move on to the next pair
  n1 = sn2;
  n2 = en2;
  if (nodesVisitedTree1.find(n1) != nodesVisitedTree1.end() ||
      nodesVisitedTree1.find(n2) != nodesVisitedTree1.end()) {
    return false;
  }
  while (n1 != n2) {
    // we alternatively move up, depend on which one is smaller
    int nodeNew;
    if (n1 < n2) {
      // move n1
      n1 = GetParent(n1);
      nodeNew = n1;
    } else {
      // move n2
      n2 = GetParent(n2);
      nodeNew = n2;
    }

    //
    if (nodesVisitedTree1.find(nodeNew) != nodesVisitedTree1.end()) {
      return false;
    }
  }

  return true;
}

int MarginalTree ::GetPath(int sn, int en, set<int> &edgesOnPath) const {
  // find edges on the path, and return the MRCA
  int n1 = sn, n2 = en;
  edgesOnPath.insert(n1);
  edgesOnPath.insert(n2);
  while (n1 != n2) {
    // we alternatively move up, depend on which one is smaller
    int nodeNew;
    if (n1 < n2) {
      // move n1
      n1 = GetParent(n1);
      nodeNew = n1;
    } else {
      // move n2
      n2 = GetParent(n2);
      nodeNew = n2;
    }

    //
    edgesOnPath.insert(nodeNew);
  }
  // remove MRCA from result
  YW_ASSERT_INFO(edgesOnPath.find(n1) != edgesOnPath.end(), "wrong2");
  edgesOnPath.erase(n1);

  return n1;
}

double MarginalTree ::GetPathLen(int sn, int en) {
  // get the branch lenggth on the path
  double res = 0.0;

  set<int> edgesOnPath;
  int mrca = GetPath(sn, en, edgesOnPath);
  YW_ASSERT_INFO(edgesOnPath.find(mrca) == edgesOnPath.end(), "Fail to find");
  for (set<int>::iterator it = edgesOnPath.begin(); it != edgesOnPath.end();
       ++it) {
    res += GetEdgeLen(*it);
  }
  return res;
}

void MarginalTree ::OutputGML(const char *fileName) const {
  // Now output a file in GML format
  // First create a new name
  string name = fileName;
  // cout << "num edges = " << listEdges.size() << endl;

  DEBUG("FileName=");
  DEBUG(name);
  DEBUG("\n");
  // Now open file to write out
  ofstream outFile(name.c_str());

  // First output some header info
  outFile << "graph [\n";
  outFile << "comment ";
  OutputQuotedString(outFile, "Automatically generated by Graphing tool");
  outFile << "\ndirected  1\n";
  outFile << "id  1\n";
  outFile << "label ";
  OutputQuotedString(outFile, "Marginal Tree....\n");

  // Now output all the vertices
  //	int i;

  // cout << "a.1.1\n";
  for (int i = 0; i < (int)listNodeLabels.size(); ++i) {
    outFile << "node [\n";

    outFile << "id " << i << endl;
    outFile << "label ";
    char buf[80];
    //        sprintf(buf, "n%d",  listNodeLabels[i]  );
    sprintf(buf, "n%d", i);

    OutputQuotedString(outFile, buf);
    outFile << endl;

    // See if we need special shape here
    outFile << "defaultAtrribute   1\n";

    outFile << "]\n";
  }
  // cout << "a.1.3\n";

  // Now output all the edges, by again starting from root and output all nodes
  for (int i = 0; i < (int)listParentNodePos.size(); ++i) {
    int parpos = listParentNodePos[i];

    // cout << "Output an edge \n";
    outFile << "edge [\n";
    outFile << "source " << parpos << endl;
    outFile << "target  " << i << endl;
    outFile << "label ";
    OutputQuotedString(outFile, "");
    outFile << "\n";
    outFile << "]\n";
  }

  // Finally quite after closing file
  outFile << "\n]\n";
  outFile.close();
}

string MarginalTree ::GetNewick() const {
  // return the newick format of the tree (with length)
  // method: just get the newick at the root node
  return GetNewickAt(GetTotNodesNum() - 1);
}
string MarginalTree ::GetNewickSorted(bool fLen) const {
  //
  return GetNewickAt(GetTotNodesNum() - 1, true, fLen);
}

string MarginalTree ::GetNewickAt(int node, bool fSort, bool fLen) const {
  // find its descendents
  string res;
  int childLeft = GetLeftDescendant(node);
  int childRight = GetRightDescendant(node);
  if (childLeft < 0) {
    // must be leaf
    YW_ASSERT_INFO(IsLeaf(node) == true, "Wrong node in MT");
    // for leaf, only ouput its label together with its length
    char buf[100];
    if (fLen == true) {
      sprintf(buf, "%d:%f", GetLabel(node), GetEdgeLen(node));
    } else {
      sprintf(buf, "%d", GetLabel(node));
    }
    res = buf;
  } else {
    // append two children's
    if (childRight < 0) {
      Dump();
    }
    YW_ASSERT_INFO(childRight >= 0, "Left/right mismatch");
    res = "(";
    // res += GetNewickAt(childLeft);
    // res +=",";
    // res += GetNewickAt(childRight);
    string strPart1 = GetNewickAt(childLeft, fSort, fLen);
    string strPart2 = GetNewickAt(childRight, fSort, fLen);
    string strToAdd;
    if (fSort == false || strPart1 <= strPart2) {
      res += strPart1;
      res += ",";
      res += strPart2;
    } else {
      res += strPart2;
      res += ",";
      res += strPart1;
    }
    res += strToAdd;
    res += ")";
    if (fLen == true && node < GetTotNodesNum() - 1) {
      char buf[100];
      sprintf(buf, ":%f", GetEdgeLen(node));
      res += buf;
    }
  }
  return res;
}

void MarginalTree ::GetLeavesUnder(int nn, set<int> &leavesUnder) const {
  //
  if (IsLeaf(nn) == true) {
    leavesUnder.insert(nn);
  } else {
    set<int> listChildren;
    GetChildren(nn, listChildren);
    for (set<int>::iterator it = listChildren.begin(); it != listChildren.end();
         ++it) {
      GetLeavesUnder(*it, leavesUnder);
    }
  }
}

void MarginalTree ::GetlabelsFor(const set<int> &setPos,
                                 set<int> &setLbls) const {
  //
  setLbls.clear();
  for (set<int>::const_iterator it = setPos.begin(); it != setPos.end(); ++it) {
    setLbls.insert(GetLabel(*it));
  }
}

void MarginalTree ::GetLeafSetsForCuts(const vector<int> &listCuts,
                                       vector<set<int> > &listLeafSets) const {
  // this function finds the cutted subtrees' leaf sets for the given set of cut
  // edges
  listLeafSets.clear();

  // we first create a map of whether an edge mutate or not
  vector<bool> mapEdgeMutFlags;
  for (int i = 0; i < this->GetTotNodesNum(); ++i) {
    mapEdgeMutFlags.push_back(false);
  }
  for (int i = 0; i < (int)listCuts.size(); ++i) {
    mapEdgeMutFlags[listCuts[i]] = true;
  }

  // we start by bottom up way to traversal all nodes
  vector<set<int> > nodesLeaves(this->GetTotNodesNum());
  for (int i = 0; i < this->GetNumLeaves(); ++i) {
    // all leave nodes are trivial
    nodesLeaves[i].insert(i);
  }
  // test for all nodes
  for (int i = 0; i < this->GetTotNodesNum(); ++i) {
    // if the edge is cut, we have found an partition or it is a root
    if (mapEdgeMutFlags[i] == true || i == this->GetTotNodesNum() - 1) {
      if (nodesLeaves[i].size() > 0) {
        // cout << "Found one partition: ";
        // DumpIntSet( nodesLeaves[i] );
        listLeafSets.push_back(nodesLeaves[i]);
      }
    } else {
      // otherwise propagate to above
      UnionSets(nodesLeaves[this->GetParent(i)], nodesLeaves[i]);
    }
  }
}

int MarginalTree ::GetMRCAForNodes(const set<int> &listNodes) const {
  // find mrca of a list of nodes
  // we use a priority queue, each time, we try to find
  priority_queue<int> queueNodesToCheck;
  set<int> nodesVisited;

  for (set<int>::iterator it = listNodes.begin(); it != listNodes.end(); ++it) {
    queueNodesToCheck.push((*it) * (-1));
  }
  while (queueNodesToCheck.size() > 1) {
    int curn = -queueNodesToCheck.top();
    queueNodesToCheck.pop();

    // in case there are duplicate ones, remove these duplicate copies
    // this can happen if one node is another node's parent
    if (-queueNodesToCheck.top() == curn) {
      // don't work on this, wait for the next one
      continue;
    }

    // is this visited
    int pp = this->GetParent(curn);
    // cout << "Processing curn: " << curn << ", parent: " << pp << endl;
    if (nodesVisited.find(pp) == nodesVisited.end()) {
      // new node
      nodesVisited.insert(pp);
      // push to queue
      queueNodesToCheck.push(-1 * pp);
    }
  }
  int res = -1 * queueNodesToCheck.top();
  return res;
}

bool MarginalTree ::IsNodeUnder(int nn, int ancesNode) const {
  //
  if (nn > ancesNode) {
    return false;
  }
  int curn = nn;
  while (curn < ancesNode && curn >= 0) {
    curn = this->GetParent(curn);
  }
  if (curn == ancesNode) {
    return true;
  } else {
    return false;
  }
}

void MarginalTree ::RandPermuateLeaves() {
  // randomly permuate the leaves of the tree
  // we do this by shuffeling the parent of the leaves
  vector<int> parentsNewIndices;
  GetRandVector(parentsNewIndices, 0, GetNumLeaves() - 1);
  // cout << "Dump Random vector: ";
  // DumpIntVec( parentsNewIndices );
  // now shuffling it
  vector<int> leavesParNew;
  for (int i = 0; i < (int)parentsNewIndices.size(); ++i) {
    leavesParNew.push_back(GetParent(parentsNewIndices[i]));
  }
  // now assign it
  for (int i = 0; i < (int)parentsNewIndices.size(); ++i) {
    SetParent(i, leavesParNew[i]);
  }
}

int MarginalTree ::GetTriple(int i, int j, int k) const {
  // ensure order of a,b,c first
  OrderInt(i, j);
  OrderInt(i, k);
  OrderInt(j, k);

  //
  // is these have different triples on T1 and T2?
  // we do this by getting MRCA for all pairs of MRCAs
  int mrcaij1 = GetMRCA(i, j);
  int mrcajk1 = GetMRCA(j, k);
  int mrcaik1 = GetMRCA(i, k);

  // now just test exhustively
  if (mrcaij1 == mrcajk1) {
    return 3;
  } else if (mrcaij1 == mrcaik1) {
    return 2;
  } else {
    return 1;
  }
}

int MarginalTree ::GetSibling(int a) const {
  // get sibling of the node (leaf or non-leaf)
  int par = GetParent(a);
  int lc = GetLeftDescendant(par);
  int rc = GetRightDescendant(par);

  YW_ASSERT_INFO(a == lc || a == rc, "Very wrong");
  if (a == lc) {
    return rc;
  } else {
    return lc;
  }
}

bool MarginalTree ::AreNodesSibling(int a, int b) const {
  //
  return GetSibling(a) == b;
}

void MarginalTree ::SortByLeafId() {
  // sort based on leaf id. That is, leaf ids = 0,1,2,3,.. in the list
  vector<int> listNodeLabelsNew(this->listNodeLabels.size());
  vector<int> listParentNodePosNew(this->listParentNodePos.size());
  vector<double> listEdgeDistNew(this->listEdgeDist.size());

  listNodeLabelsNew = listNodeLabels;
  listParentNodePosNew = listParentNodePos;
  listEdgeDistNew = listEdgeDist;

  // now sort and swap the leaf part
  // collect leaves
  vector<int> listLeafIds;
  for (int i = 0; i < GetNumLeaves(); ++i) {
    listLeafIds.push_back(listNodeLabels[i]);
  }
  // vector<int> listLeafIdsOld = listLeafIds;
  SortIntVec(listLeafIds);
  // cout << "listLeafIds = ";
  // DumpIntVec( listLeafIds );
  // create a map
  // map<int,int> mapLeafIdToOldPos;
  // for( int i=0; i<(int)listLeafIdsOld.size(); ++i )
  //{
  //    mapLeafIdToOldPos.insert(map<int,int> :: value_type(listLeafIdsOld[i],i)
  //    );
  //}
  map<int, int> mapLeafIdToNewPos;
  for (int i = 0; i < (int)listLeafIds.size(); ++i) {
    // cout << "Set map from id " << listLeafIds[i] << " to position " << i <<
    // endl;
    mapLeafIdToNewPos.insert(map<int, int>::value_type(listLeafIds[i], i));
  }
  // now swap the info in each in the old list
  for (int i = 0; i < (int)GetNumLeaves(); ++i) {
    int vid = listNodeLabels[i];
    YW_ASSERT_INFO(mapLeafIdToNewPos.find(vid) != mapLeafIdToNewPos.end(),
                   "FAIL to find");
    int posNew = mapLeafIdToNewPos[vid];
    // cout << "vid = " << vid << ", Set " << posNew << " to position " << i <<
    // endl;
    listNodeLabelsNew[posNew] = vid;
    listParentNodePosNew[posNew] = listParentNodePos[i];
    listEdgeDistNew[posNew] = listEdgeDist[i];
  }

#if 0 // there is some issues with this piece of code: namely, it can not deal
      // with non-distinct id in trees properly. Although ids are expected to be
      // distinct but sometime they are not; so change it on 8/13/13
      // list leaf in order
	vector<int> listLeaves = this->listNodeLabels;
	SortIntVec(listLeaves);
cout << "after sorting, leaf list = ";
DumpIntVec( listLeaves);
	for(int i=0; i<(int)this->listNodeLabels.size(); ++i)
	{
		// keep a sorted list
		int lvid = this->listNodeLabels[i];
		// find which leaf is it
		int lvpos = GetItemIndexInVec(listLeaves, lvid);
		YW_ASSERT_INFO(lvpos >= 0, "FATAL ERROR10");
		listNodeLabelsNew[ lvpos ] = lvid;
		listParentNodePosNew[ lvpos ] = listParentNodePos[i];
		if( lvpos < (int) listEdgeDistNew.size() )
		{
			listEdgeDistNew[lvpos] = listEdgeDist[i];
		}
	}
#endif
  // now write back
  this->listNodeLabels = listNodeLabelsNew;
  this->listParentNodePos = listParentNodePosNew;
  this->listEdgeDist = listEdgeDistNew;

  // redo the descendents
  BuildDescendantInfo();
}

void MarginalTree ::FixDupIds() {
  // remove redundent ids with something new
  // sort based on leaf id. That is, leaf ids = 0,1,2,3,.. in the list
  // also, keep the leaf and internal nodes id separated
  vector<int> listNodeLabelsNew(this->listNodeLabels.size());
  int numLeaves = GetNumLeaves();

  // list leaf in order
  set<int> setNids;
  PopulateSetByVec(setNids, this->listNodeLabels);
  int idNext = *(setNids.rbegin()) + 1;

  set<int> idsSeenBefore;

  for (int i = 0; i < (int)this->listNodeLabels.size(); ++i) {
    // keep a sorted list
    int lvid = this->listNodeLabels[i];
    if (idsSeenBefore.find(lvid) != idsSeenBefore.end()) {
      lvid = idNext++;
    }
    listNodeLabelsNew[i] = lvid;
    idsSeenBefore.insert(lvid);
  }

  // now inc the id of the internal nodes
  for (int i = numLeaves; i < (int)listNodeLabelsNew.size(); ++i) {
    listNodeLabelsNew[i] += 3 * numLeaves;
  }

  // now write back
  this->listNodeLabels = listNodeLabelsNew;
}

void MarginalTree ::RearrangeParIncOrder() {
  // cout << "--RearrangeParIncOrder:\n";
  // sometimes the parent position is out of order, say 1,3,3,2,2,1,...
  // we can rearrange the internal node so that it becomes 1,2,2,3,3,1...
  // check the order of the appreance of the parent node
  // CAUTION: after this, need to perform descendent list rebuilt

  //#if 0
  int curParOrderIndex = GetNumLeaves();
  map<int, int> mapCurParPosToNewParPos;
  set<int> setSeePars;
  queue<int> nodesToProc;
  // add in the leaves first
  // vector<int> parposListNew(listParentNodePos.size() );
  // // the new par and dist list parposListNew[ parposListNew.size()-1 ] = -1;
  // // the last one is always -1 vector<double> distListNew(listEdgeDist.size()
  // );
  for (int i = 0; i < GetNumLeaves(); ++i) {
    nodesToProc.push(i);
  }
  vector<int> listNewParsInOrder;
  while (nodesToProc.empty() == false) {
    int nodeCur = nodesToProc.front();
    nodesToProc.pop();
    int parpos = GetParent(nodeCur);
    if (parpos < 0) {
      // root, do nothing
      continue;
    }
    // cout << "nodecur = " << nodeCur << ", parpos = " << parpos << endl;
    if (setSeePars.find(parpos) == setSeePars.end()) {
      // add to set and move to next
      setSeePars.insert(parpos);
    } else {
      // have seen before, so record the mapping
      YW_ASSERT_INFO(mapCurParPosToNewParPos.find(parpos) ==
                         mapCurParPosToNewParPos.end(),
                     "Should not be here");
      ;
      mapCurParPosToNewParPos.insert(
          map<int, int>::value_type(parpos, curParOrderIndex));
      // cout << "map old pos " << parpos << " to " << curParOrderIndex << endl;
      ++curParOrderIndex;
      // when a parent node is done, process it
      nodesToProc.push(parpos);

      listNewParsInOrder.push_back(parpos);
    }
  }

  // now swap the par positions
  vector<int> parposListNew;
  vector<double> distListNew = listEdgeDist;
  for (int ii = 0; ii < (int)GetNumLeaves(); ++ii) {
    YW_ASSERT_INFO(mapCurParPosToNewParPos.find(listParentNodePos[ii]) !=
                       mapCurParPosToNewParPos.end(),
                   "False");
    parposListNew.push_back(mapCurParPosToNewParPos[listParentNodePos[ii]]);
  }
  // then output the internal node in the given order
  for (int ii = 0; ii < (int)listNewParsInOrder.size(); ++ii) {
    int nindex = listNewParsInOrder[ii];
    if (mapCurParPosToNewParPos.find(listParentNodePos[nindex]) !=
        mapCurParPosToNewParPos.end()) {
      parposListNew.push_back(
          mapCurParPosToNewParPos[listParentNodePos[nindex]]);
      // cout << "set dist of edge " << nindex << " (old dist " <<
      // distListNew[nindex] << " to node " << mapCurParPosToNewParPos[ nindex ]
      // << " w/ dist "; cout << listEdgeDist[ mapCurParPosToNewParPos[ nindex ]
      // ] << endl;
      distListNew[nindex] = listEdgeDist[mapCurParPosToNewParPos[nindex]];
    } else {
      parposListNew.push_back(-1);
    }
  }

  // finally set up the new lists
  this->listParentNodePos = parposListNew;
  this->listEdgeDist = distListNew;
  //#endif
}

string MarginalTree ::GetNewickNoBrLen() const {
  // get the newick format w/o branch length
  string strCurr = this->GetNewick();
  PhylogenyTreeBasic trPhy;
  trPhy.ConsOnNewick(strCurr);
  trPhy.Order();
  string res;
  trPhy.ConsNewick(res);
  return res;
}

string MarginalTree ::GetNewickNoBrLen2() const {
  //
  return GetNewickAt(GetTotNodesNum() - 1, true, false);
}

void MarginalTree ::RemapLeafLabels(const map<int, int> &mapLeafLblsToNew) {
#if 0
cout << "RemapLeafLabels: ";
this->Dump();
cout << "mapLeafLblsToNew: ";
for(map<int,int> :: const_iterator it = mapLeafLblsToNew.begin(); it != mapLeafLblsToNew.end(); ++it)
{
cout << "[" << it->first << "," << it->second << "]  ";
}
cout << endl;
this->Dump();
#endif
  // convert each existing labels to consecutive labels e.g. 0, 1, 2, ...
  for (int i = 0; i < (int)listNodeLabels.size(); ++i) {
    int lblCur = listNodeLabels[i];
    // cout << "lblCur: " << lblCur << endl;
    YW_ASSERT_INFO(lblCur < 0 ||
                       mapLeafLblsToNew.find(lblCur) != mapLeafLblsToNew.end(),
                   "Fail to find123");
    if (lblCur >= 0) {
      listNodeLabels[i] = (*(mapLeafLblsToNew.find(lblCur))).second;
    }
  }
  // rebuild descendent info
  BuildDescendantInfo();
}

void MarginalTree ::MapLeafLblConsecutiveOrder(vector<int> &listLeafLblsOld) {
  listLeafLblsOld.clear();
  int idNext = 0;
  MapLeafLblConsecutiveOrderAt(this->GetRoot(), idNext, listLeafLblsOld);
  // adding the remaining internal nodes
  for (int i = GetNumLeaves(); i < GetTotNodesNum(); ++i) {
    listLeafLblsOld.push_back(GetLabel(i));
    SetLabel(i, idNext);
    ++idNext;
  }
}

void MarginalTree ::MapLeafLblConsecutiveOrderAt(int rootST, int &idNext,
                                                 vector<int> &listLeafLblsOld) {
  if (IsLeaf(rootST)) {
    listLeafLblsOld.push_back(GetLabel(rootST));
    SetLabel(rootST, idNext);
    ++idNext;
  } else {
    MapLeafLblConsecutiveOrderAt(GetLeftDescendant(rootST), idNext,
                                 listLeafLblsOld);
    MapLeafLblConsecutiveOrderAt(GetRightDescendant(rootST), idNext,
                                 listLeafLblsOld);
  }
}

void MarginalTree ::ResetIncLabel() {
  //
  for (int i = 0; i < GetNumLeaves(); ++i) {
    listNodeLabels[i] = i;
  }
}

void MarginalTree ::IncLabels() {
  for (int i = 0; i < GetNumLeaves(); ++i) {
    ++listNodeLabels[i];
  }
}

void MarginalTree ::FindSibLeafPairs(
    vector<pair<int, int> > &listSibPairs) const {
  // cout << "FindSibLeafPairs:\n";
  // Dump();
  // find leaves that are siblings (return the index (note not label) of the sib
  // pairs)
  for (int i = GetNumLeaves(); i < GetTotNodesNum(); ++i) {
    //
    int nvleft = GetLeftDescendant(i);
    int nvRight = GetRightDescendant(i);
    if (IsLeaf(nvleft) == true && IsLeaf(nvRight) == true) {
      pair<int, int> pp(nvleft, nvRight);
      listSibPairs.push_back(pp);
    }
  }
  YW_ASSERT_INFO(listSibPairs.size() > 0, "Must have at least one pair");
}

void MarginalTree ::MakeLeafSubtreeOfTwo(int posLeaf, int lblChild1,
                                         int lblChild2, double len1,
                                         double len2) {
  // cout << "MakeLeafSubtreeOfTwo: posLeaf: " << posLeaf << ", child1:" <<
  // lblChild1 << ", child2:" << lblChild2 << ", len1:" << len1 << ", len2:" <<
  // len2 << endl;
  // add two new leaves below a leaf (here, the two new leaves are located at
  // the end of leaves; and the new internal (original elaf) is right next to
  // these new leaves) also clean up the tree a bit (set labels of internal
  // nodes to be -1)
  vector<int> listNodeLabelsNew;
  vector<int> listParentNodePosNew;
  vector<double> listEdgeDistNew;

  //
  for (int i = 0; i < GetNumLeaves(); ++i) {
    if (i != posLeaf) {
      listNodeLabelsNew.push_back(GetLabel(i));
      listParentNodePosNew.push_back(GetParent(i) + 2);
      listEdgeDistNew.push_back(GetEdgeLen(i));
    }
  }
  // add the two new leaves
  listNodeLabelsNew.push_back(lblChild1);
  listNodeLabelsNew.push_back(lblChild2);
  int posCur = (int)listNodeLabelsNew.size();
  listNodeLabelsNew.push_back(-1);
  listParentNodePosNew.push_back(posCur);
  listParentNodePosNew.push_back(posCur);
  listParentNodePosNew.push_back(GetParent(posLeaf) + 2);
  listEdgeDistNew.push_back(len1);
  listEdgeDistNew.push_back(len2);
  listEdgeDistNew.push_back(GetEdgeLen(posLeaf));
  for (int i = GetNumLeaves(); i < GetTotNodesNum(); ++i) {
    listNodeLabelsNew.push_back(-1);
    if (GetParent(i) >= 0) {
      listParentNodePosNew.push_back(GetParent(i) + 2);
    } else {
      listParentNodePosNew.push_back(-1);
    }
    listEdgeDistNew.push_back(GetEdgeLen(i));
  }

  // now update the info
  ++this->numLeaves;
  this->listNodeLabels = listNodeLabelsNew;
  this->listParentNodePos = listParentNodePosNew;
  this->listEdgeDist = listEdgeDistNew;
  listLeftDescs.clear();
  listRightDescs.clear();
  BuildDescendantInfo();
}

void MarginalTree ::GetLabelListForLeaf(vector<int> &listLbls) const {
  //
  listLbls.clear();
  for (int i = 0; i < GetNumLeaves(); ++i) {
    listLbls.push_back(GetLabel(i));
  }
}

void MarginalTree ::FindDiffSubtreesFrom(const MarginalTree &mtreeRef,
                                         set<int> &setDiffBrs,
                                         set<int> &setDiffRefMissed) const {
  // find all branches (subtrees below them) that are not in the reference tree
  // setDiffBrs: in this tree but not in reference tree
  // setDiffRefMissed: in reference tree but not in this tree
  vector<set<int> > listSubtreesRef;
  mtreeRef.ConsDecedentLeavesInfoLabels(listSubtreesRef);
  vector<set<int> > listSubtreesThis;
  ConsDecedentLeavesInfoLabels(listSubtreesThis);
  set<set<int> > setSubtreesRef;
  PopulateSetByVecGen(setSubtreesRef, listSubtreesRef);
  set<set<int> > setSubtreesThis;
  PopulateSetByVecGen(setSubtreesThis, listSubtreesThis);
  //
  setDiffBrs.clear();
  for (int i = 0; i < (int)listSubtreesThis.size(); ++i) {
    if (setSubtreesRef.find(listSubtreesThis[i]) == setSubtreesRef.end()) {
      setDiffBrs.insert(i);

#if 0
            // alsoinsert any ancestral edge into it
            set<int> nodesAnces;
            GetPath( i, GetRoot(), nodesAnces );
            UnionSets( setDiffBrs, nodesAnces );

            // and also all the silbings
            for( set<int> :: iterator itg = nodesAnces.begin(); itg != nodesAnces.end(); ++itg )
            {
                if( IsLeaf(*itg) == false )
                {
                    setDiffBrs.insert( GetLeftDescendant(*itg) );
                    setDiffBrs.insert( GetRightDescendant(*itg) );
                }
            }
#endif
    }
  }
  setDiffRefMissed.clear();
  for (int i = 0; i < (int)listSubtreesRef.size(); ++i) {
    if (setSubtreesThis.find(listSubtreesRef[i]) == setSubtreesThis.end()) {
      setDiffRefMissed.insert(i);
    }
  }
}

bool MarginalTree ::IsOutgroup(int lvid) const {
  // cout << "IsOutgroup: lvid = " << lvid << ", tree is: ";
  // Dump();
  int rtn = GetRoot();
  // check two children of root
  int lc = GetLeftDescendant(rtn);
  if (IsLeaf(lc)) {
    if (GetLabel(lc) == lvid) {
      // cout << "good OG\n";
      return true;
    }
  }
  int rc = GetRightDescendant(rtn);
  if (IsLeaf(rc)) {
    if (GetLabel(rc) == lvid) {
      // cout << "good OG\n";
      return true;
    }
  }
  // cout << "BAD OG\n";
  return false;
}
