#include "PhylogenyTreeBasic.h"
#include "Utils3.h"
#include "Utils4.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stack>

// ***************************************************************************
// The following code is largely based on Gusfield's 1991 Paper
// ***************************************************************************
extern void OutputQuotedString(ofstream &outFile, const char *buf);

string GetStringFromId(int id) {
  char buf[100];
  sprintf(buf, "%d", id);
  return buf;
}

int GetNewickNumLeaves(const string &strNewick, char chSepLeft, char chSepRight,
                       char midSep) {
  // the number of leaves of newick is equal to the number of separator char
  // plus one
  int res = 0;
  bool fCount = false; // only count when seeing left sep (
  for (int i = 0; i < (int)strNewick.length(); ++i) {
    if (strNewick[i] == chSepLeft) {
      fCount = true;
    } else if (strNewick[i] == chSepRight) {
      if (fCount == true) {
        // add the last one
        res++;
      }

      fCount = false;
    } else if (strNewick[i] == midSep) {
      if (fCount == true) {
        res++;
      } else {
        fCount = true;
      }
    }
  }
  return res;
}

bool GetTripleType(TreeNode *pn1, TreeNode *pn2, TreeNode *pn3,
                   pair<pair<TreeNode *, TreeNode *>, TreeNode *> &triple) {
  TreeNode *pmrca12 = pn1->GetMRCA(pn2);
  TreeNode *pmrca13 = pn1->GetMRCA(pn3);
  TreeNode *pmrca23 = pn2->GetMRCA(pn3);
  //
  int dummy;
  if (pmrca13 != pmrca12) {
    if (pmrca13->IsAncesterOf(pmrca12, dummy) == true) {
      triple.first.first = pn1;
      triple.first.second = pn2;
      triple.second = pn3;
      return true;
    } else if (pmrca12->IsAncesterOf(pmrca13, dummy) == true) {
      triple.first.first = pn1;
      triple.first.second = pn3;
      triple.second = pn2;
      return true;
    } else {
      YW_ASSERT_INFO(false, "Impossible");
    }
  }
  // if( pmrca23 != pmrca12 &&  pmrca12->IsAncesterOf(pmrca23, dummy) == true )
  else if (pmrca23 != pmrca12) {
    triple.first.first = pn1;
    triple.first.second = pn2;
    triple.second = pn3;
    return true;
  }
  // triple not found
  return false;
}

// different from Marginal tree, we allow mulfurcating trees here
// this can be convenient in some cases
bool ReadinPhyloTreesNewick(ifstream &inFile, int numLeaves,
                            vector<PhylogenyTreeBasic *> &treePtrList,
                            TaxaMapper *pTMapper) {
  // NOTE: RETURN TRUE IF NO LABEL ADJUSTMENT IS DONE
  // RETURN FALSE IF WE SWITCHED LABEL BY DECREASING BY ONE
  // figure out leave num
  bool fNoChange = true;
  int nLvs = numLeaves;

  // read marginal trees in newick format
  // here there is no preamble, one line per tree
  while (inFile.eof() == false) {
    // ensure the first char is '('; otherwise stop
    char ch;
    inFile >> ch;
    inFile.putback(ch);
    if (ch != '(') {
      break;
    }

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
#if 0
for(multiset<string> :: iterator it22 = setLabels.begin(); it22 != setLabels.end(); ++it22)
{
cout << "Label found: " << *it22 << endl;
}
#endif
    nLvs = setLabels.size();
    //#endif
    //
    PhylogenyTreeBasic *pphTree = new PhylogenyTreeBasic;
    // if( fDup == false )
    //{
    pphTree->ConsOnNewick(treeNewick, -1, false, pTMapper);
    // cout << "Done phylogenetic tree construction...\n";
    // pphTree->OutputGML("tmp.gml");
    //}
    // else
    //{
    //	phTree.ConsOnNewickDupLabels(treeNewick, pTMapper);
    //}

    if (pTMapper != NULL) {
      pTMapper->SetInitialized(true);
    }
    // string strTr;
    // pphTree->ConsNewick(strTr);
    // cout << "After reconstruction: strTr = " << strTr << endl;
    // see if zero is in, if not, must have 1 and decrease by 1
    set<int> lvids;
    pphTree->GetLeaveIds(lvids);
    // cout << "lvids : ";
    // DumpIntSet( lvids );
    int idInternal = lvids.size();
    YW_ASSERT_INFO(lvids.find(0) != lvids.end(),
                   "Must adjust leaf label first (to start with 0)");

    //	YW_ASSERT_INFO( lvids.find(1) != lvids.end(), "Wrong" );

    // decrease by one
    PhylogenyTreeIterator itorTree(*pphTree);
    itorTree.Init();
    // pphTree->InitPostorderWalk();
    while (itorTree.IsDone() == false) {
      //				TreeNode *pn =
      // pphTree->NextPostorderWalk( ) ;
      TreeNode *pn = itorTree.GetCurrNode();
      itorTree.Next();
      if (pn == NULL) {
        break; // done with all nodes
      }
      if (pn->IsLeaf() == false) {
        pn->SetID(idInternal++);
      }
    }

    // mark the change
    //	fNoChange = false;

    vector<int> nidsList, nparsList;
    pphTree->GetNodeParInfo(nidsList, nparsList);
    // phTree.GetNodeParInfoNew(nidsList, nparsList);
    // phTree.GetNodeParInfo(nidsList, nparsList);
    // if( nLvs <= 0 )
    //{
    // string strTrNW;
    // pphTree->ConsNewick(strTrNW);
    // cout << "strTrNW: " << strTrNW << endl;
    treePtrList.push_back(pphTree);
    // cout << "Newick format of this marginal tree: ";
    // cout << tree.GetNewick() << endl;
  }
  return fNoChange;
}

// create a random tree
void InitRandomTree(PhylogenyTreeBasic &treeToInit, int numTaxa, int rndSeed) {
  //
  if (rndSeed >= 0) {
    InitRandom(rndSeed);
  }
  // create leaves first
  int idToUseNext = 0;
  vector<TreeNode *> listActiveNodes;
  for (int i = 0; i < numTaxa; ++i) {
    //
    TreeNode *pLeaf = new TreeNode(idToUseNext++);
    // label it
    pLeaf->SetLabel(GetStringFromId(i));
    listActiveNodes.push_back(pLeaf);
  }
  // now create random coalescence
  while (listActiveNodes.size() > 1) {
    // get two random nodes and coalesce them
    int rndpos1 = (int)(listActiveNodes.size() * GetRandFraction());
    YW_ASSERT_INFO(rndpos1 < (int)listActiveNodes.size(), "overflow");
    TreeNode *node1 = listActiveNodes[rndpos1];
    RemoveVecElementAt(listActiveNodes, rndpos1);
    int rndpos2 = (int)(listActiveNodes.size() * GetRandFraction());
    YW_ASSERT_INFO(rndpos2 < (int)listActiveNodes.size(), "overflow");
    TreeNode *node2 = listActiveNodes[rndpos2];
    RemoveVecElementAt(listActiveNodes, rndpos2);
    //
    TreeNode *pnodeNew = new TreeNode(idToUseNext++);
    vector<int> listEmpty;
    pnodeNew->AddChild(node1, listEmpty);
    pnodeNew->AddChild(node2, listEmpty);
    // add this node to list of active nodes
    listActiveNodes.push_back(pnodeNew);
  }
  // now here is the root
  YW_ASSERT_INFO(listActiveNodes.size() == 1, "Only one root");
  treeToInit.SetRoot(listActiveNodes[0]);
}

void CreatePhyTreeWithRootedSplits(PhylogenyTreeBasic &treeToProc, int numTaxa,
                                   const set<set<int> > &setGivenSplits) {
  // create a phy tree with the given rooted splits
  // ASSUME: taxa starts from 0 to numTaxa-1
  // result can be a non-binary tree
  // first order them
  vector<set<set<int> > > listGivenSplits(numTaxa + 1);
  for (set<set<int> >::const_iterator it = setGivenSplits.begin();
       it != setGivenSplits.end(); ++it) {
    int sz = it->size();
    listGivenSplits[sz].insert(*it);
  }
  // if the whole set is not in, add it so that we have a single lin in the end
  if (listGivenSplits[numTaxa].size() == 0) {
    //
    set<int> sall;
    PopulateSetWithInterval(sall, 0, numTaxa - 1);
    listGivenSplits[numTaxa].insert(sall);
  }
#if 0
cout << "Set of given splits: ";
for(int i=0; i<(int)listGivenSplits.size(); ++i)
{
if( listGivenSplits[i].size() > 0 )
{
for( set<set<int> > :: iterator it = listGivenSplits[i].begin(); it != listGivenSplits[i].end(); ++it)
{
DumpIntSet( *it );
}
}
}
#endif

  // active list of lineages indexed by their set
  map<set<int>, TreeNode *> mapActiveLins;
  // initially all the leaf lins
  int idToUse = 0;
  for (int i = 0; i < numTaxa; ++i) {
    TreeNode *pLeaf = new TreeNode(idToUse++);
    set<int> sint;
    sint.insert(i);
    string strLbl = GetStringFromId(i);
    pLeaf->SetLabel(strLbl);
    mapActiveLins.insert(map<set<int>, TreeNode *>::value_type(sint, pLeaf));
  }
  // now scan through the entire list
  for (int k = 2; k < (int)listGivenSplits.size(); ++k) {
    // start from 2 so that avoid trivial sets
    if (listGivenSplits[k].size() == 0) {
      continue;
    }
    // for each input list, find those lins that is contained within the
    // clusters
    for (set<set<int> >::iterator it2 = listGivenSplits[k].begin();
         it2 != listGivenSplits[k].end(); ++it2) {
      // each subset corresponds to a new internal node
      TreeNode *pnode = new TreeNode(idToUse++);

      // cout << "list of active lins: ";
      // for( map< set<int>, TreeNode *> :: iterator iggg =
      // mapActiveLins.begin(); iggg != mapActiveLins.end(); ++iggg )
      //{
      //    DumpIntSet( iggg->first);
      //}
      // cout << "Considering given split: ";
      // DumpIntSet( *it2 );
      // find the proper node in the previous set
      set<set<int> > setMatached;
      int szTot = 0;
      for (map<set<int>, TreeNode *>::iterator it3 = mapActiveLins.begin();
           it3 != mapActiveLins.end(); ++it3) {
        //
        // cout <<  "treat this active lineage: ";
        // DumpIntSet( it3->first );
        if (IsSetContainer(*it2, it3->first) == true) {
          // cout << "yes, continer!\n";
          //
          setMatached.insert(it3->first);
          szTot += it3->first.size();
          //
          vector<int> sempty;
          pnode->AddChild(it3->second, sempty);
        }
      }
      YW_ASSERT_INFO(szTot == (int)it2->size(), "Size: mismatch1");
      // remove the old ones and add the newly created one
      for (set<set<int> >::iterator it4 = setMatached.begin();
           it4 != setMatached.end(); ++it4) {
        mapActiveLins.erase(*it4);
      }
      mapActiveLins.insert(map<set<int>, TreeNode *>::value_type(*it2, pnode));
    }
  }
  YW_ASSERT_INFO(mapActiveLins.size() == 1,
                 "Wrong: must have only a single lineage left");
  treeToProc.SetRoot(mapActiveLins.begin()->second);

  // string strNW;
  // treeToProc.ConsNewick(strNW);
  // cout << "Result of createtreebyrootedplits: " << strNW << endl;
  // cout << "SetGivenSplits: \n";
  // for(set<set<int> > :: iterator it = setGivenSplits.begin(); it !=
  // setGivenSplits.end(); ++it)
  //{
  // DumpIntSet( *it);
  //}
  // cout << "numTaxa: " << numTaxa << endl;
}

void DumpAllSubtreesWithTaxaSize(
    const vector<PhylogenyTreeBasic *> &listPtrGTrees, int numTaxonSubtree,
    const char *fileNameOut) {
  ofstream outfile(fileNameOut);

  // dump out subtrees with certain number of taxa (if the tree contains fewer
  // than this number, just dump out the entire tree)
  for (int tr = 0; tr < (int)listPtrGTrees.size(); ++tr) {
    //
    set<string> listLeafLabelsSet;
    vector<string> listLeafLabels, listLeafLabelsSetDistinct;
    listPtrGTrees[tr]->GetAllLeafLabeles(listLeafLabels);
    PopulateSetByVecGen(listLeafLabelsSet, listLeafLabels);
    PopulateVecBySetGen(listLeafLabelsSetDistinct, listLeafLabelsSet);

    //
    int numSubsetSz = numTaxonSubtree;
    if (numSubsetSz > (int)listLeafLabelsSetDistinct.size()) {
      numSubsetSz = listLeafLabelsSetDistinct.size();
    }

    // find all subsets
    vector<int> posvec;
    GetFirstCombo(numSubsetSz, (int)listLeafLabelsSetDistinct.size(), posvec);
    while (true) {
      set<string> setTaxaStep;
      for (int i = 0; i < (int)posvec.size(); ++i) {
        setTaxaStep.insert(listLeafLabelsSetDistinct[posvec[i]]);
      }

      //
      PhylogenyTreeBasic *ptreeNew = new PhylogenyTreeBasic;
      listPtrGTrees[tr]->CreatePhyTreeFromLeavesWithLabels(setTaxaStep,
                                                           *ptreeNew, true);
      string nwTree;
      ptreeNew->ConsNewick(nwTree);
      outfile << nwTree << endl;
      delete ptreeNew;

      if (GetNextCombo(numSubsetSz, (int)listLeafLabelsSetDistinct.size(),
                       posvec) == false) {
        break;
      }
    }
  }

  outfile.close();
}

void DumpAllSubtreesWithBoundedSize(
    const vector<PhylogenyTreeBasic *> &listPtrGTrees, int maxSzSubtree,
    int maxIdentSubtreeSz, const char *fileNameOut) {
  //
  // cout << "DumpAllSubtreesWithBoundedSize: maxSzSubtree: " << maxSzSubtree <<
  // ", maxIdentSubtreeSz: " << maxIdentSubtreeSz << ", filenameOut: " <<
  // fileNameOut << endl;
  // dump all subtrees with at most maxSzSubtree leaves into a file (that is,
  // breaking trees into pieces) in order to avoid issues that large subtrees
  // with identical labels, we first shrink such subtree within the size (if
  // exists) e.g. maxIdentSubtreeSz = 5 and maxSzSubtree = 10 YW: 12/09/15: in
  // case of a non-binary tree, we may have multiple subtrees as siblings; if
  // this is the case, output each pair of subtrees YW: 12/10/15: don't output
  // trees with only two siblings
  ofstream outfile(fileNameOut);

  // dump out subtrees with certain number of taxa (if the tree contains fewer
  // than this number, just dump out the entire tree)
  bool fTreeOut = false;
  for (int tr = 0; tr < (int)listPtrGTrees.size(); ++tr) {
    // cout << "Processing tree: " << tr << endl;
    // create a new tree where identical subtrees match what we want
    PhylogenyTreeBasic *ptreeWork =
        ConsPhyTreeShrinkIdentSubtrees(listPtrGTrees[tr], maxIdentSubtreeSz);
    // cout << "tree working: ";
    // ptreeWork->Dump();
    // find all subtrees that are no bigger than the desired ones
    set<TreeNode *> setSTRoots;
    ptreeWork->GetSubtreesWithMaxSize(setSTRoots, maxSzSubtree);
    // cout << "Number of subtrees: " << setSTRoots.size() << endl;

    // find any missing
    // set<string> setLabelsPresent;
    // PhylogenyTreeBasic :: FindAllLabelsInSubtrees(setSTRoots,
    // setLabelsPresent); set<string> setLabelsMiss;
    // ptreeWork->GetRoot()->GetAllDistinctLeafLabeles(setLabelsMiss);
    // SubtractSetsGen(setLabelsMiss, setLabelsPresent);

    // list of all subtrees that are uniform
    set<TreeNode *> setSTUniform;
    for (set<TreeNode *>::iterator it = setSTRoots.begin();
         it != setSTRoots.end(); ++it) {
      //
      set<string> strLblsStep;
      (*it)->GetAllDistinctLeafLabeles(strLblsStep);

      if (strLblsStep.size() == 1) {
        setSTUniform.insert(*it);
      }
    }
    // cout << "Number of uniform subtrees: " << setSTUniform.size() << endl;

    // output each subtree one by one
    PhylogenyTreeBasic *ptreeNew = new PhylogenyTreeBasic;

    while (setSTRoots.size() >= 1) {
#if 0
cout << "Start of each iteration: tree is: ";
ptreeWork->Dump();
cout << "Set of subtrees: \n";
for(set<TreeNode *> :: iterator itt = setSTRoots.begin(); itt != setSTRoots.end(); ++itt)
{
////
(*itt)->Dump();
cout << endl;
}
#endif

      TreeNode *pnSTRootCurr = NULL;
      set<TreeNode *> setSTToRemove;

      // rule: if there is a non-uniform subtree, output it
      for (set<TreeNode *>::iterator itg = setSTRoots.begin();
           itg != setSTRoots.end(); ++itg) {
        //
        if (setSTUniform.find(*itg) == setSTUniform.end()) {
          //
          pnSTRootCurr = *itg;
          break;
        }
      }
      if (pnSTRootCurr == NULL) {
        // if no non-uniform subtrees are found, find a sibling pairs of
        // subtrees and take the whole subtree to output YW: need to be careful;
        // I don't want to have left-over
        set<TreeNode *> ppSibs;
        bool fres =
            PhylogenyTreeBasic ::GetSiblingsNodesFrom(setSTRoots, ppSibs);
        YW_ASSERT_INFO(fres == true, "Fail to find silblings");
        pnSTRootCurr = (*ppSibs.begin())->GetParent();
        while (true) {
          // find out how many subtrees covered if taking this
          set<TreeNode *> setSTCoveredStep;
          PhylogenyTreeBasic ::FindDescendentsOfNodeWithin(
              pnSTRootCurr, setSTRoots, setSTCoveredStep);
          // if there are at least two left, use it or we have reached the root
          if ((int)setSTCoveredStep.size() + 1 < (int)setSTRoots.size() ||
              pnSTRootCurr == ptreeWork->GetRoot()) {
            //
            break;
          } else {
            // move up
            pnSTRootCurr = pnSTRootCurr->GetParent();
          }
        }

      } else {
        // if there are only one leftover and it is uniform one, output all the
        // tree
        if (setSTRoots.size() == 2 && setSTUniform.size() > 0) {
          //
          pnSTRootCurr = ptreeWork->GetRoot();
        }
      }

      // remove any subtrees that are descendent of the output subtree
      // just output it
      YW_ASSERT_INFO(pnSTRootCurr != NULL, "Cannot be NULL");
      PhylogenyTreeBasic ::FindDescendentsOfNodeWithin(pnSTRootCurr, setSTRoots,
                                                       setSTToRemove);

      // cout << "******** pnSTRootCurr: ";
      // pnSTRootCurr->Dump();

      // if this is a single node or degree of this node is two, just output it
      if (setSTRoots.find(pnSTRootCurr) != setSTRoots.end() ||
          pnSTRootCurr->GetChildrenNum() == 2) {
        // cout << "******** outputing subtree rooted at: ";
        // pnSTRootCurr->Dump();

        ptreeNew->SetRootPlain(pnSTRootCurr);
        // if the tree has at least one intermediate node, output it
        // if( ptreeNew->GetNumInternalNodes() >= 2 )
        {
          string nwTree;
          ptreeNew->ConsNewick(nwTree);
          outfile << nwTree << endl;
          fTreeOut = true;
        }
      } else {
        YW_ASSERT_INFO(pnSTRootCurr->GetChildrenNum() >= 3,
                       "Must be a mulfurcating node");
        // now enumerate all pairs of children of this node
        TreeNode *pnRootNew = new TreeNode;
        ptreeNew->SetRootPlain(pnRootNew);
        vector<TreeNode *> listChildren;
        set<TreeNode *> listChildrenSet;
        pnSTRootCurr->GetAllChildren(listChildrenSet);
        PopulateVecBySetGen(listChildren, listChildrenSet);
        vector<int> posvec;
        GetFirstCombo(2, (int)listChildren.size(), posvec);
        while (true) {
          vector<int> vecdummy;
          pnRootNew->AddChild(listChildren[posvec[0]], vecdummy);
          pnRootNew->AddChild(listChildren[posvec[1]], vecdummy);
          // if( ptreeNew->GetNumInternalNodes() >= 2 )
          {
            string nwTree;
            ptreeNew->ConsNewick(nwTree);
            outfile << nwTree << endl;
            fTreeOut = true;
          }
          pnRootNew->DetachAllChildren();
          listChildren[posvec[0]]->SetParent(pnSTRootCurr);
          listChildren[posvec[1]]->SetParent(pnSTRootCurr);
          if (GetNextCombo(2, (int)listChildren.size(), posvec) == false) {
            break;
          }
        }
        //
        delete pnRootNew;
      }

      // now detach this node from the rest of tree
      TreeNode *pnparcurr = pnSTRootCurr->GetParent();
      pnSTRootCurr->DetachSelf();

      if (pnSTRootCurr == ptreeWork->GetRoot()) {
        break;
      }

      delete pnSTRootCurr;
      pnSTRootCurr = NULL;

      // cout << "Before degree-one cleainup, tree is: ";
      // ptreeWork->Dump();
      // exit(1);

      // if( pnparcurr != NULL && pnparcurr != ptreeWork->GetRoot() )
      if (pnparcurr != NULL) {
        ptreeWork->RemoveDegreeOneNodeAt(pnparcurr);
      }
      for (set<TreeNode *>::iterator it = setSTToRemove.begin();
           it != setSTToRemove.end(); ++it) {
        setSTRoots.erase(*it);
      }
    }
    ptreeNew->SetRootPlain(NULL);
    // delete ptreeNew;
    // cout << "output tree deleted\n";

#if 0
        // create a psuedo tree just for outputing
        PhylogenyTreeBasic *ptreeNew = new PhylogenyTreeBasic;
        for( set<TreeNode *> :: iterator it = setSTRoots.begin(); it != setSTRoots.end(); ++it )
        {
            //
            ptreeNew->SetRootPlain( *it );

            // if there are at least two taxa in the subtree, output it
            //set<string> listLeafLabelsSet;
            //vector<string> listLeafLabels;
            //ptreeNew->GetAllLeafLabeles(listLeafLabels);
            //PopulateSetByVecGen( listLeafLabelsSet, listLeafLabels );
            //if(  listLeafLabelsSet.size() >= 2 )
            //{
                string nwTree;
                ptreeNew->ConsNewick(nwTree, false, 1.0, true);
                outfile << nwTree << endl;

                fTreeOut = true;
            //}
        }
        ptreeNew->SetRootPlain(NULL);
        delete ptreeNew;
cout << "output tree deleted\n";
#endif

    delete ptreeWork;
    // cout << "Shrunk tree deleted.\n";
  }

  outfile.close();
  // cout << "Tree outoupt finished.\n";
  YW_ASSERT_INFO(fTreeOut == true,
                 "ERROR: no subtrees output. Your trees appear to be either "
                 "very clustered into uniform subtrees or the parameters (size "
                 "of subtree and identical trees size upper bounds are wrong.");
}

PhylogenyTreeBasic *ConsPhyTreeShrinkIdentSubtrees(PhylogenyTreeBasic *ptreeIn2,
                                                   int maxIdentSubtreeSz,
                                                   bool fIdConsecutive) {
  // create a new tree
  PhylogenyTreeBasic *ptreeRes = new PhylogenyTreeBasic;
  // construct according to Newick format
  string strNW;
  ptreeIn2->ConsNewick(strNW, false, 1.0, true);
  ptreeRes->ConsOnNewick(strNW);

  // cout << "ConsPhyTreeShrinkIdentSubtrees: tree in: " << strNW << endl;

  // create a tree with identical subtree that is no greater than the given size
  // (i.e. if a subtree is of the same label, shrink it if needed) first obtain
  // the max identity subtrees
  set<TreeNode *> setSTRootsIdents;
  ptreeRes->GetMaxSubtrees(setSTRootsIdents);
  // cout << "Number of maximiaml subtrees: " << setSTRootsIdents.size() <<
  // endl;

  // find all the leaves
  vector<set<TreeNode *> > listMaxSubtreesLeaves;
  for (set<TreeNode *>::iterator it = setSTRootsIdents.begin();
       it != setSTRootsIdents.end(); ++it) {
    set<TreeNode *> setLeavesUnder;
    (*it)->GetAllLeavesUnder(setLeavesUnder);
    listMaxSubtreesLeaves.push_back(setLeavesUnder);
    // cout << "Sz of subtree found: " << setLeavesUnder.size() << endl;
  }

  // now remove nodes until the subtrees is no longer too large
  for (int i = 0; i < (int)listMaxSubtreesLeaves.size(); ++i) {
    //
    if ((int)listMaxSubtreesLeaves[i].size() > maxIdentSubtreeSz) {
      vector<TreeNode *> listNodes;
      PopulateVecBySetGen(listNodes, listMaxSubtreesLeaves[i]);

      // remove some leaves
      for (int j = maxIdentSubtreeSz; j < (int)listMaxSubtreesLeaves[i].size();
           ++j) {
        ptreeRes->RemoveNodeKeepChildren(listNodes[j]);

        // cout << "After removing a leaf: current tree: ";
        // string tr;
        // ptreeRes->ConsNewick(tr);
        // cout << tr << endl;
      }
    }
  }

  // cout << "ConsPhyTreeShrinkIdentSubtrees: resulting tree: ";
  // string tr;
  // ptreeRes->ConsNewick(tr);
  // cout << tr << endl;

  // YW: set consecutive id?
  if (fIdConsecutive == true) {
    AssignConsecutiveIdsForTree(*ptreeRes);
  }

  return ptreeRes;
}

void ChangebackLeafLabelForTreeWithZeroBaseId(PhylogenyTreeBasic *ptree,
                                              TaxaMapper *pTMapper) {
  // cout << "Before ChangebackLeafLabelForTreeWithZeroBaseId: ";
  // ptree->Dump();
  //
  YW_ASSERT_INFO(pTMapper != NULL, "Must have a mapper");
  vector<TreeNode *> listLeafNodes;
  ptree->GetAllLeafNodes(listLeafNodes);
  for (int i = 0; i < (int)listLeafNodes.size(); ++i) {
    // get the int id
    int lbl = listLeafNodes[i]->GetIntLabel();
    string lblOrig = pTMapper->GetString(lbl);
    // cout << "lbl:" << lbl << ", lblOrig: " << lblOrig << endl;;
    listLeafNodes[i]->SetLabel(lblOrig);
  }
  // cout << "After ChangebackLeafLabelForTreeWithZeroBaseId: ";
  // ptree->Dump();
}

bool ConvPhyloTreesToZeroBasedId(vector<PhylogenyTreeBasic *> &treePtrList,
                                 TaxaMapper *pTMapper) {
  // the given trees are not zero-based; so convert them to be; pTMMapeer: not
  // initialied upon entry; then store the mapping between id to string
  for (int i = 0; i < (int)treePtrList.size(); ++i) {
    vector<TreeNode *> listLeafNodes;
    treePtrList[i]->GetAllLeafNodes(listLeafNodes);
    if (pTMapper->IsInitialized() == false) {
      //
      for (int i = 0; i < (int)listLeafNodes.size(); ++i) {
        // get the int id
        string lbl = listLeafNodes[i]->GetLabel();
        int idTouse = pTMapper->AddTaxaString(lbl);
        // cout << "lbl:" << lbl << ", lblOrig: " << lblOrig << endl;;
        listLeafNodes[i]->SetIntLabel(idTouse);
      }
    } else {
      //
      for (int i = 0; i < (int)listLeafNodes.size(); ++i) {
        // get the int id
        int lbl = listLeafNodes[i]->GetIntLabel();
        string lblOrig = pTMapper->GetString(lbl);
        // cout << "lbl:" << lbl << ", lblOrig: " << lblOrig << endl;;
        listLeafNodes[i]->SetLabel(lblOrig);
      }
    }
  }
  return true; // for now, just true
}

void ChangeLeafIntLabelOfTree(PhylogenyTreeBasic &treeToChange,
                              const map<int, int> &mapOldIntLblToNewIntLbl,
                              bool fSetUserLblToo) {
#if 0
cout << "Before ChangeLeafIntLabelOfTree: ";
treeToChange.Dump();
cout << "mapOldIntLblToNewIntLbl: ";
for(map<int,int> :: const_iterator it = mapOldIntLblToNewIntLbl.begin(); it != mapOldIntLblToNewIntLbl.end(); ++it)
{
cout << "[" << it->first << "," << it->second << "] ";
}
cout << endl;
#endif
  //
  vector<TreeNode *> listLeafNodes;
  treeToChange.GetAllLeafNodes(listLeafNodes);
  for (int i = 0; i < (int)listLeafNodes.size(); ++i) {
    // get the int id
    int lbl = listLeafNodes[i]->GetIntLabel();

    if (mapOldIntLblToNewIntLbl.find(lbl) == mapOldIntLblToNewIntLbl.end()) {
      treeToChange.Dump();
      cout << "lbl: " << lbl << endl;
      cout << "mapOldIntLblToNewIntLbl: ";
      for (map<int, int>::const_iterator it = mapOldIntLblToNewIntLbl.begin();
           it != mapOldIntLblToNewIntLbl.end(); ++it) {
        cout << "[" << it->first << ", " << it->second << "]   ";
      }
      cout << endl;
    }

    YW_ASSERT_INFO(mapOldIntLblToNewIntLbl.find(lbl) !=
                       mapOldIntLblToNewIntLbl.end(),
                   "Fail to find the orignal label");
    int lblIntNew = (*(mapOldIntLblToNewIntLbl.find(lbl))).second;
    // cout << "lbl:" << lbl << ", lblIntNew: " << lblIntNew << endl;;
    listLeafNodes[i]->SetIntLabel(lblIntNew);
    if (fSetUserLblToo) {
      char buf[100];
      sprintf(buf, "%d", lblIntNew);
      string strbuf(buf);
      listLeafNodes[i]->SetUserLabel(strbuf);
    }
#if 0
// for now, also set user label as well
char buf[100];
sprintf(buf, "%d",lblIntNew);
string strbuf(buf);
listLeafNodes[i]->SetUserLabel(strbuf);
#endif
  }
#if 0
cout << "After ChangeLeafIntLabelOfTree: ";
treeToChange.Dump();
#endif
}

void AssignConsecutiveIdsForTree(PhylogenyTreeBasic &treeToChange) {
  //
  vector<TreeNode *> listAllNodes;
  treeToChange.GetAllNodes(listAllNodes);
  int idToUse = 0;
  for (int i = 0; i < (int)listAllNodes.size(); ++i) {
    // leaves assigned to a distinct id first
    if (listAllNodes[i]->IsLeaf() == true) {
      listAllNodes[i]->SetID(idToUse++);
    }
  }
  for (int i = 0; i < (int)listAllNodes.size(); ++i) {
    // leaves assigned to a distinct id first
    if (listAllNodes[i]->IsLeaf() == false) {
      listAllNodes[i]->SetID(idToUse++);
    }
  }
}

void RandTrimLeavesFromTree(PhylogenyTreeBasic *ptreeToTrim,
                            int numLeavesRemain) {
  // do nothing if the gene trees are small
  if (ptreeToTrim->GetNumLeaves() <= numLeavesRemain) {
    return;
  }

  // cout << "RandTrimLeavesFromTree: before trimming: tree is: ";
  // string strNW;
  // ptreeToTrim->ConsNewick( strNW, false, 1.0, true );
  // cout << strNW << endl;

  // for a large tree, we want to randomly trim some leaves to make the tree
  // smaller rule: never completely delete some leaf label; prefer to deleting
  // leaves that appear more frequently
  map<int, set<TreeNode *> > mapLeafLblToNodes;
  vector<TreeNode *> listLeafNodes;
  ptreeToTrim->GetAllLeafNodes(listLeafNodes);
  for (int i = 0; i < (int)listLeafNodes.size(); ++i) {
    int lbl = listLeafNodes[i]->GetIntLabel();
    if (mapLeafLblToNodes.find(lbl) == mapLeafLblToNodes.end()) {
      set<TreeNode *> ss;
      mapLeafLblToNodes.insert(map<int, set<TreeNode *> >::value_type(lbl, ss));
    }
    mapLeafLblToNodes[lbl].insert(listLeafNodes[i]);
  }
  // create a list of nodes to remove
  vector<set<TreeNode *> > listNodesToRemove;
  vector<double> listNodesToRemoveSz;
  for (map<int, set<TreeNode *> >::iterator it = mapLeafLblToNodes.begin();
       it != mapLeafLblToNodes.end(); ++it) {
    //
    listNodesToRemove.push_back(it->second);
    listNodesToRemoveSz.push_back(it->second.size());
  }
  // now start removing
  int numLeavesCurr = ptreeToTrim->GetNumLeaves();
  while (numLeavesCurr > numLeavesRemain) {
    int indexChosen = GetWeightedRandItemIndex(listNodesToRemoveSz);

    if (listNodesToRemoveSz[indexChosen] < 1.01) {
      // cannot delete the one with only one copy left
      continue;
    }
    YW_ASSERT_INFO(listNodesToRemove[indexChosen].size() >= 2, "Wrong");
    TreeNode *pnToRm = *(listNodesToRemove[indexChosen].begin());
    listNodesToRemove[indexChosen].erase(pnToRm);
    --numLeavesCurr;
    TreeNode *pnPar = pnToRm->GetParent();
    ptreeToTrim->RemoveNode(pnToRm);
    ptreeToTrim->RemoveDegreeOneNodeAt(pnPar);
    listNodesToRemoveSz[indexChosen] -= 1.0;
  }
  AssignConsecutiveIdsForTree(*ptreeToTrim);
  // cout << "RandTrimLeavesFromTree: After trimming: tree is: ";
  // string strNW2;
  // ptreeToTrim->ConsNewick( strNW2, false, 1.0, true );
  // cout << strNW2 << endl;
}

// ***************************************************************************
void NewickUtils ::RetrieveLabelSet(const string &strNW,
                                    multiset<string> &setLabels) {
  // cout << "RetrieveLabelSet: strNW = " << strNW << endl;
  //
  setLabels.clear();

  string strIdDirect = strNW;
  int curpos = 0;
  int lastposOut = 0;
  // char *strIdBuf = (char *)strIdDirect.c_str();
  while (curpos < (int)strNW.length()) {
    // cout << "curpos = " << curpos << endl;
    bool fIdentifier = false;
    if ((strNW[curpos] == '(' || strNW[curpos] == ',') &&
        (curpos == (int)strNW.length() - 1 || strNW[curpos + 1] != '(')) {
      fIdentifier = true;
    }
    // cout << "Adding it: " << strId[curpos] << endl;
    lastposOut++;
    curpos++;

    // should we search for id
    if (fIdentifier == true) {
      // cout << "Now searching for identifier\n";
      // now scan to the right to find the position to read the identifier
      while (curpos < (int)strNW.length()) {
        if (strNW[curpos] != ')' && strNW[curpos] != ':' &&
            strNW[curpos] != ',') {
          curpos++;
        } else {
          break;
        }
      }
      //
      // curpos--;
      string strFoundId;
      // cout << "lastposOut = " << lastposOut << ", curpos = " << curpos <<
      // endl;
      strFoundId = strNW.substr(lastposOut, curpos - lastposOut);
      setLabels.insert(strFoundId);
      lastposOut = curpos;
      // cout << "One identifier found: " << strFoundId << endl;
    }
  }
}

bool NewickUtils ::FindSplitIn(const string &strNW, string &strPart1,
                               string &strPart2) {
  // break up the NW into two parts by the center ,
  // return false if atomic
  int posSplit = -1;
  int level = 0;
  for (int i = 0; i < (int)strNW.length(); ++i) {
    if (strNW[i] == '(') {
      level++;
    } else if (strNW[i] == ')') {
      level--;
    } else if (strNW[i] == ',') {
      if (level == 1) {
        posSplit = i;
        break;
      }
    }
  }

  if (posSplit < 0) {
    return false;
  }
  //
  int posLeft = strNW.find('(');
  int posRight = strNW.rfind(')');
  strPart1 = strNW.substr(posLeft + 1, posSplit - posLeft - 1);
  strPart2 = strNW.substr(posSplit + 1, posRight - posSplit - 1);

  return true;
}

void NewickUtils ::UpdateLabells(string &strNW,
                                 const map<string, string> &mapOldLabelToNew) {
  // change the taxa name in the old newick format to the new ones as recorded
  // in the map
  string strNWNew;
  string strIdDirect = strNW;
  int curpos = 0;
  int lastposOut = 0;
  map<string, string> &mapOldLabelToNewRef =
      const_cast<map<string, string> &>(mapOldLabelToNew);
  // bool fOutputCurChar = true;
  // char *strIdBuf = (char *)strIdDirect.c_str();
  while (curpos < (int)strNW.length()) {
    // cout << "curpos = " << curpos << endl;
    bool fIdentifier = false;
    if ((strNW[curpos] == '(' || strNW[curpos] == ',') &&
        (curpos == (int)strNW.length() - 1 || strNW[curpos + 1] != '(')) {
      fIdentifier = true;
    }

    // add it always since this is deliminator
    strNWNew += strNW[curpos];

    // cout << "Adding it: " << strId[curpos] << endl;
    lastposOut++;
    curpos++;

    // should we search for id
    if (fIdentifier == true) {
      // cout << "Now searching for identifier\n";
      // now scan to the right to find the position to read the identifier
      while (curpos < (int)strNW.length()) {
        if (strNW[curpos] != ')' && strNW[curpos] != ':' &&
            strNW[curpos] != ',') {
          curpos++;
        } else {
          break;
        }
      }
      //
      // curpos--;
      string strFoundId;
      // cout << "lastposOut = " << lastposOut << ", curpos = " << curpos <<
      // endl;
      strFoundId = strNW.substr(lastposOut, curpos - lastposOut);

      //
      YW_ASSERT_INFO(mapOldLabelToNew.find(strFoundId) !=
                         mapOldLabelToNew.end(),
                     "Fail to find the id in the map");
      strNWNew.append(mapOldLabelToNewRef[strFoundId]);

      lastposOut = curpos;
      // cout << "One identifier found: " << strFoundId << endl;

      // now move back by one letter
      //--curpos;
    }
  }

  // cout << "UpdateLabells: before update, newick = " << strNW << ", after
  // update: " << strNWNew << endl;
  strNW = strNWNew;
}

string NewickUtils ::RemoveBrLenFromTree(string &strNW) {
  //
  int curpos = 0;
  bool fSkip = false;
  string strNWNew;
  // char *strIdBuf = (char *)strIdDirect.c_str();
  while (curpos < (int)strNW.length()) {
    // cout << "curpos = " << curpos << endl;
    if (strNW[curpos] == ':') {
      fSkip = true;
    } else if (strNW[curpos] == ',' || strNW[curpos] == ')' ||
               strNW[curpos] == ';') {
      // continue skipping until reaching the separate: , or )
      fSkip = false;
    }

    if (fSkip == false) {
      strNWNew += strNW[curpos];
    }

    curpos++;
  }
  return strNWNew;
}

void NewickUtils ::ConsolidateSinglChildChain(string &strNW) {
  if (strNW[0] != '(') {
    // nothing needs to be done
    return;
  }

  // cout << "conslidate: " << strNW << endl;
  // sometime there may be a nested chain of enclosed parenthesis
  // consolidate these; and maintain the proper branch length if there are
  string strRes = strNW;
  double lenTot = 0.0;
  bool fLen = false;
  // bool fParenthRemoved = false;

  while (true) {
    // cout << "current string: " << strRes << endl;
    // stop if it become automic
    string str1, str2;
    bool fNonAtom = FindSplitIn(strRes, str1, str2);

    // if( fNonAtom == false )
    //{
    // fParenthRemoved = true;

    //
    YW_ASSERT_INFO(strRes[0] == '(', "wrong");
    int posRight = strRes.rfind(')');
    YW_ASSERT_INFO(posRight > 0, "wrong1");
    // cout << "posRight: " << posRight << endl;
    if (posRight != (int)strRes.length() - 1) {
      int posLen = strRes.find(':', posRight);
      // cout << "posLen: " << posLen << endl;
      if (posLen > 0) {
        // if( lenTot > 0.0)
        //{
        // cout << "*HHHHH\n";
        //}
        fLen = true;
        lenTot += GetLenAt(strRes, posLen + 1);
        // cout << "lenTot: " << lenTot << endl;
      }
    }
    //}

    int len = posRight - 1;
    strRes = strRes.substr(1, len);

    if (fNonAtom == true) {
      break;
    }
  }
  string strRes1;
  // if( fParenthRemoved )
  //{
  strRes1 += "(";
  //}
  strRes1 += strRes;
  // if( fParenthRemoved )
  //{
  strRes1 += ")";
  //}
  if (fLen) {
    strRes1 += ":" + std::to_string(lenTot);
  }

  strNW = strRes1;
  // cout << "conslidate to " << strNW << endl;
}

double NewickUtils ::GetLenAt(const string &strNW, int posLen) {
  //
  int posLenEnd = strNW.length() - 1;
  int sepPos1 = strNW.find(',', posLen);
  int sepPos2 = strNW.find(')', posLen);
  if (sepPos1 > 0 && sepPos1 - 1 < posLenEnd) {
    posLenEnd = sepPos1 - 1;
  }
  if (sepPos2 > 0 && sepPos2 - 1 < posLenEnd) {
    posLenEnd = sepPos2 - 1;
  }
  if (posLenEnd <= posLen) {
    cout << "posLen: " << posLen << ", posLenEnd: " << posLenEnd
         << ", tree: " << strNW << endl;
  }
  YW_ASSERT_INFO(posLenEnd >= posLen, "No length found");
  string lenstr = strNW.substr(posLen, posLenEnd - posLen + 1);
  return atof(lenstr.c_str());
}

// ***************************************************************************

TaxaMapper ::TaxaMapper() {
  curId = 0;
  fInit = false;
}

// utility
bool TaxaMapper ::IsEmpty() { return mapStrToId.size() == 0; }

int TaxaMapper ::AddTaxaString(const string &str) {
  // cout << "AddTaxaString : curId = " << curId  << " for new taxa string " <<
  // str << endl;
  if (mapStrToId.find(str) == mapStrToId.end()) {
    mapStrToId.insert(map<string, int>::value_type(str, curId));
    mapIdToStr.insert(map<int, string>::value_type(curId, str));
    curId++;
  }
  // else
  //{
  return mapStrToId[str];
  //}
}

void TaxaMapper ::AddTaxaStringWithId(int tid, const string &str) {
  // caution: don't mix up with the previous auto-id mode
  mapStrToId.insert(map<string, int>::value_type(str, tid));
  mapIdToStr.insert(map<int, string>::value_type(tid, str));
}

int TaxaMapper ::GetId(const string &str) {
  // cout << "Num of entries in str mapper : " << mapStrToId.size() << endl;
  // for( map<string,int> :: iterator it =mapStrToId.begin(); it !=
  // mapStrToId.end(); ++it )
  //{
  // cout << it->first << ", " << it->second << endl;
  //}

  if (mapStrToId.find(str) == mapStrToId.end()) {
    // when the str is not pre-recorded, return negative value
    return -1;
    // cout << "This taxa: " << str << " seems to be wrong\n";
    // YW_ASSERT_INFO( false, "Fail to find the taxa" );
  }
  return mapStrToId[str];
}
bool TaxaMapper ::IsIdIn(int id) {
  return mapIdToStr.find(id) != mapIdToStr.end();
}

string TaxaMapper ::GetString(const int id) {
  if (mapIdToStr.find(id) == mapIdToStr.end()) {
    cout << "mapIdToStr: ";
    for (map<int, string>::iterator it = mapIdToStr.begin();
         it != mapIdToStr.end(); ++it) {
      cout << "[" << it->first << "," << it->second << "]  ";
    }
    cout << endl;

    cout << "This taxa id: " << id << " seems to be wrong\n";
    YW_ASSERT_INFO(false, "Fail to find the taxa");
  }
  return mapIdToStr[id];
}

string TaxaMapper ::ConvIdStringWithOrigTaxa(const string &strId) {
#if 0
cout << "strID: " << strId << ": Num of entries in str mapper : " << mapIdToStr.size() << endl;
for( map<int,string> :: iterator it =mapIdToStr.begin(); it != mapIdToStr.end(); ++it )
{
cout << it->first << ", " << it->second << endl;
}
#endif
  // convert a string with id (i.e. integer-based identifier) back
  // to user-specified format
  // Simple approach: find everything bebetween ( and , (or :),  and ) and
  // convert to
  // YW: 05/02/19: also allow '#' as seperator to support mutation tree
  string res;
  string strIdDirect = strId;
  int curpos = 0;
  int lastposOut = 0;
  //	char *strIdBuf = (char *)strIdDirect.c_str();
  while (curpos < (int)strId.length()) {
    // cout << "curpos = " << curpos << ", res = " << res << endl;
    bool fIdentifier = false;
    if ((strId[curpos] == '(' || strId[curpos] == ',' ||
         strId[curpos] == '#') &&
        (curpos == (int)strId.length() - 1 ||
         (strId[curpos + 1] != '(' && strId[curpos + 1] != '#'))) {
      fIdentifier = true;
    }
    // cout << "Adding it: " << strId[curpos] << endl;
    res += strId[curpos];
    lastposOut++;
    curpos++;

    // should we search for id
    if (fIdentifier == true) {
      // cout << "Now searching for identifier\n";
      // now scan to the right to find the position to read the identifier
      while (curpos < (int)strId.length()) {
        if (strId[curpos] != ')' && strId[curpos] != ':' &&
            strId[curpos] != ',' && strId[curpos] != '#') {
          curpos++;
        } else {
          break;
        }
      }
      // cout << "lastposOut: " << lastposOut << ", curpos = " << curpos <<
      // endl;
      //
      // curpos--;
      int idnum = -1;
      string strSub = strId.substr(lastposOut, curpos - lastposOut);
      // char buftmp[100];
      // memcpy(buftmp, &strIdBuf[lastposOut], curpos-lastposOut );
      // sscanf(buftmp, "%d", &idnum);
      sscanf(strSub.c_str(), "%d", &idnum);
      string idNew = GetString(idnum);
      ////cout << "After searching, curpos = " << curpos << ", buftmp = " <<
      /// buftmp  << ", idnum = " << idnum << ", idNew = " << idNew << endl;
      // cout << "After searching, curpos = " << curpos << ", strSub = " <<
      // strSub  << ", idnum = " << idnum << ", idNew = " << idNew << endl; char
      // buf[100]; sprintf(buf, "%d", idNew);
      res += idNew;
      lastposOut = curpos;
    }
  }
  return res;
}

string TaxaMapper ::ExtractIdPartFromStr(const string &strIdNW) {
  // extract id part of the string
  string strToUse = strIdNW;
  size_t posSeparator = strIdNW.find(':');

  if (posSeparator != string::npos) {
    strToUse = strIdNW.substr(0, (int)posSeparator);
  }
  return strToUse;
}

int TaxaMapper ::GetIdFromStr(const string &strPart, TaxaMapper *pTMapper) {
  // cout << "GetIdFromStr: " << strPart << endl;

  string strToUse = strPart;
  size_t posSeparator = strPart.find(':');

  if (posSeparator != string::npos) {
    strToUse = strPart.substr(0, (int)posSeparator);
  }

  // 05/07/15: it is also possible user add gene index (in # sign)
  size_t posSeparator2 = strToUse.find('#');
  if (posSeparator2 != string::npos) {
    strToUse = strToUse.substr(0, (int)posSeparator2);
  }
  // cout << "strPart: " << strPart << ",strUse: " << strToUse << endl;

  // get rid of
  int res = -1;
  if (pTMapper == NULL) {
    sscanf(strToUse.c_str(), "%d", &res);
    // cout << "Empty mapper\n";
  } else {
    // are we reading in the first tree or not
    res = pTMapper->GetId(strToUse);
    // if( pTMapper->IsInitialized() == true )
    //{
    //	res  = pTMapper->GetId(strToUse);
    // cout << "GetIdFromStr: GetId: " << strToUse << ": " << res << endl;
    //}
    // else
    if (res < 0) {
      // this label is not seen before, so we add a new record
      // this is new
      res = pTMapper->AddTaxaString(strToUse);
      // cout << "GetIdFromStr: New id: " << strToUse << ": " << res << endl;
    }
  }
  return res;
}

void TaxaMapper ::GetAllTaxaIds(set<int> &taxaIndices) const {
  //
  taxaIndices.clear();
  for (map<int, string>::const_iterator it = mapIdToStr.begin();
       it != mapIdToStr.end(); ++it) {
    taxaIndices.insert(it->first);
  }
}

void TaxaMapper ::GetAllTaxaStrs(set<string> &setStrs) const {
  //
  setStrs.clear();
  for (map<int, string>::const_iterator it = mapIdToStr.begin();
       it != mapIdToStr.end(); ++it) {
    setStrs.insert(it->second);
  }
}

void TaxaMapper ::InitToDec1Mode(int numTaxa) {
  // assume taxa is in the format as 1, 2, 3 and so on
  // init as follows: 1 ==> 0, 2 ==> 1 and so on
  for (int taxa = 1; taxa <= numTaxa; ++taxa) {
    char buf[100];
    sprintf(buf, "%d", taxa);
    string strid = buf;
    AddTaxaString(strid);
  }
  SetInitialized(true);
}

void TaxaMapper ::Dump() const {
  //
  cout << "curId = " << curId;
  if (fInit == true) {
    cout << "initialized. ";
  } else {
    cout << "not initialized yet. ";
  }
  for (map<string, int>::const_iterator it = mapStrToId.begin();
       it != mapStrToId.end(); ++it) {
    //
    cout << "Mapping taxa " << it->first << " to id: " << it->second << "  ";
  }
  cout << endl;
}

// ***************************************************************************
// Tree class functions
// ***************************************************************************
TreeNode ::TreeNode()
    : parent(NULL), id(-1), label("-"), shape(PHY_TN_DEFAULT_SHAPE),
      lenBranchAbove(-1.0) {}

TreeNode ::TreeNode(int iid)
    : parent(NULL), id(iid), label("-"), shape(PHY_TN_DEFAULT_SHAPE),
      lenBranchAbove(-1.0) {
  //    id = iid;
  //    cout << "Creating tree node " << iid << endl;
}

TreeNode ::~TreeNode() {
  // cout << "Deleting tree node " << id << ", number of children: " <<
  // GetChildrenNum() << endl; cout << "Dump: "; Dump();
  // We recursively delete all its children here
  for (int i = 0; i < (int)listChildren.size(); ++i) {
    delete listChildren[i];
  }
  listChildren.clear();
}

void TreeNode::Dump() const {
  //
  cout << "<node: " << GetLabel() << ", id=" << GetID();
  if (lenBranchAbove >= 0.0) {
    cout << ", length = " << lenBranchAbove;
  }
  cout << ", num of child = " << GetChildrenNum() << ">   ";
}

TreeNode *TreeNode ::Copy() {
  // make a copy (and its descendents)
  TreeNode *pCopy = new TreeNode(GetID());
  pCopy->SetLabel(this->GetLabel());
  pCopy->SetUserLabel(this->GetUserLabel());
  pCopy->lenBranchAbove = this->lenBranchAbove;
  pCopy->nodeValues = this->nodeValues;
  for (int i = 0; i < GetChildrenNum(); ++i) {
    TreeNode *pccopy = GetChild(i)->Copy();
    vector<int> listLbelsCopy;
    if ((int)this->listEdgeLabels.size() >= i + 1) {
      listLbelsCopy = this->listEdgeLabels[i];
    }
    pCopy->AddChild(pccopy, listLbelsCopy);
  }
  return pCopy;
}

void TreeNode ::AddChild(TreeNode *pChild, const vector<int> &labels) {
  // This function add an edge. The edge can be labeled with a set of labels
  // (for now, only integers)
  YW_ASSERT(pChild != NULL);

  // make sure this child is not already a children
  // not sure if really need it

  pChild->parent = this;
  listChildren.push_back(pChild);
  listEdgeLabels.push_back(labels);
}

void TreeNode ::AddEdgeLabelToChild(int cIndex, int lbl) {
  YW_ASSERT_INFO(cIndex < GetChildrenNum(), "Overflow");
  this->listEdgeLabels[cIndex].push_back(lbl);
}

void TreeNode ::RemoveChild(TreeNode *pChild) {
  YW_ASSERT_INFO(pChild != NULL, "RemoveChild: wrong");
  pChild->parent = NULL;
  vector<TreeNode *> listChildrenNew;
  vector<vector<int> > listEdgeLabelsNew;
  YW_ASSERT_INFO(listChildrenNew.size() == listEdgeLabelsNew.size(),
                 "must be same size");
  for (int i = 0; i < (int)listChildren.size(); ++i) {
    if (listChildren[i] != pChild) {
      listChildrenNew.push_back(listChildren[i]);
      listEdgeLabelsNew.push_back(listEdgeLabels[i]);
    }
  }
  // update
  listChildren = listChildrenNew;
  listEdgeLabels = listEdgeLabelsNew;
}

void TreeNode ::RemoveAllChildren() {
  // remove all children of this node
  // listChildren.clear();
  // listEdgeLabels.clear();
  while (GetChildrenNum() > 0) {
    TreeNode *pc = GetChild(0);
    // cout << "Removing pc = ";
    // pc->Dump();
    // cout << endl;
    RemoveChild(pc);
  }
  // cout << "Done with removeallchildren\n";
}

void TreeNode ::DetachAllChildren() {
  // diff from RemoveAllChildren, simply detach the children from the parent
  // (i.e. parent no longer has record for these children)
  this->listChildren.clear();
  this->listEdgeLabels.clear();
}

void TreeNode ::DetachSelf() {
  // detach this node from parent (but don't perform any memory release)
  TreeNode *pp = GetParent();

  if (pp != NULL) {
    //
    pp->RemoveChild(this);
  }
}

void TreeNode ::GetDescendentLabelSet(set<int> &labelSet) {
  // This function accumulate the set of descendents in the label sets
  // CAUTION: assume labelset is EMPTY!!!!
  // if( IsLeaf() == true)
  //{
  string lbl = GetLabel();
  // cout << "lbl = " << lbl << endl;

  if (lbl != "-" && lbl != "?" && lbl != "()" && lbl != "(?)") {
    const char *buf = lbl.c_str();
    int rowIndex;
    if (buf[0] < '0' || buf[0] > '9') {
      sscanf(buf + 1, "%d", &rowIndex);
    } else {
      // This is a plain label, use it
      sscanf(buf, "%d", &rowIndex);
    }
    // cout << "rowIndex = " << rowIndex << endl;
    labelSet.insert(rowIndex);
  } else if (nodeValues.size() >= 1) {
    // simply insert a single value here
    // labelSet.insert( nodeValues[0] );
  }

#if 0
        // set every label into the set
        for(int i=0; i<nodeValues.size(); ++i)
        {
            if( nodeValues[i] >= 0 )
            {
                labelSet.insert( nodeValues[i] );
            }
        }
#endif
  //}
  // else
  if (IsLeaf() == false) {
    for (int i = 0; i < GetChildrenNum(); ++i) {
      GetChild(i)->GetDescendentLabelSet(labelSet);
    }
  }
}

bool TreeNode ::IsAncesterOf(TreeNode *pAssumedDescend, int &branchIndex) {
  // This function check to see if pAssumedDescend is descedent of the current
  // node If so, we also find the branch index that comes to this node
  if (pAssumedDescend == NULL) {
    return false;
  }
  if (pAssumedDescend == this) {
    branchIndex = -1;
    return true;
  }

  TreeNode *pCurrent = pAssumedDescend;
  TreeNode *pParent = pAssumedDescend->parent;

  while (pParent != NULL) {
    if (pParent == this) {
      // Find out which branch leads to it
      branchIndex = -1;
      for (int i = 0; i < (int)listChildren.size(); ++i) {
        if (listChildren[i] == pCurrent) {
          branchIndex = i;
        }
      }
      YW_ASSERT(branchIndex >= 0);
      // Tell the good news
      return true;
    }
    pCurrent = pParent;
    pParent = pParent->parent;
  }

  return false;
}

void TreeNode ::GetAllDescendents(set<TreeNode *> &setDescendents) {
  // Note: include itself
  setDescendents.insert(this);
  for (int i = 0; i < (int)listChildren.size(); ++i) {
    listChildren[i]->GetAllDescendents(setDescendents);
  }
}

void TreeNode ::GetAllLeavesUnder(set<TreeNode *> &setDescendents) {
  // Note: include itself
  if (this->IsLeaf() == true) {
    setDescendents.insert(this);
  }
  for (int i = 0; i < (int)listChildren.size(); ++i) {
    listChildren[i]->GetAllLeavesUnder(setDescendents);
  }
}

void TreeNode ::GetAllLeavesIdUnder(set<int> &setDescendents) {
  set<TreeNode *> ss;
  GetAllLeavesUnder(ss);
  setDescendents.clear();
  for (set<TreeNode *>::iterator it = ss.begin(); it != ss.end(); ++it) {
    setDescendents.insert((*it)->GetID());
  }
}

void TreeNode ::GetAllDescendIntLbls(set<int> &setIntLbs) {
  //
  if (this->IsLeaf() == true) {
    setIntLbs.insert(this->GetIntLabel());
  } else {
    for (int i = 0; i < (int)listChildren.size(); ++i) {
      listChildren[i]->GetAllDescendIntLbls(setIntLbs);
    }
  }
}

void TreeNode ::GetAllLeafLabeles(vector<string> &listLeafLabels) {
  //
  if (IsLeaf() == true) {
    listLeafLabels.push_back(GetLabel());
  } else {
    for (int i = 0; i < (int)listChildren.size(); ++i) {
      listChildren[i]->GetAllLeafLabeles(listLeafLabels);
    }
  }
}
void TreeNode ::GetAllLeafIntLabeles(vector<int> &listLeafLabels) {
  //
  if (IsLeaf() == true) {
    listLeafLabels.push_back(GetIntLabel());
  } else {
    for (int i = 0; i < (int)listChildren.size(); ++i) {
      listChildren[i]->GetAllLeafIntLabeles(listLeafLabels);
    }
  }
}

void TreeNode ::GetAllDistinctLeafLabeles(set<string> &setLeafLabels) {
  //
  vector<string> listLeafLabels;
  GetAllLeafLabeles(listLeafLabels);
  PopulateSetByVecGen(setLeafLabels, listLeafLabels);
}

string TreeNode ::GetShapeLabel(const set<int> &idTerms,
                                map<int, int> &mapNodeLabel) const {
  // cout << "idTerms = ";
  // DumpIntSet( idTerms );
  string res;

  // return a shape label:
  // at present, shape label is like ((),(())). That is, no leaf labels
  // just the type of topology. Note if we have (S1,S2), then S1 <= S2
  if (idTerms.find(GetID()) != idTerms.end()) {
    int idNum = 1;
    if (mapNodeLabel.find(GetID()) != mapNodeLabel.end()) {
      idNum = mapNodeLabel[GetID()];
    }
    char buf[100];
    sprintf(buf, "%d", idNum);
    res = buf;
    // string str1 = "A";
    // return str1;
  }
  // else
  //	{
  //		string strEmpty;
  //		res = strEmpty;
  //	}
  //}
  else {
    // otherwise get its descendent
    vector<string> listLabels;
    for (int i = 0; i < (int)listChildren.size(); ++i) {
      listLabels.push_back(
          listChildren[i]->GetShapeLabel(idTerms, mapNodeLabel));
    }
    // now sort it
    for (int i = 0; i < (int)listLabels.size(); ++i) {
      for (int j = i + 1; j < (int)listLabels.size(); ++j) {
        // swap if needed
        if (listLabels[i] > listLabels[j]) {
          string tmp = listLabels[i];
          listLabels[i] = listLabels[j];
          listLabels[j] = tmp;
        }
      }
    }

    // how many are not empty?
    int numNonEmpty = 0;
    for (int i = 0; i < (int)listLabels.size(); ++i) {
      if (listLabels[i].length() > 0) {
        numNonEmpty++;
      }
    }

    // add it
    bool fStart = false;
    for (vector<string>::iterator it = listLabels.begin();
         it != listLabels.end(); ++it) {
      if (it->length() > 0) {
        if (fStart == false) {
          if (numNonEmpty > 1) {
            // add a header
            res = "(";
          }
        } else {
          res += ",";
        }
        res += *it;

        fStart = true;
      }
    }
    if (fStart == true && numNonEmpty > 1)
    // if( fStart == true  )
    {
      res += ")";
    }
  }
  // cout << "res label for this node: " << res << endl;
  return res;
}

// differeent from above, this one will apply label to the string label
string TreeNode ::GetShapeLabel(const set<int> &idTerms, bool fSort) const {
  // cout << "idTerms = ";
  // DumpIntSet( idTerms );
  string res;

  // return a shape label:
  // at present, shape label is like ((),(())). That is, no leaf labels
  // just the type of topology. Note if we have (S1,S2), then S1 <= S2
  if (idTerms.find(GetID()) != idTerms.end()) {
    // int idNum = 1;
    if (fSort == true) {
      res = "1";
    } else {
      char buf[100];
      sprintf(buf, "%d", GetID());
      res = buf;
    }
  }

  else {
    // otherwise get its descendent
    vector<string> listLabels;
    for (int i = 0; i < (int)listChildren.size(); ++i) {
      listLabels.push_back(listChildren[i]->GetShapeLabel(idTerms, fSort));
    }
    // now sort it
    if (fSort == true) {
      for (int i = 0; i < (int)listLabels.size(); ++i) {
        for (int j = i + 1; j < (int)listLabels.size(); ++j) {
          // swap if needed
          if (listLabels[i] > listLabels[j]) {
            string tmp = listLabels[i];
            listLabels[i] = listLabels[j];
            listLabels[j] = tmp;
          }
        }
      }
    }

    // how many are not empty?
    int numNonEmpty = 0, numEmpty = 0;
    for (int i = 0; i < (int)listLabels.size(); ++i) {
      if (listLabels[i].length() > 0) {
        numNonEmpty++;
      } else {
        numEmpty++;
      }
    }

    // add it
    bool fStart = false;
    // bool fFirst = true;
    bool fParenth = false;
    // bool fSpaceAdded = false;
    for (vector<string>::iterator it = listLabels.begin();
         it != listLabels.end(); ++it) {
      // YW: only add "(" if there are more than 1 non-empty below
      if (fStart == false && it->length() > 0) {
        // YW: just add a "("
        // if(  (numNonEmpty >= 1 && numEmpty > 0 ) || numNonEmpty >= 2  )
        //{
        // add a header
        if (numNonEmpty > 1) {
          res = "(";
          fParenth = true;
        }
        res += *it;
        fStart = true;
        //}
      } else if (fStart == true) {
        // YW: only add "," if there is something
        if (it->length() > 0) {
          res += ",";
        }
        // fFirst = false;
        if (it->length() > 0) {

          res += *it;
          // fStart = true;
        }
        // YW: donot add anything if the branch is empty
#if 0
				else
				{
					// for empty branches, put a mark to it
					// when there is something under it (that is shrink the entire subtree of unknown to a symbol -
					if(numNonEmpty >= 1 && fSpaceAdded == false)
					{
						//
						res += ",-";
						fSpaceAdded = true;
					}
				}
#endif
      }
    }
    // if( fStart == true  && numNonEmpty >= 1)
    if (fParenth == true) {
      res += ")";
    }
  }
  // cout << "res label for this node: " << res << endl;
  return res;
}

string TreeNode::GetShapeLabelNodeBrNum(
    map<TreeNode *, pair<int, int> > &mapNodeNumBrannches,
    vector<int> &listOrderedLeaves) {
  // format: <num of underlying branches, event id>, negative for internal nodes
  // the ordered leaves: correspond to their order of appearing in the output
  // newick shape string this can be useful when you want to know how to match
  // the leaves when some sort of comparision is needed get shape label.
  // Different from above, the input is: <treenode, #ofbranches out of this
  // node> convention: if #br < 0, it means all branches have descendents
  listOrderedLeaves.clear();
  if (this->IsLeaf() == true) {
    YW_ASSERT_INFO(mapNodeNumBrannches.find(this) != mapNodeNumBrannches.end(),
                   "Leaf: not in map");
    // cout << "Find one leaf: " << mapNodeNumBrannches[this].second << endl;
    listOrderedLeaves.push_back(mapNodeNumBrannches[this].second);
    return string("()");
  } else {
    YW_ASSERT_INFO(mapNodeNumBrannches.find(this) != mapNodeNumBrannches.end(),
                   "Fail to find222");
    // const TreeNode *pn = const_cast<const TreeNode *>( this );
    int numBrWOChildRecur = mapNodeNumBrannches[this].first;
    // cout << "numBrWOChildRecur = " << numBrWOChildRecur << endl;
    multiset<string> setDescStrings;
    map<string, set<vector<int> > > mapStringToVecLeaves;
    for (int i = 0; i < (int)GetChildrenNum(); ++i) {
      //
      TreeNode *pnchild = GetChild(i);
      //
      if (mapNodeNumBrannches.find(pnchild) != mapNodeNumBrannches.end()) {
        //
        vector<int> listOrderedLeavesStep;
        string str = pnchild->GetShapeLabelNodeBrNum(mapNodeNumBrannches,
                                                     listOrderedLeavesStep);
        setDescStrings.insert(str);
        if (mapStringToVecLeaves.find(str) == mapStringToVecLeaves.end()) {
          //
          set<vector<int> > ssint;
          mapStringToVecLeaves.insert(
              map<string, set<vector<int> > >::value_type(str, ssint));
        }
        mapStringToVecLeaves[str].insert(listOrderedLeavesStep);

        //
        --numBrWOChildRecur;
      }
    }
    // add the remaiing by just filling the item
    // vector<int> listLvIds;
    for (int i = 0; i < numBrWOChildRecur; ++i) {
      string strLv = "()";
      setDescStrings.insert(strLv);

      //
      if (mapStringToVecLeaves.find(strLv) == mapStringToVecLeaves.end()) {
        //
        set<vector<int> > ssint;
        mapStringToVecLeaves.insert(
            map<string, set<vector<int> > >::value_type(strLv, ssint));
      }
      vector<int> vec1;
      vec1.push_back(mapNodeNumBrannches[this].second);
      mapStringToVecLeaves[strLv].insert(vec1);
    }
    // cout << "setdescstrings: ";
    // for(multiset<string> :: iterator itgg = setDescStrings.begin(); itgg !=
    // setDescStrings.end(); ++itgg)
    //{
    // cout << *itgg << "   ";
    //}
    // cout << endl;
    // now creat the contacation
    YW_ASSERT_INFO(setDescStrings.size() > 1, "Can not be empty2");
    string res = "(";
    for (multiset<string>::iterator it = setDescStrings.begin();
         it != setDescStrings.end(); ++it) {
      if (it != setDescStrings.begin()) {
        res += ",";
      }
      res += *it;
    }
    res += ")";

    // now assemble the list of ordered nodes
    for (map<string, set<vector<int> > >::iterator itg =
             mapStringToVecLeaves.begin();
         itg != mapStringToVecLeaves.end(); ++itg) {
      for (set<vector<int> >::iterator itg2 = itg->second.begin();
           itg2 != itg->second.end(); ++itg2) {
        // cout << "In GetShapeLabelNodeBrNum: find a vector of sites: ";
        // DumpIntVec(*itg2);
        ConcatIntVec(listOrderedLeaves, *itg2);
      }
    }

    return res;
  }
}

int TreeNode ::GetLevel() const {
  // choose a not efficient but simple coding
  int res = 0;
  for (int i = 0; i < (int)listChildren.size(); ++i) {
    int lvDesc = listChildren[i]->GetLevel();
    if (lvDesc + 1 > res) {
      res = lvDesc + 1;
    }
  }
  return res;
}

void TreeNode ::GetEdgeLabelsToChild(TreeNode *pChild, vector<int> &lbls) {
  YW_ASSERT_INFO(listChildren.size() == listEdgeLabels.size(),
                 "Child num and edge label num do not match");
  lbls.clear();
  for (int i = 0; i < (int)listChildren.size(); ++i) {
    if (listChildren[i] == pChild) {
      GetEdgeLabelsAtBranch(i, lbls);
    }
  }
  // YW_ASSERT_INFO(false, "GetEdgeLabelsToChild :: Fail to find such child");
}

TreeNode *TreeNode ::GetMRCA(TreeNode *pOther) {
  TreeNode *pRes = this;
  int dummy;
  while (pRes != NULL && pRes->IsAncesterOf(pOther, dummy) == false) {
    pRes = pRes->GetParent();
  }
  YW_ASSERT_INFO(pRes != NULL, "Fail to find MRCA");
  return pRes;
}

int TreeNode ::GetNumEdgesToAncestor(TreeNode *pAssumedAncestor) {
  // get # of edges betwene this node to its ancestor
  // return -1 if the ancestor is not true ancestor
  int res = 0;
  TreeNode *pRes = this;
  while (pRes != NULL && pRes != pAssumedAncestor) {
    ++res;
    pRes = pRes->GetParent();
  }
  if (pRes == NULL) {
    res = -1;
  }

  return res;
}

void TreeNode ::GetSiblings(vector<TreeNode *> &listSibs) {
  // siblings are parent's children (except itself)
  listSibs.clear();
  if (this->GetParent() != NULL) {
    //
    for (int i = 0; i < this->GetParent()->GetChildrenNum(); ++i) {
      TreeNode *pn = this->GetParent()->GetChild(i);
      if (pn != this) {
        listSibs.push_back(pn);
      }
    }
  }
}

void TreeNode ::Order() {
  // do nothing if leaf
  if (IsLeaf() == true) {
    return;
  }
  // first order the leaves
  for (int i = 0; i < (int)listChildren.size(); ++i) {
    listChildren[i]->Order();
  }

  //
  vector<multiset<string> > listDescLeaves;
  for (int i = 0; i < (int)listChildren.size(); ++i) {
    vector<string> vecLeafStrings;
    listChildren[i]->GetAllLeafLabeles(vecLeafStrings);
    multiset<string> setLeafStrings;
    for (int j = 0; j < (int)vecLeafStrings.size(); ++j) {
      setLeafStrings.insert(vecLeafStrings[j]);
    }
    listDescLeaves.push_back(setLeafStrings);
  }
  //
  YW_ASSERT_INFO(listEdgeLabels.size() == listChildren.size(),
                 "Same size must be");
  for (int i = 0; i < (int)listChildren.size(); ++i) {
    for (int j = i + 1; j < (int)listChildren.size(); ++j) {
      //
      if (listDescLeaves[i] > listDescLeaves[j]) {
        // exhcnage everything
        TreeNode *ptmp = listChildren[i];
        listChildren[i] = listChildren[j];
        listChildren[j] = ptmp;

        vector<int> vtmp = listEdgeLabels[i];
        listEdgeLabels[i] = listEdgeLabels[j];
        listEdgeLabels[j] = vtmp;

        //
        multiset<string> stmp = listDescLeaves[i];
        listDescLeaves[i] = listDescLeaves[j];
        listDescLeaves[j] = stmp;
      }
    }
  }
}

int TreeNode ::GetIntLabel() const {
  int res = -1;
  sscanf(label.c_str(), "%d", &res);
  return res;
}

void TreeNode ::SetIntLabel(int lbl) {
  //
  char buf[1024];
  sprintf(buf, "%d", lbl);
  label = buf;
}

bool TreeNode ::IsMulfurcate() {
  if (IsLeaf() == true) {
    return false;
  } else {
    if (GetChildrenNum() > 2) {
      return true;
    }
    for (int ii = 0; ii < GetChildrenNum(); ++ii) {
      if (GetChild(ii)->IsMulfurcate() == true) {
        return true;
      }
    }

    return false;
  }
}

TreeNode *TreeNode ::GetRoot() const {
  TreeNode *pself = const_cast<TreeNode *>(this);
  TreeNode *proot = pself;
  while (proot->GetParent() != NULL) {
    proot = proot->GetParent();
  }
  YW_ASSERT_INFO(proot != NULL, "Root is null");
  return proot;
}

void TreeNode ::GetAllAncestors(set<TreeNode *> &listAncestors) {
  if (GetParent() != NULL) {
    listAncestors.insert(GetParent());
    GetParent()->GetAllAncestors(listAncestors);
  }
}

void TreeNode ::GetAllChildren(set<TreeNode *> &setChildren) const {
  //
  // TreeNode *pthis = const_cast<TreeNode *>(this);
  // PopulateSetByVecGen( setChildren, pthis->listChildren );
  setChildren.clear();
  for (int i = 0; i < GetChildrenNum(); ++i) {
    setChildren.insert(listChildren[i]);
  }
}

int TreeNode ::GetChildIndex(TreeNode *pchild) const {
  // get the index of this particular child; if not found, the error
  TreeNode *pself = const_cast<TreeNode *>(this);
  int res = -1;
  for (int i = 0; i < (int)listChildren.size(); ++i) {
    if (pself->GetChild(i) == pchild) {
      res = i;
      break;
    }
  }
  YW_ASSERT_INFO(res >= 0, "Fail to find666");
  return res;
}

void TreeNode ::RemoveLabels() {
  // remove all edge labels (i.e. make them empty)
  int numLLs = listEdgeLabels.size();
  listEdgeLabels.clear();
  listEdgeLabels.resize(numLLs);

  // then reurrisve do it
  for (int i = 0; i < GetChildrenNum(); ++i) {
    GetChild(i)->RemoveLabels();
  }
}

void TreeNode ::RemoveLabelsPar() {
  // remove the parent to this node's label
  TreeNode *ppar = GetParent();
  if (ppar == NULL) {
    return;
  }
  int childIndex = ppar->GetChildIndex(this);
  YW_ASSERT_INFO(childIndex < (int)ppar->listEdgeLabels.size(), "Overflow");
  ppar->listEdgeLabels[childIndex].clear();
}

void TreeNode ::IncEdgeLabelsBy(int offset, bool fSub) {
  //
  for (int i = 0; i < (int)listEdgeLabels.size(); ++i) {
    for (int j = 0; j < listEdgeLabels[i].size(); ++j) {
      listEdgeLabels[i][j] += offset;
    }
  }
  if (fSub) {
    for (int i = 0; i < (int)listChildren.size(); ++i) {
      listChildren[i]->IncEdgeLabelsBy(offset, fSub);
    }
  }
}

void TreeNode ::Binarize(int &idToUseNext) {
  // recursively make the tree binary
  // if this node has more than 2 children, create a new internal node
  if (GetChildrenNum() > 2) {
    //
    TreeNode *pnode = new TreeNode(idToUseNext++);
    for (int i = 1; i < GetChildrenNum(); ++i) {
      vector<int> ss;
      pnode->AddChild(GetChild(i), ss);
    }
    TreeNode *pn1 = GetChild(0);
    this->listChildren.clear();
    this->listChildren.push_back(pn1);
    vector<int> ss;
    AddChild(pnode, ss);
  }

  for (int i = 0; i < GetChildrenNum(); ++i) {
    //
    GetChild(i)->Binarize(idToUseNext);
  }
}

int TreeNode ::GetMaxIdWithinSubtree() const {
  //
  int res = GetID();
  TreeNode *pthis = const_cast<TreeNode *>(this);
  for (int i = 0; i < GetChildrenNum(); ++i) {
    TreeNode *pnc = pthis->GetChild(i);
    int nc = pnc->GetMaxIdWithinSubtree();
    if (nc > res) {
      //
      res = nc;
    }
  }
  return res;
}

int TreeNode ::GetNumNodesUnder(bool fInternalOnly, bool fAddNonBinary) const {
  // fInternalOnly: true if only count internal node
  // include itself if this is an internal node
  // fAddNonBinary: true if an internal node is considered to have multiple
  // (hidden) nodes
  int res = 0;
  if (fInternalOnly == false || IsLeaf() == false) {
    res = 1;
  }
  // recursively check all children
  TreeNode *pn = const_cast<TreeNode *>(this);
  for (int i = 0; i < GetChildrenNum(); ++i) {
    res += pn->GetChild(i)->GetNumNodesUnder(fInternalOnly, fAddNonBinary);
  }
  return res;
}

// ***************************************************************************
// Utilites functions
// ***************************************************************************

void PhylogenyTreeIteratorBacktrack ::Init() {
  while (stackNodesToExplore.empty() == false) {
    stackNodesToExplore.pop();
  }
  // cout << "Nnow stack empty.\n";
  // Now recurisvely store the order of the walk
  TreeNode *rootNode = phyTree.GetRoot();
  if (rootNode != NULL) {
    stackNodesToExplore.push(rootNode);
  }
}

void PhylogenyTreeIteratorBacktrack ::Next() {
  if (stackNodesToExplore.empty() == true) {
    return;
  }
  TreeNode *pn = stackNodesToExplore.top();
  // push its descendent in
  stackNodesToExplore.pop();
  for (int i = 0; i < (int)pn->GetChildrenNum(); ++i) {
    //
    stackNodesToExplore.push(pn->GetChild(i));
  }
}
void PhylogenyTreeIteratorBacktrack ::Back() {
  if (stackNodesToExplore.empty() == true) {
    return;
  }
  // simply get rid of the current node
  stackNodesToExplore.pop();
}

bool PhylogenyTreeIteratorBacktrack ::IsDone() {
  return stackNodesToExplore.empty();
}

TreeNode *PhylogenyTreeIteratorBacktrack ::GetCurrNode() {
  if (IsDone() == false) {
    return stackNodesToExplore.top();
  } else {
    return NULL;
  }
}

///////////////////////////////////////////////////////////////////
void PhylogenyTreeIterator ::Init() {
  while (stackPostorder.empty() == false) {
    stackPostorder.pop();
  }
  // cout << "Nnow stack empty.\n";
  // Now recurisvely store the order of the walk
  TreeNode *rootNode = phyTree.GetRoot();
  if (rootNode != NULL) {
    phyTree.PostOrderPushStack(rootNode, stackPostorder);
  }
}

void PhylogenyTreeIterator ::Next() {
  if (stackPostorder.empty() == true) {
    return;
  }
  // TreeNode *pn = stackPostorder.top();
  stackPostorder.pop();
}

bool PhylogenyTreeIterator ::IsDone() { return stackPostorder.empty(); }

TreeNode *PhylogenyTreeIterator ::GetCurrNode() {
  if (IsDone() == false) {
    return stackPostorder.top();
  } else {
    return NULL;
  }
}

// ***************************************************************************
// Main functions
// ***************************************************************************

PhylogenyTreeBasic ::PhylogenyTreeBasic() : rootNode(NULL), numLeaves(-1) {}

PhylogenyTreeBasic ::~PhylogenyTreeBasic() {
  // cout << "Deleting tree: ";
  // Dump();

  // Should delete the tree
  if (rootNode != NULL) {
    delete rootNode;
    rootNode = NULL;
  }
}

PhylogenyTreeBasic *PhylogenyTreeBasic ::Copy() {
  PhylogenyTreeBasic *pCopy = new PhylogenyTreeBasic;
  pCopy->numLeaves = pCopy->numLeaves;
  pCopy->SetRoot(this->GetRoot()->Copy());
  return pCopy;
}

void PhylogenyTreeBasic ::PostOrderPushStack(
    TreeNode *treeNode, stack<TreeNode *> &stackPostorder) {
  stackPostorder.push(treeNode);
  // cout << "Pusing node " << treeNode->GetLabel() << endl;

  for (int i = 0; i < (int)treeNode->listChildren.size(); ++i) {
    PostOrderPushStack(treeNode->listChildren[i], stackPostorder);
  }
}

void PhylogenyTreeBasic ::ConsOnNewick(const string &nwString, int numLeaves,
                                       bool fBottomUp, TaxaMapper *pTMapper) {
  // Here we try to reconstruct from a newick string here
  // This function creates the tree by creating and linking tree nodes
  // Make sure the tree is empty
  if (rootNode != NULL) {
    delete rootNode;
    rootNode = NULL;
  }

  // we perform this by recursively
  int invId = 1000000;
  if (numLeaves > 0) {
    // here we assume leaf id starts from 0, will check it
    invId = numLeaves;
  }
  int leafId = 0;
  rootNode = ConsOnNewickSubtree(nwString, leafId, invId, numLeaves, fBottomUp,
                                 pTMapper);
}

void PhylogenyTreeBasic ::ConsOnNewickDupLabels(const string &nwString,
                                                TaxaMapper *pTMapper) {
  // Here we try to reconstruct from a newick string here
  // This function creates the tree by creating and linking tree nodes
  // Make sure the tree is empty
  if (rootNode != NULL) {
    delete rootNode;
    rootNode = NULL;
  }

  // we perform this by recursively
  int numLeaves = GetNewickNumLeaves(nwString);
  // we start counting leaves from 0
  int invId = numLeaves;
  int leafId = 0;
  // cout << "Num of leaves = " << numLeaves << endl;
  rootNode = ConsOnNewickSubtreeDupLabels(nwString, invId, leafId, pTMapper);
}

// ********************************************************************************
// Utitlieis for construcing edge label trees

static int GetEdgeLabelPosFrom(const string &strMutTreeCur, int posCur) {
  //
  int posCurGNTPF = posCur;
  while (posCurGNTPF < (int)strMutTreeCur.length()) {
    // printf "getNextTaxaPosFrom: %d: curr ch: %s\n", posCurGNTPF,
    // substr(strMutTreeCur,posCurGNTPF,1);
    if (strMutTreeCur[posCurGNTPF] == '#') {
      break;
    }
    ++posCurGNTPF;
  }
  if (posCurGNTPF >= (int)strMutTreeCur.length()) {
    posCurGNTPF = -1;
  }
  return posCurGNTPF;
}

static int getNextTaxaPosFromLevelUp(const string &strMutTreeCur, int posCur) {
  int posCurGNTPF = posCur;
  int level = 0;
  bool fUpperOnly = false;
  while (posCurGNTPF < (int)strMutTreeCur.length()) {
    char chGNTPF = strMutTreeCur[posCurGNTPF];
    if (chGNTPF == '#' && ((level >= 0 && fUpperOnly == false) || level > 0)) {
      break;
    }
    if (chGNTPF == '(') {
      --level;
    } else if (chGNTPF == ')') {
      ++level;
    } else if (chGNTPF == ',') {
      fUpperOnly = true;
    }

    ++posCurGNTPF;
  }
  if (posCurGNTPF >= (int)strMutTreeCur.length()) {
    posCurGNTPF = -1;
  }
  return posCurGNTPF;
}

static string getTaxaAt(const string &strMutTreeCur, int posCur) {
  int posGTA = posCur;
  if (strMutTreeCur[posCur] == '#') {
    posGTA = posCur + 1;
  }
  //  now find where it ends
  int posGTA2 = posGTA;
  while (posGTA2 < (int)strMutTreeCur.length()) {
    char chGTA = strMutTreeCur[posGTA2];
    if (chGTA == '#' || chGTA == ',' || chGTA == ')') {
      break;
    }
    ++posGTA2;
  }
  if (posGTA2 > (int)strMutTreeCur.length()) {
    posGTA2 = (int)strMutTreeCur.length() - 1;
  }
  return strMutTreeCur.substr(posGTA, posGTA2 - posGTA);
}

void PhylogenyTreeBasic ::ConsOnNewickEdgeLabelTree(const string &nwString) {
  // view each edge label as taxon; a stand-alone edge label is the leaf;
  // edge label may or may not have a leading seperator (# in this
  // implementation); e.g. ((#1,#2#3)#4)  this give four node, one for each edge
  // label
  if (rootNode != NULL) {
    delete rootNode;
    rootNode = NULL;
  }
  // find all edge labels and how they are related
  map<string, string> mapEdgeLabelPar;
  int posEdgeLbl = 0;
  while (posEdgeLbl < (int)nwString.length()) {
    //
    posEdgeLbl = GetEdgeLabelPosFrom(nwString, posEdgeLbl);
    if (posEdgeLbl < 0) {
      break;
    }
    string strTaxon = getTaxaAt(nwString, posEdgeLbl);
    // find its parent
    int posEdgeLblPar = getNextTaxaPosFromLevelUp(nwString, posEdgeLbl + 1);
    string strPar;
    if (posEdgeLblPar >= 0) {
      //
      strPar = getTaxaAt(nwString, posEdgeLblPar);
    }
    mapEdgeLabelPar[strTaxon] = strPar;
    // cout << "Taxon: " << strTaxon << " is child of " << strPar << endl;
    ++posEdgeLbl;
  }
  // now create nodes
  int nidNext = 1;
  this->rootNode = new TreeNode(nidNext++);
  string strLblRoot = "-";
  int posRootLbl = -1;
  std::size_t pos1 = nwString.find_last_of(')');
  std::size_t pos2 = nwString.find_last_of('#');
  if (pos1 != string::npos && pos2 != string::npos) {
    posRootLbl = max(pos1, pos2);
  } else if (pos1 != string::npos) {
    posRootLbl = pos1;
  } else if (pos2 != string::npos) {
    posRootLbl = pos2;
  }
  if (posRootLbl >= 0) {
    strLblRoot = getTaxaAt(nwString, posRootLbl);
  }

  // cout << "root label: " << strLblRoot << endl;
  // now create all descendents
  map<string, TreeNode *> mapNodes;
  mapNodes[strLblRoot] = this->rootNode;
  while (true) {
    // find direct descendents
    TreeNode *pnPar = NULL;
    string strChildUse;
    for (map<string, string>::iterator it = mapEdgeLabelPar.begin();
         it != mapEdgeLabelPar.end(); ++it) {
      string strChild = it->first;
      string strPar = it->second;
      if (mapNodes.find(strChild) == mapNodes.end() &&
          mapNodes.find(strPar) != mapNodes.end()) {
        pnPar = mapNodes[strPar];
        strChildUse = strChild;
      }
    }
    if (pnPar == NULL) {
      break;
    }
    TreeNode *pnode = new TreeNode(nidNext++);
    pnode->SetLabel(strChildUse);
    vector<int> listLblsDummy;
    pnPar->AddChild(pnode, listLblsDummy);

    mapNodes[strChildUse] = pnode;
  }

  if (strLblRoot.length() == 0) {
    strLblRoot = "-";
  }
  this->rootNode->SetLabel(strLblRoot);
}

void PhylogenyTreeBasic ::InitPostorderWalk() {
  // cout << "InitPostorderWalk() entry\n";
  // when walk, return the value of the node if any
  // Clearup the previous storage if any
  while (stackPostorder.empty() == false) {
    stackPostorder.pop();
  }
  // cout << "Nnow stack empty.\n";
  // Now recurisvely store the order of the walk
  if (rootNode != NULL) {
    PostOrderPushStack(rootNode, stackPostorder);
  }
}

TreeNode *PhylogenyTreeBasic ::NextPostorderWalk() {
  // Return false, when nothing to go any more
  if (stackPostorder.empty() == true) {
    return NULL;
  }
  TreeNode *pn = stackPostorder.top();
  stackPostorder.pop();

//    node = pn;
#if 0
    if( pn->nodeValues.size() > 0 )
    {
        // There is valid node value stored here
        nodeValue = pn->nodeValues[0];
    }
    else
    {
        nodeValue = -1;     // no node value is stored here
    }
#endif
  return pn;
}

void PhylogenyTreeBasic ::OutputGML(const char *inFileName) {
  // Now output a file in GML format
  // First create a new name
  string name = inFileName;
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
  OutputQuotedString(outFile, "Phylogeny Tree....\n");

  // Now output all the vertices
  //	int i;
  stack<TreeNode *> nodesStack;
  if (rootNode != NULL) {
    nodesStack.push(rootNode);
  }
  // cout << "a.1.1\n";
  while (nodesStack.empty() == false) {
    TreeNode *pn = nodesStack.top();
    nodesStack.pop();

    outFile << "node [\n";

    outFile << "id " << pn->id << endl;
    outFile << "label ";
    string nameToUse = " ";
    if (pn->GetLabel() != "-") {
      nameToUse = pn->GetLabel();
    }
#if 0
        else
        {
            // we take the nonde value here
            char buf[100];
            if( pn->nodeValues.size() > 0 )
            {
                sprintf(buf, "(%d)", pn->nodeValues[0] );        // CAUTION, here we assume each leaf has exactly 1 label
                nameToUse = buf;
            }
            else
            {
                // if no nodes value is set, still use label
         //       nameToUse = pn->GetLabel();

                // YW: TBD change
                nameToUse.empty();
            }
        }
#endif
    const char *name = nameToUse.c_str();

    // 		char name[100];
    //       if( pn->IsLeaf() == false)
    //        {
    //		    name[0] = 'v';
    //		    sprintf(&name[1], "%d", pn->id);
    //        }
    //        else
    //        {
    // For leaf, we simply output their value (row number)
    //            sprintf(name, "%d", pn->nodeValues[0] );        // CAUTION,
    //            here we assume each leaf has exactly 1 label
    //        }
    OutputQuotedString(outFile, name);
    outFile << endl;

    // See if we need special shape here
    if (pn->GetShape() == PHY_TN_RECTANGLE) {
      outFile << "vgj [ \n shape  ";
      OutputQuotedString(outFile, "Rectangle");
      outFile << "\n]\n";
    } else {
      outFile << "defaultAtrribute   1\n";
    }

    outFile << "]\n";

    // Now try to get more nodes
    for (int i = 0; i < (int)pn->listChildren.size(); ++i) {
      nodesStack.push(pn->listChildren[i]);
    }
    // cout << "a.1.2\n";
  }
  // cout << "a.1.3\n";

  // Now output all the edges, by again starting from root and output all nodes
  YW_ASSERT(nodesStack.empty() == true);
  if (rootNode != NULL) {
    nodesStack.push(rootNode);
  }
  while (nodesStack.empty() == false) {
    TreeNode *pn = nodesStack.top();
    nodesStack.pop();

    for (int i = 0; i < (int)pn->listChildren.size(); ++i) {

      // cout << "Output an edge \n";
      outFile << "edge [\n";
      outFile << "source " << pn->id << endl;
      outFile << "target  " << pn->listChildren[i]->id << endl;
      outFile << "label ";
      if (pn->listEdgeLabels[i].size() > 0) {
        string lblName;
        char name[100];
        //		        name[0] = 'e';
        for (int iel = 0; iel < (int)pn->listEdgeLabels[i].size(); ++iel) {
          sprintf(name, "e%d  ", pn->listEdgeLabels[i][iel]);
          lblName += name;
        }
        OutputQuotedString(outFile, lblName.c_str());
      } else {
        OutputQuotedString(outFile, "");
      }
      outFile << "\n";
      outFile << "]\n";

      // Store next one to stack
      nodesStack.push(pn->listChildren[i]);
    }
  }

  // Finally quite after closing file
  outFile << "\n]\n";
  outFile.close();
}

void PhylogenyTreeBasic ::OutputGMLNoLabel(const char *inFileName) {
  //
  // Now output a file in GML format
  // First create a new name
  string name = inFileName;
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
  OutputQuotedString(outFile, "Phylogeny Tree....\n");

  // Now output all the vertices
  //	int i;
  stack<TreeNode *> nodesStack;
  if (rootNode != NULL) {
    nodesStack.push(rootNode);
  }
  // cout << "a.1.1\n";
  while (nodesStack.empty() == false) {
    TreeNode *pn = nodesStack.top();
    nodesStack.pop();

    outFile << "node [\n";

    outFile << "id " << pn->id << endl;
    outFile << "label ";
    string nameToUse = " ";
    const char *name = nameToUse.c_str();

    // 		char name[100];
    //       if( pn->IsLeaf() == false)
    //        {
    //		    name[0] = 'v';
    //		    sprintf(&name[1], "%d", pn->id);
    //        }
    //        else
    //        {
    // For leaf, we simply output their value (row number)
    //            sprintf(name, "%d", pn->nodeValues[0] );        // CAUTION,
    //            here we assume each leaf has exactly 1 label
    //        }
    OutputQuotedString(outFile, name);
    outFile << endl;

    // See if we need special shape here
    if (pn->GetShape() == PHY_TN_RECTANGLE) {
      outFile << "vgj [ \n shape  ";
      OutputQuotedString(outFile, "Rectangle");
      outFile << "\n]\n";
    } else {
      outFile << "defaultAtrribute   1\n";
    }

    outFile << "]\n";

    // Now try to get more nodes
    for (int i = 0; i < (int)pn->listChildren.size(); ++i) {
      nodesStack.push(pn->listChildren[i]);
    }
    // cout << "a.1.2\n";
  }
  // cout << "a.1.3\n";

  // Now output all the edges, by again starting from root and output all nodes
  YW_ASSERT(nodesStack.empty() == true);
  if (rootNode != NULL) {
    nodesStack.push(rootNode);
  }
  while (nodesStack.empty() == false) {
    TreeNode *pn = nodesStack.top();
    nodesStack.pop();

    for (int i = 0; i < (int)pn->listChildren.size(); ++i) {

      // cout << "Output an edge \n";
      outFile << "edge [\n";
      outFile << "source " << pn->id << endl;
      outFile << "target  " << pn->listChildren[i]->id << endl;
      outFile << "label ";
      if (pn->listEdgeLabels[i].size() > 0) {
        string lblName;
        char name[100];
        //		        name[0] = 'e';
        for (int iel = 0; iel < (int)pn->listEdgeLabels[i].size(); ++iel) {
          sprintf(name, "e%d  ", pn->listEdgeLabels[i][iel]);
          lblName += name;
        }
        OutputQuotedString(outFile, lblName.c_str());
      } else {
        OutputQuotedString(outFile, "");
      }
      outFile << "\n";
      outFile << "]\n";

      // Store next one to stack
      nodesStack.push(pn->listChildren[i]);
    }
  }

  // Finally quite after closing file
  outFile << "\n]\n";
  outFile.close();
}

// construct a newick string for this tree
void PhylogenyTreeBasic ::ConsNewick(string &strNewick, bool wGridLen,
                                     double gridWidth, bool fUseCurLbl) {
  strNewick.empty();

  // work from this node
  YW_ASSERT_INFO(rootNode != NULL, "Root is not set");
  strNewick =
      ConsNewickTreeNode(rootNode, wGridLen, gridWidth, fUseCurLbl, false);
}

void PhylogenyTreeBasic ::ConsNewickSorted(string &strNewick, bool wGridLen,
                                           double gridWidth, bool fUseCurLbl) {
  strNewick.empty();

  // work from this node
  YW_ASSERT_INFO(rootNode != NULL, "Root is not set");
  strNewick =
      ConsNewickTreeNode(rootNode, wGridLen, gridWidth, fUseCurLbl, true);
}

void PhylogenyTreeBasic ::ConsNewickEdgeLabel(string &strNewick) {
  strNewick.empty();

  // work from this node
  YW_ASSERT_INFO(rootNode != NULL, "Root is not set");
  strNewick = ConsNewickTreeNode(rootNode, false, 1.0, true, true, true);
}

string PhylogenyTreeBasic ::ConsNewickTreeNode(TreeNode *pNode, bool wGridLen,
                                               double gridWidth,
                                               bool fUseCurLbl, bool fSort,
                                               bool fEdgeLbel) {
  // cout << "--------------------------------In ConsNewickTreeNode: I am
  // here\n";
  string resNodeStr;
  // Is this node a leaf? If so, we output the label of it
  if (pNode->IsLeaf() == true) {
    // Add this label if this label is not there
    string tmpstr = pNode->GetUserLabel();
    if (fUseCurLbl == true) {
      tmpstr = pNode->GetLabel();
    }
    resNodeStr = tmpstr;
  } else {
    string tmpstr = pNode->GetLabel();
    YW_ASSERT_INFO(pNode->listChildren.size() >= 1,
                   "Must have some children here.");

    // When there is only one child and no self-label
    if (tmpstr.size() <= 2 && pNode->listChildren.size() == 1) {
      resNodeStr = ConsNewickTreeNode(pNode->listChildren[0], wGridLen,
                                      gridWidth, fUseCurLbl, fSort, fEdgeLbel);
    } else {

      // Otherwise, we simply collect all sub strings here, and sepearate by a ,
      string comboStrName = "(";

      bool fAddSep = false;
      // does this node has a label by itself? if so, output it
      if (tmpstr.size() > 2) {
        comboStrName += tmpstr.substr(1, tmpstr.size() - 2);
        // comboStrName += ",";

        // all others should be added sep.
        fAddSep = true;
      }

      // handle its children
      if (fSort == false) {
        for (unsigned int i = 0; i < pNode->listChildren.size(); ++i) {
          string stepRes =
              ConsNewickTreeNode(pNode->listChildren[i], wGridLen, gridWidth,
                                 fUseCurLbl, fSort, fEdgeLbel);

          if (stepRes.size() > 0) {
            if (fAddSep == true) {
              comboStrName += ",";
            }

            comboStrName += stepRes;

            // from now on, add sep
            fAddSep = true;

            // if( i+1 < pNode->listChildren.size() )
            //{
            //    comboStrName += ",";
            //}
          }
        }
      } else {
        // sort the labels from children
        multiset<string> strsChildren;
        for (unsigned int i = 0; i < pNode->listChildren.size(); ++i) {
          string stepRes =
              ConsNewickTreeNode(pNode->listChildren[i], wGridLen, gridWidth,
                                 fUseCurLbl, fSort, fEdgeLbel);
          if (stepRes.size() > 0) {
            strsChildren.insert(stepRes);
          }
        }
        for (multiset<string>::iterator it = strsChildren.begin();
             it != strsChildren.end(); ++it) {
          //
          if (fAddSep == true) {
            comboStrName += ",";
          }

          comboStrName += *it;

          // from now on, add sep
          fAddSep = true;
        }
      }
      comboStrName += ")";
      // cout << "comboStrName = " << comboStrName << endl;
      resNodeStr = comboStrName;
    }
  }

  // now see if we need to add length info
  //
  if (wGridLen == true) {
    //
    TreeNode *pNodePar = pNode->GetParent();
    if (pNodePar != NULL) {
      double len = gridWidth * (pNodePar->GetLevel() - pNode->GetLevel());
      // cout << "**************************PhylogenyTreeBasic::len = " << len
      // << endl;
      char buf[100];
      sprintf(buf, ":%f", len);
      resNodeStr += buf;
    }
  } else if (pNode->GetLength() >= 0.0) {
#if 0
        // if length is set, add it
        resNodeStr += ":";
        resNodeStr += ConvToString(pNode->GetLength() );
#endif
  }

  if (fEdgeLbel) {
    TreeNode *pParNode = pNode->GetParent();
    if (pParNode != NULL) {
      int cIndex = pParNode->GetChildIndex(pNode);

      // add edge label in the format: s1s2s3....
      string strEdgeLbel;
      vector<int> listEdgeLabels;
      pParNode->GetEdgeLabelsAtBranch(cIndex, listEdgeLabels);

      // cout << "cIndex: " << cIndex <<", listEdgeLabels: ";
      // DumpIntVec(listEdgeLabels);

      for (int i = 0; i < (int)listEdgeLabels.size(); ++i) {
        char buf[10000];
        sprintf(buf, "#%d", listEdgeLabels[i]);
        strEdgeLbel += buf;
      }
      if (strEdgeLbel.length() > 0) {
        resNodeStr += ":";
        resNodeStr += strEdgeLbel;
      }
    }
  }

  return resNodeStr;
}

// This function adds a new tree node, and return it. Also set the parent node
// to the pareamter
TreeNode *PhylogenyTreeBasic ::AddTreeNode(TreeNode *parNode, int id) {
  if (id < 0) {
    id = GetNumVertices();
  }

  TreeNode *pnode = new TreeNode(id);
  pnode->AddNodeValue(id);

  // Should delete the tree
  if (parNode == NULL) {
    YW_ASSERT_INFO(
        rootNode == NULL,
        "Can not add a node with no parent if the tree is not empty");
    rootNode = pnode;
    return pnode;
  }

  // Otherwise, set the parent
  SEQUENCE emptySeq;
  parNode->AddChild(pnode, emptySeq);
  return pnode;
}

int PhylogenyTreeBasic ::GetNumVertices() const {
  int res = 0;
  stack<TreeNode *> stackNodes;
  if (rootNode != NULL) {
    stackNodes.push(rootNode);
  }
  while (stackNodes.empty() == false) {
    TreeNode *pcurr = stackNodes.top();
    stackNodes.pop();
    ++res;
    // Now enque its children
    for (int i = 0; i < (int)pcurr->listChildren.size(); ++i) {
      stackNodes.push(pcurr->listChildren[i]);
    }
  }
  return res;
}

// int PhylogenyTreeBasic :: GetIdFromStr( const string &strPart, TaxaMapper
// *pTMapper )
//{
// cout << "GetIdFromStr: " << strPart << endl;
//	string strToUse = strPart;
//	size_t posSeparator = strPart.find( ':' );
//	if( posSeparator != string::npos )
//	{
//		strToUse = strPart.substr(0, (int)posSeparator  );
//	}
//	// get rid of
//	int res = -1;
//	if( pTMapper == NULL)
//	{
//		sscanf( strToUse.c_str(), "%d", &res  );
// cout << "Empty mapper\n";
//	}
//	else
//	{
//		// are we reading in the first tree or not
//		if( pTMapper->IsInitialized() == true )
//		{
//			res  = pTMapper->GetId(strToUse);
// cout << "GetIdFromStr: GetId: " << strToUse << ": " << res << endl;
//		}
//		else
//		{
//			// this is new
//			res = pTMapper->AddTaxaString( strToUse );
// cout << "GetIdFromStr: New id: " << strToUse << ": " << res << endl;
//		}
//	}
//	return res;
//}

TreeNode *PhylogenyTreeBasic ::ConsOnNewickSubtree(const string &nwStringPart,
                                                   int &leafId, int &invId,
                                                   int numLeaves,
                                                   bool fBottomUp,
                                                   TaxaMapper *pTMapper) {
  // cout << "Entry nwStringPart = "<< nwStringPart << endl;

  TreeNode *pres = NULL;
  int posLenBegin = -1;

  // this function builds recursively subtrees for this part of string
  // First, is this string a leaf or not
  if (nwStringPart[0] != '(') {
    // TreeNode *pLeaf = new TreeNode( nodeId  );
    //// also set its label this way
    // pLeaf->AddNodeValue( nodeId );

    // 7/27/10 YW: for now, we take this convention:
    // tree node id = label  if no mapper is passed
    // Why? This case is by default for internal use only
    // while mapper is used for external (user) specified
    // Yes, this is a leaf
    int nodeId = TaxaMapper ::GetIdFromStr(nwStringPart, pTMapper);
    //	sscanf( nwStringPart.c_str(), "%d", &nodeId  );

    if (numLeaves > 0) {
      if (nodeId >= numLeaves) {
        cout << "Wrong: nodeId = " << nodeId << ", numLeaves = " << numLeaves
             << endl;
      }
      YW_ASSERT_INFO(nodeId < numLeaves,
                     "We assume in phylogeny tree, leaf id starts from 0");
    }
    // cout << "node id = " << nodeId << endl;

    int idtouse = leafId;
    if (pTMapper == NULL) {
      // in this case take the same as node id
      idtouse = nodeId;
    } else {
      // update leafid since we are using it
      leafId++;
    }

    TreeNode *pLeaf = new TreeNode(idtouse);
    // also set its label this way
    pLeaf->AddNodeValue(idtouse);
    // leafId ++;

    // get rid of any part after : if there is length info
    // string strLeafLabel = nwStringPart;
    // if( strLa )
    //{
    //}
    string strLbl = GetStringFromId(nodeId);
    pLeaf->SetLabel(strLbl);

    string strLblUser = TaxaMapper ::ExtractIdPartFromStr(nwStringPart);
    pLeaf->SetUserLabel(strLblUser);

    // cout << "ConsOnNewickSubtree: set leaf label: " << strLbl << endl;
    // return pLeaf;
    pres = pLeaf;

    size_t posLenSep = nwStringPart.find(':');
    if (posLenSep != string::npos) {
      //
      posLenBegin = posLenSep + 1;
    }
  } else {
    // This is not a leaf
    // so we create underlying level for it
    int idToUse = 1000;
    if (fBottomUp == false) {
      idToUse = invId++;
    }
    TreeNode *pInternal = new TreeNode(idToUse);
    int lastpos = 1;
    int curpos = 0;
    int parnet = 0; // (: +1, ) -1
    while (true) {
      // cout << "curpos = " << curpos << endl;

      if (curpos >= (int)nwStringPart.size()) {
        // we are done
        break;
      }

      // keep balance
      if (nwStringPart[curpos] == '(') {
        parnet++;
      } else if (nwStringPart[curpos] == ')') {
        parnet--;

        // when parnet = 0, we know we end
        if (parnet == 0) {
          // now adding the last piece
          // create a new node
          int strl = curpos - lastpos;
          string subs = nwStringPart.substr(lastpos, strl);
          //    cout << "last subs = " << subs << endl;
          TreeNode *pChild = ConsOnNewickSubtree(subs, leafId, invId, numLeaves,
                                                 fBottomUp, pTMapper);

          // also append it as child
          vector<int> empytLabels;
          pInternal->AddChild(pChild, empytLabels);

          // aslo update lastpos
          lastpos = curpos + 1;
        }

      } else if (nwStringPart[curpos] == ',') {
        // Yes, this is a sepeartor, but we only start to process it when the
        // balance of parenetnis is right
        if (parnet == 1) {
          // create a new node
          int strl = curpos - lastpos;
          string subs = nwStringPart.substr(lastpos, strl);
          //    cout << "subs = " << subs << endl;
          TreeNode *pChild = ConsOnNewickSubtree(subs, leafId, invId, numLeaves,
                                                 fBottomUp, pTMapper);

          // also append it as child
          vector<int> empytLabels;
          pInternal->AddChild(pChild, empytLabels);

          // aslo update lastpos
          lastpos = curpos + 1;
        }
      } else if (nwStringPart[curpos] == ':') {
        // keep track of length
        if (parnet == 0) {
          posLenBegin = curpos + 1;
        }
      }

      // now move to next pos
      curpos++;
    }

    // if we go bottom up labeling the node, we should re-label the node here
    if (fBottomUp == true) {
      pInternal->SetID(invId++);
    }
    // return pInternal;
    pres = pInternal;
  }

  //
  if (posLenBegin >= 0) {
    // also read in length
    size_t posRightExt = nwStringPart.find(')', posLenBegin);
    int rightPos = (int)nwStringPart.size() - 1;
    if (posRightExt != string::npos) {
      rightPos = posRightExt - 1;
    }
    string subs =
        nwStringPart.substr(posLenBegin, posRightExt - posLenBegin + 1);
    double len = StrToDouble(subs);
    pres->SetLength(len);
  }
  return pres;
}

TreeNode *PhylogenyTreeBasic ::ConsOnNewickSubtreeDupLabels(
    const string &nwStringPart, int &invId, int &leafId, TaxaMapper *pTMapper) {
  // cout << "Entry nwStringPart = "<< nwStringPart << endl;

  // this function builds recursively subtrees for this part of string
  // First, is this string a leaf or not
  if (nwStringPart[0] != '(') {
    // ensure no internal has every been set yet
    // YW_ASSERT_INFO( invId < 0, "invId should not be set when leaf is being
    // processed" );

    // Yes, this is a leaf
    int nodeId = leafId;
    leafId++;
    int leafLabel = TaxaMapper ::GetIdFromStr(nwStringPart, pTMapper);
    // sscanf( nwStringPart.c_str(), "%d", &leafLabel  );

    // cout << "leaf id = " << nodeId << endl;
    TreeNode *pLeaf = new TreeNode(nodeId);
    // also set its label this way
    pLeaf->AddNodeValue(nodeId);

    // get rid of any part after : if there is length info
    // string strLeafLabel = nwStringPart;
    // if( strLa )
    //{
    //}
    char buf[1000];
    sprintf(buf, "%d", leafLabel);
    string strLabel = buf;
    pLeaf->SetLabel(strLabel);

    string strLabelUser = TaxaMapper ::ExtractIdPartFromStr(nwStringPart);
    pLeaf->SetUserLabel(strLabelUser);

    // cout << "ConsOnNewickSubtree: set leaf label: " << strLabel << endl;
    return pLeaf;
  } else {

    // This is not a leaf
    // so we create underlying level for it
    int idToUse = invId;
    TreeNode *pInternal = new TreeNode(idToUse);
    int lastpos = 1;
    int curpos = 0;
    int parnet = 0; // (: +1, ) -1
    while (true) {
      // cout << "curpos = " << curpos << endl;

      if (curpos >= (int)nwStringPart.size()) {
        // we are done
        break;
      }

      // keep balance
      if (nwStringPart[curpos] == '(') {
        parnet++;
      } else if (nwStringPart[curpos] == ')') {
        parnet--;

        // when parnet = 0, we know we end
        if (parnet == 0) {
          // now adding the last piece
          // create a new node
          int strl = curpos - lastpos;
          string subs = nwStringPart.substr(lastpos, strl);
          //    cout << "last subs = " << subs << endl;
          TreeNode *pChild =
              ConsOnNewickSubtreeDupLabels(subs, invId, leafId, pTMapper);

          // also append it as child
          vector<int> empytLabels;
          pInternal->AddChild(pChild, empytLabels);

          // aslo update lastpos
          lastpos = curpos + 1;
        }

      } else if (nwStringPart[curpos] == ',') {
        // Yes, this is a sepeartor, but we only start to process it when the
        // balance of parenetnis is right
        if (parnet == 1) {
          // create a new node
          int strl = curpos - lastpos;
          string subs = nwStringPart.substr(lastpos, strl);
          //    cout << "subs = " << subs << endl;
          TreeNode *pChild =
              ConsOnNewickSubtreeDupLabels(subs, invId, leafId, pTMapper);

          // also append it as child
          vector<int> empytLabels;
          pInternal->AddChild(pChild, empytLabels);

          // aslo update lastpos
          lastpos = curpos + 1;
        }
      }

      // now move to next pos
      curpos++;
    }

    // if we go bottom up labeling the node, we should re-label the node here
    // if(invId < 0 )
    //{
    //	invId = leafId;
    //}

    pInternal->SetID(invId++);
    // cout << "Set internal node to " << pInternal->GetID() << endl;
    return pInternal;
  }
}

// Get nodes info
// 7/27/10: we want to get node label (NOT id!)
void PhylogenyTreeBasic ::GetNodeParInfo(vector<int> &nodeIds,
                                         vector<int> &parPos) {
  // cout << "GetNodeParInfo: \n";
  // simply put consecutive node ids but keep track of node parent positions
  // ensure we get the correct node mapping between id and pointer to node
  map<TreeNode *, int> mapNodeIds;

  // id is simply consecutive
  int numTotVerts = GetNumVertices();
  nodeIds.resize(numTotVerts);
  for (int i = 0; i < numTotVerts; ++i) {
    nodeIds[i] = i;
  }
  parPos.resize(numTotVerts);
  for (int i = 0; i < numTotVerts; ++i) {
    parPos[i] = -1;
  }

  // IMPORTANT: assume binary tree, otherwise all bets are off!!!!
  // int numLeaves = ( numTotVerts+1 )/2;
  int numLeaves = GetNumLeaves();
  // cout << "numLeaves: " << numLeaves << endl;
  // do traversal
  int curNodeNum = 0;
  // InitPostorderWalk();
  PhylogenyTreeIterator itorTree(*this);
  itorTree.Init();
  while (itorTree.IsDone() == false) {
    TreeNode *pn = itorTree.GetCurrNode();
    itorTree.Next();
    // TreeNode *pn = NextPostorderWalk( ) ;
    if (pn == NULL) {
      // cout << "No node here. Stop.\n";
      break; // done with all nodes
    }

    //
    if (pn->IsLeaf() == true) {
      // skip it for now
      continue;
    }

    //
    int nonleafInd = numLeaves + curNodeNum;
    curNodeNum++;
    // remember it
    mapNodeIds.insert(map<TreeNode *, int>::value_type(pn, nonleafInd));
    // now set its descendents to this index, either leaf or non-leaf
    // if it is non-leaf, do a lookup of the stored id. Leaf: just go by its id
    for (int jj = 0; jj < pn->GetChildrenNum(); ++jj) {
      TreeNode *pnjj = pn->GetChild(jj);
      int pnjjid;
      int pnjjlabel = -1;
      if (pnjj->IsLeaf() == true) {
        pnjjid = pnjj->GetID();
        // assume id is distinct, while label can be duplicate
        pnjjlabel = pnjj->GetIntLabel();
        // cout << "pnjjid = " << pnjjid << ", pnjjlabel: " << pnjjlabel << ",
        // numLeaves: " << numLeaves << endl;
        YW_ASSERT_INFO(pnjjid >= 0 && pnjjid < numLeaves,
                       "Leaf id: out of range");
      } else {
        YW_ASSERT_INFO(mapNodeIds.find(pnjj) != mapNodeIds.end(),
                       "Fail to find the node");
        pnjjid = mapNodeIds[pnjj];
      }
      parPos[pnjjid] = nonleafInd;
      // this says whether we change the label of the node
      // this is needed when there are duplicate labels in the tree
      if (pnjjlabel >= 0) {
        nodeIds[pnjjid] = pnjjlabel;
      }
    }
  }

  // print out
  // cout << "original tree:  ";
  // string strTree;
  // ConsNewick(strTree);
  // cout << strTree << endl;
  // cout << "Parent position : ";
  // DumpIntVec( parPos );
}

void PhylogenyTreeBasic ::GetNodeParInfoNew(vector<int> &nodeIds,
                                            vector<int> &parPos) {
  // cout << "In GetNodeParInfoNew: tree is: ";
  // this->Dump();
  // the previous version has various of problems, but it is being used by some
  // programs so I decide to add a new function Note this one assume all nodes
  // are labeled consecutively simply put consecutive node ids but keep track of
  // node parent positions ensure we get the correct node mapping between id and
  // pointer to node
  // map<TreeNode *,int> mapNodeIds;

  // id is simply consecutive
  int numTotVerts = GetNumVertices();
  // nodeIds.resize(numTotVerts);
  // for(int i=0; i<numTotVerts; ++i)
  //{
  //	nodeIds[i] = i;
  //}
  // parPos.resize(numTotVerts);
  // for(int i=0; i<numTotVerts; ++i)
  //{
  //	parPos[i] = -1;
  //}

  // IMPORTANT: assume binary tree, otherwise all bets are off!!!!
  // int numLeaves = ( numTotVerts+1 )/2;
  int numLeaves = GetNumLeaves();
  // cout << "Numleaves = " << numLeaves << endl;
  // do traversal
  // int curNodeNum = 0;
  // InitPostorderWalk();
  PhylogenyTreeIterator itorTree(*this);
  itorTree.Init();
  while (itorTree.IsDone() == false) {
    TreeNode *pn = itorTree.GetCurrNode();
    itorTree.Next();
    // TreeNode *pn = NextPostorderWalk( ) ;
    if (pn == NULL) {
      // cout << "No node here. Stop.\n";
      break; // done with all nodes
    }

    //
    int curNodeId = pn->GetID();
    // cout << "curNodeId: " << curNodeId << endl;
    YW_ASSERT_INFO(
        curNodeId < numTotVerts,
        "curNodeId exceeds limit (the node ids must be consecutive from 0)");
    if (pn->IsLeaf() == true) {
      // skip it for now
      YW_ASSERT_INFO(
          curNodeId < numLeaves,
          "The tree violates assumption that tree leaf id start from 0");
    }

    // add a record
    nodeIds.push_back(pn->GetID());
    TreeNode *pnPar = pn->GetParent();
    if (pnPar == NULL) {
      parPos.push_back(-1);
    } else {
      // simply its id
      parPos.push_back(pnPar->GetID());
    }

    //	continue;
    //}
#if 0
		//
		//int nonleafInd = numLeaves + curNodeNum;
		int nonleafInd = curNodeId;
		//curNodeNum++;
		// remember it
		mapNodeIds.insert( map<TreeNode *,int> :: value_type( pn, curNodeId ) );
		// now set its descendents to this index, either leaf or non-leaf
		// if it is non-leaf, do a lookup of the stored id. Leaf: just go by its id
		for(int jj=0; jj<pn->GetChildrenNum(); ++jj)
		{
			TreeNode *pnjj = pn->GetChild(jj);
			int pnjjid;
			if( pnjj->IsLeaf() == true )
			{
				pnjjid = pnjj->GetID();
				YW_ASSERT_INFO( pnjjid >=0 && pnjjid < numLeaves, "Leaf id: out of range" );
			}
			else
			{
				YW_ASSERT_INFO( mapNodeIds.find( pnjj ) != mapNodeIds.end(), "Fail to find the node"  );
				pnjjid = mapNodeIds[pnjj];
			}
			parPos[pnjjid] = nonleafInd;
#endif
    //}
  }

  // print out
  // cout << "original tree:  ";
  // string strTree;
  // ConsNewick(strTree);
  // cout << strTree << endl;
  // cout << "Parent position : ";
  // DumpIntVec( parPos );
}

//
bool PhylogenyTreeBasic ::ConsOnParPosList(const vector<int> &parPos,
                                           int numLeaves, bool fBottupUpLabel) {
  //
  string strNewick;
  if (ConvParPosToNewick(parPos, strNewick) == false) {
    return false;
  }
  // cout << "Newick string = " << strNewick << endl;
  ConsOnNewick(strNewick, numLeaves, fBottupUpLabel);
  return true;
}

bool PhylogenyTreeBasic ::ConvParPosToNewick(const vector<int> &parPos,
                                             string &strNewick) {
  // convert par position representation to newick
  // we always assume the last item is -1
  YW_ASSERT_INFO(parPos[parPos.size() - 1] == -1,
                 "Must be -1 for the last value in parPos");
  ConvParPosToNewickSubtree(parPos.size() - 1, parPos, strNewick);
  return true;
}

void PhylogenyTreeBasic ::ConvParPosToNewickSubtree(int nodeInd,
                                                    const vector<int> &parPos,
                                                    string &strNewick) {
  // this function generate under a single node (leaf or non-leaf), the newick
  // under the subtree
  vector<int> listUnderNodeInds;
  for (int i = 0; i < (int)parPos.size(); ++i) {
    if (parPos[i] == nodeInd) {
      listUnderNodeInds.push_back(i);
    }
  }
  // leaf if empty
  if (listUnderNodeInds.size() == 0) {
    char buf[100];
    sprintf(buf, "%d", nodeInd);
    strNewick = buf;
    return;
  }
  YW_ASSERT_INFO(listUnderNodeInds.size() == 2,
                 "Only binary trees are supported for now");

  // now get newick for the two part and merge it
  string strFirst, strSecond;
  ConvParPosToNewickSubtree(listUnderNodeInds[0], parPos, strFirst);
  ConvParPosToNewickSubtree(listUnderNodeInds[1], parPos, strSecond);
  strNewick = "(";
  strNewick += strFirst;
  strNewick += ",";
  strNewick += strSecond;
  strNewick += ")";
}

void PhylogenyTreeBasic ::GetLeaveIds(set<int> &lvids) {
  lvids.clear();

  PhylogenyTreeIterator itorTree(*this);
  itorTree.Init();
  while (itorTree.IsDone() == false) {
    TreeNode *pn = itorTree.GetCurrNode();
    itorTree.Next();
    if (pn == NULL) {
      break; // done with all nodes
    }
    if (pn->IsLeaf() == true) {
      lvids.insert(pn->GetID());
    }
  }
}
void PhylogenyTreeBasic ::GetLeafIntLabels(set<int> &setIntLabels) {
  vector<TreeNode *> listLeafNodes;
  GetAllLeafNodes(listLeafNodes);
  setIntLabels.clear();
  for (int i = 0; i < (int)listLeafNodes.size(); ++i) {
    setIntLabels.insert(listLeafNodes[i]->GetIntLabel());
  }
}

void PhylogenyTreeBasic::GetLeavesIdsWithLabel(const string &label,
                                               set<int> &lvids) {
  lvids.clear();
  PhylogenyTreeIterator itorTree(*this);
  itorTree.Init();
  while (itorTree.IsDone() == false) {
    TreeNode *pn = itorTree.GetCurrNode();
    // cout << "GetLeavesIdsWithLabel: ";
    // cout << pn->GetLabel() << endl;
    itorTree.Next();
    if (pn == NULL) {
      break; // done with all nodes
    }
    if (pn->GetLabel() == label) {
      lvids.insert(pn->GetID());
    }
  }
}

void PhylogenyTreeBasic ::GetLeavesWithLabels(const set<string> &setLabels,
                                              set<TreeNode *> &setLvNodes) {
  //
  setLvNodes.clear();
  PhylogenyTreeIterator itorTree(*this);
  itorTree.Init();
  while (itorTree.IsDone() == false) {
    TreeNode *pn = itorTree.GetCurrNode();
    // cout << "GetLeavesIdsWithLabel: ";
    // cout << pn->GetLabel() << endl;
    itorTree.Next();
    if (pn == NULL) {
      break; // done with all nodes
    }
    if (setLabels.find(pn->GetLabel()) != setLabels.end()) {
      setLvNodes.insert(pn);
    }
  }
}

void PhylogenyTreeBasic ::UpdateIntLabel(const vector<int> &listLabels) {
  // by assumption, id is from 0 to the following
  PhylogenyTreeIterator itorTree(*this);
  itorTree.Init();
  while (itorTree.IsDone() == false) {
    TreeNode *pn = itorTree.GetCurrNode();
    itorTree.Next();
    if (pn == NULL) {
      break; // done with all nodes
    }
    // cout << "node id = " << pn->GetID() << endl;

    YW_ASSERT_INFO(pn->GetID() < (int)listLabels.size(), "Tree id: over limit");
    int lblInt = listLabels[pn->GetID()];
    char strbuf[100];
    sprintf(strbuf, "%d", lblInt);
    string lblNew = strbuf;
    pn->SetLabel(lblNew);
  }
}

void PhylogenyTreeBasic ::Reroot(TreeNode *pRootDesc) {
  YW_ASSERT_INFO(pRootDesc != NULL, "Can not take NULL pointer");
  // if the node is set ot be root, nothing to be done
  if (pRootDesc == rootNode) {
    return;
  }
  // cout << "pass1\n";
  // create a new node
  // vector<int> dummyLbls;
  TreeNode *pRootNew = new TreeNode(rootNode->GetID());
  TreeNode *pRootOtherDesc = pRootDesc->GetParent();
  YW_ASSERT_INFO(pRootOtherDesc != NULL, "TBD");
  vector<int> lblsNew;
  // for now, concerntrate the labels without SPLITTING
  pRootOtherDesc->GetEdgeLabelsToChild(pRootDesc, lblsNew);
  pRootOtherDesc->RemoveChild(pRootDesc);
  pRootNew->AddChild(pRootDesc, lblsNew);
  // cout << "pass2\n";
  //
  TreeNode *pCurNode = pRootOtherDesc;
  TreeNode *pCurNodePar = pRootNew;
  while (true) {
    // setup the ancestral relationship
    YW_ASSERT_INFO(pCurNode != NULL && pCurNodePar != NULL, "Something wrong");
    // cout << "BEFORE CHANGING...\n";
    // cout << "pCurNode: label =" << pCurNode->GetLabel() << ", ID = " <<
    // pCurNode->GetID() << ", num of children " << pCurNode->GetChildrenNum()
    // << endl; for( int pp=0; pp< pCurNode->GetChildrenNum(); ++pp )
    //{
    // cout << "** Child: " << pCurNode->GetChild(pp)->GetID() << endl;
    //}
    // cout << "pCurNodePar: label =" << pCurNodePar->GetLabel() << ", ID = " <<
    // pCurNodePar->GetID()  << ", num of children " <<
    // pCurNodePar->GetChildrenNum()  << endl; for( int pp=0; pp<
    // pCurNodePar->GetChildrenNum(); ++pp )
    //{
    // cout << "** Child: " << pCurNodePar->GetChild(pp)->GetID() << endl;
    //}
    vector<int> lblsNew;
    pCurNode->GetEdgeLabelsToChild(pCurNodePar, lblsNew);
    TreeNode *pNodeNext = pCurNode->GetParent();
    pCurNode->RemoveChild(pCurNodePar);
    // pCurNode->SetParent(pCurNodePar);
    pCurNodePar->AddChild(pCurNode, lblsNew);

#if 0
		vector<TreeNode *> listParChildren;
		for(int c=0; c<(int)pCurNode->GetChildrenNum(); ++c  )
		{
			//if( pCurNode->GetChild(c) != pCurNode )
			//{
			listParChildren.push_back( pCurNode->GetChild(c) ) ;
			//}
		}
		for(int c=0; c<(int)listParChildren.size(); ++c  )
		{
			//if( pCurNode->GetChild(c) != pCurNode )
			//{
			pCurNode->RemoveChild( listParChildren[c] ) ;
			//}
		}
		// add these to the descendent of the new par
		for( int c=0; c<(int)listParChildren.size(); ++c )
		{
			vector<int> emptyLbls;
			pCurNodePar->AddChild(listParChildren[c], emptyLbls);
		}
#endif

    // cout << "AFTER CHANGING...\n";
    // cout << "pCurNode: label =" << pCurNode->GetLabel() << ", ID = " <<
    // pCurNode->GetID() << ", num of children " << pCurNode->GetChildrenNum()
    // << endl; for( int pp=0; pp< pCurNode->GetChildrenNum(); ++pp )
    //{
    // cout << "** Child: " << pCurNode->GetChild(pp)->GetID() << endl;
    //}
    // cout << "pCurNodePar: label =" << pCurNodePar->GetLabel() << ", ID = " <<
    // pCurNodePar->GetID()  << ", num of children " <<
    // pCurNodePar->GetChildrenNum()  << endl; for( int pp=0; pp<
    // pCurNodePar->GetChildrenNum(); ++pp )
    //{
    // cout << "** Child: " << pCurNodePar->GetChild(pp)->GetID() << endl;
    //}

    // find the other descendents of the par
    if (pNodeNext == NULL) {
      vector<TreeNode *> listParChildren;
      for (int c = 0; c < (int)pCurNode->GetChildrenNum(); ++c) {
        // if( pCurNode->GetChild(c) != pCurNode )
        //{
        listParChildren.push_back(pCurNode->GetChild(c));
        //}
      }
      for (int c = 0; c < (int)listParChildren.size(); ++c) {
        // if( pCurNode->GetChild(c) != pCurNode )
        //{
        pCurNode->RemoveChild(listParChildren[c]);
        //}
      }
      // add these to the descendent of the new par
      for (int c = 0; c < (int)listParChildren.size(); ++c) {
        vector<int> lblsNew;
        pCurNode->GetEdgeLabelsToChild(listParChildren[c], lblsNew);

        // vector<int> emptyLbls;
        pCurNodePar->AddChild(listParChildren[c], lblsNew);
      }
      pCurNodePar->RemoveChild(pCurNode);

      // cout << "FINALLY...\n";
      // cout << "pCurNode: label =" << pCurNode->GetLabel() << ", ID = " <<
      // pCurNode->GetID() << ", num of children " << pCurNode->GetChildrenNum()
      // << endl; for( int pp=0; pp< pCurNode->GetChildrenNum(); ++pp )
      //{
      // cout << "** Child: " << pCurNode->GetChild(pp)->GetID() << endl;
      //}
      // cout << "pCurNodePar: label =" << pCurNodePar->GetLabel() << ", ID = "
      // << pCurNodePar->GetID()  << ", num of children " <<
      // pCurNodePar->GetChildrenNum()  << endl; for( int pp=0; pp<
      // pCurNodePar->GetChildrenNum(); ++pp )
      //{
      // cout << "** Child: " << pCurNodePar->GetChild(pp)->GetID() << endl;
      //}
      // done. pCurNode is the root, we should by-pass this node and assign
      // their children to pCurNodePar
      break;
    }
    //
    pCurNodePar = pCurNode;
    pCurNode = pNodeNext;
  }

  // finally get rid of the original root
  delete rootNode;
  rootNode = pRootNew;
}

int PhylogenyTreeBasic ::GetNumLeaves() {
  if (numLeaves > 0) {
    return numLeaves;
  }
  set<int> lvids;
  GetLeaveIds(lvids);
  numLeaves = lvids.size();
  return numLeaves;
}

int PhylogenyTreeBasic ::GetNumInternalNodes() {
  //
  vector<TreeNode *> listAllNodes;
  GetAllNodes(listAllNodes);
  int res = 0;
  for (int i = 0; i < (int)listAllNodes.size(); ++i) {
    if (listAllNodes[i]->IsLeaf() == false) {
      //
      ++res;
    }
  }
  return res;
}

void PhylogenyTreeBasic ::GetAllLeafNodes(
    vector<TreeNode *> &listLeafNodes) const {
  listLeafNodes.clear();

  PhylogenyTreeBasic &refSelf = const_cast<PhylogenyTreeBasic &>(*this);
  PhylogenyTreeIterator itorTree(refSelf);
  itorTree.Init();
  while (itorTree.IsDone() == false) {
    TreeNode *pn = itorTree.GetCurrNode();
    itorTree.Next();
    if (pn == NULL) {
      break; // done with all nodes
    }
    if (pn->IsLeaf() == true) {
      listLeafNodes.push_back(pn);
    }
  }
}

void PhylogenyTreeBasic ::GetAllNodes(vector<TreeNode *> &listLeafNodes) const {
  listLeafNodes.clear();

  PhylogenyTreeBasic &refSelf = const_cast<PhylogenyTreeBasic &>(*this);
  PhylogenyTreeIterator itorTree(refSelf);
  itorTree.Init();
  while (itorTree.IsDone() == false) {
    TreeNode *pn = itorTree.GetCurrNode();
    itorTree.Next();
    if (pn == NULL) {
      break; // done with all nodes
    }
    listLeafNodes.push_back(pn);
  }
}

// remove all leaf nodes without taxa ids
void PhylogenyTreeBasic ::CleanNonLabeledLeaves() {
  // cout << "CleanNonLabeledLeaves:\n";
  // mark all nodes that are on the path from a labeled leaf node to root
  set<TreeNode *> setNodesNonredundent;

  vector<TreeNode *> listLeafNodes;
  GetAllLeafNodes(listLeafNodes);
  for (int ii = 0; ii < (int)listLeafNodes.size(); ++ii) {
    // cout << "Leaflabel: " << listLeafNodes[ii]->GetLabel() << endl;
    if (listLeafNodes[ii]->GetLabel().empty() == true ||
        listLeafNodes[ii]->GetLabel() == "-") {
      //
      // cout << "This leaf is REDUNDENT\n";
      continue;
    }

    TreeNode *pncurr = listLeafNodes[ii];
    while (pncurr != NULL &&
           setNodesNonredundent.find(pncurr) == setNodesNonredundent.end()) {

      //
      setNodesNonredundent.insert(pncurr);

      //
      pncurr = pncurr->GetParent();
    }
  }

  // now clean it by removing each node that does not appear in that
  PhylogenyTreeIterator itorTree(*this);
  itorTree.Init();
  vector<TreeNode *> listNodesToClean;
  while (itorTree.IsDone() == false) {
    TreeNode *pn = itorTree.GetCurrNode();
    itorTree.Next();
    if (pn == NULL) {
      break; // done with all nodes
    }
    // cout << "node id = " << pn->GetID() << endl;

    //
    if (setNodesNonredundent.find(pn) == setNodesNonredundent.end()) {
      // remove it
      listNodesToClean.push_back(pn);
    }
  }
  // now clean
  for (int ii = 0; ii < (int)listNodesToClean.size(); ++ii) {
    // cout << "Remove one node\n";
    RemoveNode(listNodesToClean[ii]);
  }
}

void PhylogenyTreeBasic ::RemoveNode(TreeNode *pn) {
  // remove the node (but does not do anything to its descendent if it has; that
  // is, we assume the node has no children)
  YW_ASSERT_INFO(pn->IsLeaf() == true, "Wrong: it still have children");
  TreeNode *pnpar = pn->GetParent();
  if (pnpar != NULL) {
    pnpar->RemoveChild(pn);
  }
  delete pn;
}

void PhylogenyTreeBasic ::RemoveNodeKeepChildren(TreeNode *pn) {
  YW_ASSERT_INFO(pn != NULL, "null");
  // cout << "RemoveNodeKeepChildren: pn: ";
  // pn->Dump();

  // remove node (and move all its children to be the nodes of the grand par
  // YW: cannot remove the root this way
  YW_ASSERT_INFO(pn != GetRoot(), "Cannot remove root this way");
  TreeNode *pnpar = pn->GetParent();
  YW_ASSERT_INFO(pnpar != NULL, "Wrong3");
  pnpar->RemoveChild(pn);

  for (int i = 0; i < pn->GetChildrenNum(); ++i) {
    vector<int> emptyLbls;
    pnpar->AddChild(pn->GetChild(i), emptyLbls);
  }
  pn->DetachAllChildren();
  delete pn;

  // remove newly created degree one node
  RemoveDegreeOneNodeAt(pnpar);
}
void PhylogenyTreeBasic ::RemoveDegreeOneNodeAt(TreeNode *pn) {
  // return;
  // cout << "removing degree one node: ";
  // pn->Dump();
  // cout << "Current tree: ";
  // this->Dump();
  // exit(1);
  // remove this node if it is a degree-1 node
  int numChildren = pn->GetChildrenNum();
  YW_ASSERT_INFO(numChildren >= 1, "Num of children: at least 1");
  if (numChildren == 1) {
    // if root, then delete it and re-set the root
    if (pn == GetRoot()) {
      // cout << "The degree one node is root!\n";
      TreeNode *pnchild = pn->GetChild(0);
      YW_ASSERT_INFO(pnchild != NULL, "pnchild: null");
      // cout << "pnchild: ";
      // pnchild->Dump();
      pnchild->DetachSelf();
      // cout << "After detach: root: ";
      // pn->Dump();
      // pn->DetachAllChildren();
      pnchild->SetParent(NULL);
      delete pn;
      SetRootPlain(pnchild);
    } else {
      // then invoke the removekeepchild
      RemoveNodeKeepChildren(pn);
    }
  }
  // cout << "Done: RemoveDegreeOneNodeAt. Tree is now: ";
  // this->Dump();
}

void PhylogenyTreeBasic ::RemoveDegreeOneNodes() {
  //
  vector<TreeNode *> listNodesAll;
  this->GetAllNodes(listNodesAll);
  for (int i = 0; i < (int)listNodesAll.size(); ++i) {
    if (listNodesAll[i]->IsLeaf() == false) {
      RemoveDegreeOneNodeAt(listNodesAll[i]);
    }
  }
}

void PhylogenyTreeBasic ::RemoveDescendentsFrom(set<TreeNode *> &setTreeNodes) {
  // only keep those whose ancestor is ot in the set given
  set<TreeNode *> setTreeNodeNew;
  for (set<TreeNode *>::iterator it = setTreeNodes.begin();
       it != setTreeNodes.end(); ++it) {
    // check whether any of its parent is in the list
    bool fKeep = true;
    TreeNode *ppar = (*it)->GetParent();
    while (ppar != NULL) {
      if (setTreeNodes.find(ppar) != setTreeNodes.end()) {
        fKeep = false;
        break;
      }
      ppar = ppar->GetParent();
    }
    if (fKeep == true) {
      setTreeNodeNew.insert(*it);
    }
  }
  setTreeNodes = setTreeNodeNew;
}

// given a set of clusters (subsets of tree taxa), construct the corresponding
// phylo trees YW: need to allow mulfurcating trees
void PhylogenyTreeBasic ::ConsPhyTreeFromClusters(
    const set<set<int> > &setClusters) {
  // cout << "ConsPhyTreeFromClusters :: Cluseters: \n";
  // for( set< set<int> > :: const_iterator it = setClusters.begin(); it !=
  // setClusters.end(); ++it )
  //{
  // DumpIntSet( *it );
  //}
  // assume all leaves are given as singleton taxon. So first collect those
  // singleton subsets
  set<set<int> > setSubsetsActive;
  TreeNode *nodeLast = NULL;
  map<set<int>, TreeNode *> mapClusterToNode;
  for (set<set<int> >::const_iterator it = setClusters.begin();
       it != setClusters.end(); ++it) {
    if (it->size() == 1) {
      // add in setClusters
      setSubsetsActive.insert(*it);
      // also create nodes
      TreeNode *pnode = new TreeNode(*(it->begin()));
      char buf[100];
      sprintf(buf, "%d", *(it->begin()));
      string sbuf = buf;
      pnode->SetLabel(sbuf);
      nodeLast = pnode;
      mapClusterToNode.insert(
          map<set<int>, TreeNode *>::value_type(*it, pnode));
    }
  }
  // setup num of leaves now
  this->numLeaves = mapClusterToNode.size();

  // need to allow mulfurcating trees
  // approach: for each cluster, maintain a pointer that points to the cluster
  // that is its parent then, each time, loop through to find all parents
  map<set<int>, set<int> > mapClustrToPar;
  // try to see whether we can create new nodes
  for (set<set<int> >::iterator it1 = setClusters.begin();
       it1 != setClusters.end(); ++it1) {
    set<set<int> >::iterator it2 = setClusters.begin();
    ++it2;
    for (; it2 != setClusters.end(); ++it2) {
      //
      set<int> sLarger = *it1;
      set<int> sSmaller = *it2;
      if (sLarger.size() < sSmaller.size()) {
        sLarger = *it2;
        sSmaller = *it1;
      }
      // can these two coalesce into a single cluster known
      if (sLarger.size() > sSmaller.size() &&
          IsSetContainer(sLarger, sSmaller) == true) {
        if (mapClustrToPar.find(sSmaller) == mapClustrToPar.end() ||
            mapClustrToPar[sSmaller].size() > sLarger.size()) {
          mapClustrToPar.erase(sSmaller);
          mapClustrToPar.insert(
              map<set<int>, set<int> >::value_type(sSmaller, sLarger));
        }
      }
    }
  }

  // loop until there is only a single subset
  while (setSubsetsActive.size() > 1) {
    set<set<int> > setSubsetsActiveNext = setSubsetsActive;
    // cout << "Current active sets: \n";
    // for( set< set<int> > :: const_iterator it = setSubsetsActiveNext.begin();
    // it != setSubsetsActiveNext.end(); ++it )
    //{
    // DumpIntSet( *it );
    //}
    // try to find several clusters that have the same parent cluster
    // try to see whether we can create new nodes
    map<set<int>, set<set<int> > > mapClusterCoal;
    for (set<set<int> >::iterator it1 = setSubsetsActive.begin();
         it1 != setSubsetsActive.end(); ++it1) {
      // get parent
      YW_ASSERT_INFO(mapClustrToPar.find(*it1) != mapClustrToPar.end(),
                     "Cluster: not found");
      if (mapClusterCoal.find(mapClustrToPar[*it1]) == mapClusterCoal.end()) {
        set<set<int> > sempty;
        mapClusterCoal.insert(map<set<int>, set<set<int> > >::value_type(
            mapClustrToPar[*it1], sempty));
      }
      // cout << "Having child cluster: ";
      // DumpIntSet( mapClustrToPar[*it1] );
      // cout << ", for child ";
      // DumpIntSet(*it1);
      mapClusterCoal[mapClustrToPar[*it1]].insert(*it1);
    }

    // now process each record
    for (map<set<int>, set<set<int> > >::iterator it2 = mapClusterCoal.begin();
         it2 != mapClusterCoal.end(); ++it2) {
      // YW_ASSERT_INFO( it2->second.size() > 1, "Must have at least two
      // coalescing" );
      // cout << "Set parent: ";
      // DumpIntSet(it2->first);
      set<int> sunion;
      for (set<set<int> >::iterator it3 = it2->second.begin();
           it3 != it2->second.end(); ++it3) {
        // cout << "Set child: ";
        // DumpIntSet(*it3);
        // can these two coalesce into a single cluster known
        UnionSets(sunion, *it3);
      }
      // cout << "sunion = ";
      // DumpIntSet( sunion );
      // ensure these do coal into some meaningful cluster
      if (setClusters.find(sunion) == setClusters.end()) {
        // cout << "This set not complete\n";
        // this cluster not done yet
        continue;
      }

      // create this new node
      TreeNode *pnode = new TreeNode;
      nodeLast = pnode;
      for (set<set<int> >::iterator it3 = it2->second.begin();
           it3 != it2->second.end(); ++it3) {
        // cout << "Processing first subset: ";
        // DumpIntSet( *it1 );
        // cout << "Processing second subset: ";
        // DumpIntSet( *it2 );
        // these two add up to an input cluster and so create a new node for it
        YW_ASSERT_INFO(mapClusterToNode.find(*it3) != mapClusterToNode.end(),
                       "Fail1");
        vector<int> emptyLabels;
        pnode->AddChild(mapClusterToNode[*it3], emptyLabels);
        setSubsetsActiveNext.erase(*it3);
      }
      mapClusterToNode.insert(
          map<set<int>, TreeNode *>::value_type(sunion, pnode));
      setSubsetsActiveNext.insert(sunion);
      // cout << "Creating node: " << endl;
    }
    // must make progress
    YW_ASSERT_INFO(setSubsetsActive != setSubsetsActiveNext,
                   "Did not make progress");
    setSubsetsActive = setSubsetsActiveNext;
  }
  YW_ASSERT_INFO(nodeLast != NULL, "nodeLast: NULL");
  SetRoot(nodeLast);
}

// find the set of clades in the subtree specified by the given leaf nodes
void PhylogenyTreeBasic ::FindCladeOfSubsetLeaves(
    const set<TreeNode *> &setLeaves, set<set<TreeNode *> > &setSubtreeClades) {
  // caution: do not check whether these are true leaves
  TreeNode *pRoot = this->GetRoot();
  set<TreeNode *> setAllNodes;
  pRoot->GetAllDescendents(setAllNodes);

  //
  for (set<TreeNode *>::iterator it = setAllNodes.begin();
       it != setAllNodes.end(); ++it) {
    //
    set<TreeNode *> setLeavesUnder;
    (*it)->GetAllLeavesUnder(setLeavesUnder);
    set<TreeNode *> setLeavesSS;
    JoinSetsGen(setLeavesUnder, setLeaves, setLeavesSS);
    if (setLeavesSS.size() > 0) {
      setSubtreeClades.insert(setLeavesSS);
    }
  }
}

// find the set of clades in the subtree specified by the given leaf nodes
void PhylogenyTreeBasic ::FindCladeOfSubsetLeavesExact(
    const set<TreeNode *> &setLeaves, set<set<TreeNode *> > &setSubtreeClades) {
  // caution: do not check whether these are true leaves
  TreeNode *pRoot = this->GetRoot();
  set<TreeNode *> setAllNodes;
  pRoot->GetAllDescendents(setAllNodes);

  //
  for (set<TreeNode *>::iterator it = setAllNodes.begin();
       it != setAllNodes.end(); ++it) {
    //
    set<TreeNode *> setLeavesUnder;
    (*it)->GetAllLeavesUnder(setLeavesUnder);
    set<TreeNode *> setLeavesSS;
    JoinSetsGen(setLeavesUnder, setLeaves, setLeavesSS);
    if (setLeavesSS == setLeavesUnder) {
      setSubtreeClades.insert(setLeavesSS);
    }
  }
}

void PhylogenyTreeBasic ::GroupLeavesToSubtrees(
    const set<TreeNode *> &setLeaves,
    const set<set<TreeNode *> > &cladeNodesToProc,
    set<set<TreeNode *> > &setSubtreeClades) {
  // group the leaves into subtrees (i.e. the subtrees contains exactly those
  // appear in the leaves YW: note this is not the most realistic way (say you
  // have one noisy leaf sepearting two otherwise fully connected catepillar
  // tree, then the result willl be a lot more trees to use). But this servers
  // as a starting point YW: here, we are given some subset out of some
  // pre-specified leaf set, and some subsets (clades) over these leaves; we
  // want to find the set of maximal clades containing partition these leaves
  // TreeNode *pRoot = this->GetRoot();
  // set<TreeNode *> setAllNodes;
  // pRoot->GetAllDescendents(setAllNodes);

  // order based on the size
  map<int, set<set<TreeNode *> > > mapSubtreeSz;
  // for( set<TreeNode *> :: iterator it = setAllNodes.begin(); it !=
  // setAllNodes.end(); ++it)
  for (set<set<TreeNode *> >::const_iterator it = cladeNodesToProc.begin();
       it != cladeNodesToProc.end(); ++it) {
    //
    // set<TreeNode *> setLeavesUnder;
    //(*it)->GetAllLeavesUnder( setLeavesUnder );
    if (mapSubtreeSz.find(it->size()) == mapSubtreeSz.end()) {
      set<set<TreeNode *> > ss;
      mapSubtreeSz.insert(
          map<int, set<set<TreeNode *> > >::value_type(it->size(), ss));
    }
    mapSubtreeSz[it->size()].insert(*it);
  }

  // reverse order
  set<TreeNode *> setNodesProc = setLeaves;
  for (map<int, set<set<TreeNode *> > >::reverse_iterator rit =
           mapSubtreeSz.rbegin();
       rit != mapSubtreeSz.rend(); ++rit) {
    //
    for (set<set<TreeNode *> >::iterator itg = rit->second.begin();
         itg != rit->second.end(); ++itg) {
      //
      set<TreeNode *> setLeavesSS;
      JoinSetsGen(*itg, setNodesProc, setLeavesSS);
      if (setLeavesSS.size() == itg->size()) {
        // find a good match here, use it
        setSubtreeClades.insert(*itg);
        SubtractSetsGen(setNodesProc, *itg);
      }
    }
    if (setNodesProc.size() == 0) {
      break;
    }
  }
  YW_ASSERT_INFO(setNodesProc.size() == 0, "Fail to classify all subtrees");
}

void PhylogenyTreeBasic ::GroupLeavesToSubtreesSamePar(
    const set<TreeNode *> &setLeaves,
    const set<set<TreeNode *> > &cladeNodesToProc,
    set<set<TreeNode *> > &setSubtreeClades) {
  // group leaves that form subtrees w/ same parents. Difference from above: for
  // two subtrees that share the same parent but could be other branches, put
  // the together
  GroupLeavesToSubtrees(setLeaves, cladeNodesToProc, setSubtreeClades);
  // now see whether we can combine subtrees s.t. the combined one is still
  // contined in some parent
  map<set<TreeNode *>, set<TreeNode *> > mapSubtreesToPar;
  for (set<set<TreeNode *> >::iterator it = setSubtreeClades.begin();
       it != setSubtreeClades.end(); ++it) {
    for (set<set<TreeNode *> >::iterator itg = cladeNodesToProc.begin();
         itg != cladeNodesToProc.end(); ++itg) {
      //
      if (*itg != *it && itg->size() > it->size() &&
          (mapSubtreesToPar.find(*it) == mapSubtreesToPar.end() ||
           mapSubtreesToPar[*it].size() > itg->size())) {
        //
        set<TreeNode *> sint;
        JoinSetsGen(*itg, *it, sint);
        if (sint.size() == it->size()) {
          //
          if (mapSubtreesToPar.find(*it) == mapSubtreesToPar.end()) {
            mapSubtreesToPar.insert(
                map<set<TreeNode *>, set<TreeNode *> >::value_type(*it, *itg));
          } else {
            mapSubtreesToPar[*it] = *itg;
          }
        }
      }
    }
  }
  map<set<TreeNode *>, set<TreeNode *> > mapRevParToSubtrees;
  for (map<set<TreeNode *>, set<TreeNode *> >::iterator it =
           mapSubtreesToPar.begin();
       it != mapSubtreesToPar.end(); ++it) {
    //
    if (mapRevParToSubtrees.find(it->second) == mapRevParToSubtrees.end()) {
      mapRevParToSubtrees.insert(
          map<set<TreeNode *>, set<TreeNode *> >::value_type(it->second,
                                                             it->first));
    } else {
      UnionSetsGen(mapRevParToSubtrees[it->second], it->first);
    }
  }
  setSubtreeClades.clear();
  for (map<set<TreeNode *>, set<TreeNode *> >::iterator it =
           mapRevParToSubtrees.begin();
       it != mapRevParToSubtrees.end(); ++it) {
    setSubtreeClades.insert(it->second);
  }
}

void PhylogenyTreeBasic ::GetAllClades(set<set<int> > &setClades) {
  //
  setClades.clear();
  // now clean it by removing each node that does not appear in that
  PhylogenyTreeIterator itorTree(*this);
  itorTree.Init();
  while (itorTree.IsDone() == false) {
    TreeNode *pn = itorTree.GetCurrNode();
    itorTree.Next();
    if (pn == NULL) {
      break; // done with all nodes
    }
    // cout << "node id = " << pn->GetID() << endl;
    set<TreeNode *> setDescendents;
    pn->GetAllLeavesUnder(setDescendents);
    set<int> sint;
    for (set<TreeNode *>::iterator itg = setDescendents.begin();
         itg != setDescendents.end(); ++itg) {
      sint.insert((*itg)->GetIntLabel());
    }
    setClades.insert(sint);
  }
}

void PhylogenyTreeBasic ::GetAllCladesList(vector<set<int> > &listClades) {
  listClades.clear();
  // now clean it by removing each node that does not appear in that
  PhylogenyTreeIterator itorTree(*this);
  itorTree.Init();
  while (itorTree.IsDone() == false) {
    TreeNode *pn = itorTree.GetCurrNode();
    itorTree.Next();
    if (pn == NULL) {
      break; // done with all nodes
    }
    // cout << "node id = " << pn->GetID() << endl;
    set<TreeNode *> setDescendents;
    pn->GetAllLeavesUnder(setDescendents);
    set<int> sint;
    for (set<TreeNode *>::iterator itg = setDescendents.begin();
         itg != setDescendents.end(); ++itg) {
      sint.insert((*itg)->GetIntLabel());
    }
    listClades.push_back(sint);
  }
}

// different from the above, (1) we allow duplicate int-labels (and thus
// multiset) (2) group clades by common parents
void PhylogenyTreeBasic ::GetAllCladeGroupsIntLabel(
    multiset<multiset<multiset<int> > > &setCladeGroupsDupLabels,
    multiset<int> &rootClade) {
  // group all clades by parent nodes (i.e. clades with same parent are in one
  // class) root clade: the one with all leaves
  map<TreeNode *, multiset<multiset<int> > > mapCladeGroupsForNode;

  //
  setCladeGroupsDupLabels.clear();
  rootClade.clear();
  // now clean it by removing each node that does not appear in that
  PhylogenyTreeIterator itorTree(*this);
  itorTree.Init();
  while (itorTree.IsDone() == false) {
    TreeNode *pn = itorTree.GetCurrNode();
    itorTree.Next();
    if (pn == NULL) {
      break; // done with all nodes
    }
    // cout << "node id = " << pn->GetID() << endl;
    set<TreeNode *> setDescendents;
    pn->GetAllLeavesUnder(setDescendents);
    multiset<int> sint;
    for (set<TreeNode *>::iterator itg = setDescendents.begin();
         itg != setDescendents.end(); ++itg) {
      sint.insert((*itg)->GetIntLabel());
    }
    TreeNode *pnPar = pn->GetParent();
    if (pnPar == NULL) {
      // this is the root clade
      rootClade = sint;
    } else {
      if (mapCladeGroupsForNode.find(pnPar) == mapCladeGroupsForNode.end()) {
        multiset<multiset<int> > mms;
        mapCladeGroupsForNode.insert(
            map<TreeNode *, multiset<multiset<int> > >::value_type(pnPar, mms));
      }
      mapCladeGroupsForNode[pnPar].insert(sint);
    }
  }
  YW_ASSERT_INFO(rootClade.size() > 0, "Fail to collect root clade");
  for (map<TreeNode *, multiset<multiset<int> > >::iterator it =
           mapCladeGroupsForNode.begin();
       it != mapCladeGroupsForNode.end(); ++it) {
    //
    setCladeGroupsDupLabels.insert(it->second);
  }
}

void PhylogenyTreeBasic ::GetAllCladesById(set<set<int> > &setClades) {
  //
  setClades.clear();
  // now clean it by removing each node that does not appear in that
  PhylogenyTreeIterator itorTree(*this);
  itorTree.Init();
  while (itorTree.IsDone() == false) {
    TreeNode *pn = itorTree.GetCurrNode();
    itorTree.Next();
    if (pn == NULL) {
      break; // done with all nodes
    }
    // cout << "node id = " << pn->GetID() << endl;
    set<TreeNode *> setDescendents;
    pn->GetAllLeavesUnder(setDescendents);
    set<int> sint;
    for (set<TreeNode *>::iterator itg = setDescendents.begin();
         itg != setDescendents.end(); ++itg) {
      sint.insert((*itg)->GetID());
    }
    setClades.insert(sint);
  }
}

void PhylogenyTreeBasic ::GetAllCladeNodess(set<set<TreeNode *> > &setClades) {
  //
  setClades.clear();
  // now clean it by removing each node that does not appear in that
  PhylogenyTreeIterator itorTree(*this);
  itorTree.Init();
  while (itorTree.IsDone() == false) {
    TreeNode *pn = itorTree.GetCurrNode();
    itorTree.Next();
    if (pn == NULL) {
      break; // done with all nodes
    }
    // cout << "node id = " << pn->GetID() << endl;
    set<TreeNode *> setDescendents;
    pn->GetAllLeavesUnder(setDescendents);

    setClades.insert(setDescendents);
  }
}

TreeNode *PhylogenyTreeBasic ::GetSubtreeRootForLeaves(
    const set<TreeNode *> &setLvNodes) {
  PhylogenyTreeIterator itorTree(*this);
  itorTree.Init();
  while (itorTree.IsDone() == false) {
    TreeNode *pn = itorTree.GetCurrNode();
    itorTree.Next();
    if (pn == NULL) {
      break; // done with all nodes
    }
    // cout << "node id = " << pn->GetID() << endl;
    set<TreeNode *> setDescendents;
    pn->GetAllLeavesUnder(setDescendents);

    if (setLvNodes == setDescendents) {
      return pn;
    }
  }
  return NULL;
}

void PhylogenyTreeBasic ::GroupNodesWithCommonPars(
    const set<TreeNode *> &setNodes,
    map<TreeNode *, set<TreeNode *> > &mapNodesWithSamePar) {
  //
  mapNodesWithSamePar.clear();
  for (set<TreeNode *>::const_iterator it = setNodes.begin();
       it != setNodes.end(); ++it) {
    //
    TreeNode *ppar = (*it)->GetParent();
    if (mapNodesWithSamePar.find(ppar) == mapNodesWithSamePar.end()) {
      set<TreeNode *> ss;
      mapNodesWithSamePar.insert(
          map<TreeNode *, set<TreeNode *> >::value_type(ppar, ss));
    }
    mapNodesWithSamePar[ppar].insert(*it);
  }
}

void PhylogenyTreeBasic ::RemoveEdgeLabels() {
  //
  this->rootNode->RemoveLabels();
}

void PhylogenyTreeBasic ::RemoveEdgeLabelsToLeaves() {
  // get all leaves
  vector<TreeNode *> vecLeaves;
  GetAllLeafNodes(vecLeaves);
  for (int i = 0; i < (int)vecLeaves.size(); ++i) {
    vecLeaves[i]->RemoveLabelsPar();
  }
}

void PhylogenyTreeBasic ::IncEdgeLabelsBy(int offset) {
  // inc edge label of this node (and subtree if needed)
  this->rootNode->IncEdgeLabelsBy(offset, true);
}

string PhylogenyTreeBasic ::GetShapeLabelNodeBrNum(
    map<TreeNode *, pair<int, int> > &mapNodeNumBrannches,
    vector<int> &listOrderedLeaves) {
  // format: <num of underlying branches, event id>, negative for internal nodes
  map<TreeNode *, pair<int, int> > mapNodeNumBrannchesUse = mapNodeNumBrannches;
  // given: num of branches at each node,
  // return shape label as empty Newick format
  // for this, first need to find out all nodes that all descendents have
  // appeared in the tree
  set<TreeNode *> setAncesNotGiven;
  for (map<TreeNode *, pair<int, int> >::iterator it =
           mapNodeNumBrannches.begin();
       it != mapNodeNumBrannches.end(); ++it) {
    set<TreeNode *> setAllAnces;
    it->first->GetAllAncestors(setAllAnces);
    for (set<TreeNode *>::iterator itg = setAllAnces.begin();
         itg != setAllAnces.end(); ++itg) {
      if (mapNodeNumBrannches.find(*itg) == mapNodeNumBrannches.end()) {
        //
        pair<int, int> pp(-1, -1);
        mapNodeNumBrannchesUse.insert(
            map<TreeNode *, pair<int, int> >::value_type(*itg, pp));
      }
    }
  }
  // now call the root to find the label
  return this->rootNode->GetShapeLabelNodeBrNum(mapNodeNumBrannchesUse,
                                                listOrderedLeaves);
}

void PhylogenyTreeBasic ::MakeSubtreeUnrefined(TreeNode *pSubtree) {
  // make this subtree unrefined (i.e. each leaf points to the root
  // CAUTION: all edge labels are LOST!!!!
  set<TreeNode *> setAllLeavesUnder;
  pSubtree->GetAllLeavesUnder(setAllLeavesUnder);
  // cout << "setAllLeavesUnder: ";
  // for( set<TreeNode *> :: iterator it = setAllLeavesUnder.begin(); it !=
  // setAllLeavesUnder.end(); ++it)
  //{
  //(*it)->Dump();
  //}
  // cout << endl;
  set<TreeNode *> setAllDescUnder;
  pSubtree->GetAllDescendents(setAllDescUnder);
  // cout << "setAllDescUnder: ";
  // for( set<TreeNode *> :: iterator it = setAllDescUnder.begin(); it !=
  // setAllDescUnder.end(); ++it)
  //{
  //(*it)->Dump();
  //}
  // cout << endl;

  // detach all leaves from their parent
  for (set<TreeNode *>::iterator it = setAllLeavesUnder.begin();
       it != setAllLeavesUnder.end(); ++it) {
    //
    TreeNode *ppar = (*it)->GetParent();
    ppar->RemoveChild(*it);
  }

  pSubtree->RemoveAllChildren();

  // remove all descendent except the leaves
  for (set<TreeNode *>::iterator it = setAllDescUnder.begin();
       it != setAllDescUnder.end(); ++it) {
    // need to be careful b/c node deletion is recurisvely
    if (setAllLeavesUnder.find(*it) == setAllLeavesUnder.end() &&
        (*it) != pSubtree && ((*it)->GetParent() == pSubtree)) {
      // cout << "Delete this node: ";
      //(*it)->Dump();
      delete *it;
    }
  }
  // then add the leaves directly under the subtree root
  for (set<TreeNode *>::iterator it = setAllLeavesUnder.begin();
       it != setAllLeavesUnder.end(); ++it) {
    vector<int> lblEmpty;
    pSubtree->AddChild(*it, lblEmpty);
  }
  // string strTree;
  // ConsNewick(strTree);
  // cout << "After MakeSubtreeUnrefiined: tree is " << strTree << endl;
}

void PhylogenyTreeBasic ::Binarize() {
  // make the tree binary
  int idToUseNext = this->rootNode->GetMaxIdWithinSubtree() + 1;
  this->rootNode->Binarize(idToUseNext);
  // string strTree;
  // ConsNewick(strTree);
  // cout << "After binarization: tree is " << strTree << endl;
}

void PhylogenyTreeBasic ::CreatePhyTreeFromLeavesWithLabels(
    const set<string> &setLeafLabels, PhylogenyTreeBasic &treeSubsetLeaves,
    bool fUseOldTaxonName) {
  // given a set of leaf labels, construct another phylogenetic tree that is
  // extracted from the current tree by only taking those leaves with one of the
  // given labels YW: caution: all taxa names are mapped to 0,1,2,... according
  // to their order in list if fUseOldTaxonName=false otherwise, keep the
  // original flag
  set<int> setSubsetLeaves;
  map<int, string> mapOrigIdToOrigStrLbl;
  int idToUseFirst = 0;
  for (set<string>::const_iterator it = setLeafLabels.begin();
       it != setLeafLabels.end(); ++it) {
    string lblcur = *it;
    set<int> setSubsetLeavesStep;
    GetLeavesIdsWithLabel(lblcur, setSubsetLeavesStep);
    // cout << "CreatePhyTreeFromLeavesWithLabels: lblcur: " << lblcur <<"
    // setSubsetLeavesStep: "; DumpIntSet(setSubsetLeavesStep);
    UnionSets(setSubsetLeaves, setSubsetLeavesStep);

    string lblToUse = lblcur;
    if (fUseOldTaxonName == false) {
      char buf[100];
      sprintf(buf, "%d", idToUseFirst++);
      lblToUse = buf;
    }
    for (set<int>::iterator it2 = setSubsetLeavesStep.begin();
         it2 != setSubsetLeavesStep.end(); ++it2) {
      // mapOrigIdToOrigStrLbl.insert( map<int,string> :: value_type(*it2,
      // lblcur) );
      mapOrigIdToOrigStrLbl.insert(
          map<int, string>::value_type(*it2, lblToUse));
      // cout << "mapOrigIdToOrigStrLbl: " << *it2 << ", lblToUse: " << lblToUse
      // << endl;
    }
  }

  // get all clades first
  set<set<int> > setClades;
  GetAllCladesById(setClades);
  // cout << "All clades: \n";
  // for(set<set<int> > :: iterator it = setClades.begin(); it !=
  // setClades.end(); ++it)
  //{
  // DumpIntSet(*it);
  //}

  // map the remaining id to 0,1,2....
  map<int, int> mapIdToContinue;
  map<int, string> mapContIdToOrigStr;
  int idToUse = 0;
  for (set<int>::iterator it = setSubsetLeaves.begin();
       it != setSubsetLeaves.end(); ++it) {
    YW_ASSERT_INFO(
        mapOrigIdToOrigStrLbl.find(*it) != mapOrigIdToOrigStrLbl.end(), "Fail");
    mapContIdToOrigStr.insert(
        map<int, string>::value_type(idToUse, mapOrigIdToOrigStrLbl[*it]));
    // cout << "mapContIdToOrigStr: idtouse: " << idToUse << ", string orig: "
    // << mapOrigIdToOrigStrLbl[*it] << endl;
    mapIdToContinue.insert(map<int, int>::value_type(*it, idToUse++));
  }

  set<set<int> > setCladesSub;
  // now extract those with only those given
  for (set<set<int> >::iterator it = setClades.begin(); it != setClades.end();
       ++it) {
    set<int> sintstep;
    JoinSets(*it, setSubsetLeaves, sintstep);
    if (sintstep.size() > 0) {
      // convert to continuios id first
      set<int> sintstep2;
      MapIntSetTo(sintstep, mapIdToContinue, sintstep2);

      setCladesSub.insert(sintstep2);

      // cout << "Adding a clade: ";
      // DumpIntSet( sintstep2);
      // cout << "for orig clade: ";
      // DumpIntSet(sintstep);
    }
  }

  // now build a tree with these labels
  CreatePhyTreeWithRootedSplits(treeSubsetLeaves, setSubsetLeaves.size(),
                                setCladesSub);

  // now map the leaves of the new tree to the original ids
  treeSubsetLeaves.AssignLeafLabels(mapContIdToOrigStr);

  // cout << "This is the phylogenetic tree constructed from subset of leaves: "
  // this->OutputGML("tree1.gml");
  // treeSubsetLeaves.OutputGML("t1.gml");
  // exit(1);
}

void PhylogenyTreeBasic ::AssignLeafLabels(
    const map<int, string> &mapLeafLbls) {
  // assign labels stored in the map (format: node id to lbl)
  vector<TreeNode *> listLeafNodes;
  GetAllLeafNodes(listLeafNodes);
  for (int i = 0; i < (int)listLeafNodes.size(); ++i) {
    int idn = listLeafNodes[i]->GetID();
    map<int, string>::const_iterator itg = mapLeafLbls.find(idn);
    YW_ASSERT_INFO(itg != mapLeafLbls.end(), "Fail");
    string strLblNew = itg->second;
    listLeafNodes[i]->SetLabel(strLblNew);
    listLeafNodes[i]->SetUserLabel(strLblNew);
  }
}
void PhylogenyTreeBasic ::ReassignLeafLabels(
    const map<string, string> &mapLeafLbls) {
  vector<TreeNode *> listLeafNodes;
  GetAllLeafNodes(listLeafNodes);
  for (int i = 0; i < (int)listLeafNodes.size(); ++i) {
    string str = listLeafNodes[i]->GetLabel();
    // cout << "leaf label curr: " << str << endl;
    map<string, string>::const_iterator itg = mapLeafLbls.find(str);

    if (itg == mapLeafLbls.end()) {
      // TBD. YW: for now. Need to look at later...
      continue;
    }

    YW_ASSERT_INFO(itg != mapLeafLbls.end(), "Fail");
    string strLblNew = itg->second;
    listLeafNodes[i]->SetLabel(strLblNew);
    listLeafNodes[i]->SetUserLabel(strLblNew);
  }
}

void PhylogenyTreeBasic ::SetUserLabelToCurrLabels() {
  vector<TreeNode *> listLeafNodes;
  GetAllLeafNodes(listLeafNodes);
  for (int i = 0; i < (int)listLeafNodes.size(); ++i) {
    listLeafNodes[i]->SetUserLabel(listLeafNodes[i]->GetLabel());
  }
}

void PhylogenyTreeBasic ::SetLabelsToCurrUserLabels() {
  vector<TreeNode *> listLeafNodes;
  GetAllLeafNodes(listLeafNodes);
  for (int i = 0; i < (int)listLeafNodes.size(); ++i) {
    listLeafNodes[i]->SetLabel(listLeafNodes[i]->GetUserLabel());
  }
}

int PhylogenyTreeBasic ::GetMaxDegree() const {
  int res = 0;

  PhylogenyTreeBasic &thisTree = const_cast<PhylogenyTreeBasic &>(*this);
  PhylogenyTreeIterator itor(thisTree);
  itor.Init();
  while (itor.IsDone() == false) {
    TreeNode *pn = itor.GetCurrNode();

    int degThis = pn->GetChildrenNum();
    if (degThis > res) {
      res = degThis;
    }

    itor.Next();
  }
  return res;
}

void PhylogenyTreeBasic ::Dump() const {
  // dump all nodes
  PhylogenyTreeBasic &thisTree = const_cast<PhylogenyTreeBasic &>(*this);
  PhylogenyTreeIterator itor(thisTree);
  itor.Init();
  while (itor.IsDone() == false) {
    TreeNode *pn = itor.GetCurrNode();
    pn->Dump();
    cout << endl;
    itor.Next();
  }
}

void PhylogenyTreeBasic ::GetSubtreesWithMaxSize(set<TreeNode *> &setSTRoots,
                                                 int maxSzSubtree) const {
#if 0
    // YW: this piece of code is not used
    set<TreeNode *> setSTRootsStep;
    set<string> setLblsEmpty;
    GetSubtreesWithMaxSizeExcludeTaxa( setSTRootsStep, maxSzSubtree, setLblsEmpty );

    // find out any missing taxa
    set<string> setLblsCovered;
    for(set<TreeNode *> :: iterator it = setSTRootsStep.begin(); it != setSTRootsStep.end(); ++it )
    {
        set<string> setLblsCoveredStep;
        (*it)->GetAllDistinctLeafLabeles(setLblsCoveredStep);
        UnionSetsGen( setLblsCovered, setLblsCoveredStep );

cout << "setLblsCoveredStep: ";
for(set<string> :: iterator it = setLblsCoveredStep.begin(); it != setLblsCoveredStep.end(); ++it)
{
cout << *it << "  ";
}
cout << endl;
    }
    set<string> setLblsExc;
    this->GetRoot()->GetAllDistinctLeafLabeles( setLblsExc );
    SubtractSetsGen( setLblsExc, setLblsCovered );
cout << "setLblsExc: ";
for(set<string> :: iterator it = setLblsExc.begin(); it != setLblsExc.end(); ++it)
{
cout << *it << "  ";
}
cout << endl;

    // now do it again so to allow those with these "free" taxa
    GetSubtreesWithMaxSizeExcludeTaxa( setSTRoots, maxSzSubtree, setLblsExc);
#endif

  //#if 0
  // retrieve roots of subtrees that are no biggere than the specified size
  setSTRoots.clear();
  stack<TreeNode *> stackTrNodes;
  stackTrNodes.push(this->GetRoot());

  while (stackTrNodes.size() > 0) {
    //
    TreeNode *pncurr = stackTrNodes.top();
    stackTrNodes.pop();

    // save it if this subtree size is not too big
    set<TreeNode *> setDescendents;
    pncurr->GetAllLeavesUnder(setDescendents);
    // cout << "pncur: number of descendents: " << setDescendents.size() << " ";
    // pncurr->Dump();
    if ((int)setDescendents.size() <= maxSzSubtree) {
      // cout << "Adding tis node.\n";
      setSTRoots.insert(pncurr);
    } else {
      // cout << "Process each of its descendents.\n";
      // then check all its descendents
      for (int i = 0; i < pncurr->GetChildrenNum(); ++i) {
        TreeNode *pnc = pncurr->GetChild(i);
        stackTrNodes.push(pnc);
        // cout << "pushing child: ";
        // pnc->Dump();
      }
    }
  }
  //#endif
}

void PhylogenyTreeBasic ::GetMaxSubtrees(set<TreeNode *> &setSTRootsIdents) {
  // obtain max subtrees with identical leaf labels
  setSTRootsIdents.clear();
  stack<TreeNode *> stackNodes;
  stackNodes.push(GetRoot());
  while (stackNodes.empty() == false) {
    //
    TreeNode *pncurr = stackNodes.top();
    stackNodes.pop();

    vector<string> strLblLeaves;
    pncurr->GetAllLeafLabeles(strLblLeaves);
    set<string> strLblLeavesSet;
    PopulateSetByVecGen(strLblLeavesSet, strLblLeaves);
    YW_ASSERT_INFO(strLblLeavesSet.size() >= 1, "Must have at least one label");
    if (strLblLeavesSet.size() == 1) {
      //
      setSTRootsIdents.insert(pncurr);
    } else {
      // consider all children
      for (int i = 0; i < pncurr->GetChildrenNum(); ++i) {
        stackNodes.push(pncurr->GetChild(i));
      }
    }
  }
}

bool PhylogenyTreeBasic ::GetSiblingsPairFrom(
    const set<TreeNode *> &setNodesToChoose,
    pair<TreeNode *, TreeNode *> &pairSibs) {
  // find which pairs of given nodes have the same paent
  bool fres = false;

  //
  map<TreeNode *, TreeNode *> mapParToOrigNode;
  for (set<TreeNode *>::const_iterator it = setNodesToChoose.begin();
       it != setNodesToChoose.end(); ++it) {
    TreeNode *pp = (*it)->GetParent();
    if (mapParToOrigNode.find(pp) == mapParToOrigNode.end()) {
      mapParToOrigNode.insert(map<TreeNode *, TreeNode *>::value_type(pp, *it));
    } else {
      pairSibs.first = mapParToOrigNode[pp];
      pairSibs.second = *it;
      fres = true;
      break;
    }
  }
#if 0
cout << "GetSiblingsPairFrom: \n";
if( fres == true )
{
cout << "sib1: ";
pairSibs.first->Dump();
cout << "sib2: ";
pairSibs.second->Dump();
}
#endif

  return fres;
}

bool PhylogenyTreeBasic ::GetSiblingsNodesFrom(
    const set<TreeNode *> &setNodesToChoose, set<TreeNode *> &setSibs) {
  // find which nodes from given nodes have the same paent
  // YW: we prefer the lower if there are multiple choices
  bool fres = false;

  //
  map<TreeNode *, set<TreeNode *> > mapParToOrigNode;
  for (set<TreeNode *>::const_iterator it = setNodesToChoose.begin();
       it != setNodesToChoose.end(); ++it) {
    TreeNode *pp = (*it)->GetParent();
    if (mapParToOrigNode.find(pp) == mapParToOrigNode.end()) {
      set<TreeNode *> ss;
      mapParToOrigNode.insert(
          map<TreeNode *, set<TreeNode *> >::value_type(pp, ss));
    }
    mapParToOrigNode[pp].insert(*it);
  }
  // assign one with at least two nodes
  for (map<TreeNode *, set<TreeNode *> >::iterator it =
           mapParToOrigNode.begin();
       it != mapParToOrigNode.end(); ++it) {
    if (it->second.size() > 1) {
      bool fGood = true;

      for (map<TreeNode *, set<TreeNode *> >::iterator it2 =
               mapParToOrigNode.begin();
           it2 != mapParToOrigNode.end(); ++it2) {
        //
        int dummy;
        if (it->first != it2->first &&
            (it->first)->IsAncesterOf(it2->first, dummy) == true) {
          // this one is not lowest
          fGood = false;
          break;
        }
      }

      if (fGood) {
        setSibs = it->second;
        fres = true;
        break;
      }
    }
  }

#if 0
    cout << "GetSiblingsPairFrom: \n";
    if( fres == true )
    {
        cout << "sib1: ";
        pairSibs.first->Dump();
        cout << "sib2: ";
        pairSibs.second->Dump();
    }
#endif

  return fres;
}

void PhylogenyTreeBasic ::FindAllLabelsInSubtrees(
    const set<TreeNode *> &setSTRoots, set<string> &setLabels) {
  // get all labels
  setLabels.clear();
  for (set<TreeNode *>::const_iterator it = setSTRoots.begin();
       it != setSTRoots.end(); ++it) {
    set<string> setLblsCoveredStep;
    (*it)->GetAllDistinctLeafLabeles(setLblsCoveredStep);
    UnionSetsGen(setLabels, setLblsCoveredStep);
  }
}

void PhylogenyTreeBasic ::FindDescendentsOfNodeWithin(
    TreeNode *pAnc, const set<TreeNode *> &setNodesToChoose,
    set<TreeNode *> &setDescendents) {
  //
  setDescendents.clear();
  for (set<TreeNode *>::const_iterator itg = setNodesToChoose.begin();
       itg != setNodesToChoose.end(); ++itg) {
    int dummy;
    if (pAnc->IsAncesterOf(*itg, dummy) == true) {
      setDescendents.insert(*itg);
    }
  }
}

bool PhylogenyTreeBasic ::TestIsomorphic(
    PhylogenyTreeBasic &treeOther,
    map<TreeNode *, TreeNode *> &mapOldNodeToNew) const {
#if 0
cout << "TestIsomorphic: current tree: ";
this->Dump();
cout << "treeOther: ";
treeOther.Dump();
#endif
  // return true if isomorphic (and set the mapping between the leaf nodes
  // collect shape label of two trees. Here, we map each current tree node to
  // the corresponding one of the other
  PhylogenyTreeBasic *pthis = const_cast<PhylogenyTreeBasic *>(this);
  set<int> lvidsThis, lvidsOther;
  pthis->GetLeaveIds(lvidsThis);
  treeOther.GetLeaveIds(lvidsOther);
#if 0
cout << "lvidsThis:";
DumpIntSet(lvidsThis);
cout << "lvidsOther:";
DumpIntSet(lvidsOther);
#endif

  map<TreeNode *, string> mapNodeShapeThis, mapNodeShapeOther;
  vector<TreeNode *> listNodesThis, listNodesOther;
  GetAllNodes(listNodesThis);
  treeOther.GetAllNodes(listNodesOther);
  for (int i = 0; i < (int)listNodesThis.size(); ++i) {
    string strShape = listNodesThis[i]->GetShapeLabel(lvidsThis, true);
    mapNodeShapeThis.insert(
        map<TreeNode *, string>::value_type(listNodesThis[i], strShape));
    // cout << "Find a shape (this):" << strShape << endl;
  }
  for (int i = 0; i < (int)listNodesOther.size(); ++i) {
    string strShape = listNodesOther[i]->GetShapeLabel(lvidsOther, true);
    mapNodeShapeOther.insert(
        map<TreeNode *, string>::value_type(listNodesOther[i], strShape));
    // cout << "Find a shape (other):" << strShape << endl;
  }
  if (mapNodeShapeThis[GetRoot()] != mapNodeShapeOther[treeOther.GetRoot()]) {
    // cout << "Root label mismatch: " << mapNodeShapeThis[ GetRoot()]  << " vs
    // " << mapNodeShapeOther[ treeOther.GetRoot() ]  << endl;
    return false; // not isomorphic if the root symbol is not isomorhphic
  }
  // we also list the matching nodes of each node (incl. internal)
  mapOldNodeToNew.clear();
  mapOldNodeToNew.insert(
      map<TreeNode *, TreeNode *>::value_type(GetRoot(), treeOther.GetRoot()));
  stack<TreeNode *> stackNodesToProc;
  stackNodesToProc.push(GetRoot());
  while (stackNodesToProc.empty() == false) {
    TreeNode *pnCurrOld = stackNodesToProc.top();
    stackNodesToProc.pop();
#if 0
cout << "Processing node: ";
pnCurrOld->Dump();
cout << endl;
#endif
    // get all children
    set<TreeNode *> setChildren;
    pnCurrOld->GetAllChildren(setChildren);
    map<string, set<TreeNode *> > setChildrenShape;
    for (set<TreeNode *>::iterator it = setChildren.begin();
         it != setChildren.end(); ++it) {
      TreeNode *pchild = *it;
      string strchild = mapNodeShapeThis[pchild];
      if (setChildrenShape.find(strchild) == setChildrenShape.end()) {
        set<TreeNode *> ss;
        setChildrenShape.insert(
            map<string, set<TreeNode *> >::value_type(strchild, ss));
      }
      setChildrenShape[strchild].insert(pchild);
#if 0
cout << "Adding a string:node pair: " << strchild << ": node: ";
pchild->Dump();
cout << endl;
#endif
      // also save for more processing
      stackNodesToProc.push(pchild);
    }
    // now find the matching one
    set<TreeNode *> setChildOther;
    YW_ASSERT_INFO(mapOldNodeToNew.find(pnCurrOld) != mapOldNodeToNew.end(),
                   "Fai to find1");
    TreeNode *pnCurrOther = mapOldNodeToNew[pnCurrOld];
    pnCurrOther->GetAllChildren(setChildOther);
#if 0
cout << "Now check the matching other: ";
pnCurrOther->Dump();
cout << endl;
#endif
    for (set<TreeNode *>::iterator it = setChildOther.begin();
         it != setChildOther.end(); ++it) {
      TreeNode *pchildother = *it;
      string strchildother = mapNodeShapeOther[pchildother];
#if 0
cout << "child(other): ";
pchildother->Dump();
cout << ": stringshape: " << strchildother << endl;
#endif
      YW_ASSERT_INFO(setChildrenShape.find(strchildother) !=
                             setChildrenShape.end() &&
                         setChildrenShape[strchildother].size() > 0,
                     "Fail to find2");
      // assign to the first one in the list
      TreeNode *pnmatch = *(setChildrenShape[strchildother].begin());
      setChildrenShape[strchildother].erase(pnmatch);
#if 0
cout << "Matching: pncurold: ";
pnCurrOld->Dump();
cout << " to pnmatch:";
pnmatch->Dump();
cout << endl;
#endif
      // remember the matching
      mapOldNodeToNew.insert(
          map<TreeNode *, TreeNode *>::value_type(pnmatch, pchildother));
    }
  }
#if 0
cout << "mapOldNodeToNew: \n";
for( map<TreeNode*, TreeNode*> :: iterator it = mapOldNodeToNew.begin(); it != mapOldNodeToNew.end(); ++it)
{
cout << "oldnode ";
it->first->Dump();
cout << "   to   ";
it->second->Dump();
cout << endl;
}
#endif

  return true;
}

PhylogenyTreeBasic *ConsPhyTreeSubsetTaxa(PhylogenyTreeBasic *ptreeIn,
                                          const set<int> &setTaxaKept) {
  // construct a phylogeny tree by keeping subset of taxa
  PhylogenyTreeBasic *pCopy = ptreeIn->Copy();
  vector<TreeNode *> listLeafNodes;
  pCopy->GetAllLeafNodes(listLeafNodes);

  for (int i = 0; i < (int)listLeafNodes.size(); ++i) {
    int lbl = listLeafNodes[i]->GetIntLabel();
    if (setTaxaKept.find(lbl) == setTaxaKept.end()) {
      // remove this node
      TreeNode *pParOrig = listLeafNodes[i]->GetParent();
      pCopy->RemoveNode(listLeafNodes[i]);
      pCopy->RemoveDegreeOneNodeAt(pParOrig);
    }
  }
  // pCopy->RemoveDegreeOneNodes();

  return pCopy;
}

// implement needed
string ConsEdgeLabeTreeSeg(const string &strNWWithLabels, int regBeg,
                           int regEnd) {
  // cout << "ConsEdgeLabeTreeSeg: [" << regBeg << "," << regEnd << "]: \n";
  // if there is edge outside any parenthesis, keep it
  int posRightParenths = regEnd;
  while (posRightParenths > 0 && strNWWithLabels[posRightParenths] != ')') {
    --posRightParenths;
  }
  string strChild;
  if (posRightParenths > 0) {
    // search for children, perform search for each segment between separator ,
    // (on the same level)
    vector<string> listChildStrs;
    int level = 0;
    int regChildStart = regBeg + 1;
    for (int p = regBeg + 1; p <= posRightParenths - 1; ++p) {
      if ((strNWWithLabels[p] == ',' || p == posRightParenths - 1) &&
          level == 0) {
        int regChildEnd = p - 1;
        if (p == posRightParenths - 1) {
          regChildEnd = p;
        }
        string strChildStep =
            ConsEdgeLabeTreeSeg(strNWWithLabels, regChildStart, regChildEnd);
        if (strChildStep.length() > 0) {
          listChildStrs.push_back(strChildStep);
        }
        regChildStart = p + 1;
      } else if (strNWWithLabels[p] == '(') {
        ++level;
      } else if (strNWWithLabels[p] == ')') {
        --level;
      }
    }
    if (listChildStrs.size() > 0) {
      strChild = "(";
      for (int i = 0; i < (int)listChildStrs.size(); ++i) {
        strChild += listChildStrs[i];
        if (i < (int)listChildStrs.size() - 1) {
          strChild += ",";
        }
      }
      strChild += ")";
    }
  }

  string strEdgeLbelCur;
  if (regEnd != posRightParenths) {
    // search for :
    int pos = regEnd;
    while (pos >= regBeg && strNWWithLabels[pos] != ':') {
      --pos;
    }
    if (pos >= regBeg) {
      strEdgeLbelCur = strNWWithLabels.substr(pos + 1, regEnd - pos);
    }
  }
  string strRes = strChild + strEdgeLbelCur;
  // cout << "strRes: " << strRes << endl;
  return strRes;
}

string ConsEdgeLabeTree(const string &strNWWithLabels) {
  // construct newick format of edge label tree; that is,
  // delete all taxa, only leave edge label
  // e.g. ((2,4:#4):#3,(3:#5,5):#2,1):#1  ==> ((#4)#3,(#5)#2)#1
  return ConsEdgeLabeTreeSeg(strNWWithLabels, 0, strNWWithLabels.length() - 1);
}
