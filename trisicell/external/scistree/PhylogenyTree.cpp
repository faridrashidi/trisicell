#include "PhylogenyTree.h"
#include <fstream>
#include <iostream>
#include <stack>

// ***************************************************************************
// The following code is largely based on Gusfield's 1991 Paper
// ***************************************************************************
extern void OutputQuotedString(ofstream &outFile, const char *buf);

// ***************************************************************************
// Utilites functions
// ***************************************************************************

int PhylogenyTree ::GetIntLabelFromParenthStr(const string &strLabelWParenth) {
  //
  YW_ASSERT_INFO(strLabelWParenth[0] == '(' &&
                     strLabelWParenth[strLabelWParenth.length() - 1] == ')',
                 "String does not come with ()");
  string strPrune = strLabelWParenth.substr(1, strLabelWParenth.length() - 2);
  int res = -1;
  sscanf(strPrune.c_str(), "%d", &res);
  return res;
}

void PhylogenyTree ::GetARoot(const BinaryMatrix &mat, vector<int> &root) {
  if (knownRoot.size() > 0) {
    root = knownRoot;
    return;
  }

  // We take the majority sequence as root. Refer to the paper for details
  root.clear();
  for (int c = 0; c < mat.GetColNum(); ++c) {
    int rc = 0;
    int numOne = 0;
    for (int r = 0; r < mat.GetRowNum(); ++r) {
      if (mat(r, c) == 1) {
        numOne++;
      }
    }
    // 12/08/07: fixed. Must consider the case say 6 0 and 5 1,
    // has to plus one to ensure correctness
    if (numOne >= (mat.GetRowNum() + 1) / 2) {
      rc = 1;
    }
    root.push_back(rc);
  }
  //    cout << "Root = ";
  //    DumpIntVec ( root );
}

void PhylogenyTree ::RadixSortByCol(const BinaryMatrix &mat,
                                    const vector<int> &root,
                                    vector<int> &sortList) {
  // cout << "root = ";
  // DumpIntVec( root );
  // This is the step 1 of Gusfield tree building algorithm
  // We treat each column as a number, encoded by the binary vecgtor stored in
  // the column row 1 contains the MSB of the number The result is stored in a
  // sorted list, with LARGEST number comes in first For details of radix sort,
  // refer CLR
  sortList.clear();
  for (int i = 0; i < mat.GetColNum(); ++i) {
    sortList.push_back(i);
  }

  // Now sort from LSB of the number, i.e. last row first
  for (int i = mat.GetRowNum() - 1; i >= 0; --i) {
    SortByOneBit(i, mat, root, sortList);
  }
}

void PhylogenyTree ::SortByOneBit(int bitPosRow, const BinaryMatrix &mat,
                                  const vector<int> &root,
                                  vector<int> &sortList) {
  // cout << "bitPosRow = " << bitPosRow << endl;
  // cout << "root here = ";
  // DumpIntVec( root );
  // cout << "entry sortList = ";
  // DumpIntVec( sortList );
  // Sort the list by one bit (the ith row)
  // Initailize a pre-list, holding the last sorted list. Simply initailize to
  // original order
  vector<int> preList = sortList;
  sortList.clear();

  // We do two path, first to find 1 cells in that row and next one cell (since
  // we want the LARGEST first) This is in fact counting sort, with k (the
  // limit) == 1
  for (int i = 0; i < preList.size(); ++i) {
    // Note that we 1 = NON-ROOT-VALUE
    // cout << "mat(bitPosRow, preList[i] ) = " << mat( bitPosRow, preList[i]  )
    // << endl;
    if (mat(bitPosRow, preList[i]) != root[preList[i]]) {
      sortList.push_back(preList[i]);
    }
  }
  // cout << "parital sortList = ";
  // DumpIntVec( sortList );

  for (int i = 0; i < preList.size(); ++i) {
    if (mat(bitPosRow, preList[i]) == root[preList[i]]) {
      sortList.push_back(preList[i]);
    }
  }
  // cout << "exit sortList = ";
  // DumpIntVec( sortList );
}

void PhylogenyTree ::RemoveDupSites(const BinaryMatrix &mat,
                                    vector<int> &sortedPosList,
                                    vector<vector<int> > &duplicates) {
  // This function takes the sorted list, and then remove the duplicate sites
  // by comparing one site to its left row, if duplicate, do not put into new
  // list
  vector<int> noDupList;
  if (sortedPosList.size() > 0) {
    noDupList.push_back(sortedPosList[0]);
  }
  vector<int> dupList; // store which sites are duplicates to this one
  for (int i = 1; i < sortedPosList.size(); ++i) {
    bool match = true;
    // Check to see if this column is the same as its immediate left one
    for (int r = 0; r < mat.GetRowNum(); ++r) {
      if (mat(r, sortedPosList[i]) != mat(r, sortedPosList[i - 1])) {
        match = false;
        break;
      }
    }
    if (match == false) {
      noDupList.push_back(sortedPosList[i]);

      // Now we maintian the duplicate list
      // cout << "for site " << noDupList[noDupList.size() - 2] << ", duplicate
      // sites are: "; DumpIntVec( dupList );

      duplicates.push_back(dupList);
      dupList.clear();
    } else {
      // This site is the same as its immediate left one
      dupList.push_back(sortedPosList[i]);
    }
  }

  // Finally, add the final list to it
  duplicates.push_back(dupList);
  // cout << "for site " << noDupList[noDupList.size() - 1] << ", duplicate
  // sites are: "; DumpIntVec( dupList );
  dupList.clear();

  // Now set the noDupList to result
  sortedPosList.clear();
  sortedPosList = noDupList;
}

void PhylogenyTree ::ComputeLijLj(const BinaryMatrix &mat,
                                  const vector<int> &root,
                                  const vector<int> &sortedPosList,
                                  vector<int *> &Lij, vector<int> &Lj) {
  //    cout << "sortedPosList = ";
  //    DumpIntVec( sortedPosList );

  // Build Lij and Lj according to the algorithm
  // CAUTION: you have to keep in mind that Lij, Lj are all based on M', not M
  // so do a conversion before use
  for (int i = 0; i < mat.GetRowNum(); ++i) {
    int last1Pos = -1;
    for (int j = 0; j < sortedPosList.size(); ++j) {
      if (mat(i, sortedPosList[j]) != root[sortedPosList[j]]) {
        // We find a one here, good
        Lij[i][j] = last1Pos;

        // cout << "at (" << i << ", " << j << "), Lij = " << last1Pos << endl;

        // Remember it
        last1Pos = j;
      }
    }
  }

  // Now we computes the Lj vector
  Lj.clear();
  for (int j = 0; j < sortedPosList.size(); ++j) {
    int max = -1;
    for (int r = 0; r < mat.GetRowNum(); ++r) {
      if (mat(r, sortedPosList[j]) != root[sortedPosList[j]] &&
          Lij[r][j] > max) {
        max = Lij[r][j];
      }
    }
    // Now set Lj
    Lj.push_back(max);
    // cout << "At j = " << j << ", Lj = " << max << endl;
  }
}

bool PhylogenyTree ::ExamineLijLj(const BinaryMatrix &mat,
                                  const vector<int> &root,
                                  const vector<int> &sortedPosList,
                                  const vector<int *> &Lij,
                                  const vector<int> &Lj) {
  // cout << "Examine here...\n";
  for (int i = 0; i < mat.GetRowNum(); ++i) {
    for (int j = 0; j < sortedPosList.size(); ++j) {
      if (mat(i, sortedPosList[j]) != root[sortedPosList[j]] &&
          Lj[j] != Lij[i][j]) {
        // cout << "At (" << i << ", " << j << "), Lij = " << Lij[i][j] << ",
        // but Lj = " << Lj[j] << endl;
        return false;
      }
    }
  }
  // cout << "done here.\n";
  return true; // yes, there is a tree
}

void PhylogenyTree ::BuildTree(const BinaryMatrix &mat, const vector<int> &root,
                               const vector<int> &sortedPosList,
                               const vector<vector<int> > &duplicates,
                               const vector<int> &Lj) {
  // This function creates the tree by creating and linking tree nodes
  // Make sure the tree is empty
  if (rootNode != NULL) {
    delete rootNode;
    rootNode = NULL;
  }

  // root is labeled as -1, since all other (column) nodes are labeled by a site
  rootNode = new TreeNode(-1);

  // Create a node for each site
  vector<TreeNode *> colNodes;
  for (int i = 0; i < sortedPosList.size(); ++i) {
    TreeNode *pNode =
        new TreeNode(sortedPosList[i]); // for now, use original labels to do it
    colNodes.push_back(pNode);
  }

  // Link each node Nj (where L(j) >= 0) to that L(j) node
  for (int j = 0; j < Lj.size(); ++j) {
    // Figure out the labels
    vector<int> labels;
    labels.push_back(sortedPosList[j]);
    // Add those in the duplicates
    for (int dup = 0; dup < duplicates[j].size(); ++dup) {
      labels.push_back(duplicates[j][dup]);
    }
    if (Lj[j] >= 0) {
      // Link it
      TreeNode *nodeLj = colNodes[Lj[j]];

      // Add it
      nodeLj->AddChild(colNodes[j], labels);
      // cout << "Add col node " << sortedPosList[j]  << " under node " <<
      // sortedPosList[ Lj[j] ] << ".\n";
    } else {
      // For this node, we link it from the root
      rootNode->AddChild(colNodes[j], labels);
      // cout << "Add col node " << sortedPosList[j]  << " under root.\n";
    }
  }

  // Now add rows into this tree
  for (int i = 0; i < mat.GetRowNum(); ++i) {
    int ci = -1;
    // Find ci that is the largest cell has one in row i
    for (int j = sortedPosList.size() - 1; j >= 0; j--) {
      if (mat(i, sortedPosList[j]) != root[sortedPosList[j]]) {
        ci = j;
        break;
      }
    }
    if (ci < 0) {
      //    cout << "trouble here.\n";
      //    YW_ASSERT(false);
      // This is the same as the root sequence
      TreeNode *pLeaf =
          new TreeNode(mat.GetColNum() + i); // Use id=row index + colNum
      pLeaf->AddNodeValue(i);
      // also set its label
      char buf[100], buf1[100];
      sprintf(buf, "(%d)", i);
      sprintf(buf1, "%d", i);
      pLeaf->SetLabel(buf);
      pLeaf->SetUserLabel(buf1);

      vector<int> emptyLabel;
      rootNode->AddChild(pLeaf, emptyLabel);
      // cout << "Add row " << i << " under root node.\n";

    } else {
      // Here we always add a node as children. CAUTION: here we may create
      // degree-2 nodes, we need to cleanup after this 06/05/05: actually I
      // decided to go another way: put to leaf first, then splits the multiple
      // labels into different leaves if needed
      TreeNode *pn = colNodes[ci];
      if (pn->IsLeaf() == true) {
        // also set its label
        char buf[100], buf1[100];
        sprintf(buf, "(%d)", i);
        sprintf(buf1, "%d", i);
        pn->SetLabel(buf);
        pn->SetUserLabel(buf1);

        // Now attach this row to the existing leaf, HOW?
        pn->AddNodeValue(i);
        // cout << "Add row " << i << " to a leaf (col node) " <<
        // sortedPosList[ci]  << ".\n";

      } else {
        TreeNode *pLeaf =
            new TreeNode(mat.GetColNum() + i); // Use id=row index + colNum
        pLeaf->AddNodeValue(i);
        // also set its label
        char buf[100], buf1[100];
        sprintf(buf, "(%d)", i);
        sprintf(buf1, "%d", i);
        pLeaf->SetLabel(buf);
        pLeaf->SetUserLabel(buf1);

        vector<int> emptyLabel;
        pn->AddChild(pLeaf, emptyLabel);
        // cout << "Add row " << i << " to a non-leaf (col node) " <<
        // sortedPosList[ci]  << ".\n";
      }
    }
  }
}

void PhylogenyTree ::CleanupTree(const BinaryMatrix &mat) {
  // 06/05/05: take another route, breakup multiple labels
  TreeNode *curTN = NULL;
  stack<TreeNode *> stackNodes;
  if (rootNode != NULL) {
    stackNodes.push(rootNode);
  }

  while (stackNodes.empty() == false) {
    // Move to next node in stack
    curTN = stackNodes.top();
    stackNodes.pop();

    // For a leaf, we try to split it
    if (curTN->IsLeaf() == true && curTN->nodeValues.size() > 1) {
      for (int i = 0; i < curTN->nodeValues.size(); ++i) {
        // Find one to split
        TreeNode *pLeaf =
            new TreeNode(mat.GetColNum() +
                         curTN->nodeValues[i]); // Use id=row index + colNum
        pLeaf->AddNodeValue(curTN->nodeValues[i]);
        vector<int> emptyLabel;
        curTN->AddChild(pLeaf, emptyLabel);
        // cout << "Spliting row " << curTN->nodeValues[i] << " from leaf " <<
        // curTN->id  << ".\n";

        // Set the label to the individual nodes values
        char buf[100], buf1[100];
        sprintf(buf, "(%d)", curTN->nodeValues[i]);
        sprintf(buf1, "%d", curTN->nodeValues[i]);
        pLeaf->SetLabel(buf);
        pLeaf->SetUserLabel(buf1);
      }

      // Finally, clear the labels at parent node
      curTN->nodeValues.clear();

      // We also clear the old label
      curTN->SetLabel("-");
      curTN->SetUserLabel("-");
    }

    // push children into stack
    for (int i = 0; i < curTN->listChildren.size(); ++i) {
      stackNodes.push(curTN->listChildren[i]);
    }
  }
}

void PhylogenyTree ::RemoveDegreeTwoNodes() {
  // This function removes all degree-2 nodes
  // we start from the root and remove any node with degree 2
  TreeNode *curTN = NULL;
  stack<TreeNode *> stackNodes;
  if (rootNode != NULL) {
    stackNodes.push(rootNode);
  }

  while (stackNodes.empty() == false) {
    // Move to next node in stack
    curTN = stackNodes.top();
    stackNodes.pop();

    // push children into stack
    for (int i = 0; i < curTN->listChildren.size(); ++i) {
      stackNodes.push(curTN->listChildren[i]);
    }

    // any node, if it has only a single child, remove the current node
    if (curTN->IsLeaf() == false && curTN->GetChildrenNum() == 1) {
      // remove it
      TreeNode *pcnode = curTN->listChildren[0];
      TreeNode *ppar = curTN->GetParent();

      vector<int> listLblpn;
      curTN->GetEdgeLabelsAtBranch(0, listLblpn);

      // change cur's par if exist
      if (ppar != NULL) {
        // construct the concatnated label list
        int pindex = ppar->GetChildIndex(curTN);
        vector<int> listLblpn2;
        ppar->GetEdgeLabelsAtBranch(pindex, listLblpn2);
        AppendIntVec(listLblpn, listLblpn2);

        // here need to maintian the edge labesl
        ppar->RemoveChild(curTN);
        // vector<int> labelsEmpty;
        ppar->AddChild(pcnode, listLblpn);
      } else {
        // cur node is root, then change the root
        YW_ASSERT_INFO(curTN == rootNode, "Must be root");
        rootNode = pcnode;
      }

      // set new parent
      pcnode->SetParent(ppar);
    }
  }
}

// ***************************************************************************
// Main functions
// ***************************************************************************

PhylogenyTree ::PhylogenyTree() {}

PhylogenyTree ::~PhylogenyTree() {}

bool PhylogenyTree ::ConsOnBinMatrix(const BinaryMatrix &mat) {
  // Build tree from binary matrix
  vector<int> sortedPosList;

  // We first find a good root from data
  vector<int> root;
  GetARoot(mat, root);

  // We first sort columns (treated as binary number) by putting the largest
  // first
  RadixSortByCol(mat, root, sortedPosList);

  // cout << "the sorted column list is: \n";
  // DumpIntVec( sortedPosList);

  // Remove Duplicate columns
  vector<vector<int> > listDuplicates; // used to save for each one in
                                       // sortedPosList the sites to its right
                                       // that is duplicate as it, in ORIGINAL
                                       // numbering
  RemoveDupSites(mat, sortedPosList, listDuplicates);
  // cout << "the no duplicate sorted column list is: \n";
  // DumpIntVec( sortedPosList);

  // Now we compute the Lij and Lj values, from Gusfield's algorithm
  vector<int *> Lij;
  for (int i = 0; i < mat.GetRowNum(); ++i) {
    int *pbuf = new int[sortedPosList.size()];
    Lij.push_back(pbuf);
  }
  vector<int> Lj;
  ComputeLijLj(mat, root, sortedPosList, Lij, Lj);
  if (ExamineLijLj(mat, root, sortedPosList, Lij, Lj) == false) {
    cout << "No tree.\n";
    return false; // no tree
  }
  // Now we start to build tree here
  BuildTree(mat, root, sortedPosList, listDuplicates, Lj);
  // cout << "Yes, there is a tree here.\n";

  // Finally, we cleanup
  CleanupTree(mat);

  // Now we have to do cleanup
  for (int i = 0; i < Lij.size(); ++i) {
    delete[] Lij[i];
  }

  return true;
}

void PhylogenyTree ::GetLeavesWithMatRowIndices(const set<int> &setMatRows,
                                                set<TreeNode *> &setLeaves) {
  // cout << "GetLeavesWithMatRowIndices: setMatRows = ";
  // DumpIntSet( setMatRows );
  // given a set of row indices in mat (assume this is one where phylogeny is
  // constructed)
  set<string> setLabel;
  for (set<int>::iterator it = setMatRows.begin(); it != setMatRows.end();
       ++it) {
    // use the same naming convention
    char buf[100];
    // sprintf(buf, "%d", *it);
    sprintf(buf, "(%d)", *it);
    string lbl(buf);
    setLabel.insert(lbl);
  }
  GetLeavesWithLabels(setLabel, setLeaves);
}

// ***************************************************************************

string ConsRootedPerfectPhylogenyFromMat(const BinaryMatrix &matInput,
                                         bool fEdgeLabel, bool fOneBase) {
  // constructed tree assuming zero-rooted tree
  // collect rooted splits
  set<set<int> > setRootedSplits;
  map<set<int>, set<int> > mapSplitSites;
  set<int> setAll1sSites;
  for (int s = 0; s < matInput.GetColNum(); ++s) {
    set<int> split;
    matInput.GetRowsWithAllele(s, 1, split);
    mapSplitSites[split].insert(s + 1); // let site start from index 1
    setRootedSplits.insert(split);

    if (split.size() == matInput.GetRowNum()) {
      setAll1sSites.insert(s);
    }

    // cout << "Site " << s << " split: ";
    // DumpIntSet(split);
  }

  // cout << "Set of all-1 sites: ";
  // DumpIntSet(setAll1sSites);

#if 0
    vector<string> listSiteNames;
    if(pListSiteNames != NULL )
    {
        listSiteNames = *pListSiteNames;
    }
    else
    {
        //  use simple site name, starting from 1
        for(int s=0; s<matInput.GetColNum(); ++s)
        {
            string str = std::to_string(s+1);
            listSiteNames.push_back(str);
        }
    }
#endif

  //
  PhylogenyTreeBasic tree;
  CreatePhyTreeWithRootedSplits(tree, matInput.GetRowNum(), setRootedSplits);

  // setup edge labels if needed
  if (fEdgeLabel) {
    tree.RemoveEdgeLabels();
    //
    vector<TreeNode *> listNodes;
    tree.GetAllNodes(listNodes);
    for (int i = 0; i < (int)listNodes.size(); ++i) {
      if (listNodes[i]->IsLeaf()) {
        continue;
      }
      // check all children
      for (int j = 0; j < listNodes[i]->GetChildrenNum(); ++j) {
        TreeNode *pChild = listNodes[i]->GetChild(j);
        set<int> setLeavesUnder;
        pChild->GetAllLeavesIdUnder(setLeavesUnder);
        // cout << "The " << j << " th child: leaves under: ";
        // DumpIntSet(setLeavesUnder);
        if (mapSplitSites.find(setLeavesUnder) != mapSplitSites.end()) {
          set<int> setEdgeLbels = mapSplitSites[setLeavesUnder];
          for (set<int>::iterator it = setEdgeLbels.begin();
               it != setEdgeLbels.end(); ++it) {
            listNodes[i]->AddEdgeLabelToChild(j, *it);
          }
        }
      }
    }
  }

  if (fOneBase) {
    map<int, int> mapIncLeafLbls;
    for (int i = 0; i < matInput.GetRowNum(); ++i) {
      mapIncLeafLbls[i] = i + 1;
    }
    ChangeLeafIntLabelOfTree(tree, mapIncLeafLbls);
  }

  string res;
  if (fEdgeLabel == false) {
    tree.ConsNewick(res);
  } else {
    tree.ConsNewickEdgeLabel(res);
    if (setAll1sSites.size() > 0) {
      res += ":";
      // add all-1 labels at the top
      for (set<int>::iterator it = setAll1sSites.begin();
           it != setAll1sSites.end(); ++it) {
        int ss = *it;
        string strId = std::to_string(ss + 1);
        res += "#" + strId;
      }
    }
  }
  return res;
}
