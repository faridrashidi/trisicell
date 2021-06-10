#ifndef PHYLOGENY_TREE_H
#define PHYLOGENY_TREE_H

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

#include "BinaryMatrix.h"
#include "PhylogenyTreeBasic.h"
#include "Utils.h"

using namespace std;

// ***************************************************************************
// Define phylogeny tree class, which extends from the basic phylogeny class
// the main purpose is to support building from matrix (perfect phylogeny)
// ***************************************************************************

class PhylogenyTree : public PhylogenyTreeBasic {
public:
  PhylogenyTree(); // Empty tree
  virtual ~PhylogenyTree();
  bool ConsOnBinMatrix(const BinaryMatrix &mat); // Build tree from binary
                                                 // matrix
  void SetRoot(const vector<int> &rootToSet) { knownRoot = rootToSet; }
  void RemoveDegreeTwoNodes();
  static int GetIntLabelFromParenthStr(const string &strLabelWParenth);
  void GetLeavesWithMatRowIndices(const set<int> &setMatRows,
                                  set<TreeNode *> &setLeaves);

private:
  void GetARoot(const BinaryMatrix &mat, vector<int> &root);
  void RadixSortByCol(const BinaryMatrix &mat, const vector<int> &root,
                      vector<int> &sortList);
  void SortByOneBit(int bitPosRow, const BinaryMatrix &mat,
                    const vector<int> &root, vector<int> &sortList);
  void RemoveDupSites(const BinaryMatrix &mat, vector<int> &sortedPosList,
                      vector<vector<int> > &duplicates);
  void ComputeLijLj(const BinaryMatrix &mat, const vector<int> &root,
                    const vector<int> &sortedPosList, vector<int *> &Lij,
                    vector<int> &Lj);
  bool ExamineLijLj(const BinaryMatrix &mat, const vector<int> &root,
                    const vector<int> &sortedPosList, const vector<int *> &Lij,
                    const vector<int> &Lj);
  void BuildTree(const BinaryMatrix &mat, const vector<int> &root,
                 const vector<int> &sortedPosList,
                 const vector<vector<int> > &duplicates, const vector<int> &Lj);
  void CleanupTree(const BinaryMatrix &mat);

  vector<int> knownRoot;
};

// ***************************************************************************

string ConsRootedPerfectPhylogenyFromMat(const BinaryMatrix &matInput,
                                         bool fEdgeLabel,
                                         bool fOneBase = false);

#endif // PHYLOGENY_TREE_H
