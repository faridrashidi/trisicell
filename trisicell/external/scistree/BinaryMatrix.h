#ifndef BINARY_MATRIX_H
#define BINARY_MATRIX_H

#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>
using namespace std;

#include "BioSequenceMatrix.h"
#include "UnWeightedGraph.h"
#include "Utils3.h"

typedef vector<set<int> > COLUMN_EQUIV_CLASS;

// ***************************************************************************
// Define a reusable binary matrix class
// ***************************************************************************
class BinaryMatrix : public BioSequenceMatrix {
public:
  BinaryMatrix();
  ~BinaryMatrix();
  BinaryMatrix(int nr, int nc);

  // Support assignment/copy constructor
  BinaryMatrix(const BinaryMatrix &rhs);
  BinaryMatrix &operator=(const BinaryMatrix &rhs);

  // Important interface functions we need
  virtual bool IsDataValid(int val); // check to see if this data is good for
                                     // this class e.g. for genotype data, 0, 1,
                                     // 2

  // Matrix editing functions specific to Binary (i.e. haplotype) Matrix
  void TrimDupSites(set<int> *pRemovedSites = NULL, bool fTrimSubsumed = false);
  int FindDupRow();
  void FindNonInformativeSites(set<int> &sitesNoinfo);
  bool TrimNonInformativeSites(set<int> *pRemovedSet = NULL);
  void TrimUniformSites(set<int> *pRemovedSet = NULL);
  void FindUniformSites(set<int> &sitesUniform) const;
  void TrimFullyCompatibleSites(set<int> *pRemovedSet = NULL);
  virtual void TrimNgbrDupCompSites(set<int> *pRemovedSet = NULL);
  void TrimSubsumedRows();
  bool IsRowSubsumedBy(int r1, int r2);
  bool IsColSubsumedBy(int c1, int c2);
  void FindSubsumedSites(set<int> &ssSites);

  // Matrix property checking
  bool IsColNonInformative(int c, int *singletonState);
  bool IsColNonInformative(int c);
  bool IsColTrivial(int c);
  void GetTrivialSites(vector<int> &trivSites);
  bool IsCompatible(int c1, int c2);
  bool IsCompatibleRooted(int c1, int c2, int rallele1, int rallele2);
  bool IsSiteCompatibleWithRegion(int s, int rc1, int rc2);
  bool IsRegionFullyCompatible(int rc1, int rc2);
  void GetGamates(int c1, int c2, bool &f00, bool &f01, bool &f10, bool &f11);
  virtual bool IsColComplement(int c1, int c2);
  virtual bool IsColDuplicate(int c1, int c2);
  bool IsPerfectPhylogeny();
  bool IsZeroColumn(int c);
  bool IsAllColumnsUnique();
  int GetZeroColNum();
  void GetAllIncompatiblePairs(set<pair<int, int> > &incompatibles);
  virtual int GetMajorityState(int site);
  int GetMinorStateNum(int site, int &minorState) const;
  void GetMinorStateRows(int site, int &minorState,
                         set<int> &listRowsWMinor) const;
  void GetRowsWithAllele(int site, int alleleState, set<int> &setRows) const;
  static int GetTheOtherAllele(int allele);

  // Construct interval-speceific equivalance row classes,
  // i.e. sets of row indexes that are same
  void BuildColEquivClasses();
  void GetUniqueColsInRange(int c1, int c2, set<int> &setUniques);
  bool IsSequencesMatch(int r1, int r2, vector<int> &seqColPos);
  void GetSequencesDiffSites(int r1, int r2, set<int> &seqColDiffs) const;

  // Ohter utilities
  void ConstructConflictGraph(UnWeightedGraph &graph);
  void ConflictGraphComponents(vector<BinaryMatrix> &listSubMatrix);
  void ConfigZeroMajSeq(); // make majority elem all-0 for each position
  void ConfigZeroAncesSeq(const vector<int> &seqAnces); // make the matrix s.t.
                                                        // the ancestral state
                                                        // is always 0 in matrix
  void DumpConvGenotypes();
  void GreedyRemoveIncompatSites(
      BinaryMatrix &matReduced); // greedily remove incompatible sites (i.e.
                                 // first remove site that is incompatible w/
                                 // most sites and continue)
  void CalcSFS(vector<double> &listSFSFrac) const;
  int GetDiffSitesForTwoRows(int r1, int r2) const;
  double CalcAvePairRowsDiff() const;
  double CalcAvePairRowsDiffBetween(const set<int> &rowsSet1,
                                    const set<int> &rowsSet2,
                                    double &valMindiffOut) const;
  void CollectAllPairwiseDiffs(const set<int> &rowsSet1,
                               const set<int> &rowsSet2,
                               vector<double> &listRowPairsDiff) const;

  // Missing data utilities
  bool IsColumnBinary(int c) const;
  bool IsRowBinary(int r) const;
  void TrimNonBinaryRows();
  bool IsRowRangeBinary(int r, int left, int right);

  // Lower/upper recombination bound utilities
  int ComputeHKBound();
  int ComputeFastHapBound();
  int ComputeFastRecombUpperBound();
  int ComputeMinRecombWeight(int rowIndex);

private:
  // Interval-based equivlance classes
  COLUMN_EQUIV_CLASS setColEquiv;
};

// some other useful functions
// this structure defines what rows to keep and what not to, and for each
// removed row, which row it comes from (i.e. duplicate) NOTE: we are dealing
// with the current rows only. THat is, the removal may be in stages in each
// stage, we only consider what we have so far
typedef struct {
  set<int> rowsRemoved;
  vector<pair<int, int> > pairsRmKeepRows;
} REMOVED_ROWS_INFO;

void GetNoninformativeRowsInMat(const BinaryMatrix &mat, set<int> &trimedRows,
                                vector<REMOVED_ROWS_INFO> &trimedRowInfo,
                                set<int> &trimedCols, BinaryMatrix &matUpdated,
                                bool fRmDup = false);
void SplitMatrixIntoMaximalFullyCompatRegs(
    const BinaryMatrix &mat, vector<pair<int, int> > &listFullyCompatRegs);

void ReadSitePosFromFirstRowInFile(const char *filename, int numSites,
                                   vector<double> &listSitePos);

#endif // BINARY_MATRIX_H
