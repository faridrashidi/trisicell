#ifndef BIO_SEQUENCE_MATRIX_H
#define BIO_SEQUENCE_MATRIX_H

#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>
using namespace std;

#include "Utils.h"

#define MAX_SITE_NUM 1024000

// ***************************************************************************
// Common Utilities

// ***************************************************************************
// Define a reusable binary matrix class
// ***************************************************************************
class BioSequenceMatrix {
public:
  // BioSequenceMatrix();
  virtual ~BioSequenceMatrix() = 0;

  // Important interface functions we need
  virtual bool IsDataValid(int val) = 0; // check to see if this data is good
                                         // for this class e.g. for genotype
                                         // data, 0, 1, 2
  void SetSize(int nr, int nc);

  // Matrix editing functions
  void AppendRow(const vector<int> &row);
  void AppendSetOfRows(const set<SEQUENCE> &rows);
  void AppendRows(const vector<SEQUENCE> &rows);
  void InsertColumns(const vector<SEQUENCE> &sitesValue,
                     const vector<int> &sitesPos);
  void SetRow(int i, const vector<int> &valNew);
  void SetCol(int i, const vector<int> &valNew);
  void Clear();
  void Copy(const BioSequenceMatrix &rhs);
  virtual bool ReadFromFile(ifstream &inFile, bool fSkipFirstLine = true);
  // virtual bool ReadFromFilePartial( ifstream &inFile, bool fSkipFirstLine );
  void Dump() const;
  void OutputToFile(const char *fileName) const;
  void OutputToFile(ofstream &outFile) const;
  void RemoveRow(int rowIndex);
  void ExchangeColumns(int r1, int r2);
  void RemoveColumns(set<int> &duplicateSites);
  void RemoveRows(set<int> &setRows);
  void TrimDupRows(set<int> *pTrimedRows = NULL,
                   vector<pair<int, int> > *pTrimRowInfo = NULL);
  virtual void FindNgbrDupCompSites(set<int> *pRemovedSet = NULL);
  virtual bool IsColComplement(int c1, int c2) = 0;
  virtual bool IsColDuplicate(int c1, int c2) = 0;
  void AppendMatrixByCol(const BioSequenceMatrix &appendedMat);
  void AppendMatrixByRow(const BioSequenceMatrix &appendedMat);

  // Overload operator for [], like a[1, 2]
  const int &operator()(int r, int c) const;
  int &operator()(int r, int c);
  const int &GetValAt(int r, int c) const;
  void SetValAt(int r, int c, int val);

  // Access matrix
  bool IsEmpty() const { return GetColNum() == 0 || GetRowNum() == 0; }
  int GetColNum() const { return nCols; }
  int GetRowNum() const { return rowsArray.size(); }
  void GetRow(int r, vector<int> &row) const;
  void GetCol(int c, vector<int> &col) const;
  int FindRow(const SEQUENCE &seq) const;
  int FindColumn(const SEQUENCE &seq) const;
  void SubMatrix(int rt, int rb, int cl, int cr,
                 BioSequenceMatrix &submat) const;
  void SubMatrixSelectedSites(const vector<int> &sites,
                              BioSequenceMatrix &submat) const;
  void SubMatrixSelectedRows(const vector<int> &rows,
                             BioSequenceMatrix &submat) const;
  void GetAllSequences(vector<SEQUENCE> &seqs) const;
  void GetSeqsFeqs(map<SEQUENCE, int> &mapSeqFreqs);
  void GetSeqsOccurrence(map<SEQUENCE, set<int> > &mapSeqOccurs);
  virtual int GetMajorityState(int site) = 0;
  void DumpRowMultiplicity() const;
  int GetMultiplictyForRow(int r) const;
  int GetMultiplictyForRow(const SEQUENCE &seq) const;
  int GetMultiplictyForRow(const SEQUENCE &seq, set<int> &identRows) const;
  int GetMultiplictyForRowIV(int r, int left, int right) const;
  void GetColMultiplicityMap(vector<int> &listColMulti) const;
  bool IsIntervalConsistent(int r1, int left1, int right1, int r2, int left2,
                            int right2) const;
  bool IsMissingValue();
  bool IsMissingValueInSite(int c);
  bool IsMissingValueInRow(int r);
  int GetMissingValueNumInRow(int r);
  void MapDupToNodup(map<int, int> &mapDupToNodup) const;
  int GetNodupRowsNum(vector<int> *pListUniqeRowIndex) const;

protected:
  // Some functions
  bool CmpColumns(int c1, int c2);

  // Internal data
  // we represent a binary matrix as bool type
  vector<int *> rowsArray; // array of rows
  int nCols;               // number of sites (columns)

private:
  // Disable certain operations
  BioSequenceMatrix &operator=(const BioSequenceMatrix &rhs) { return *this; }
  // bool fMissingValue;
};

#endif // BIO_SEQUENCE_MATRIX_H
