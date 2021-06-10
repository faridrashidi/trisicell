#ifndef GENOTYPE_MATRIX_H
#define GENOTYPE_MATRIX_H

#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>
using namespace std;

#include "BioSequenceMatrix.h"
#include "Utils.h"

typedef pair<int, int> COLUMN_PAIR;

// ***************************************************************************
// Define a reusable binary matrix class
// ***************************************************************************
class GenotypeMatrix : public BioSequenceMatrix {
public:
  GenotypeMatrix();
  ~GenotypeMatrix();
  GenotypeMatrix(int nr, int nc);

  // Support assignment/copy constructor
  GenotypeMatrix(const GenotypeMatrix &rhs);
  GenotypeMatrix &operator=(const GenotypeMatrix &rhs);

  // Important interface functions we need
  virtual bool IsDataValid(int val); // check to see if this data is good for
                                     // this class e.g. for genotype data, 0, 1,
                                     // 2
  virtual bool IsColComplement(int c1, int c2);
  virtual bool IsColDuplicate(int c1, int c2);
  virtual int GetMajorityState(int site);

  // DPPH needs these functions
  void PreSolve(); // perform neccessary preprocessing. For now, assume DPPH
  bool AreColumnsCompanion(int c1, int c2);
  bool AreColumnsForcedInPhase(int c1, int c2);
  bool AreColumnsForcedOutPhase(int c1, int c2);
  bool AreColumnsComplete(int c1, int c2);
  int GetNumTwosInRow(int r);
  bool IsSiteTrival(int site);

private:
  // Internal functions
  void SetupCompanionColumns(); // Initialize companion rows

  // Private data structures
  typedef map<COLUMN_PAIR, set<int> > COMPANION_ROW_MAP;
  COMPANION_ROW_MAP companionRows;
  typedef map<COLUMN_PAIR, int> FORCED_COL_MAP;
  FORCED_COL_MAP forcedColumnPairs; // value = 1 if out of phase, = 0 if in
                                    // phase, otherwise no entry
  vector<COLUMN_PAIR> completePairs;
};

#endif // GENOTYPE_MATRIX_H
