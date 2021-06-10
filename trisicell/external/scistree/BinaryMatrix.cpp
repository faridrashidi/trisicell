#include "BinaryMatrix.h"
#include "Utils2.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

// ***************************************************************************
// Define a reusable binary matrix class
// ***************************************************************************

BinaryMatrix ::BinaryMatrix() { nCols = 0; }

BinaryMatrix ::~BinaryMatrix() {
  // Need to free up data if needed
  Clear();
}

BinaryMatrix ::BinaryMatrix(int nr, int nc) { SetSize(nr, nc); }

BinaryMatrix ::BinaryMatrix(const BinaryMatrix &rhs) { Copy(rhs); }

BinaryMatrix &BinaryMatrix ::operator=(const BinaryMatrix &rhs) {
  Clear();

  Copy(rhs);

  return *this;
}

bool BinaryMatrix ::IsDataValid(int val) {
  if (val == 0 || val == 1) {
    return true;
  } else {
    return false;
  }
}

//#if 0
void BinaryMatrix ::TrimNgbrDupCompSites(set<int> *pRemovedSet) {
  set<int> setOfRemovals; // contains sites to be removed
  int cleft = 0;
  while (cleft < nCols - 1) {
    // Check to see if the next row  immediately is complement or not
    if (IsColComplement(cleft, cleft + 1) == true ||
        IsColDuplicate(cleft, cleft + 1) == true) {
      setOfRemovals.insert(cleft + 1);
      // cout << "Site " << cleft+1 << " is same/complement." << endl;
    }
    // Consider  next site
    cleft++;
  }
  if (pRemovedSet != NULL) {
    pRemovedSet->clear();
    *pRemovedSet = setOfRemovals;
  }
  // Finally, remove columns
  RemoveColumns(setOfRemovals);
}
//#endif

// Consolidate columns in matrix
void BinaryMatrix::TrimDupSites(set<int> *pRemovedSites, bool fTrimSubsumbed) {
  int c1, c2;
  unsigned int r;
  set<int> setOfDuplicates; // contains sites to be removed

  for (c1 = 0; c1 < nCols; ++c1) {
    for (c2 = c1 + 1; c2 < nCols; ++c2) {
      // now we compare these two cols: c1, c2
      // if they match, we put c2 into set
      bool f = false;
      for (r = 0; r < rowsArray.size(); ++r) {
        // compare each cell
        if (rowsArray[r][c1] != rowsArray[r][c2]) {
          f = true;
          break;
        }
      }

      // Check against size
      if (r == rowsArray.size()) {
        // we find a duplicate
        if (setOfDuplicates.find(c2) == setOfDuplicates.end()) {
          //					cout <<  "Site " << c2 << " is duplicate of
          //site
          //"; 					cout << c1 << endl;
        }
        setOfDuplicates.insert(c2);
      }
    }
  }
  if (fTrimSubsumbed == true) {
    // cout << "Now start to find subsumbed sites...\n";
    FindSubsumedSites(setOfDuplicates);
  }

  // Now save the trimed sites info, if needed
  if (pRemovedSites != NULL) {
    *pRemovedSites = setOfDuplicates;
  }

  // Finally, remove columns
  RemoveColumns(setOfDuplicates);
}

void BinaryMatrix ::TrimSubsumedRows() {
  // Dump();
  set<int> ssRows;
  for (int r1 = 0; r1 < GetRowNum(); ++r1) {
    for (int r2 = 0; r2 < GetRowNum(); ++r2) {
      if (r1 == r2) {
        continue;
      }
      if (IsRowSubsumedBy(r1, r2) == true) {
        ssRows.insert(r1);
      }
    }
  }
  // cout << "ssRows = ";
  // DumpIntSet( ssRows );
  // if( ssRows.size() > 0 )
  //{
  //    exit(1);
  //}
  RemoveRows(ssRows);
}

bool BinaryMatrix ::IsRowSubsumedBy(int r1, int r2) {
  // Test whether a row is subsumed by another row
  bool fRes = true;
  bool fEqual = true;

  for (int c = 0; c < nCols; ++c) {
    if (rowsArray[r1][c] != rowsArray[r2][c]) {
      fEqual = false;
      if (IsMissingValueBit(rowsArray[r1][c]) == false) {
        fRes = false;
        break;
      }
    }
  }

  if (fEqual == true) {
    // do not consider two identical rows are subsumbed by another
    return false;
  }
  return fRes;
}

bool BinaryMatrix ::IsColSubsumedBy(int c1, int c2) {
  // Test whether a row is subsumed by another row
  bool fRes = true;
  bool fEqual = true;

  for (int r = 0; r < GetRowNum(); ++r) {
    if (rowsArray[r][c1] != rowsArray[r][c2]) {
      fEqual = false;
      if (IsMissingValueBit(rowsArray[r][c1]) == false) {
        fRes = false;
        break;
      }
    }
  }

  if (fEqual == true) {
    // do not consider two identical rows are subsumbed by another
    return false;
  }
  return fRes;
}

void BinaryMatrix ::FindSubsumedSites(set<int> &ssSites) {
  // Dump();
  for (int c1 = 0; c1 < GetColNum(); ++c1) {
    for (int c2 = 0; c2 < GetColNum(); ++c2) {
      if (c1 == c2) {
        continue;
      }
      if (IsColSubsumedBy(c1, c2) == true) {
        // cout << "site c1 = " << c1 << " is subsumed by c2 = " << c2 << endl;
        ssSites.insert(c1);
        break;
      }
    }
  }
  // cout << "ssSites = ";
  // DumpIntSet( ssSites );
  // if( ssSites.size() > 0 )
  //{
  //    exit(1);
  //}
}

int BinaryMatrix ::FindDupRow() {
  // This function tracking any removal of rows, but
  // in addition to it, we track which rows remains
  unsigned int r1, r2;
  int c;

  for (r1 = 0; r1 < rowsArray.size(); ++r1) {
    for (r2 = r1 + 1; r2 < rowsArray.size(); ++r2) {
      /*
              Now test whether row 1 and row 2 are the same
      */
      bool fSame = true;
      for (c = 0; c < nCols; ++c) {
        if (rowsArray[r1][c] != rowsArray[r2][c]) {
          fSame = false;
          break;
        }
      }
      if (fSame) {
        // cout << "row " << r2 << " is duplicate." << endl;
        return r2;
      }
    }
  }

  return -1;
}

void BinaryMatrix ::FindNonInformativeSites(set<int> &sitesNoinfo) {
  sitesNoinfo.clear();

  // find set of non-informative sites
  int c1;
  unsigned int r;

  for (c1 = 0; c1 < nCols; ++c1) {
    int numZeros = 0, numOnes = 0;
    // now we compare these two cols: c1, c2
    // if they match, we put c2 into set
    for (r = 0; r < rowsArray.size(); ++r) {
      if (rowsArray[r][c1] == 0) {
        numZeros++;
#if 0
				if(numZeros >=2 && numOnes >=2)
				{
					break;
				}
#endif
      } else if (rowsArray[r][c1] == 1) {
        numOnes++;
#if 0
				if(numZeros >=2 &&  numOnes >= 2)
				{
					break;
				}
#endif
      }
    }
    // Check to see if this is non-informative
    if (numZeros <= 1 || numOnes <= 1) {
      // we find a duplicate
      //			cout << "Site  " << c1+1 << "is non-informative"
      //<< endl;
      sitesNoinfo.insert(c1);
    }
  }
}

void BinaryMatrix ::FindUniformSites(set<int> &sitesUniform) const {
  //
  sitesUniform.clear();

  // find set of non-informative sites
  int c1;
  unsigned int r;

  for (c1 = 0; c1 < nCols; ++c1) {
    int numZeros = 0, numOnes = 0;
    // now we compare these two cols: c1, c2
    // if they match, we put c2 into set
    for (r = 0; r < rowsArray.size(); ++r) {
      if (rowsArray[r][c1] == 0) {
        numZeros++;
      } else if (rowsArray[r][c1] == 1) {
        numOnes++;
      }
    }
    // Check to see if this is non-informative
    if (numZeros == 0 || numOnes == 0) {
      // we find a duplicate
      //			cout << "Site  " << c1+1 << "is non-informative"
      //<< endl;
      sitesUniform.insert(c1);
    }
  }
}

/*
        Remove all non-informative sites
        A site is non-informative if it is all 0 (1), or has only single 0(1)
*/
bool BinaryMatrix ::TrimNonInformativeSites(set<int> *pRemovedSet) {
  set<int> setOfDuplicates;
  FindNonInformativeSites(setOfDuplicates);
  if (pRemovedSet != NULL) {
    *pRemovedSet = setOfDuplicates;
  }

  // Finally, remove columns
  bool res = false;
  if (setOfDuplicates.size() > 0) {
    res = true;
    RemoveColumns(setOfDuplicates);
  }
  return res;
}

void BinaryMatrix ::TrimUniformSites(set<int> *pRemovedSet) {
  set<int> setOfDuplicates;
  FindUniformSites(setOfDuplicates);
  if (pRemovedSet != NULL) {
    *pRemovedSet = setOfDuplicates;
  }

  // Finally, remove columns
  if (setOfDuplicates.size() > 0) {
    RemoveColumns(setOfDuplicates);
  }
}

void BinaryMatrix ::TrimFullyCompatibleSites(set<int> *pRemovedSet) {
  int c1, c2;
  set<int> setOfDuplicates; // contains sites to be removed
  for (c1 = 0; c1 < nCols; ++c1) {
    // now we compare these two cols: c1, c2
    // if they match, we put c2 into set
    bool f = true; // by default, we say f is fully-compatible
    // Now we test whether sites c1 is compatible with c2
    for (c2 = 0; c2 < nCols; ++c2) {
      if (IsCompatible(c1, c2) == false) {
        f = false;
        break;
      }
    }
    if (f == true && IsColumnBinary(c1) == true) {
      // cout << "Site " << c1+1 << " is fully compatible" << endl;
      setOfDuplicates.insert(c1);
    }
  }

  // Now remember the set if needed
  if (pRemovedSet != NULL) {
    pRemovedSet->clear();
    *pRemovedSet = setOfDuplicates;
  }

  // Finally, remove columns
  RemoveColumns(setOfDuplicates);
}

bool BinaryMatrix ::IsAllColumnsUnique() {
  bool res = true;

  for (int i = 0; i < nCols - 1; ++i) {
    for (int j = i + 1; j < nCols; ++j) {
      // check to see if column i, j are duplicate
      if (CmpColumns(i, j) == true) {
        return false;
      }
    }
  }

  return res;
}

bool BinaryMatrix ::IsColNonInformative(int c) {
  int numZeros = 0, numOnes = 0, numMissing = 0;
  // now we compare these two cols: c1, c2
  // if they match, we put c2 into set
  for (unsigned int r = 0; r < rowsArray.size(); ++r) {
    if (rowsArray[r][c] == 0) {
      numZeros++;
#if 0
			if(numZeros >=2 && numOnes >=2)
			{
				break;
			}
#endif
    } else if (rowsArray[r][c] == 1) {
      numOnes++;
#if 0
			if(numZeros >=2 &&  numOnes >= 2)
			{
				break;
			}
#endif
    } else if (IsMissingValueBit(rowsArray[r][c]) == true) {
      numMissing++;
    }
  }
  // Check to see if this is non-informative
  if ((numZeros == 1 || numOnes == 1) && numMissing == 0) {
    // we find a duplicate
    //			cout << "Site  " << c1+1 << "is  non-informative" <<
    // endl;
    return true;
  } else {
    return false;
  }
}

bool BinaryMatrix ::IsColNonInformative(int c, int *singletonState) {
  int numZeros = 0, numOnes = 0;
  // now we compare these two cols: c1, c2
  // if they match, we put c2 into set
  for (unsigned int r = 0; r < rowsArray.size(); ++r) {
    if (rowsArray[r][c] == 0) {
      numZeros++;
#if 0
			if(numZeros >=2 && numOnes >=2)
			{
				break;
			}
#endif
    } else if (rowsArray[r][c] == 1) {
      numOnes++;
#if 0
			if(numZeros >=2 &&  numOnes >= 2)
			{
				break;
			}
#endif
    }
  }
  // Check to see if this is non-informative
  if (numZeros == 1 || numOnes == 1) {
    if (singletonState != NULL) {
      if (numZeros == 1) {
        *singletonState = 0;
      } else {
        *singletonState = 1;
      }
    }
    // we find a duplicate
    //			cout << "Site  " << c1+1 << "is  non-informative" <<
    // endl;
    return true;
  } else {
    return false;
  }
}

bool BinaryMatrix ::IsColTrivial(int c) {
  // check whether column c is trivial or not
  // a column is trivial if the column is all 0 or all 1
  bool hasZero = false;
  bool hasOne = false;
  for (int i = 0; i < GetRowNum(); ++i) {
    if (rowsArray[i][c] == 0) {
      hasZero = true;
    } else {
      hasOne = true;
    }
  }
  if (hasZero && hasOne) {
    return false;
  } else {
    return true;
  }
}

void BinaryMatrix ::GetTrivialSites(vector<int> &trivSites) {
  trivSites.clear();
  for (int c = 0; c < GetColNum(); ++c) {
    if (IsColTrivial(c) == true) {
      trivSites.push_back(c);
    }
  }
}

bool BinaryMatrix ::IsSequencesMatch(int r1, int r2, vector<int> &seqColPos) {
  bool res = true;
  // cout << "r1 = " << r1 << ", r2 = " << r2 << ", seqeucne lpocation are ";
  // DumpIntVec( seqColPos );

  // This function test whether (non-continuous) sequences for two rows match or
  // not
  for (unsigned int i = 0; i < seqColPos.size(); ++i) {
    if (rowsArray[r1][seqColPos[i]] != rowsArray[r2][seqColPos[i]]) {
      res = false;
      break;
    }
  }
  return res;
}

void BinaryMatrix ::GetSequencesDiffSites(int r1, int r2,
                                          set<int> &seqColDiffs) const {
  // colect the set of sites that the two rows are different
  seqColDiffs.clear();
  for (int c = 0; c < GetColNum(); ++c) {
    if (rowsArray[r1][c] != rowsArray[r2][c]) {
      seqColDiffs.insert(c);
    }
  }
}

bool BinaryMatrix ::IsZeroColumn(int c) {
  bool res = true;
  for (unsigned int i = 0; i < rowsArray.size(); ++i) {
    if (rowsArray[i][c] == 1) {
      res = false;
      break;
    }
  }
  return res;
}

int BinaryMatrix ::GetZeroColNum() {
  int res = 0;
  for (int i = 0; i < nCols; ++i) {
    if (IsZeroColumn(i)) {
      res++;
    }
  }
  return res;
}

void BinaryMatrix ::BuildColEquivClasses() {
  for (int i = 0; i < nCols; ++i) {
    bool f = false;
    for (COLUMN_EQUIV_CLASS::iterator it = setColEquiv.begin();
         it != setColEquiv.end(); ++it) {
      set<int> &s = *it;

      // check to see if column i/j are the same
      if (CmpColumns(i, *(s.begin())) == true) {
        // remember this fact in the map
        f = true;
        s.insert(i);
        break;
      }
    }

    if (f == false) {
      // Create a new set
      set<int> s1;
      s1.insert(i);
      setColEquiv.push_back(s1);
    }
  }
}

void BinaryMatrix ::GetUniqueColsInRange(int c1, int c2, set<int> &setUniques) {
  // make sure equiv class is pre-processed
  if (setColEquiv.empty()) {
    BuildColEquivClasses();
  }

  // exam each column equivlance classes
  // Put the mostly diesired (for now, it is the one near the center)
  // into result set (which must be in range)
  int center = (c1 + c2) / 2;
  for (unsigned int i = 0; i < setColEquiv.size(); ++i) {
    set<int> &s = setColEquiv[i];
    int cand = -100;
    for (set<int>::iterator it = s.begin(); it != s.end(); ++it) {
      int c = *it;
      if (c >= c1 && c <= c2 && abs(cand - center) > abs(c - center)) {
        cand = c;
      }
    }
    if (cand >= 0) {
      setUniques.insert(cand);
    }
  }
}

bool BinaryMatrix ::IsPerfectPhylogeny() {
  for (int i = 0; i < nCols - 1; ++i) {
    for (int j = i + 1; j < nCols; ++j) {
      if (IsCompatible(i, j) == false) {
        // cout << "Site i=" << i << ", j=" << j << " are incompatible.\n";
        return false;
      }
    }
  }
  return true;
}

void BinaryMatrix ::ConstructConflictGraph(UnWeightedGraph &graph) {
  // Conflict graph vertex num = # of columns
  // Edge is whether col i conflict with col j
  LIST_VERTEX vertList;
  LIST_EDGE edgeList;

  for (int i = 0; i < nCols; ++i) {
    char buf[100];
    buf[0] = 'c';
    sprintf(&buf[1], "%d", i + 1);
    BGVertex v(buf);
    vertList.push_back(v);
  }
  graph.SetVertices(vertList);

  // Now check for all pair of columns for conflict
  for (int i = 0; i < nCols - 1; ++i) {
    for (int j = i + 1; j < nCols; ++j) {
      if (IsCompatible(i, j) == false) {
        // cout << "Add one edge (" << i << " , " << j << ")" << endl;
        BGEdge eg("e", i, j, graph.GetListVerts());
        edgeList.push_back(eg);
      }
    }
  }

  // Finally, setup the vertex\edge lists
  graph.SetEdges(edgeList);
}

bool BinaryMatrix ::IsColumnBinary(int c) const {
  for (int i = 0; i < GetRowNum(); ++i) {
    if (rowsArray[i][c] != 0 && rowsArray[i][c] != 1) {
      return false;
    }
  }
  return true;
}

bool BinaryMatrix ::IsRowBinary(int r) const {
  for (int i = 0; i < nCols; ++i) {
    if (rowsArray[r][i] != 0 && rowsArray[r][i] != 1) {
      return false;
    }
  }
  return true;
}

void BinaryMatrix ::TrimNonBinaryRows() {
  set<int> setOfDuplicates;
  setOfDuplicates.clear();
  unsigned int r1;
  // int c;

  for (r1 = 0; r1 < rowsArray.size(); ++r1) {
    if (IsRowBinary(r1) == false) {
      // The row with duplicated rows are treated the same
      // cout << "row " << r2 << " is not binary." << endl;
      setOfDuplicates.insert(r1);
    }
  }
  /*
          Now we remove all duplicate rows
  */
  RemoveRows(setOfDuplicates);
}

bool BinaryMatrix ::IsRowRangeBinary(int r, int left, int right) {
  for (int i = left; i <= right; ++i) {
    if (rowsArray[r][i] == 2) {
      return false;
    }
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////
//		Inernal utility functions
////////////////////////////////////////////////////////////////////////////////

bool BinaryMatrix ::IsCompatible(int c1, int c2) {
  bool f00 = false;
  bool f01 = false;
  bool f10 = false;
  bool f11 = false;

  // if c1==c2, we assume it is compatible
  if (c1 == c2) {
    return true;
  }
#if 0 // no, acutally, we need to be more cautious, unless we see evidence, we
      // put it
      // For now, if a column is not binary, we consider it is not compatible
	if( IsColumnBinary(c1) == false || IsColumnBinary(c2) == false)
	{
		return false;
	}
#endif
  // 4-gamet test
  for (unsigned int i = 0; i < rowsArray.size(); ++i) {
    if (rowsArray[i][c1] == 0 && rowsArray[i][c2] == 0) {
      f00 = true;
    }
    if (rowsArray[i][c1] == 0 && rowsArray[i][c2] == 1) {
      f01 = true;
    }
    if (rowsArray[i][c1] == 1 && rowsArray[i][c2] == 0) {
      f10 = true;
    }
    if (rowsArray[i][c1] == 1 && rowsArray[i][c2] == 1) {
      f11 = true;
    }
  }

  // Now check to see if all flags are set
  if (f00 && f01 && f10 && f11)
    return false;
  else
    return true;
}

bool BinaryMatrix ::IsCompatibleRooted(int c1, int c2, int rallele1,
                                       int rallele2) {
  bool f00 = false;
  bool f01 = false;
  bool f10 = false;
  bool f11 = false;

  // if c1==c2, we assume it is compatible
  if (c1 == c2) {
    return true;
  }

  // 3-gamet test
  for (unsigned int i = 0; i < rowsArray.size(); ++i) {
    if (rowsArray[i][c1] == rallele1 && rowsArray[i][c2] == rallele2) {
      f00 = true;
    }
    if (rowsArray[i][c1] == rallele1 && rowsArray[i][c2] != rallele2) {
      f01 = true;
    }
    if (rowsArray[i][c1] != rallele1 && rowsArray[i][c2] == rallele2) {
      f10 = true;
    }
    if (rowsArray[i][c1] != rallele1 && rowsArray[i][c2] != rallele2) {
      f11 = true;
    }
  }

  // Now check to see if all flags are set
  if (f01 && f10 && f11)
    return false;
  else
    return true;
}

bool BinaryMatrix ::IsSiteCompatibleWithRegion(int s, int rc1, int rc2) {
  bool res = true;
  for (int rci = rc1; rci <= rc2; ++rci) {
    if (IsCompatible(s, rci) == false) {
      res = false;
      break;
    }
  }
  return res;
}

bool BinaryMatrix ::IsRegionFullyCompatible(int rc1, int rc2) {
  for (int rci = rc1; rci <= rc2; ++rci) {
    for (int rcj = rci + 1; rcj <= rc2; ++rcj) {
      if (IsCompatible(rci, rcj) == false) {
        return false;
      }
    }
  }
  return true;
}

void BinaryMatrix ::GetGamates(int c1, int c2, bool &f00, bool &f01, bool &f10,
                               bool &f11) {
  // init to all false upon start
  f00 = false;
  f01 = false;
  f10 = false;
  f11 = false;

  // if c1==c2, we assume it is compatible
  if (c1 == c2) {
    return;
  }

  // 4-gamet test
  for (unsigned int i = 0; i < rowsArray.size(); ++i) {
    if (rowsArray[i][c1] == 0 && rowsArray[i][c2] == 0) {
      f00 = true;
    }
    if (rowsArray[i][c1] == 0 && rowsArray[i][c2] == 1) {
      f01 = true;
    }
    if (rowsArray[i][c1] == 1 && rowsArray[i][c2] == 0) {
      f10 = true;
    }
    if (rowsArray[i][c1] == 1 && rowsArray[i][c2] == 1) {
      f11 = true;
    }
  }
}

bool BinaryMatrix ::IsColComplement(int c1, int c2) {
  for (unsigned int i = 0; i < rowsArray.size(); ++i) {
    // cout << "[i, c1] = " << rowsArray[i][c1] << ", rowsArray[i][c2] = " <<
    // rowsArray[i][c2] << endl;
    if (rowsArray[i][c1] == rowsArray[i][c2]) {
      return false;
    }
  }
  // cout << "col " << c1 << ", " << c2 << " are compl.\n";
  return true;
}
bool BinaryMatrix ::IsColDuplicate(int c1, int c2) {
  for (unsigned int i = 0; i < rowsArray.size(); ++i) {
    if (rowsArray[i][c1] != rowsArray[i][c2]) {
      return false;
    }
  }
  // cout << "col " << c1 << ", " << c2 << " are identical.\n";
  return true;
}

void BinaryMatrix ::GetAllIncompatiblePairs(
    set<pair<int, int> > &incompatibles) {
  incompatibles.clear();
  for (int i = 0; i < nCols; i++) {
    for (int j = i + 1; j < nCols; ++j) {
      // Test to see if site i, j are compatible
      if (IsCompatible(i, j) == false) {
        pair<int, int> p(i, j);
        incompatibles.insert(p);
      }
    }
  }
}

int BinaryMatrix ::ComputeHKBound() {
  // The idea is to test for incompatible between each column
  // Then create an incompatibility map, and compute the bound
  map<INTERVAL, int> bounds;

  int nCols = GetColNum();
  int nRows = GetRowNum();
  if (nCols <= 1 || nRows <= 3) {
    return 0;
  }

  for (int i = 0; i < nCols - 1; ++i) {
    for (int j = i + 1; j < nCols; ++j) {
      // Check if site i, j conflict
      int val = 0;
      if (IsCompatible(i, j) == false) {
        val = 1;
      }
      INTERVAL iv(i, j);
      bounds.insert(map<INTERVAL, int>::value_type(iv, val));
    }
  }
  vector<int> locBreakpoints; // do not really need this, but...
  return CalcCompositeBound(bounds, 0, nCols - 1, locBreakpoints);
}

int BinaryMatrix ::ComputeFastHapBound() {
  // Simply test for each submatrix for a rough haplotype bound
  // Then create an incompatibility map, and compute the bound
  // To speed things up, we do not perform optimal RecMin
  // Rather simply no-subset

  map<INTERVAL, int> bounds;

  int nc = GetColNum();
  int nr = GetRowNum();
  if (nc <= 1 || nr <= 3) {
    return 0;
  }

  for (int i = 0; i < nc - 1; ++i) {
    for (int j = i + 1; j < nc; ++j) {
      // Check if site i, j conflict
      int val = 0;

      BinaryMatrix submat;
      SubMatrix(0, GetRowNum() - 1, i, j, submat);
      submat.TrimFullyCompatibleSites();
      submat.TrimDupRows();

      val = submat.GetRowNum() - submat.GetColNum() - 1;
      if (val < 0) {
        val = 0;
      }

      INTERVAL iv(i, j);
      bounds.insert(map<INTERVAL, int>::value_type(iv, val));
      // cout << "interval " << i << ", " << j  << " quick bd = " << val <<
      // endl;
    }
  }
  vector<int> locBreakpoints; // do not really need this, but...
  return CalcCompositeBound(bounds, 0, nc - 1, locBreakpoints);
}

// This function computes a fast recombination upper bound, which can be useful
// in applications like branch and bound The idea is to remove a sequence from
// inputmat a time, and take the min to recombine them
int BinaryMatrix ::ComputeFastRecombUpperBound() {
  // Create a new sequence for operation
  BinaryMatrix matToOp = *this;

  int res = 0;
  // Whenver the matrix is too small, we stop
  while (true) {
    // First perform cleanup: drop non-informatives rows, collapse identical
    // rows
    set<int> setOfRemoved;
    matToOp.TrimFullyCompatibleSites(&setOfRemoved);
    matToOp.FindNgbrDupCompSites(&setOfRemoved);
    matToOp.RemoveColumns(setOfRemoved);
    matToOp.TrimDupRows();

    if (matToOp.GetRowNum() <= 3) {
      break;
    }

    // Find the smallest cost row
    int minRmCost = HAP_MAX_INT;
    int minRow = -1;
    // Try every leftover row in matToOp
    for (int r = 0; r < matToOp.GetRowNum(); ++r) {
      //            SEQUENCE row;
      //            matToOp.GetRow( r, row );
      int recCost = matToOp.ComputeMinRecombWeight(r);
      if (recCost < minRmCost) {
        minRmCost = recCost;
        minRow = r;
      }
    }
    YW_ASSERT_INFO(minRow >= 0, "Error: minRow must be updated at least once.");
    // cout << "minRmCost = " << minRmCost << ", minRow = " << minRow << endl;
    // Now we remove this sequence
    res += minRmCost;
    set<int> seqsToRemove;
    seqsToRemove.insert(minRow);
    matToOp.RemoveRows(seqsToRemove);
  }
  // cout << "A fast recomb. upper bound = " << res << endl;
  return res;
}

int BinaryMatrix ::ComputeMinRecombWeight(int rowIndex) {
  // This function computes a recombination number given the rows in matrix
  // that are ancesters of rowIndex
  // This function computes the minimum recombination weight for the given
  // hapRow when restricted to interval [left, right] in mat
  int res = 0;
  // cout << "ComputeMinRecombWeight :: rowIndex = " << rowIndex << endl;
  // cout <<"matrix here is: ";
  // Dump();
  set<int> lastTrackRows; // set of rows that matching the hapRow

  // Initially every row is a match
  for (int i = 0; i < GetRowNum(); ++i) {
    if (i != rowIndex) {
      lastTrackRows.insert(i);
    }
  }

  for (int curpos = 0; curpos < GetColNum(); ++curpos) {
    // Each time, we intersect the set with the sets matching the current bit
    set<int> trackRows;
    for (int i = 0; i < GetRowNum(); ++i) {
      if (i == rowIndex) {
        continue;
      }

      if (GetValAt(i, curpos) == GetValAt(rowIndex, curpos)) {
        // Yes, this row matches
        trackRows.insert(i);
      }
    }
    YW_ASSERT_INFO(trackRows.size() > 0, "trackRows must contain some rows.");

    // Now we test if there is intersection, if non-empty, we contiinue
    set<int> sint;
    JoinSets(trackRows, lastTrackRows, sint);
    if (sint.size() == 0) {
      // No intersection, so we have to increase the result (we know there must
      // be one recomb here, from the right-maximal proof)
      ++res;

      // Re-initialize lastTrackRows here
      lastTrackRows = trackRows;
      //            PopulateSetWithInterval( lastTrackRows, 0, mat.size() - 1 );
    } else {
      // In this case, we still continue
      lastTrackRows = sint;
    }
  }

  // cout << "Min recomb = " << res << endl;
  return res;
}

int BinaryMatrix ::GetMajorityState(int site) {
  int res = 0;
  for (int r = 0; r < GetRowNum(); ++r) {
    if (GetValAt(r, site) == 0) {
      res++;
    }
  }
  if (res >= (GetRowNum() + 1) / 2) {
    return 0;
  } else {
    return 1;
  }
}

int BinaryMatrix ::GetMinorStateNum(int site, int &minorState) const {
  int res = 0;
  for (int r = 0; r < GetRowNum(); ++r) {
    if (GetValAt(r, site) == 0) {
      res++;
    }
  }
  if (res >= (GetRowNum() + 1) / 2) {
    minorState = 1;
    return GetRowNum() - res;
  } else {
    minorState = 0;
    return res;
  }
}

void BinaryMatrix ::GetMinorStateRows(int site, int &minorState,
                                      set<int> &listRowsWMinor) const {
  GetMinorStateNum(site, minorState);
  for (int r = 0; r < GetRowNum(); ++r) {
    if (GetValAt(r, site) == minorState) {
      listRowsWMinor.insert(r);
    }
  }
}

void BinaryMatrix ::GetRowsWithAllele(int site, int alleleState,
                                      set<int> &setRows) const {
  //
  setRows.clear();
  for (int r = 0; r < GetRowNum(); ++r) {
    if (GetValAt(r, site) == alleleState) {
      setRows.insert(r);
    }
  }
}

int BinaryMatrix ::GetTheOtherAllele(int allele) {
  //
  if (allele == 0) {
    return 1;
  } else {
    return 0;
  }
}

void BinaryMatrix ::ConfigZeroMajSeq() {
  // make majority elem all-0 for each position
  //
  for (int c = 0; c < GetColNum(); ++c) {
    int mc = GetMajorityState(c);
    if (mc == 1) {
      // switch it
      for (int r = 0; r < GetRowNum(); ++r) {
        //
        if (GetValAt(r, c) == 0) {
          rowsArray[r][c] = 1;
        } else {
          rowsArray[r][c] = 0;
        }
      }
    }
  }
}

void BinaryMatrix ::ConfigZeroAncesSeq(const vector<int> &seqAnces) {
  // if seqAnces[i] = 1, then swap 0/1 in the matrix
  YW_ASSERT_INFO((int)seqAnces.size() == GetColNum(), "Size: mismatch2");
  for (int c = 0; c < GetColNum(); ++c) {
    int mc = seqAnces[c];
    if (mc == 1) {
      // switch it
      for (int r = 0; r < GetRowNum(); ++r) {
        //
        if (GetValAt(r, c) == 0) {
          rowsArray[r][c] = 1;
        } else {
          rowsArray[r][c] = 0;
        }
      }
    }
  }
}

void BinaryMatrix ::DumpConvGenotypes() {
  // for 00: 1
  YW_ASSERT_INFO((GetRowNum() % 2) == 0,
                 "To get genotypes, must have EVEN number of rows");

  cout << "Converted genotype: " << GetRowNum() / 2 << " by " << GetColNum()
       << " sites\n";

  for (int i = 0; i < GetRowNum(); i += 2) {
    for (int c = 0; c < GetColNum(); ++c) {
      if (GetValAt(i, c) == 0 && GetValAt(i + 1, c) == 0) {
        cout << "0";
      } else if (GetValAt(i, c) == 1 && GetValAt(i + 1, c) == 1) {
        cout << "1";
      } else {
        cout << "2";
      }
    }
    cout << endl;
  }
}

void BinaryMatrix ::GreedyRemoveIncompatSites(BinaryMatrix &matReduced) {
  // greedily remove incompatible sites (i.e. first remove site that is
  // incompatible w/ most sites and continue) approach: try to find some subset
  // of columns that fits the perfect phylogeny; and use that to estimate the
  // number of migrations hopefully this works reasonably well for low
  // reombinaiton rates
  vector<vector<bool> > listPairCompatibles;

  //
  listPairCompatibles.resize(this->GetColNum());
  for (int s1 = 0; s1 < this->GetColNum(); ++s1) {
    listPairCompatibles[s1].resize(this->GetColNum());
    for (int s2 = s1 + 1; s2 < this->GetColNum(); ++s2) {
      listPairCompatibles[s1][s2] = IsCompatible(s1, s2);
    }
  }
  // keep track of which sites are incompaiblw with which
  vector<set<int> > listIncompatSitesPerSite(this->GetColNum());
  for (int s1 = 0; s1 < this->GetColNum(); ++s1) {
    listPairCompatibles[s1].resize(this->GetColNum());
    for (int s2 = s1 + 1; s2 < this->GetColNum(); ++s2) {
      if (listPairCompatibles[s1][s2] == false) {
        //
        listIncompatSitesPerSite[s1].insert(s2);
        listIncompatSitesPerSite[s2].insert(s1);
      }
    }
  }
  // cout << "List of incompatible sites: \n";
  // for( int jj=0; jj<(int)listIncompatSitesPerSite.size(); ++jj )
  //{
  // cout << "site: " << jj << ": ";
  // DumpIntSet(listIncompatSitesPerSite[jj]);
  //}

  // remove the matrix sites by dropping the one w/ largest incompatible pairs
  // until all sites become compatible w/ each other
  set<int> setChosenRemoveSites;
  while (true) {
    // find the site w/ largest incompat sites
    vector<int> listIncSize;
    for (int ii = 0; ii < (int)listIncompatSitesPerSite.size(); ++ii) {
      listIncSize.push_back(listIncompatSitesPerSite[ii].size());
    }
    int sChosen = std::max_element(listIncSize.begin(), listIncSize.end()) -
                  listIncSize.begin();
    int siteChosen = sChosen;
    if (listIncSize[siteChosen] == 0) {
      // all remaining sites are compatible. Stop
      break;
    }
    // cout << "List of inompat size: ";
    // DumpIntVec(listIncSize);
    // cout << "Choosen site: " << siteChosen << endl;

    // add this site; then remove this site from each incomp site list
    setChosenRemoveSites.insert(siteChosen);
    listIncompatSitesPerSite[siteChosen].clear();
    for (int jj = 0; jj < (int)listIncompatSitesPerSite.size(); ++jj) {
      listIncompatSitesPerSite[jj].erase(siteChosen);
    }
  }
  // cout << "List of sites to remove: ";
  // DumpIntSet(setChosenRemoveSites);
  //
  vector<int> listKeptSites;
  for (int s1 = 0; s1 < (int)this->GetColNum(); ++s1) {
    //
    if (setChosenRemoveSites.find(s1) == setChosenRemoveSites.end()) {
      listKeptSites.push_back(s1);
    }
  }
  YW_ASSERT_INFO(listKeptSites.size() > 0, "ListKeptSites: wrong");
  SubMatrixSelectedSites(listKeptSites, matReduced);
  // cout << "GreedyRemoveIncompatSites: original mat = ";
  // this->Dump();
  // cout << "After removing incompatible sites greedyly, matrix = ";
  // matReduced.Dump();
}

void BinaryMatrix ::CalcSFS(vector<double> &listSFSFrac) const {
  // compute SFS; that is, list[i] = frac of sites with minor allele (assumed to
  // be 1) appears i times note: assume 0 is ancestral
  listSFSFrac.clear();
  int numRows = GetRowNum();
  for (int r = 0; r <= numRows; ++r) {
    listSFSFrac.push_back(0.0);
  }
  int numCols = GetColNum();
  for (int s = 0; s < GetColNum(); ++s) {
    //
    int minorState;
    int numTimes = GetMinorStateNum(s, minorState);
    if (minorState == 0) {
      //
      numTimes = numRows - numTimes;
    }
    YW_ASSERT_INFO(numTimes >= 0 && numTimes <= numRows, "Wrong");
    listSFSFrac[numTimes] += 1.0 / numCols;
  }
}

int BinaryMatrix ::GetDiffSitesForTwoRows(int r1, int r2) const {
  //
  int res = 0;
  for (int c = 0; c < GetColNum(); ++c) {
    if (GetValAt(r1, c) != GetValAt(r2, c)) {
      ++res;
    }
  }
  return res;
}

double BinaryMatrix ::CalcAvePairRowsDiff() const {
  // average pairwise diff (normalized by row length)
#if 0
    double res = 0.0;
    int numPairs = 0;
    for(int r1=0; r1<GetRowNum(); ++r1)
    {
        for(int r2=r1+1; r2<GetRowNum(); ++r2)
        {
            ++numPairs;
            res += GetDiffSitesForTwoRows(r1,r2);
        }
    }
    return res/( GetColNum()* numPairs);
#endif

  // use a faster approach
  // first accumlate the num of 1s in the first i rows at each site
  vector<vector<int> > vecNum1sAtSites(GetRowNum());
  for (int r = 0; r < (int)vecNum1sAtSites.size(); ++r) {
    // accumlate for each col
    for (int c = 0; c < GetColNum(); ++c) {
      int num1s = 0;
      if (r > 0) {
        num1s = vecNum1sAtSites[r - 1][c];
      }
      if (GetValAt(r, c) == 1) {
        ++num1s;
      }
      vecNum1sAtSites[r].push_back(num1s);
    }
  }
  // now accumate diffs
  double totDiffs = 0.0;
  for (int r = 1; r < GetRowNum(); ++r) {
    // calc tot diffs here
    for (int c = 0; c < GetColNum(); ++c) {
      int stepVal = 0;
      if (GetValAt(r, c) == 0) {
        stepVal = vecNum1sAtSites[r - 1][c];
      } else {
        stepVal = r - vecNum1sAtSites[r - 1][c];
      }
      YW_ASSERT_INFO(stepVal >= 0, "Cannot be negative");
      totDiffs += stepVal;
    }
  }
  int numPairs = GetRowNum() * (GetRowNum() - 1) / 2;
  return totDiffs / (GetColNum() * numPairs);
}

double BinaryMatrix ::CalcAvePairRowsDiffBetween(const set<int> &rowsSet1,
                                                 const set<int> &rowsSet2,
                                                 double &resMinDiffOut) const {
  //
  double res = 0.0;
  int numPairs = 0;
  double resMaxDiff = GetColNum();
  for (set<int>::iterator it1 = rowsSet1.begin(); it1 != rowsSet1.end();
       ++it1) {
    int r1 = *it1;
    for (set<int>::iterator it2 = rowsSet2.begin(); it2 != rowsSet2.end();
         ++it2) {
      int r2 = *it2;
      ++numPairs;
      int valdiff = GetDiffSitesForTwoRows(r1, r2);
      res += valdiff;
      if (resMaxDiff > valdiff) {
        resMaxDiff = valdiff;
      }
    }
  }
  resMinDiffOut = resMaxDiff / GetColNum();
  return res / (GetColNum() * numPairs);
}

void BinaryMatrix ::CollectAllPairwiseDiffs(
    const set<int> &rowsSet1, const set<int> &rowsSet2,
    vector<double> &listRowPairsDiff) const {
  //
  listRowPairsDiff.clear();
  for (set<int>::iterator it1 = rowsSet1.begin(); it1 != rowsSet1.end();
       ++it1) {
    int r1 = *it1;
    for (set<int>::iterator it2 = rowsSet2.begin(); it2 != rowsSet2.end();
         ++it2) {
      int r2 = *it2;
      int valdiff = GetDiffSitesForTwoRows(r1, r2);
      listRowPairsDiff.push_back(((double)valdiff) / GetColNum());
    }
  }
  // sort results
  SortDoubleVec(listRowPairsDiff);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void GetNoninformativeRowsInMat(const BinaryMatrix &mat, set<int> &trimedRows,
                                vector<REMOVED_ROWS_INFO> &trimedRowInfo,
                                set<int> &trimedCols, BinaryMatrix &matUpdated,
                                bool fRmDup) {
  //
  BinaryMatrix matUse = mat;

  //
  trimedRows.clear();
  trimedRowInfo.clear();

  // we perform trimming dup rows and then noninformative columns, repeatively.
  // Stop when there is no more work to do
  vector<int> curRowsRemoved;
  vector<int> curColsRemoved;
  while (true) {
    // cout << "cur mat = ";
    // matUse.Dump();
    // simply check to see if anything can be done
    // first remove non-inform sites
    //

    set<int> removedCols;
    matUse.FindNonInformativeSites(removedCols);
    // cout << "Removed cols = ";
    // DumpIntSet( removedCols );

    // now also remove dup sites
    if (fRmDup == true) {
      set<int> sitesDupRm;
      matUse.FindNgbrDupCompSites(&sitesDupRm);
      // cout << "Dup sites removed: ";
      // DumpIntSet( sitesDupRm );
      // also remember sites being trimmed
      UnionSets(removedCols, sitesDupRm);
    }

    // now removed stuff
    matUse.RemoveColumns(removedCols);

    if (removedCols.size() > 0) {
      // do the same for cols
      vector<int> remColsVec;
      PopulateVecBySet(remColsVec, removedCols);
      vector<int> posOrigCol;
      RecoverOrigIndicesAfterDeletion(curColsRemoved, remColsVec, posOrigCol);
      AppendIntVec(curColsRemoved, posOrigCol);
    }

    // now see if anything becomes identical
    set<int> removedRows;
    vector<pair<int, int> > listRowRemInfo;
    matUse.TrimDupRows(&removedRows, &listRowRemInfo);
    // cout << "Trimmed rows = ";
    // DumpIntSet(removedRows);
    // for(int jjj=0; jjj<(int)listRowRemInfo.size(); ++jjj)
    //{
    // cout << "Deleting row " << listRowRemInfo[jjj].first << ", since exists a
    // duplicate " << listRowRemInfo[jjj].second << endl;
    //}
    // remember which rows are rmeoved and which row it gets its value from
    if (removedRows.size() > 0) {
      REMOVED_ROWS_INFO rri;
      rri.rowsRemoved = removedRows;
      rri.pairsRmKeepRows = listRowRemInfo;
      trimedRowInfo.push_back(rri);
    }

    // stop if found nothing
    if (removedRows.size() == 0) {
      break;
    }
    // cout << "Removed these rows: ";
    // DumpIntSet(removedRows);
    // save it
    vector<int> remRowsVec;
    PopulateVecBySet(remRowsVec, removedRows);

    vector<int> posOrig;
    RecoverOrigIndicesAfterDeletion(curRowsRemoved, remRowsVec, posOrig);

    // append finally
    AppendIntVec(curRowsRemoved, posOrig);
  }
  // cout << "Finally, removed rows are: ";
  // DumpIntVec(curRowsRemoved);
  // cout << "Finally, removed cols are: ";
  // DumpIntVec(curColsRemoved);
  // after trimming redundent rows
  // cout << "After trimming, matrix rows = ";
  // matUse.Dump();

  // conver to set
  PopulateSetByVec(trimedRows, curRowsRemoved);

  // also other output
  matUpdated = matUse;
  PopulateSetByVec(trimedCols, curColsRemoved);
}

void SplitMatrixIntoMaximalFullyCompatRegs(
    const BinaryMatrix &mat, vector<pair<int, int> > &listFullyCompatRegs) {
  BinaryMatrix &matInst = const_cast<BinaryMatrix &>(mat);

  // divide a (potentially very large) matrix into maximal fully compatible
  // regions
  int posLeft = 0;
  int posCur = posLeft + 1;
  while (posCur < mat.GetColNum()) {
    // check compaibility for all previous ones
    bool fFullyCompat = true;
    for (int c = posLeft; c < posCur; ++c) {
      if (matInst.IsCompatible(c, posCur) == false) {
        fFullyCompat = false;
        break;
      }
    }
    if (fFullyCompat == false) {
      pair<int, int> pp(posLeft, posCur - 1);
      listFullyCompatRegs.push_back(pp);
      posLeft = posCur;
    }
    ++posCur;
  }
  // add last segment if remain a lst one
  pair<int, int> pp(posLeft, mat.GetColNum() - 1);
  listFullyCompatRegs.push_back(pp);
}

void ReadSitePosFromFirstRowInFile(const char *filename, int numSites,
                                   vector<double> &listSitePos) {
  //
  ifstream inFile(filename);
  if (!inFile) {
    YW_ASSERT_INFO(false, "Fatal error: cannot open the file");
  }
  string whitespace = " ";
  int MAX_NUM_SITES = 102400;
  const int BUF_SZ = MAX_NUM_SITES * sizeof(int);
  char buf[BUF_SZ];
  inFile.getline(buf, BUF_SZ);
  string strbuf(buf);
  size_t strEnd = strbuf.find_last_not_of(whitespace);
  strbuf = strbuf.substr(0, strEnd);
  std::istringstream is(strbuf);
  listSitePos.clear();
  while (is.eof() == false) {
    double pos;
    is >> pos;
    listSitePos.push_back(pos);
  }
  // cout << "numSites: " << numSites << endl;
  // cout << "ListSitePos: " << listSitePos.size() << " ";
  // DumpDoubleVec(listSitePos);
  YW_ASSERT_INFO((int)listSitePos.size() == numSites, "Wrong");
}
