#include "GenotypeMatrix.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>

// ***************************************************************************
// Define a reusable binary matrix class
// ***************************************************************************

GenotypeMatrix ::GenotypeMatrix() { nCols = 0; }

GenotypeMatrix ::~GenotypeMatrix() {
  // Need to free up data if needed
  Clear();
}

GenotypeMatrix ::GenotypeMatrix(int nr, int nc) { SetSize(nr, nc); }

GenotypeMatrix ::GenotypeMatrix(const GenotypeMatrix &rhs) { Copy(rhs); }

GenotypeMatrix &GenotypeMatrix ::operator=(const GenotypeMatrix &rhs) {
  Clear();

  Copy(rhs);

  return *this;
}

bool GenotypeMatrix ::IsDataValid(int val) {
  if (val == 0 || val == 1 || val == 2) {
    return true;
  } else {
    return false;
  }
}

void GenotypeMatrix ::PreSolve() {
  // Generate the companion rows
  SetupCompanionColumns();
}

bool GenotypeMatrix ::AreColumnsCompanion(int c1, int c2) {
  if (c1 == c2) {
    return false;
  }
  if (c1 > c2) {
    int tmp = c1;
    c1 = c2;
    c2 = tmp;
  }
  COLUMN_PAIR cp(c1, c2);
  if (companionRows.find(cp) == companionRows.end() ||
      companionRows[cp].size() == 0) {
    return false;
  } else {
    return true;
  }
}

bool GenotypeMatrix ::AreColumnsForcedInPhase(int c1, int c2) {
  if (c1 == c2) {
    return false;
  }
  if (c1 > c2) {
    int tmp = c1;
    c1 = c2;
    c2 = tmp;
  }
  COLUMN_PAIR cp(c1, c2);
  if (forcedColumnPairs.find(cp) == forcedColumnPairs.end() ||
      forcedColumnPairs[cp] == 1) {
    return false;
  } else {
    return true;
  }
}

bool GenotypeMatrix ::AreColumnsForcedOutPhase(int c1, int c2) {
  if (c1 == c2) {
    return false;
  }
  if (c1 > c2) {
    int tmp = c1;
    c1 = c2;
    c2 = tmp;
  }
  COLUMN_PAIR cp(c1, c2);
  if (forcedColumnPairs.find(cp) == forcedColumnPairs.end() ||
      forcedColumnPairs[cp] == 0) {
    return false;
  } else {
    return true;
  }
}

bool GenotypeMatrix ::AreColumnsComplete(int c1, int c2) {
  if (c1 == c2) {
    return false;
  }
  if (c1 > c2) {
    int tmp = c1;
    c1 = c2;
    c2 = tmp;
  }
  COLUMN_PAIR cp(c1, c2);
  for (int i = 0; i < completePairs.size(); ++i) {
    if (completePairs[i] == cp) {
      return true;
    }
  }
  return false;
}

int GenotypeMatrix ::GetNumTwosInRow(int r) {
  // For now, it is not optimized yet
  // we simply count the number of twos
  // later we can rely on preprocessing
  int res = 0;
  for (int i = 0; i < GetColNum(); ++i) {
    if (rowsArray[r][i] == 2) {
      ++res;
    }
  }
  return res;
}

bool GenotypeMatrix ::IsSiteTrival(int site) {
  int numTwos = 0;
  int numZeros = 0;
  int numOnes = 0;
  for (int i = 0; i < GetRowNum(); ++i) {
    if (rowsArray[i][site] == 0) {
      numZeros++;
    } else if (rowsArray[i][site] == 1) {
      numOnes++;
    } else if (rowsArray[i][site] == 2) {
      numTwos++;
    } else {
      YW_ASSERT(false);
    }
  }
  if (numTwos <= 1 && (numZeros == 0 || numOnes == 0)) {
    return true;
  } else {
    return false;
  }
}

bool GenotypeMatrix ::IsColComplement(int c1, int c2) {
  YW_ASSERT_INFO(false, "Not implemented");
  return false;
}

bool GenotypeMatrix ::IsColDuplicate(int c1, int c2) {
  for (int i = 0; i < GetRowNum(); ++i) {
    if (rowsArray[i][c1] != rowsArray[i][c2]) {
      return false;
    }
  }
  return true;
}

int GenotypeMatrix ::GetMajorityState(int site) {
  int numTwos = 0;
  int numZeros = 0;
  int numOnes = 0;
  for (int i = 0; i < GetRowNum(); ++i) {
    if (rowsArray[i][site] == 0) {
      numZeros++;
    } else if (rowsArray[i][site] == 1) {
      numOnes++;
    } else if (rowsArray[i][site] == 2) {
      numTwos++;
    } else {
      YW_ASSERT(false);
    }
  }
  int ma = 0;
  int max = numZeros;
  if (max < numOnes) {
    max = numOnes;
    ma = 1;
  }
  if (max < numTwos) {
    max = numTwos;
    ma = 2;
  }
  return ma;
}

/////////////////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION DETAILS
void GenotypeMatrix ::SetupCompanionColumns() {
  // This function checks the data and fill in the companion rows for every pair
  // of columns
  for (int i = 0; i < nCols; ++i) {
    for (int j = i + 1; j < nCols; ++j) {
      COLUMN_PAIR cp(i, j);
      set<int> cmpnRows;
      // The following 4 variables shows what are know already for (i,j)
      bool found00 = false, found01 = false;
      bool found10 = false, found11 = false;
      ;
      for (int k = 0; k < GetRowNum(); ++k) {
        if (rowsArray[k][i] == 2 && rowsArray[k][j] == 2) {
          cmpnRows.insert(k);
        } else if (rowsArray[k][i] == 0 &&
                   rowsArray[k][j] == 2) // Now check for forced pattern
        {
          found00 = true;
          found01 = true;
        } else if (rowsArray[k][i] == 1 && rowsArray[k][j] == 2) {
          found10 = true;
          found11 = true;
        } else if (rowsArray[k][i] == 2 && rowsArray[k][j] == 0) {
          found00 = true;
          found10 = true;
        } else if (rowsArray[k][i] == 2 && rowsArray[k][j] == 1) {
          found01 = true;
          found11 = true;
        } else if (rowsArray[k][i] == 0 && rowsArray[k][j] == 0) {
          found00 = true;
        } else if (rowsArray[k][i] == 0 && rowsArray[k][j] == 1) {
          found01 = true;
        } else if (rowsArray[k][i] == 1 && rowsArray[k][j] == 0) {
          found10 = true;
        } else if (rowsArray[k][i] == 1 && rowsArray[k][j] == 1) {
          found11 = true;
        }
      }
      // Now we add this to our map
      if (cmpnRows.size() > 0) {
        companionRows.insert(COMPANION_ROW_MAP::value_type(cp, cmpnRows));
      }

      // We also record the forced pattern
      if (found00 == true && found11 == true && found01 == true &&
          found10 == true) {
        // In this case we already have a complete pair
        completePairs.push_back(cp);
      } else {
        // This is not a complete pair
        if (found00 == true && found11 == true) {
          forcedColumnPairs.insert(
              FORCED_COL_MAP::value_type(cp, 0)); // forced in phase
        } else if (found01 == true && found10 == true) {
          forcedColumnPairs.insert(
              FORCED_COL_MAP::value_type(cp, 1)); // forced in phase
        }
      }
    }
  }
}
