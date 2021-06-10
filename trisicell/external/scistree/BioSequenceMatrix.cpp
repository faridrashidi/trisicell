#include "BioSequenceMatrix.h"
#include "Utils2.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>

// ***************************************************************************
// Define a reusable binary matrix class
// ***************************************************************************
#if 0
BioSequenceMatrix :: BioSequenceMatrix()
{
    fMissingValue = false;
}
#endif

BioSequenceMatrix ::~BioSequenceMatrix() { Clear(); }

void BioSequenceMatrix ::AppendRow(const vector<int> &row) {
  // Check to see if this is the first row, if so OK
  if (nCols == 0 && GetRowNum() == 0) {
    nCols = row.size();
  }

  if (row.size() != (unsigned int)nCols) {
    DEBUG("WRONG row width in AddRow");
    return;
  }

  int *buf = new int[nCols];
  for (int i = 0; i < nCols; ++i) {
    buf[i] = row[i];
  }
  rowsArray.push_back(buf);
}

void BioSequenceMatrix ::AppendSetOfRows(const set<SEQUENCE> &rows) {
  for (set<SEQUENCE>::iterator it = rows.begin(); it != rows.end(); ++it) {
    AppendRow(*it);
  }
}

void BioSequenceMatrix ::AppendRows(const vector<SEQUENCE> &rows) {
  for (unsigned int i = 0; i < rows.size(); ++i) {
    AppendRow(rows[i]);
  }
}

// This function removes a set of columns that are specified in the set as
// duplicateSites
void BioSequenceMatrix ::InsertColumns(const vector<SEQUENCE> &sitesValue,
                                       const vector<int> &sitesPos) {
  // we require the site contains the same number of values as rows
  YW_ASSERT_INFO(sitesPos.size() == (unsigned int)sitesValue.size(),
                 "Wrong vector size.");
  YW_ASSERT_INFO(sitesValue.size() > 0, "Can not be empty.");
  YW_ASSERT_INFO(sitesValue[0].size() == (unsigned int)GetRowNum(),
                 "Size mismatch.");

  int totalLen = GetColNum() + sitesPos.size();

  // First we need to calculate where to put these new sites
  // remember the passed-in values are BEFORE insertion. For example, when we
  // say we want to insert sites a,b at location 0, 2, when mean we want to put
  // a at 0 (this pushes the value forward
  vector<int> poses;
  int offset = 0;
  for (unsigned int i = 0; i < sitesPos.size(); ++i) {
    // Treat out of ranges as close to real
    int realPos = sitesPos[i];
    if (realPos < 0) {
      realPos = 0;
    } else if (realPos > GetColNum()) {
      realPos = GetColNum();
    }
    poses.push_back(realPos + offset);
    offset++;
  }
  // cout << "poses = ";
  // DumpIntVec( poses );
  // now we create  a new matrix with different size
  for (unsigned int i = 0; i < rowsArray.size(); ++i) {
    int *buf = new int[totalLen];
    int origPos = 0;
    int pos = 0;
    for (unsigned int j = 0; j < poses.size(); ++j) {
      // cout << "1. origPos = " << origPos << ", pos = " << pos << endl;
      for (; pos < poses[j]; ++pos) {
        buf[pos] = rowsArray[i][origPos++];
      }
      // Now add poses[j]
      // cout << "Before assign site,  origPos = " << origPos << ", pos = " <<
      // pos << endl;
      buf[pos++] = sitesValue[j][i];
    }
    // We finish any leftover
    for (; pos < totalLen; ++pos, ++origPos) {
      // cout << "2. origPos = " << origPos << ", pos = " << pos << endl;
      buf[pos] = rowsArray[i][origPos];
    }

    // now we free the memory of old buffer
    delete[] rowsArray[i];
    rowsArray[i] = buf;
  }

  nCols = totalLen;
}

void BioSequenceMatrix ::AppendMatrixByCol(
    const BioSequenceMatrix &appendedMat) {
  // Append the matrix by putting the matrix to the right
  // Make sure the row matches
  YW_ASSERT_INFO(appendedMat.IsEmpty() == false,
                 "For now, do not allow appending empty matrix.");
  YW_ASSERT_INFO(IsEmpty() || GetRowNum() == appendedMat.GetRowNum(),
                 "Can not append such matrix");

  // Figure out the size
  vector<int *> rowsArrayNew; // array of rows
  int rowNum, colNum;
  if (IsEmpty() == false) {
    rowNum = GetRowNum();
    colNum = GetColNum();
  } else {
    // Use the new matrix's value
    rowNum = appendedMat.GetRowNum();
    colNum = 0;
  }
  int numSitesNew = colNum + appendedMat.GetColNum();
  // Allocate space
  for (int r = 0; r < rowNum; ++r) {
    int *buf = new int[numSitesNew];
    rowsArrayNew.push_back(buf);
  }
  // Now copy the stuff in
  for (int r = 0; r < rowNum; ++r) {
    for (int c = 0; c < colNum; ++c) {
      rowsArrayNew[r][c] = rowsArray[r][c];
    }
    for (int c = 0; c < appendedMat.GetColNum(); ++c) {
      rowsArrayNew[r][c + colNum] = appendedMat(r, c);
    }
  }

  // Remove the old ones
  Clear();
  // Set to the new one
  nCols = numSitesNew;
  rowsArray = rowsArrayNew;
}

void BioSequenceMatrix ::AppendMatrixByRow(
    const BioSequenceMatrix &appendedMat) {
  // Append the matrix by putting the matrix to the right
  // Make sure the row matches
  YW_ASSERT_INFO(appendedMat.IsEmpty() == false,
                 "For now, do not allow appending empty matrix.");
  YW_ASSERT_INFO(IsEmpty() || GetColNum() == appendedMat.GetColNum(),
                 "Can not append such matrix");

  // Now copy the stuff in
  for (int r = 0; r < appendedMat.GetRowNum(); ++r) {
    SEQUENCE seq;
    appendedMat.GetRow(r, seq);
    this->AppendRow(seq);
  }
}

void BioSequenceMatrix ::SetRow(int r, const vector<int> &valNew) {
  if (valNew.size() != (unsigned int)nCols) {
    DEBUG("WRONG row width in SetRow");
    return;
  }
  for (int i = 0; i < nCols; ++i) {
    rowsArray[r][i] = valNew[i];
  }
}

void BioSequenceMatrix ::SetCol(int c, const vector<int> &valNew) {
  if (valNew.size() != (unsigned int)GetRowNum()) {
    DEBUG("WRONG row width in SetRow");
    return;
  }
  for (int i = 0; i < GetRowNum(); ++i) {
    rowsArray[i][c] = valNew[i];
  }
}

void BioSequenceMatrix ::Clear() {
  // Need to free up data if needed
  for (unsigned int i = 0; i < rowsArray.size(); ++i) {
    delete[] rowsArray[i];
  }
  rowsArray.clear();
  nCols = 0;
}

void BioSequenceMatrix ::Copy(const BioSequenceMatrix &rhs) {
  Clear(); // 012713: it seems we should clear this one first
  for (unsigned int i = 0; i < rhs.rowsArray.size(); ++i) {
    int *buf = new int[rhs.nCols];
    for (int j = 0; j < rhs.nCols; ++j) {
      buf[j] = rhs.rowsArray[i][j];
    }
    rowsArray.push_back(buf);
  }
  nCols = rhs.nCols;
}

void BioSequenceMatrix ::RemoveRow(int rowIndex) {
  if ((unsigned int)rowIndex >= rowsArray.size()) {
    return;
  }

  int nPos = -1;
  for (vector<int *>::iterator it = rowsArray.begin(); it != rowsArray.end();
       ++it) {
    nPos++;
    if (nPos == rowIndex) {
      delete[] * it;
      rowsArray.erase(it);
      return;
    }
  }
  DEBUG("Something very wrong inside BioSequenceMatrix :: RemoveRow");
}

// Consolidate rows in matrix
void BioSequenceMatrix::TrimDupRows(set<int> *pTrimedRows,
                                    vector<pair<int, int> > *pTrimRowInfo) {

  set<int> setOfDuplicates;
  vector<pair<int, int> > listRowsDeletedWithExistingPairs;
  setOfDuplicates.clear();
  unsigned int r1, r2;
  int c;

  bool res = false; // we stop unless we find some duplicate rows and/or
                    // non-informat site

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
        if (setOfDuplicates.find(r2) == setOfDuplicates.end()) {
          pair<int, int> pp;
          pp.first = r2; // first item is which row is removed
          pp.second =
              r1; // second item is which row is the source (to be kepted)
          listRowsDeletedWithExistingPairs.push_back(pp);
        }

        // cout << "row " << r2 << " is duplicate." << endl;
        setOfDuplicates.insert(r2);
      }
    }
  }
  /*
          Now we remove all duplicate rows
  */
  if (setOfDuplicates.size() > 0) {
    res = true;
    RemoveRows(setOfDuplicates);
  }
  if (pTrimedRows != NULL) {
    *pTrimedRows = setOfDuplicates;
  }
  if (pTrimRowInfo != NULL) {
    *pTrimRowInfo = listRowsDeletedWithExistingPairs;
  }

  return;
}

void BioSequenceMatrix ::DumpRowMultiplicity() const {
  // This function dump out duplicate row information
  map<SEQUENCE, int> mapRowMultiplicity;
  for (int r = 0; r < GetRowNum(); ++r) {
    SEQUENCE row;
    GetRow(r, row);
    if (mapRowMultiplicity.find(row) == mapRowMultiplicity.end()) {
      mapRowMultiplicity.insert(map<SEQUENCE, int>::value_type(row, 1));
    } else {
      mapRowMultiplicity[row]++;
    }
  }
  // Now dump out info
  cout << "In this matrix, the multiplicity of rows is: \n";
  for (map<SEQUENCE, int>::iterator it = mapRowMultiplicity.begin();
       it != mapRowMultiplicity.end(); ++it) {
    cout << "seq = ";
    DumpSequence(it->first);
    cout << ", multiplicty = ";
    cout << it->second << endl;
  }
}

void BioSequenceMatrix ::GetColMultiplicityMap(
    vector<int> &listColMulti) const {
  // for each col (site), find out the number of duplicate each site has (that
  // is, listMulti[i] = # of sites with the same column)
  listColMulti.clear();
  listColMulti.resize(GetColNum());
  map<SEQUENCE, set<int> > mapColMulti;
  for (int c = 0; c < GetColNum(); ++c) {
    SEQUENCE col;
    GetCol(c, col);
    mapColMulti[col].insert(c);
  }
  for (map<SEQUENCE, set<int> >::iterator it = mapColMulti.begin();
       it != mapColMulti.end(); ++it) {
    for (set<int>::iterator it2 = it->second.begin(); it2 != it->second.end();
         ++it2) {
      listColMulti[*it2] = it->second.size();
    }
  }
}

int BioSequenceMatrix ::GetMultiplictyForRow(int r) const {
  SEQUENCE seqRow;
  GetRow(r, seqRow);
  return GetMultiplictyForRow(seqRow);
}

int BioSequenceMatrix ::GetMultiplictyForRow(const SEQUENCE &seqRow) const {
  int res = 0;
  for (int i = 0; i < GetRowNum(); ++i) {
    SEQUENCE curRow;
    GetRow(i, curRow);
    if (curRow == seqRow) {
      ++res;
    }
  }
  return res;
}

int BioSequenceMatrix ::GetMultiplictyForRow(const SEQUENCE &seqRow,
                                             set<int> &identRows) const {
  identRows.clear();
  int res = 0;
  for (int i = 0; i < GetRowNum(); ++i) {
    SEQUENCE curRow;
    GetRow(i, curRow);
    if (curRow == seqRow) {
      identRows.insert(i);
      ++res;
    }
  }
  // YW_ASSERT_INFO( res > 0, "Must appear at least once." );
  return res;
}

int BioSequenceMatrix ::GetMultiplictyForRowIV(int r, int left,
                                               int right) const {
  SEQUENCE row;
  GetRow(r, row);
  SEQUENCE rowIV;
  GetSeqInterval(row, rowIV, left, right);
  int res = 0;
  for (int i = 0; i < GetRowNum(); ++i) {
    SEQUENCE curRow;
    GetRow(i, curRow);
    SEQUENCE rowIV1;
    GetSeqInterval(curRow, rowIV1, left, right);

    if (rowIV1 == rowIV) {
      ++res;
    }
  }
  return res;
}

bool BioSequenceMatrix ::ReadFromFile(ifstream &inFile, bool fSkipFirstLine) {
  bool res = true;

  // Now, we first check one row to find out how many sites
  // first read in the matrix name first
  const int BUF_SZ = MAX_SITE_NUM * sizeof(int);
  char buf[BUF_SZ]; // assume maximum sites allowed are 4096
  if (fSkipFirstLine == true) {
    inFile.getline(buf, BUF_SZ);
    //	cout << "Matrix name is " << buf << endl;
  }

  int rowLength = 0;
  while (!inFile.eof()) {
    inFile.getline(buf, BUF_SZ);
    DEBUG("strlen of buf ");
    DEBUG(strlen(buf));
    DEBUG("\n");
    DEBUG("buffer is:");
    DEBUG(buf);
    DEBUG("\n");

    // ignore any ine starting with #
    if (buf[0] == '#') {
      continue;
    }

    int curRowLen = 0;
    curRowLen = strlen(buf);
    // but we need to check to make sure there is no garbage character at the
    // end
    for (int i = curRowLen - 1;;) {
      if (i > 0 && buf[i] != '0' && buf[i] != '1' && buf[i] != '2' &&
          buf[i] != '*' && buf[i] != '?') {
        i--;
        curRowLen--;
      } else {
        break;
      }
    }

    if (rowLength == 0) {
      rowLength = curRowLen;
    }
    if (rowLength != curRowLen) {
      // for some reason, we are getting a smaller size
      // simplely terminate here
      // DEBUG("Warning: one row of fle seems to have fewer data bits.\n");
      // res = false;
      break;
    }
    int *pRow = new int[rowLength];
    for (int i = 0; i < rowLength; ++i) {
      if (buf[i] == '1') {
        pRow[i] = 1;
      } else if (buf[i] == '0') {
        pRow[i] = 0;
      } else if (buf[i] == '2') {
        pRow[i] = 2;
      } else if (buf[i] == '*' || buf[i] == '?') {
        pRow[i] = MISSING_VALUE_BIT;
      } else {
        YW_ASSERT_INFO(false, "Un-recognized characters in input.");
        exit(1);
      }
    }
    // Now put it into a list
    rowsArray.push_back(pRow);
  }

  // Now set return value
  nCols = rowLength;

  return res;
}

#if 0
bool BioSequenceMatrix :: ReadFromFilePartial( ifstream &inFile, bool fSkipFirstLine )
{
    // different from the previous version, read from the middle of the file and stop whenever needed
    bool res = true;

	// Now, we first check one row to find out how many sites
	char buf[MAX_SITE_NUM];		// assume maximum sites allowed are 4096

	// first read in the matrix name first
	const int BUF_SZ = MAX_SITE_NUM*sizeof(int);
    if( fSkipFirstLine == true )
    {
        inFile.getline(buf, BUF_SZ);
        //	cout << "Matrix name is " << buf << endl;
    }

	int rowLength = 0;
	while( !inFile.eof() )
	{
		inFile.getline(buf, BUF_SZ);
		DEBUG("strlen of buf ");
		DEBUG(strlen(buf));
		DEBUG("\n");
		DEBUG("buffer is:");
		DEBUG(buf);
		DEBUG("\n");
		int curRowLen = 0 ;
		curRowLen = strlen(buf);
		// but we need to check to make sure there is no garbage character at the end
		for(int i=curRowLen-1;;)
		{
			if(i >0 &&   buf[i]  != '0' && buf[i] != '1' && buf[i] != '2'
               && buf[i] != '*' && buf[i] != '?' )
			{
				i--;
				curRowLen --;
			}
			else
			{
				break;
			}
		}

		if(rowLength == 0)
		{
			rowLength = curRowLen;
		}
		if(rowLength != curRowLen )
		{
			// for some reason, we are getting a smaller size
			// simplely terminate here
			//DEBUG("Warning: one row of fle seems to have fewer data bits.\n");
			//res = false;
			break;
		}
		int *pRow = new int[rowLength];
		for(int i=0; i<rowLength; ++i)
		{
			if(buf[i] == '1')
			{
				pRow[i] = 1;
			}
			else if(buf[i] == '0')
			{
				pRow[i] = 0;
			}
			else if(buf[i] == '2')
			{
				pRow[i] = 2;
			}
            else if( buf[i] == '*' || buf[i] == '?' )
            {
                pRow[i] = MISSING_VALUE_BIT;
            }
            else
            {
                YW_ASSERT_INFO(false, "Un-recognized characters in input.");
                exit(1);
            }
		}
		// Now put it into a list
		rowsArray.push_back(pRow);
	}

	// Now set return value
	nCols = rowLength;

	return res;
}
#endif

// Dump the content of matrix
void BioSequenceMatrix ::Dump() const {
  cout << "positions: Matrix has  ";
  cout << nCols;
  cout << " columns and ";
  cout << rowsArray.size();
  cout << " rows.\n";
  for (unsigned int i = 0; i < rowsArray.size(); ++i) {
    for (int j = 0; j < nCols; ++j) {
      if (rowsArray[i][j] != MISSING_VALUE_BIT) {
        cout << rowsArray[i][j];
      } else {
        cout << "*";
      }
    }
    cout << endl;
  }
}

void BioSequenceMatrix ::OutputToFile(const char *fileName) const {
  //
  ofstream outFile;
  outFile.open(fileName);
  OutputToFile(outFile);
#if 0
	outFile << "Matrix has  ";
	outFile << nCols;
	outFile << " columns and ";
	outFile << rowsArray.size() ;
	outFile << " rows.\n";
	for(unsigned int i=0; i<rowsArray.size(); ++i)
	{
		for( int j=0; j< nCols; ++j)
		{
            if( rowsArray[i][j] != MISSING_VALUE_BIT)
            {
			    outFile << rowsArray[i][j];
            }
            else
            {
                outFile << "*";
            }
		}
		outFile << endl;
	}
#endif
  outFile.close();
}

void BioSequenceMatrix ::OutputToFile(ofstream &outFile) const {
  //
  // ofstream outFile;
  // outFile.open (fileName);
  outFile << "Matrix has  ";
  outFile << nCols;
  outFile << " columns and ";
  outFile << rowsArray.size();
  outFile << " rows.\n";
  for (unsigned int i = 0; i < rowsArray.size(); ++i) {
    for (int j = 0; j < nCols; ++j) {
      if (rowsArray[i][j] != MISSING_VALUE_BIT) {
        outFile << rowsArray[i][j];
      } else {
        outFile << "*";
      }
    }
    outFile << endl;
  }
  // outFile.close();
}

void BioSequenceMatrix ::ExchangeColumns(int c1, int c2) {
  // This function exchanges two columns in this matrix
  for (unsigned int i = 0; i < rowsArray.size(); ++i) {
    int tmp = rowsArray[i][c1];
    rowsArray[i][c1] = rowsArray[i][c2];
    rowsArray[i][c2] = tmp;
  }
}

// Offer direct access, but do not allow direct assignment
const int &BioSequenceMatrix ::operator()(int r, int c) const {
  return rowsArray[r][c];
}

int &BioSequenceMatrix ::operator()(int r, int c) { return rowsArray[r][c]; }

const int &BioSequenceMatrix ::GetValAt(int r, int c) const {
  return rowsArray[r][c];
}

void BioSequenceMatrix ::SetValAt(int r, int c, int val) {
  rowsArray[r][c] = val;
}

void BioSequenceMatrix ::GetAllSequences(vector<SEQUENCE> &seqs) const {
  seqs.clear();
  for (int i = 0; i < GetRowNum(); ++i) {
    SEQUENCE row;
    GetRow(i, row);
    seqs.push_back(row);
  }
}

void BioSequenceMatrix ::SubMatrix(int rt, int rb, int cl, int cr,
                                   BioSequenceMatrix &submat) const {
  // This function gets a submatrix, bounded from top row (rt), bottom row (rb)
  // left column (cl), right column cr
  submat.Clear();
  submat.SetSize(rb - rt + 1, cr - cl + 1);

  // Now we set rows
  for (int i = rt; i <= rb; ++i) {
    // get a vector of values
    vector<int> row;
    for (int j = cl; j <= cr; ++j) {
      row.push_back(rowsArray[i][j]);
    }

    // set row to submatrix
    submat.SetRow(i - rt, row);
  }
}

// This function gets a submatrix from selected sites
void BioSequenceMatrix ::SubMatrixSelectedSites(
    const vector<int> &sites, BioSequenceMatrix &submat) const {
  // This function gets a submatrix, with same number of rows but smaller number
  // of sites
  submat.Clear();
  submat.SetSize(rowsArray.size(), sites.size());

  // Now we set rows
  for (unsigned int i = 0; i < rowsArray.size(); ++i) {
    // get a vector of values
    vector<int> row;
    for (unsigned int j = 0; j < sites.size(); ++j) {
      int s = sites[j];
      YW_ASSERT_INFO(s < GetColNum(),
                     "SubMatrixSelectedSites: index out of range.");
      row.push_back(rowsArray[i][s]);
    }

    // set row to submatrix
    submat.SetRow(i, row);
  }
}

void BioSequenceMatrix ::SubMatrixSelectedRows(
    const vector<int> &rows, BioSequenceMatrix &submat) const {
  // This function gets a submatrix, with same number of rows but smaller number
  // of sites
  submat.Clear();
  submat.SetSize(rows.size(), nCols);

  // Now set rows
  for (unsigned int i = 0; i < rows.size(); ++i) {
    vector<int> r;
    GetRow(rows[i], r);
    submat.SetRow(i, r);
  }
}

void BioSequenceMatrix ::GetRow(int r, vector<int> &row) const {
  row.clear();
  for (int i = 0; i < nCols; ++i) {
    row.push_back(rowsArray[r][i]);
  }
}
void BioSequenceMatrix ::GetCol(int c, vector<int> &col) const {
  col.clear();
  for (int i = 0; i < GetRowNum(); ++i) {
    col.push_back(rowsArray[i][c]);
  }
}

int BioSequenceMatrix ::FindRow(const SEQUENCE &seq) const {
  // This function search the matrix to see if it contains this sequence
  // return -1 if not found
  YW_ASSERT_INFO(seq.size() == (unsigned int)GetColNum(),
                 "Size does not match.");

  for (int i = 0; i < GetRowNum(); ++i) {
    bool fFound = true;
    for (int j = 0; j < GetColNum(); ++j) {
      if (rowsArray[i][j] != seq[j]) {
        fFound = false;
        break;
      }
    }
    if (fFound == true) {
      return i;
    }
  }
  return -1;
}
int BioSequenceMatrix ::FindColumn(const SEQUENCE &seq) const {
  // This function search the matrix to see if it contains this sequence
  // return -1 if not found
  YW_ASSERT_INFO(seq.size() == (unsigned int)GetRowNum(),
                 "Size does not match.");

  for (int i = 0; i < GetColNum(); ++i) {
    bool fFound = true;
    for (int j = 0; j < GetRowNum(); ++j) {
      if (rowsArray[j][i] != seq[j]) {
        fFound = false;
        break;
      }
    }
    if (fFound == true) {
      // cout << "Col ";
      // DumpIntVec( seq );
      // cout << "is in this matrix: ";
      // this->Dump();
      return i;
    }
  }
  return -1;
}

// This function removes a set of columns that are specified in the set as
// duplicateSites
void BioSequenceMatrix ::RemoveColumns(set<int> &duplicateSites) {
  if (duplicateSites.size() == 0)
    return;

  // now we create  a new matrix with different size
  for (unsigned int i = 0; i < rowsArray.size(); ++i) {
    int *buf = new int[nCols - duplicateSites.size()];
    int cPos = 0;
    for (int j = 0; j < nCols; ++j) {
      if (duplicateSites.find(j) == duplicateSites.end()) {
        // j is not duplicate, so we should copy it
        buf[cPos++] = rowsArray[i][j];
      }
    }

    // now we free the memory of old buffer
    delete[] rowsArray[i];
    rowsArray[i] = buf;
  }

  nCols -= duplicateSites.size();
}

// Remove one row from matrix
void BioSequenceMatrix ::RemoveRows(set<int> &setRows) {
  vector<int *> saveMat;

  for (unsigned int i = 0; i < rowsArray.size(); ++i) {
    if (setRows.find(i) == setRows.end()) {
      // Only if row i is not inside the rows set, we will save it
      saveMat.push_back(rowsArray[i]);
    } else {
      // Ohterwise, we free it
      delete[] rowsArray[i];
    }
  }

  /*
          Now revert back
  */
  rowsArray.clear();
  rowsArray.swap(saveMat);
}

void BioSequenceMatrix ::SetSize(int nr, int nc) {
  // This function initialize a nr by nc matrix
  // and by default, fill in all 0 (false)
  nCols = nc;
  for (int i = 0; i < nr; ++i) {
    int *buf = new int[nc];
    for (int j = 0; j < nc; ++j) {
      buf[j] = 0;
    }
    rowsArray.push_back(buf);
  }
}

void BioSequenceMatrix ::FindNgbrDupCompSites(set<int> *pRemovedSet) {
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
  //	RemoveColumns( setOfRemovals );
}

void BioSequenceMatrix ::GetSeqsFeqs(map<SEQUENCE, int> &mapSeqFreqs) {
  // Insert, for each sequence, how many times they appears in the matrix
  mapSeqFreqs.clear();

  for (int r = 0; r < GetRowNum(); ++r) {
    SEQUENCE row;
    GetRow(r, row);
    if (mapSeqFreqs.find(row) == mapSeqFreqs.end()) {
      map<SEQUENCE, int>::value_type p(row, 1);
      mapSeqFreqs.insert(p);
    } else {
      mapSeqFreqs[row]++;
    }
  }
}

void BioSequenceMatrix ::GetSeqsOccurrence(
    map<SEQUENCE, set<int> > &mapSeqOccurs) {
  // For each distinct seq, find their occurance (which rows match this seq)
  mapSeqOccurs.clear();

  for (int r = 0; r < GetRowNum(); ++r) {
    SEQUENCE row;
    GetRow(r, row);
    if (mapSeqOccurs.find(row) == mapSeqOccurs.end()) {
      set<int> ss;
      map<SEQUENCE, set<int> >::value_type p(row, ss);
      mapSeqOccurs.insert(p);
    }
    mapSeqOccurs[row].insert(r);
  }
}

bool BioSequenceMatrix ::IsIntervalConsistent(int r1, int left1, int right1,
                                              int r2, int left2,
                                              int right2) const {
  // cout << "r1 = " << r1 << ", left1 = " << left1 << ", right1 = " << right1 ;
  // cout << ", r2 = " << r2 << ", left2 = " << left2 << ", right2 = " << right2
  // << endl;
  // Test if the two interval are consistent (i.e. has the same value at the
  // overlap
  INTERVAL iv1(left1, right1);
  INTERVAL iv2(left2, right2);
  INTERVAL ivInt;
  if (GetIntervalOverlap(iv1, iv2, ivInt) == false) {
    // If the interval are not overlapping, yes, they are consistent
    return true;
  }
  // cout << "intersection: left = " << ivInt.first << ", right = " <<
  // ivInt.second << endl;
  SEQUENCE row1;
  GetRow(r1, row1);
  SEQUENCE row1IV;
  GetSeqInterval(row1, row1IV, ivInt.first, ivInt.second);
  SEQUENCE row2;
  GetRow(r2, row2);
  SEQUENCE row2IV;
  GetSeqInterval(row2, row2IV, ivInt.first, ivInt.second);
  return (row1IV == row2IV);
}

bool BioSequenceMatrix ::IsMissingValue() {
#if 0
    // A rather inefficient way of doing things
    if( fMissingValue == true )
    {
        return true;
    }
#endif

  // now double check to make sure
  for (int r = 0; r < GetRowNum(); ++r) {
    for (int c = 0; c < GetColNum(); ++c) {
      if (GetValAt(r, c) == MISSING_VALUE_BIT) {
        //    fMissingValue = true;
        return true;
      }
    }
  }
  return false;
}

bool BioSequenceMatrix ::IsMissingValueInSite(int c) {
  // now double check to make sure
  for (int r = 0; r < GetRowNum(); ++r) {
    if (GetValAt(r, c) == MISSING_VALUE_BIT) {
      //    fMissingValue = true;
      return true;
    }
  }
  return false;
}

bool BioSequenceMatrix ::IsMissingValueInRow(int r) {
  return GetMissingValueNumInRow(r) > 0;
}

int BioSequenceMatrix ::GetMissingValueNumInRow(int r) {
  int res = 0;
  // now double check to make sure
  for (int c = 0; c < GetColNum(); ++c) {
    if (GetValAt(r, c) == MISSING_VALUE_BIT) {
      res++;
    }
  }
  return res;
}

void BioSequenceMatrix ::MapDupToNodup(map<int, int> &mapDupToNodup) const {
  // create a mapping from no-duplicate to duplicate indices
  set<int> rowsProcessed;

  int rowNoDup = 0;
  for (int r = 0; r < GetRowNum(); ++r) {
    if (rowsProcessed.find(r) != rowsProcessed.end()) {
      continue;
    }
    SEQUENCE seq;
    GetRow(r, seq);
    set<int> identRows;
    GetMultiplictyForRow(seq, identRows);
    // cout << "seq = ";
    // DumpSequence( seq );
    // DumpIntSet( identRows );
    // Add to map
    for (set<int>::iterator it = identRows.begin(); it != identRows.end();
         ++it) {
      int ss = *it;
      mapDupToNodup.insert(map<int, int>::value_type(ss, rowNoDup));
    }
    rowNoDup++;
    UnionSets(rowsProcessed, identRows);
  }
}

int BioSequenceMatrix ::GetNodupRowsNum(vector<int> *pListUniqeRowIndex) const {
  int res = 0;
  for (int r = 0; r < GetRowNum(); ++r) {
    SEQUENCE seq;
    GetRow(r, seq);
    // cout << "seq = ";
    // DumpSequence( seq );
    // DumpIntSet( identRows );
    // Add to map
    bool fUnique = true;
    for (int r2 = 0; r2 < r; ++r2) {
      SEQUENCE seq2;
      GetRow(r2, seq2);
      if (seq2 == seq) {
        fUnique = false;
        break;
      }
    }
    if (fUnique == true) {
      res++;
      if (pListUniqeRowIndex != NULL) {
        pListUniqeRowIndex->push_back(r);
      }
    }
  }
  return res;
}

////////////////////////////////////////////////////////////////////////////////
//		Inernal utility functions
////////////////////////////////////////////////////////////////////////////////

bool BioSequenceMatrix ::CmpColumns(int c1, int c2) {
  bool res = true;

  if (c1 == c2) {
    return true;
  }

  for (unsigned int i = 0; i < rowsArray.size(); ++i) {
    if (rowsArray[i][c1] != rowsArray[i][c2]) {
      res = false;
      break;
    }
  }

  return res;
}
