#include "Utils.h"
#include "cstdio"
#include "cstdlib"
#include "ctime"

//////////////////////////////////////////////////////////////////////////////
// Potential Junk, save it here
//////////////////////////////////////////////////////////////////////////////

static void MakeComplementSet(int start, int end, const set<int> &origSet,
                              set<int> &compSet) {
  for (int i = start; i <= end; ++i) {
    if (origSet.find(i) == origSet.end()) {
      compSet.insert(i);
    }
  }
}

#if 0
static long ComputeCombNum(int n, int k)
{
	long res = 1;
	for(int i=)
	{
	}
	return res;
}
#endif

//////////////////////////////////////////////////////////////////////////////
void SubtractSets(set<int> &s1, const set<int> &s2) {
  if (s2.size() == 0) {
    return;
  }
  set<int> res;
  // this function performs set intersection, i.e. s1=s1 ^s2
  for (set<int>::iterator it = s1.begin(); it != s1.end(); ++it) {
    if (s2.find(*it) == s2.end()) {
      res.insert(*it);
    }
  }
  s1.clear();
  s1 = res;
}

void JoinSets(const set<int> &s1, const set<int> &s2, set<int> &res) {
  res.clear();
  for (set<int>::iterator it = s1.begin(); it != s1.end(); ++it) {
    if (s2.find(*it) != s2.end()) {
      res.insert(*it);
    }
  }
}

void UnionSets(set<int> &sTotal, const set<int> &sToBeAdd) {
  for (set<int>::iterator it = sToBeAdd.begin(); it != sToBeAdd.end(); ++it) {
    sTotal.insert(*it);
  }
}

// templates
void JoinSets(const set<char> &s1, const set<char> &s2, set<char> &res) {
  res.clear();
  for (set<char>::iterator it = s1.begin(); it != s1.end(); ++it) {
    if (s2.find(*it) != s2.end()) {
      res.insert(*it);
    }
  }
}
void SubtractSets(set<char> &s1, const set<char> &s2) {
  if (s2.size() == 0) {
    return;
  }
  set<char> res;
  // this function performs set intersection, i.e. s1=s1 ^s2
  for (set<char>::iterator it = s1.begin(); it != s1.end(); ++it) {
    if (s2.find(*it) == s2.end()) {
      res.insert(*it);
    }
  }
  s1.clear();
  s1 = res;
}
void UnionSets(set<char> &sTotal, const set<char> &sToBeAdd) {
  for (set<char>::iterator it = sToBeAdd.begin(); it != sToBeAdd.end(); ++it) {
    sTotal.insert(*it);
  }
}
void DumpSet(const set<char> &s) {
  cout << "Set contains: ";
  for (set<char>::iterator it = s.begin(); it != s.end(); ++it) {
    cout << (int)*it << ",";
  }
  cout << endl;
}

void ConvIntSetToCharSet(const set<int> &si, set<char> &sc) {
  sc.clear();
  for (set<int>::iterator it = si.begin(); it != si.end(); ++it) {
    sc.insert((int)*it);
  }
}
void ConvCharSetToIntSet(const set<char> &sc, set<int> &si) {
  si.clear();
  for (set<char>::iterator it = sc.begin(); it != sc.end(); ++it) {
    si.insert(*it);
  }
}

// others
void RmIntValFromSet(set<int> &s, int v) {
  for (set<int>::iterator it = s.begin(); it != s.end(); ++it) {
    if (*it == v) {
      s.erase(it);
      return;
    }
  }
}

void DumpIntSet(const set<int> &incSet) {
  //#ifdef BG_DEBUG
  cout << "Set contains: ";
  for (set<int>::iterator it = incSet.begin(); it != incSet.end(); ++it) {
    cout << *it << ",";
  }
  cout << endl;
  //#endif
}

void DumpIntSetNoReturn(const set<int> &incSet) {
  for (set<int>::iterator it = incSet.begin(); it != incSet.end(); ++it) {
    cout << *it << ",";
  }
}

void DumpIntVec(const vector<int> &intVec) {
  cout << "Vector contains: ";
  for (int i = 0; i < intVec.size(); ++i) {
    cout << intVec[i] << ",";
  }
  cout << endl;
}

void PopulateSetByVec(set<int> &dest, const vector<int> &srcVec) {
  dest.clear();
  for (int i = 0; i < srcVec.size(); ++i) {
    dest.insert(srcVec[i]);
  }
}

void PopulateVecBySet(vector<int> &dest, const set<int> &srcSet) {
  dest.clear();
  for (set<int>::iterator it = srcSet.begin(); it != srcSet.end(); ++it) {
    dest.push_back(*it);
  }
}

void CopyIntSet(set<int> &dest, const set<int> &src) {
  dest.clear();
  for (set<int>::iterator it = src.begin(); it != src.end(); ++it) {
    dest.insert(*it);
  }
}

void CopyIntVec(vector<int> &dest, const vector<int> &src) {
  dest.clear();
  for (int i = 0; i < src.size(); ++i) {
    dest.push_back(src[i]);
  }
}

void CopySetIntVec(set<vector<int> > &dest, const set<vector<int> > &src) {
  dest.clear();
  for (set<vector<int> >::iterator it = src.begin(); it != src.end(); ++it) {
    vector<int> v;
    CopyIntVec(v, *it);
    dest.insert(v);
  }
}

bool IsVecSame(const vector<int> &v1, const vector<int> &v2) {
  if (v1.size() != v2.size()) {
    return false;
  }
  for (int i = 0; i < v1.size(); ++i) {
    if (v1[i] != v2[i]) {
      return false;
    }
  }
  return true;
}

bool IsIntVecInSet(const set<vector<int> > &s, const vector<int> &v) {
  for (set<vector<int> >::iterator it = s.begin(); it != s.end(); ++it) {
    vector<int> v1 = *it;
    if (IsVecSame(v, v1) == true) {
      return true;
    }
  }
  return false;
}

// The following two functions could be used in dynamic programing
// when, for example, we need to consider all sets
// One limitation is that it is limited to integer range
// that up to 32 bits
void ConvIntToVec(unsigned int val, vector<int> &vec, int numBits) {
  // we would store the least significant bit as vec[0]
  vec.clear();
  if (numBits <= 32) {
    for (int i = 0; i < numBits; ++i) {
      if ((val & 0x1) == 0) {
        vec.push_back(0);
        //				vec.insert( vec.begin(), 0);
      } else {
        vec.push_back(1);
        //				vec.insert( vec.begin(), 1);
      }
      val = val >> 1;
    }
  }
}

unsigned int ConvVecToInt(const vector<int> &vec) {
  // assume vec[0] is least siginicant
  unsigned int res = 0;

  for (int i = vec.size() - 1; i >= 0; --i) {
    YW_ASSERT_INFO(vec[i] == 0 || vec[i] == 1,
                   "In ConvVecToInt, vector is not binary.");
    // cout << "res = " << res << endl;
    if (vec[i] == 1) {
      res += 1;
    }
    if (i > 0) {
      res = res << 1;
    }
  }

  return res;
}

void ConvIntToVecMSB(unsigned int val, vector<int> &vec, int numBits) {
  // we would store the least significant bit as vec[0]
  YW_ASSERT_INFO(numBits <= 32, "ConvIntToVecMSB :: numBits is too large.");
  ConvIntToVec(val, vec, numBits);
  ReverseIntVec(vec);
}

unsigned int ConvVecToIntMSB(const vector<int> &vec) {
  vector<int> vecMSB = vec;
  // cout << "vec = ";
  // DumpIntVec( vec );
  ReverseIntVec(vecMSB);
  // cout << "vec = ";
  // DumpIntVec( vec );
  return ConvVecToInt(vecMSB);
}

void ReverseIntVec(vector<int> &vec) {
  // cout << "Before switching: vec = ";
  // DumpIntVec( vec );
  // This function would reverse the integer vector, i.e. vec[0] = vec[n-1] and
  // so on
  for (int i = 0; i < vec.size() / 2; ++i) {
    int tmp = vec[vec.size() - 1 - i];
    vec[vec.size() - 1 - i] = vec[i];
    vec[i] = tmp;
  }
  // cout << "After switching: vec = ";
  // DumpIntVec( vec );
}

unsigned int CalcBitInt(int pos, int width) {
  return 0x1 << (width - 1 - pos);
  //	return 0x1 << pos;
}

// This function update the eumeration index into different set size
// it works as following: suppose we have 3 sets of size 2,3,2 and
// initially the index is 0,1,0.
// Then after calling this function it becomes 0,1,1. And next time, since we
// reach the set bound, we have to change it to 0,2,0 (remember reseting to 0 at
// pos 3) This function returns false when we reach the end
bool GetNextEnumVec(vector<int> &curPos, const vector<int> &limitvec) {
  if (limitvec.size() != curPos.size()) {
    return false;
  }

  // Now we find the last position i in curPos,
  // s.t. curPos[i] < limitvec[i]
  int i = -1;
  for (i = curPos.size() - 1; i >= 0; --i) {
    if (curPos[i] >= limitvec[i]) {
      return false;
    }

    if (curPos[i] < limitvec[i] - 1) {
      break;
    }
  }

  if (i < 0) {
    // OK, we can not continue since we have run out of search space
    return false;
  }

  // Ohterwise, we increment this position and reset the following positions
  curPos[i] = curPos[i] + 1;
  for (int j = i + 1; j < curPos.size(); ++j) {
    curPos[j] = 0;
  }
  return true;
}

void YW_ASSERT(bool f) {
  if (f == false) {
    cout << "Assertion error" << endl;
    exit(1);
  }
}

void YW_ASSERT_INFO(bool f, const char *info) {
  //#if 0
  if (f == false) {
    cout << "Assertion Error: " << info << endl;
    ;
    exit(1);
  }
  //#endif
}

void RemoveFromIntSet(vector<int> &targetSet, int val) {
  for (vector<int>::iterator it = targetSet.begin(); it != targetSet.end();
       ++it) {
    if (*it == val) {
      targetSet.erase(it);
      return;
    }
  }
}

bool IsIntSetEquiv(const set<vector<int> > &s1, const set<vector<int> > &s2) {
  if (s1.size() != s2.size()) {
    return false;
  }
  // we check to see if every element in s1 is also in s2
  for (set<vector<int> >::iterator it = s1.begin(); it != s1.end(); ++it) {
    if (s2.find(*it) == s2.end()) {
      return false;
    }
  }
  return true;
}

void OrderInt(int &i1, int &i2) {
  // Exchange two number if i1 is greater than i2
  if (i1 > i2) {
    int tmp = i2;
    i2 = i1;
    i1 = tmp;
  }
}

static int QSortCompare(const void *arg1, const void *arg2) {
  /* Compare all of both strings: */
  // assume sorting in accending order
  int n1 = *((int *)arg1);
  int n2 = *((int *)arg2);
  // cout <<"arg1 = " << n1 << ", arg2 = " << n2 << endl;
  if (n1 > n2) {
    return 1;
  } else if (n1 < n2) {
    return -1;
  } else {
    return 0;
  }
}

void SortIntVec(vector<int> &vec, int start, int end) {
  //#if 0
  if (vec.size() == 0) {
    // do nothing
    return;
  }
  if (end < 0) {
    end = vec.size() - 1;
  }
  int sortLen = end - start + 1;
  int *array = new int[sortLen];
  for (int i = start; i <= end; ++i) {
    array[i - start] = vec[i];
  }
  qsort((void *)array, sortLen, sizeof(int), QSortCompare);
  // Now write back
  for (int i = start; i <= end; ++i) {
    vec[i] = array[i - start];
  }

  delete[] array;
//#endif
#if 0
	// Sort the vector, by the most obvious method
	for(int i=0; i<vec.size()-1; ++i)
	{
		for(int j=i+1; j<vec.size(); ++j)
		{
			if( vec[i] > vec[j] )
			{
				OrderInt( vec[i], vec[j] );
			}
		}
	}
#endif
#if 0
	if (end < 0 )
	{
		end = vec.size() - 1;
	}
	for(int i=start; i<end; ++i)
	{
		for(int j=i+1; j<=end; ++j)
		{
			if( vec[i] > vec[j] )
			{
				OrderInt( vec[i], vec[j] );
			}
		}
	}
#endif
}

void GetFirstCombo(int k, int n, vector<int> &posvec) {
  posvec.clear();
  for (int i = 0; i < k; ++i) {
    posvec.push_back(i);
  }
}

bool GetNextCombo(int k, int n, vector<int> &posvec) {
  // The idea is to move the rightmost (movable) value to the right
  int startpos = k - 1;
  return GetNextComboFrom(k, n, posvec, startpos);
#if 0
	while(pos >= 0)
	{
		if( posvec[pos] < pos + (n-k) )
		{
			posvec[pos] = posvec[pos] + 1;
			for(int i=pos+1; i<k; ++i)
			{
				posvec[i] = i + (posvec[pos] - pos);
			}
			return true;
		}
		else
		{
			pos --;
		}
	}
	return false;
#endif
}

bool GetNextComboFrom(int k, int n, vector<int> &posvec, int startpos) {
  // This function differs from the previous one in that it starts moving
  // forward
  // not neccessary from pos=k-1, but from the given startpos
  // this allows flexibility of bypassing searching when that space will not be
  // productive
  int pos = startpos;
  while (pos >= 0) {
    if (posvec[pos] < pos + (n - k)) {
      posvec[pos] = posvec[pos] + 1;
      for (int i = pos + 1; i < k; ++i) {
        posvec[i] = i + (posvec[pos] - pos);
      }
      return true;
    } else {
      pos--;
    }
  }
  return false;
}

#if 0
int ConvComboToIndex(int numCells, const vector<int> &posvec)
{
	// This function converts a position vector into an index
	// This is useful when performing dynamic programming
	// The idea is to check each position in vector
	// For i-th position, if it is greater than i, then
	// plus C(n-posvec[i], k-i).
	int res = 0;

	return res;
}
#endif

double GetRandFraction() {
  // Now we try random method. Flip a coin, to decide whether to take this new
  // choice
  static bool isRandInit = false;
  if (isRandInit == false) {
    srand((unsigned)time(NULL));
    isRandInit = true;
  }

  double c = (double)(rand()) / RAND_MAX;
  return c;
}

void GetBoolVec(int num, const vector<int> &posvec, vector<bool> &bvec) {
  int pos = 0;
  for (int i = 0; i < num; ++i) {

    if (pos >= posvec.size() || i < posvec[pos]) {
      bvec.push_back(false);
    } else if (i == posvec[pos]) {
      bvec.push_back(true);
      pos++;
    } else {
      YW_ASSERT_INFO(false, "GetBoolVec");
    }
  }
}

void GetIntVec(int num, const vector<int> &posvec, vector<int> &bvec) {
  bvec.clear();
  int pos = 0;
  for (int i = 0; i < num; ++i) {

    if (pos >= posvec.size() || i < posvec[pos]) {
      bvec.push_back(0);
    } else if (i == posvec[pos]) {
      bvec.push_back(1);
      pos++;
    } else {
      YW_ASSERT_INFO(false, "GetIntVec");
    }
  }
}

// Coomposite bound operations
int CalcCompositeBound(map<INTERVAL, int> &mapIntervalBds, int left, int right,
                       vector<int> &locBreakpoints) {
  // This method outputs a composite bound for the given interval
  int res = 0;

  int lenInterval = right - left + 1;
  vector<int> lbHelper;
  // Initialize our lb helper data
  for (int i = 0; i < lenInterval; ++i) {
    lbHelper.push_back(0);
  }

  // Now we scan through all the initervals in range (from 'left' to 'right')
  // we also need to make sure we start by sotring interval based on its right
  // end
  for (int re = left + 1; re <= right; ++re) {
    for (int le = left; le < re; ++le) {
      // we now consider the interval [le, re]
      INTERVAL iv(le, re);
      if (mapIntervalBds.find(iv) == mapIntervalBds.end()) {
        // nothing needs to be done, if interval is not in map
        continue;
      }
      int valInt = mapIntervalBds[iv];

      // we now figure out lbHelper value based on the value
      int lbSofar = 0;
      for (int i = le; i < re; ++i) {
        lbSofar += lbHelper[i - left];
      }
      if (lbSofar < valInt) {
        // we make up the diff in the last slot
        lbHelper[re - left - 1] += valInt - lbSofar;
      }
    }
  }

  // Finally, we tally the result
  for (int i = 0; i < lenInterval; ++i) {
    if (lbHelper[i] != 0) {
      for (int j = 0; j < lbHelper[i]; ++j) {
        locBreakpoints.push_back(i + left);
      }
#if 0
          cout << "Between site " << i+1 << " and site " << i+2 << ", there are " << lbHelper[i]  << " recombs." << endl;
#endif
    }
    res += lbHelper[i];
  }

  return res;

#if 0
	// This method outputs a composite bound for the given interval
	int res = 0;

	int lenInterval = right-left+1;
	vector<int> lbHelper;
	// Initialize our lb helper data
	for(int i=0; i<lenInterval; ++i)
	{
		lbHelper.push_back( 0 );
	}

	// Now we scan through all the initervals in range (from 'left' to 'right')
	// we also need to make sure we start by sotring interval based on its right end
	for(int re = left+1; re<=right; ++re)
	{
		for(int le = left; le <re; ++le)
		{
			// we now consider the interval [le, re]
			INTERVAL iv(le, re);
			if( mapIntervalBds.find( iv ) == mapIntervalBds.end() )
			{
				// nothing needs to be done, if interval is not in map
				continue;
			}
			int valInt = mapIntervalBds[ iv ];

			// we now figure out lbHelper value based on the value
			int lbSofar = 0;
			for(int i=le; i < re; ++i)
			{
				lbSofar += lbHelper[  i-left  ];
			}
			if( lbSofar < valInt)
			{
				// we make up the diff in the last slot
				lbHelper[ re - left -1] += valInt - lbSofar;
			}
		}
	}

	// Finally, we tally the result
	for(int i=0; i<lenInterval; ++i)
	{
        if( lbHelper[i] != 0)
        {
            for(int j=0; j<lbHelper[i]; ++j)
            {
                locBreakpoints.push_back( i+left );
            }
//        #if 0
          cout << "Between site " << i+1 << " and site " << i+2 << ", there are " << lbHelper[i]  << " recombs." << endl;
//        #endif
        }
		res += lbHelper[i];
	}

	return res;
#endif
}

void OutputBounds(char *boundsFileName, map<INTERVAL, int> &mapIntervalBds,
                  int nSites) {
  // This function outputs results (that are stored inside a map)
  // First open a file as named as passed in
  // Now open file to write out
  //    char fname[1024];
  //    strcpy( fname, boundsFileName );
  ofstream outFile(boundsFileName);
  if (outFile.is_open() == false) {
    cout << "Can not open output file: " << boundsFileName << endl;
    return;
  }

  outFile << "bounds-from-HapBound\n";
  for (int i = 0; i < nSites - 1; ++i) {
    for (int j = i + 1; j < nSites; ++j) {
      INTERVAL iv(i, j);
      if (mapIntervalBds.find(iv) != mapIntervalBds.end()) {
        outFile << i + 1 << "  " << j + 1 << "  " << mapIntervalBds[iv] << endl;
      } else {
        cout << "Warning: interval not complete. Missing (" << i << ", " << j
             << ")" << endl;
      }
    }
  }
  outFile.close();
}

// Some combinatorial tricks
void InitPermutation(vector<int> &nvec, const vector<int> &reference) {
  // We ASSUME reference is already sorted
  nvec = reference;
  //	SortIntVec( nvec );
}

bool GetNextPermutation(vector<int> &nvec, const vector<int> &reference) {
  // Now, we try to find the next position
  // The idea is to start from the right, and check if we can use it as
  // the starting location
  for (int i = nvec.size() - 1; i >= 0; --i) {
    // Make sure this number is not already maximum
    if (nvec[i] == reference[reference.size() - 1]) {
      continue;
    }
    // Now, we make sure there is at least one element to the right of it
    int minLarger = HAP_MAX_INT;
    int pos = -1;
    for (int j = i + 1; j < nvec.size(); ++j) {
      if (nvec[j] > nvec[i] && minLarger > nvec[j]) {
        pos = j;
        minLarger = nvec[j];
      }
    }

    // If no such j is found, stop
    if (pos < 0) {
      continue;
    }

    // Otherwise, we stop here by taking this position
    nvec[pos] = nvec[i];
    nvec[i] = minLarger;

    SortIntVec(nvec, i + 1, nvec.size() - 1);
    return true;
  }

  return false;
}

int ConvertToLinear(int r1, int r2, int nRows) {
  int n = nRows;
  int n1 = (r1 + 1) * (n - 1) - ((r1 + 1) * r1) / 2;
  return n1 - (n - r2);
}

int ConvertToLinearEq(int r1, int r2, int nRows) {
  // The only difference from the above is this one allow r1= r2
  int n = nRows;
  int n1 = (r1 + 1) * (n - 1) - ((r1 + 1) * r1) / 2;
  return n1 - (n - r2);
}

// void ConvertLinearToTwoIndices( int idLinear, int nRows, int &r1, int &r2 )
//{
//}

//****************************************************************************************
// Recombination/Mutation utilities
//****************************************************************************************

bool IsMissingValueBit(int bit) { return bit == MISSING_VALUE_BIT; }

int GetMissingValueBit() { return MISSING_VALUE_BIT; }

bool IsTwoStatesCompatible(int bit1, int bit2) {
  // we say two states are compatible if either they match exactly or one of
  // them is missing value
  return (bit1 == bit2) || IsMissingValueBit(bit1) || IsMissingValueBit(bit2);
}
void FillVecWithMV(SEQUENCE &seq, int len) {
  // note, do not clear up original seq
  for (int i = 0; i < len; ++i) {
    seq.push_back(MISSING_VALUE_BIT);
  }
}

bool IsSeqHasMV(const SEQUENCE &seq) {
  for (int i = 0; i < (int)seq.size(); ++i) {
    if (IsMissingValueBit(seq[i]) == true) {
      return true;
    }
  }
  return false;
}
bool AreTwoSeqsCompatible(const SEQUENCE &seq1, const SEQUENCE &seq2) {
  if (seq1.size() != seq2.size()) {
    return false; // size must match
  }
  for (int i = 0; i < (int)seq1.size(); ++i) {
    if (IsTwoStatesCompatible(seq1[i], seq2[i]) == false) {
      return false;
    }
  }
  return true;
}
void GetCompatibleSeqForTwo(const SEQUENCE &seq1, const SEQUENCE &seq2,
                            SEQUENCE &consensus) {
  YW_ASSERT_INFO(seq1.size() == seq2.size(), "Size mismatch");

  consensus.clear();
  for (int i = 0; i < (int)seq1.size(); ++i) {
    YW_ASSERT_INFO(IsTwoStatesCompatible(seq1[i], seq2[i]),
                   "Can not form compatible");
    if (IsMissingValueBit(seq1[i]) == false) {
      consensus.push_back(seq1[i]);
    } else {
      consensus.push_back(seq2[i]);
    }
  }
}

extern void DumpSequence(const SEQUENCE &seq);
void MutateSeqAtSite(SEQUENCE &seq, int site) {
  // cout << "MutateSeqAtSite: seq = ";
  // DumpSequence(seq);
  // cout<< "site = " << site << endl;
  YW_ASSERT_INFO(IsMissingValueBit(seq[site]) == false,
                 "Can not mutate a missing value");
  if (seq[site] == 0) {
    seq[site] = 1;
  } else {
    seq[site] = 0;
  }
}

void RecombSequencesAt(const SEQUENCE &s1, const SEQUENCE &s2, int brPt,
                       SEQUENCE &sr) {
  // NOTE, ordering is important here. You may need to call this function
  // twice, if you want to consider recombinatino from both direction
  // This function assume the first part of s1 is taken first
  // Another thing to note is 0 <= brPt <= s1.size()-2
  sr.clear();
  for (int i = 0; i <= brPt; ++i) {
    sr.push_back(s1[i]);
  }
  for (int i = brPt + 1; i < s2.size(); ++i) {
    sr.push_back(s2[i]);
  }
}

bool IsSeqRecombinnable(const SEQUENCE &s1, const SEQUENCE &s2,
                        const SEQUENCE &st) {
  // note here, we do not differenceitate left or right
  INTERVAL iv;
  return IsSeqRecombinnableIV(s1, s2, st, iv) ||
         IsSeqRecombinnableIV(s2, s1, st, iv);
}

bool IsSeqRecombinnableIV(const SEQUENCE &sleft, const SEQUENCE &sright,
                          const SEQUENCE &st, INTERVAL &iv) {
  // here assume s1 is left, and s2 is right. THIS IS IMPORTANT
  // Here, iv returns the interval where the breakpoiint can fall
  //    cout <<"s1.size = " << s1.size() << ", s2.size = " << s2.size() << ",
  //    st.size = " << st.size() << endl;
  YW_ASSERT((sleft.size() == sright.size()) && (sleft.size() == sright.size()));

  // Now we exam the recombination of s1 and s2, into st
  // first, we find the first location that does not match
  int pos = 0;
  while (pos < (int)sleft.size() &&
         IsTwoStatesCompatible(sleft[pos], st[pos]) == true) {
    pos++; // continue since they all match
  }
  // cout << "1. pos = " << pos << endl;
  if (pos == (int)sleft.size()) {
    iv.first = 0;
    iv.second = (int)sleft.size() - 1;
    return true;
  }

  // If there is no matching prefix, there is no solution
  if (pos == 0) {
    return false;
  }

  // cout << "2. pos = " << pos << endl;

  // If reach here, we are at the second difference, (the break point)
  // Now the iv should be set to the maximal range where two sequence match
  iv.first = pos - 1;
  iv.second = pos - 1;
  for (int i = pos - 1; i >= 0; --i) {
    if (IsTwoStatesCompatible(sleft[i], sright[i]) == true) {
      iv.first--;
    } else {
      break;
    }
  }

  // cout << "3. pos = " << pos << endl;

  while (pos < (int)sleft.size() &&
         IsTwoStatesCompatible(sright[pos], st[pos]) == true) {
    pos++;
  }
  // cout << "4. pos = " << pos << endl;

  if (pos == (int)sleft.size()) {
    return true;
  } else {
    return false;
  }
}

void AddUniqueSeqToVec(const SEQUENCE &seq, vector<SEQUENCE> &vecSeqs) {
  for (int i = 0; i < vecSeqs.size(); ++i) {
    if (vecSeqs[i] == seq) {
      return; // Duplicate here
    }
  }
  vecSeqs.push_back(seq);
}

bool IsSeqInVec(const SEQUENCE &seq, const vector<SEQUENCE> &vecSeqs) {
  for (int i = 0; i < vecSeqs.size(); ++i) {
    if (vecSeqs[i] == seq) {
      return true; // Duplicate here
    }
  }
  return false;
}

bool IsSeqInSet(const SEQUENCE &seq, const set<SEQUENCE> &setSeqs) {
  for (set<SEQUENCE>::iterator it = setSeqs.begin(); it != setSeqs.end();
       ++it) {
    if (*it == seq) {
      return true; // Duplicate here
    }
  }
  return false;
}

void GetEqualSubseq(const SEQUENCE &seq1, const SEQUENCE &seq2, int seedPos,
                    int &left, int &right) {
  // This function gets the best matching regions between seq1/se2 around the
  // seed Note, this function does not include comparision with the seedPos
  // IMPORTANT!!!! That is, it EXCLUDES seedPos
  if (seedPos < 0 || seedPos >= seq1.size()) {
    left = right = -1;
    return;
  }
  left = right = seedPos;

  // Now start checking
  int i;
  for (i = seedPos - 1; i >= 0; i--) {
    if (IsTwoStatesCompatible(seq1[i], seq2[i]) == false) {
      break;
    }
  }
  if (i >= 0) {
    left = i;
  }
  for (i = seedPos + 1; i < seq1.size(); ++i) {
    if (IsTwoStatesCompatible(seq1[i], seq2[i]) == false) {
      break;
    }
  }
  if (i < seq1.size()) {
    right = i;
  } else {
    right = seq1.size() - 1;
  }
}

// Compute the segments in the region of [left, right]
int CompareSegments(const SEQUENCE &seq, const SEQUENCE &targetSeq, int left,
                    int right) {
  int res = 0;
  // int numGap = 0;						// do we consider
  // gap here?
  for (int i = left; i <= right; ++i) {
    if (IsTwoStatesCompatible(seq[i], targetSeq[i]) == true) {
      res++;
    }
  }
  return res;
}

int IsSeqsMutPair(const SEQUENCE &seq1, const SEQUENCE &seq2) {
  // This function test if seq1/seq2 mutates at some site
  // If not, return -1. Otherwise, return the site that they differ
  // note that when there is missing data, we assume two compatible vals
  // do not form muatant pair
  int res = -1;
  for (int i = 0; i < seq1.size(); ++i) {
    if (IsTwoStatesCompatible(seq1[i], seq2[i]) == false) {
      if (res < 0) {
        res = i;
      } else {
        // We have seen one difference before, not mutation pair
        return -1;
      }
    }
  }
  return res;
}

int CalcSequencesDistance(const SEQUENCE &seq1, const SEQUENCE &seq2) {
  // similarly, when there is missing data, assume they can be fit to the best
  int res = 0;
  for (int i = 0; i < seq1.size(); ++i) {
    if (IsTwoStatesCompatible(seq1[i], seq2[i]) == false) {
      res++;
    }
  }
  return res;
}

void GetNewSequences(const set<SEQUENCE> &setNewNodes,
                     const set<SEQUENCE> &setExistingSeqs,
                     vector<SEQUENCE> &seqNews) {
  // Get the new nodes into a vector
  for (set<SEQUENCE>::iterator it = setNewNodes.begin();
       it != setNewNodes.end(); ++it) {
    if (setExistingSeqs.find(*it) == setExistingSeqs.end()) {
      seqNews.push_back(*it);
    }
  }
}
//************************************************************************************************
// Utilities for haplotyping
//************************************************************************************************

void GenHapRowsSetFromGenoRows(set<int> &hapRowsSet, int numGenoRows) {
  hapRowsSet.clear();
  for (int i = 0; i < 2 * numGenoRows; ++i) {
    hapRowsSet.insert(i);
  }
}

bool IsTwoLabelSetsCompatible(const set<int> &partition,
                              const vector<int> &genoSite, bool &fZeroOne) {
  // Two 2-label sets are compatible when we treat 2i-1, 2i heterozygotes the
  // same
  int *tblEntryOccurs = new int[genoSite.size()];
  for (int i = 0; i < genoSite.size(); ++i) {
    tblEntryOccurs[i] = 0;
  }

  for (set<int>::iterator it = partition.begin(); it != partition.end(); ++it) {
    int r = *it;
    tblEntryOccurs[r / 2]++;
  }

  // We treat tblNeeded as the 0-partition
  // We let tblNeeded to store zero-element occurance
  int *tblNeeded = new int[genoSite.size()];
  for (int i = 0; i < genoSite.size(); ++i) {
    if (genoSite[i] == 2) {
      tblNeeded[i] = 1;
    } else if (genoSite[i] == 1) {
      tblNeeded[i] = 2;
    }
  }

  // Now to see if it matches
  int res = false;

  bool earlyBreak = false;
  for (int i = 0; i < genoSite.size(); ++i) {
    if (tblEntryOccurs[i] != tblNeeded[i]) {
      earlyBreak = true;
      break;
    }
    //        if(  tblNeeded[i] ==1 &&  tblNeeded[i] != tblEntryOccurs[i]  )
    //        {
    //            earlyBreak = true;
    //            break;
    //        }
    //        if( tblNeeded[i] != 1 &&  tblNeeded[i] != tblEntryOccurs[i] )
    //        {
    //            earlyBreak = true;
    //            break;
    //        }
  }
  if (earlyBreak == false) {
    fZeroOne = true;
    res = true;
  }

  // check another possibility
  if (res == false) {
    earlyBreak = false;
    for (int i = 0; i < genoSite.size(); ++i) {
      if (tblEntryOccurs[i] != 2 - tblNeeded[i]) {
        earlyBreak = true;
        break;
      }
      //            if(  tblNeeded[i] ==1 &&  tblNeeded[i] != tblEntryOccurs[i]
      //            )
      //            {
      //                earlyBreak = true;
      //                break;
      //            }
      //            if( tblNeeded[i] != 1 &&  tblNeeded[i] == tblEntryOccurs[i]
      //            )
      //            {
      //                earlyBreak = true;
      //                break;
      //            }
    }
    if (earlyBreak == false) {
      fZeroOne = false; // indicate this is a zero set
      res = true;
    }
  }

#if 0
    for(int i=0; i<genoSite.size(); ++i)
    {
        if( genoSite[i] != 2  )
        {
            if(tblEntryOccurs[i] == 1)
            {
                return false;       // incompatible
            }
            else if(  )
            {
            }
        }
        else if(genoSite[i] == 2 && tblEntryOccurs[i] != 1)
        {
            return false;
        }
    }
#endif
  delete[] tblEntryOccurs;
  delete[] tblNeeded;
  return res;
}

void GenGenoPartitions(const vector<int> &genoSite, vector<int> &part0,
                       vector<int> &part1) {
  part0.clear();
  part1.clear();
  for (int i = 0; i < genoSite.size(); ++i) {
    if (genoSite[i] == 2) {
      // We simply put the first into part0
      part0.push_back(2 * i);
      part1.push_back(2 * i + 1);
    } else if (genoSite[i] == 0) {
      part0.push_back(2 * i);
      part0.push_back(2 * i + 1);
    } else if (genoSite[i] == 1) {
      part1.push_back(2 * i);
      part1.push_back(2 * i + 1);
    }
  }
}

bool Is2TwoLabelMatch(int lbla, int lblb) {
  // The two label matches if the belong to the same geno row
  return (lbla / 2 == lblb / 2);
}

bool IsTwoLabelSetContained(int genoLength, const vector<int> &setContainer,
                            const vector<int> &setContained) {
  if (setContained.size() > setContainer.size()) {
    // cout << "IsTwoLabelSetContained: size mismatched.\n";
    return false; // size mismatch
  }
  // cout << "setContainer: ";
  // DumpIntVec( setContainer);
  // cout << "setContained: ";
  // DumpIntVec( setContained);

  int *tblContainer = new int[genoLength];
  int *tblContained = new int[genoLength];

  // Init tbl
  for (int i = 0; i < genoLength; ++i) {
    tblContainer[i] = 0;
    tblContained[i] = 0;
  }
  for (int i = 0; i < setContainer.size(); ++i) {
    tblContainer[setContainer[i] / 2]++;
  }
  for (int i = 0; i < setContained.size(); ++i) {
    tblContained[setContained[i] / 2]++;
  }
  // Check to see if one contain another
  int res = true;
  for (int i = 0; i < genoLength; ++i) {
    if (tblContained[i] > tblContainer[i]) {
      // cout << "Not contained!.\n";
      res = false;
      break;
    }
  }

  delete[] tblContainer;
  delete[] tblContained;
  return res;
}

void CalcGenoNum(int genoLength, const vector<int> &partition,
                 vector<int> &genoNums) {
  // cout << "CalcGenoNum: partition = ";
  // DumpIntVec ( partition );
  genoNums.clear();
  for (int i = 0; i < genoLength; ++i) {
    genoNums.push_back(0);
  }
  for (int i = 0; i < partition.size(); ++i) {
    genoNums[partition[i] / 2]++;
  }
}

int Find2LabelOccNum(int lbl, const set<int> &setUniqeLables) {
  int res = 0;
  if (setUniqeLables.find(2 * lbl) != setUniqeLables.end()) {
    ++res;
  }
  if (setUniqeLables.find(2 * lbl + 1) != setUniqeLables.end()) {
    ++res;
  }
  return res;
}
