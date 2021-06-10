#include "Utils3.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

//////////////////////////////////////////////////////////////////////////////
// HashTable functions
//////////////////////////////////////////////////////////////////////////////

YWHashItem ::~YWHashItem() {}

YWHashTable ::YWHashTable(int nB) : numBuckets(nB) {}

YWHashTable ::~YWHashTable() {
  // NOTE: has to free memory here
  for (unsigned int i = 0; i < hashTable.size(); ++i) {
    delete hashTable[i];
  }
  hashTable.clear();
}

void YWHashTable ::AddItem(YWHashItem *pItem) { hashTable.push_back(pItem); }

YWHashItem *YWHashTable ::GetIdenticalItem(YWHashItem *pItem) {
  cout << "GetIdenticalItem: key = " << pItem->Key() << endl;
  for (unsigned int i = 0; i < hashTable.size(); ++i) {
    // cout << "We are here.\n";
    YW_ASSERT_INFO(hashTable[i] != NULL, "Can not be nothing here.");
    if (*hashTable[i] == *pItem) {
      cout << "find it here.\n";
      return hashTable[i];
    }
  }
  cout << "did not find.\n";
  return NULL;
}

YWHashItem *YWHashTable ::GetFirstItem() {
  cout << "GetFirstItem: size = " << hashTable.size() << endl;
  if (hashTable.size() == 0) {
    return NULL;
  }
  this->curPos = 0;
  return hashTable[0];
}

YWHashItem *YWHashTable ::GetNextItem() {
  cout << "GetNextItem: size = " << hashTable.size();
  cout << ", curPos = " << curPos << endl;
  if (this->curPos + 1 >= (int)hashTable.size()) {
    cout << "No more item.\n";
    return NULL;
  }
  this->curPos++;
  YWHashItem *pItem = hashTable[this->curPos];
  YW_ASSERT_INFO(pItem != NULL, "Can not be nothing.");
  cout << "GetNextItem.key() = " << pItem->Key() << endl;
  return pItem;
}

int YWHashTable ::GetTotalItemNum() const { return hashTable.size(); }

void YWHashTable ::Dump() const {
  for (unsigned int i = 0; i < hashTable.size(); ++i) {
    // cout << "We are here.\n";
    cout << "Key for item " << i << " = " << hashTable[i]->Key() << endl;
  }
}

//
bool SequenceCmp ::operator()(const SEQUENCE &seq1,
                              const SEQUENCE &seq2) const {
  if (seq1.size() != seq2.size()) {
    DumpSequence(seq1);
    DumpSequence(seq2);
  }

  YW_ASSERT_INFO(seq1.size() == seq2.size(),
                 "Can not compare two things with different length");

  for (int i = 0; i < (int)seq1.size(); ++i) {
    if (seq1[i] < seq2[i]) {
      return true;
    } else if (seq1[i] > seq2[i]) {
      return false;
    }
  }

  // if all are equal
  return false;
}

//////////////////////////////////////////////////////////////////////////////
// substring functions
//////////////////////////////////////////////////////////////////////////////
bool IsIntervalContained(const INTERVAL &iv1, const INTERVAL &iv2) {
  if ((iv1.first >= iv2.first && iv1.second <= iv2.second) ||
      (iv2.first >= iv1.first && iv2.second <= iv1.second)) {
    return true;
  }
  return false;
}
int GetIntervalLen(const INTERVAL &iv) { return iv.second - iv.first + 1; }

int GetRandItemInSet(const set<int> &items) {
  vector<int> itemsVec;
  PopulateVecBySet(itemsVec, items);
  return GetRandItemInVec(itemsVec);
}

int GetRandItemInVec(const vector<int> &items) {
  YW_ASSERT_INFO(items.size() > 0, "You can not sample from an empty set");

  double frac = GetRandFraction();
  return items[(int)(items.size() * frac)];
}

void GetRandVector(vector<int> &rndVec, int start, int end) {
  set<int> itemsNotUsed;
  PopulateSetWithInterval(itemsNotUsed, start, end);
  while (itemsNotUsed.size() > 0) {
    int itemRnd = GetRandItemInSet(itemsNotUsed);
    rndVec.push_back(itemRnd);
    itemsNotUsed.erase(itemRnd);
  }
}

int GetWeightedRandItemInVec(const vector<int> &items,
                             const vector<double> &itemWeights) {
  // cout << "items = ";
  // DumpIntVec( items );

  YW_ASSERT_INFO(items.size() == itemWeights.size(), "Size mismatch");
  double accum = 0.0;
  for (unsigned int i = 0; i < itemWeights.size(); ++i) {
    // cout << "one weight = " << itemWeights[i] << endl;
    accum += itemWeights[i];
  }
  YW_ASSERT_INFO(accum > 0.0000001, "2.Can not be too small");
  double frac = GetRandFraction();
  double curFract = 0.0;
  for (unsigned int i = 0; i < itemWeights.size(); ++i) {
    curFract += itemWeights[i] / accum;
    if (curFract >= frac) {
      return items[i];
    }
  }
  return -1; // should nothappen
}

// This functionreturn a weighted uniformly item index from the list
int GetWeightedRandItemIndex(const vector<double> &itemWeights) {
  double accum = 0.0;
  for (unsigned int i = 0; i < itemWeights.size(); ++i) {
    // cout << "one weight = " << itemWeights[i] << endl;
    accum += itemWeights[i];
  }
  // YW_ASSERT_INFO( accum > 0.0000001, "3. Can not be too small" );
  double frac = GetRandFraction();
  double curFract = 0.0;
  for (unsigned int i = 0; i < itemWeights.size(); ++i) {
    curFract += itemWeights[i] / accum;
    if (curFract >= frac) {
      return i;
    }
  }
  // Can not come here
  YW_ASSERT_INFO(false, "Something wrong here");
  return -1; // should nothappen
}

void GetOrigSubset(const vector<int> &origVec, const set<int> &subsetInd,
                   set<int> &subsetOrig) {
  subsetOrig.clear();
  for (set<int>::iterator it = subsetInd.begin(); it != subsetInd.end(); ++it) {
    YW_ASSERT_INFO(*it < (int)origVec.size(), "Size exceeds");
    subsetOrig.insert(origVec[*it]);
  }
}

void MutateSequenceAtSites(SEQUENCE &mutSeq, vector<int> &mutSites) {
  for (unsigned int p = 0; p < mutSites.size(); ++p) {
    MutateSeqAtSite(mutSeq, mutSites[p]);
  }
}

void DumpDoubleVec(const vector<double> &vecDoubles) {
  cout << "Double vector contains: ";
  for (unsigned int i = 0; i < vecDoubles.size(); ++i) {
    cout << vecDoubles[i] << ", ";
  }
  cout << endl;
}
void DumpDoubleVec(const vector<long double> &vecDoubles) {
  cout << "Double vector contains: ";
  for (unsigned int i = 0; i < vecDoubles.size(); ++i) {
    cout << vecDoubles[i] << ", ";
  }
  cout << endl;
}
void DumpBoolVec(const vector<bool> &vecBools) {
  cout << "Bool vector contains: ";
  for (unsigned int i = 0; i < vecBools.size(); ++i) {
    if (vecBools[i] == true) {
      cout << "1,";
    } else {
      cout << "0, ";
    }
  }
  cout << endl;
}

int GetLargestIndiceInDoubleVec(const vector<double> &vecDoubles) {
  YW_ASSERT_INFO(vecDoubles.size() > 0, "Can not have empty vec");
  double maxv = vecDoubles[0];
  int res = 0;
  for (unsigned int i = 0; i < vecDoubles.size(); ++i) {
    if (vecDoubles[i] > maxv) {
      maxv = vecDoubles[i];
      res = i;
    }
  }
  return res;
}

int GetLargestIndiceInDoubleVec(const vector<long double> &vecDoubles) {
  long double maxv = 0.0;
  int res = 0;
  for (unsigned int i = 0; i < vecDoubles.size(); ++i) {
    if (vecDoubles[i] > maxv) {
      maxv = vecDoubles[i];
      res = i;
    }
  }
  return res;
}

double FindMedian(const vector<double> &vecVals) {
  // for now, if there is nothing in the list, return 0.0
  if (vecVals.size() == 0) {
    return 0.0;
  }

  YW_ASSERT_INFO(vecVals.size() > 0, "FindMedian: Can not be empty");

  // Find median value for the vector
  // first sort the list of course
  vector<double> listToTry = vecVals;
  SortDoubleVec(listToTry);
  // now find the median one
  // int totSize = (int)listToTry.size();
  int pos = (int)(((int)listToTry.size() - 1) / 2);
  return listToTry[pos];
}

long double FindMedian(const vector<long double> &vecVals) {
  YW_ASSERT_INFO(vecVals.size() > 0, "FindMedian: Can not be empty");

  // Find median value for the vector
  // first sort the list of course
  vector<long double> listToTry = vecVals;
  SortDoubleVec(listToTry);
  // now find the median one
  // int totSize = (int)listToTry.size();
  int pos = (int)(((int)listToTry.size() - 1) / 2);
  return listToTry[pos];
}

double FindRankedItem(const vector<double> &vecVals, int rank) {
  YW_ASSERT_INFO(rank < (int)vecVals.size(), "Rank: overflow");
  vector<double> listToTry = vecVals;
  SortDoubleVec(listToTry);
  return listToTry[rank];
}

double FindMaxDouble(const vector<double> &vecVals) {
  // fnd max value of the solution here, assuming all values are non-negative
  vector<double> listToTry = vecVals;
  SortDoubleVec(listToTry);
  // cout << "vecVals = ";
  // DumpDoubleVec( vecVals );

  double res = listToTry[listToTry.size() - 1];
  // cout << "res = " << res << endl;
  return res;
}

double FindMaxDouble(const vector<long double> &vecVals) {
  // fnd max value of the solution here, assuming all values are non-negative
  vector<long double> listToTry = vecVals;
  SortDoubleVec(listToTry);
  // cout << "vecVals = ";
  // DumpDoubleVec( vecVals );

  double res = listToTry[listToTry.size() - 1];
  // cout << "res = " << res << endl;
  return res;
}

static int QSortCompareDouble(const void *arg1, const void *arg2) {
  /* Compare all of both strings: */
  // assume sorting in accending order
  double n1 = *((double *)arg1);
  double n2 = *((double *)arg2);
  // cout <<"arg1 = " << n1 << ", arg2 = " << n2 << endl;
  if (n1 > n2) {
    return 1;
  } else if (n1 < n2) {
    return -1;
  } else {
    return 0;
  }
}

void SortDoubleVec(vector<double> &vecVals, int start, int end) {
  //#if 0
  if (vecVals.size() <= 1) {
    // do nothing
    return;
  }
  // cout << "Before sort, double vec = ";
  // DumpDoubleVec( vecVals );
  if (end < 0) {
    end = vecVals.size() - 1;
  }
  int sortLen = end - start + 1;
  double *array = new double[sortLen];
  for (int i = start; i <= end; ++i) {
    array[i - start] = vecVals[i];
  }
  qsort((void *)array, sortLen, sizeof(double), QSortCompareDouble);
  // Now write back
  for (int i = start; i <= end; ++i) {
    vecVals[i] = array[i - start];
  }

  delete[] array;
  //#endif
  // cout << "After sort, double vec = ";
  // DumpDoubleVec( vecVals );
}

static int QSortCompareLongDouble(const void *arg1, const void *arg2) {
  /* Compare all of both strings: */
  // assume sorting in accending order
  long double n1 = *((long double *)arg1);
  long double n2 = *((long double *)arg2);
  // cout <<"arg1 = " << n1 << ", arg2 = " << n2 << endl;
  if (n1 > n2) {
    return 1;
  } else if (n1 < n2) {
    return -1;
  } else {
    return 0;
  }
}

void SortDoubleVec(vector<long double> &vecVals, int start, int end) {
  //#if 0
  if (vecVals.size() <= 1) {
    // do nothing
    return;
  }
  // cout << "Before sort, double vec = ";
  // DumpDoubleVec( vecVals );
  if (end < 0) {
    end = vecVals.size() - 1;
  }
  int sortLen = end - start + 1;
  long double *array = new long double[sortLen];
  for (int i = start; i <= end; ++i) {
    array[i - start] = vecVals[i];
  }
  qsort((void *)array, sortLen, sizeof(long double), QSortCompareLongDouble);
  // Now write back
  for (int i = start; i <= end; ++i) {
    vecVals[i] = array[i - start];
  }

  delete[] array;
  //#endif
  // cout << "After sort, double vec = ";
  // DumpDoubleVec( vecVals );
}

void FindUniformColumns(const vector<SEQUENCE> &listSeqs, set<int> &uniSites) {
  uniSites.clear();
  if (listSeqs.size() == 0) {
    return;
  }
  int numSites = (int)listSeqs[0].size();
  for (int i = 0; i < numSites; ++i) {
    bool f0 = false, f1 = false;
    for (int r = 0; r < (int)listSeqs.size(); ++r) {
      if (listSeqs[r][i] == 0) {
        f0 = true;
      } else if (listSeqs[r][i] == 1) {
        f1 = true;
      }
      if (f0 == true && f1 == true) {
        // not uniform,
        break;
      }
    }
    if (f0 == false || f1 == false) {
      // yes, this site is uniform
      uniSites.insert(i);
    }
  }
}

void BreakSeqAtBkpt(const SEQUENCE &seq, int bkpt, SEQUENCE &seqLeft,
                    SEQUENCE &seqRight) {
  seqLeft.clear();
  seqRight.clear();
  for (int i = 0; i < (int)seq.size(); ++i) {
    if (i <= bkpt) {
      // then the right seq get MV
      seqLeft.push_back(seq[i]);
      seqRight.push_back(MISSING_VALUE_BIT);
    } else {
      seqLeft.push_back(MISSING_VALUE_BIT);
      seqRight.push_back(seq[i]);
    }
  }
}

bool AreTwoSeqsBroken(const SEQUENCE &seqLeft, const SEQUENCE &seqRight) {
  // test whether the two sequences are broken from a single sequence
  // to avoid duplicate events mainly
  bool foundBkpt = false;
  if (seqLeft.size() != seqRight.size()) {
    return false;
  }

  for (int i = 0; i < (int)seqLeft.size(); ++i) {
    if (IsMissingValueBit(seqLeft[i]) == false &&
        IsMissingValueBit(seqRight[i]) == false) {
      return false; // no, not a broken seqs pair
    }

    if (IsMissingValueBit(seqRight[i]) == false) {
      if (foundBkpt == false) {
        foundBkpt = true;
      }
    }
    if (foundBkpt == true && IsMissingValueBit(seqLeft[i]) == false) {
      return false;
    }
  }
  return true;
}

// new stuff from treeHMM

bool GetFirstMutliChoice(int numStage, int numStageElem,
                         vector<int> &initChoice) {
  if (numStage <= 0 || numStageElem <= 0) {
    return false;
  }
  initChoice.clear();
  // Start by picking first one each time
  for (int i = 0; i < numStage; ++i) {
    initChoice.push_back(0);
  }
  return true;
}

bool GetNextMutliChoice(int numStage, int numStageElem,
                        vector<int> &indChoice) {
  // Now we move to next choice
  // bool res = false;
  // Find the last item not = numStageElem-1
  int itemToChange = -1;
  for (int i = ((int)indChoice.size()) - 1; i >= 0; --i) {
    if (indChoice[i] < numStageElem - 1) {
      itemToChange = i;
      break;
    }
  }
  if (itemToChange < 0) {
    // No solution
    return false;
  }
  // Now we clear out everything beyond it
  for (int i = itemToChange + 1; i < (int)indChoice.size(); ++i) {
    indChoice[i] = 0;
  }
  indChoice[itemToChange]++;
  return true;
}

// void DumpVecSequences( const vector<SEQUENCE> &vecSeqs )
//{
//    cout << "Vector of sequneces = \n";
//    for( unsigned int i=0; i<vecSeqs.size(); ++i )
//    {
//        DumpSequence( vecSeqs[i] );
//    }
//}

void GetVecSequencesIV(const vector<SEQUENCE> &vecSeqs, int left, int right,
                       vector<SEQUENCE> &vecSeqsIV) {
  vecSeqsIV.clear();
  for (unsigned int i = 0; i < vecSeqs.size(); ++i) {
    SEQUENCE ivRow;
    GetSeqInterval(vecSeqs[i], ivRow, left, right);
    vecSeqsIV.push_back(ivRow);
  }
}

int GetNumZerosInSeq(const SEQUENCE &seq) {
  int res = 0;
  for (unsigned int i = 0; i < seq.size(); ++i) {
    if (seq[i] == 0) {
      res++;
    }
  }
  return res;
}

void GetSeqSplit(const SEQUENCE &seq, set<int> &zeroBits, set<int> &oneBits) {
  zeroBits.clear();
  oneBits.clear();
  for (unsigned int i = 0; i < seq.size(); ++i) {
    if (seq[i] == 0) {
      zeroBits.insert(i);
    } else if (seq[i] == 1) // need to enforce this due to potential missing
                            // value, 7/4/08
    {
      oneBits.insert(i);
    }
  }
}

static int QSortCompareIntPair(const void *arg1, const void *arg2) {
  /* Compare all of both strings: */
  // assume sorting in accending order, and use the first value in the int pair
  // to sort
  int n1 = ((pair<int, int> *)arg1)->first;
  int n2 = ((pair<int, int> *)arg2)->first;
  // cout <<"arg1 = " << n1 << ", arg2 = " << n2 << endl;
  if (n1 > n2) {
    return 1;
  } else if (n1 < n2) {
    return -1;
  } else {
    return 0;
  }
}

void SortVecIntPairs(vector<pair<int, int> > &listPairs) {
  pair<int, int> *parray = new pair<int, int>[listPairs.size()];
  for (int i = 0; i < (int)listPairs.size(); ++i) {
    parray[i] = listPairs[i];
  }
  qsort((void *)parray, listPairs.size(), sizeof(pair<int, int>),
        QSortCompareIntPair);
  // Now write back
  for (int i = 0; i < (int)listPairs.size(); ++i) {
    listPairs[i] = parray[i];
  }

  delete[] parray;
}

////////////////////////////////////////////////////////////////////////////////
int GetSubstringLeftPos(const INTERVAL_SUBSTRING &substr) {
  return substr.first.first;
}
int GetSubstringRightPos(const INTERVAL_SUBSTRING &substr) {
  return substr.first.second;
}
void GetIVSubstringData(const INTERVAL_SUBSTRING &substr, SEQUENCE &seq) {
  seq = substr.second;
}
INTERVAL GetSubstringInterval(const INTERVAL_SUBSTRING &substr) {
  return substr.first;
}
bool GetSubstringSegment(const INTERVAL_SUBSTRING &substr,
                         const INTERVAL &ivToRead, SEQUENCE &segment) {
  YW_ASSERT_INFO(IsIntervalContained(ivToRead, substr.first) == true,
                 "Two intervals do not have contained");

  // remember we have to offset a little
  int startPos = GetSubstringLeftPos(substr);
  GetSeqInterval(substr.second, segment, ivToRead.first - startPos,
                 ivToRead.second - startPos);
  return true;
}

int GetSubstringValAt(const INTERVAL_SUBSTRING &substr, int pos) {
  YW_ASSERT_INFO(pos >= GetSubstringLeftPos(substr) &&
                     pos <= GetSubstringRightPos(substr),
                 "Range error.");

  int convPos = pos - GetSubstringLeftPos(substr);
  return substr.second[convPos];
}

bool IsSegmentContained(const INTERVAL_SUBSTRING &seqContained,
                        const INTERVAL_SUBSTRING &seqContainer) {
  // First the range has to match
  if (GetSubstringLeftPos(seqContained) < GetSubstringLeftPos(seqContainer) ||
      GetSubstringRightPos(seqContained) > GetSubstringRightPos(seqContainer)) {
    return false;
  }
  // Then the corresponding position must match too
  for (int p = GetSubstringLeftPos(seqContained);
       p <= GetSubstringRightPos(seqContained); p++) {
    if (GetSubstringValAt(seqContained, p) !=
        GetSubstringValAt(seqContainer, p)) {
      return false;
    }
  }
  return true;
}

bool AreSegmentsConsistent(const INTERVAL_SUBSTRING &seq1,
                           const INTERVAL_SUBSTRING &seq2) {
  // If disjoint, yes, it is consistent
  INTERVAL ivInt;
  bool fInt = GetIntervalOverlap(GetSubstringInterval(seq1),
                                 GetSubstringInterval(seq2), ivInt);
  if (fInt == false) {
    return true;
  }
  // cout << "ivInt.first = " << ivInt.first << ", ivInt.second = " <<
  // ivInt.second << endl;
  // make sure the two things matches
  SEQUENCE seqp1;
  GetSubstringSegment(seq1, ivInt, seqp1);
  // cout << "seqp1 = ";
  // DumpSequence( seqp1 );
  SEQUENCE seqp2;
  GetSubstringSegment(seq2, ivInt, seqp2);
  // cout << "seqp2 = ";
  // DumpSequence( seqp2 );

  if (seqp1 == seqp2) {
    return true;
  } else {
    return false;
  }
}

int GetSegmentsIntersection(const INTERVAL_SUBSTRING &seq1,
                            const INTERVAL_SUBSTRING &seq2, INTERVAL &iv) {
  // we simply get how larget the intersection from the interval ONLY

  bool fInt = GetIntervalOverlap(GetSubstringInterval(seq1),
                                 GetSubstringInterval(seq2), iv);
  if (fInt == false) {
    return 0;
  }
  return iv.second - iv.first + 1;
}

bool AreSegmentsNextto(const INTERVAL_SUBSTRING &seq1,
                       const INTERVAL_SUBSTRING &seq2) {
  // cout << "seq1.left = " <<  GetSubstringLeftPos(seq1) << ", right = " <<
  // GetSubstringRightPos(seq1) << endl; cout << "seq2.left = " <<
  // GetSubstringLeftPos(seq2) << ", right = " <<  GetSubstringRightPos(seq2) <<
  // endl;
  // Two segments are next to each other if the can form a single bigger
  // ungapped piece
  if (GetSubstringLeftPos(seq1) == GetSubstringRightPos(seq2) + 1 ||
      GetSubstringLeftPos(seq2) == GetSubstringRightPos(seq1) + 1) {
    // cout << "Yes, neighbours.\n";
    return true;
  } else {
    // cout << "No, not neighbours.\n";
    return false;
  }
}

void DumpSubstring(const INTERVAL_SUBSTRING &substr) {
  cout << "[" << GetSubstringLeftPos(substr) << ",";
  cout << GetSubstringRightPos(substr) << "], ";
  DumpSequence(substr.second);
}

// ***************************************************************************
// Numerical utilities
// ***************************************************************************
double GetLogSumOfLogs(const vector<double> &listLogs) {
  if (listLogs.size() == 0) {
    // nothing to process
    return 0.0;
  }
  // given a list of log terms, compute the sum of prob (need to take exp)
  // and express the sum in the log again
  // first get the largest term and use it as a base
  int posmax = GetLargestIndiceInDoubleVec(listLogs);
  double valmax = listLogs[posmax];
  double asum = 0.0;
  for (int i = 0; i < (int)listLogs.size(); ++i) {
    asum += exp(listLogs[i] - valmax);
  }
  double res = valmax + log(asum);
  // cout << "res = " << res << ", valmax = " << valmax << ", in list: ";
  // DumpDoubleVec(listLogs);
  // cout << "Direct evaluation = " << GetLogSumOfLogsDirect(listLogs) << endl;
  return res;
}

double GetLogSumOfLogsDirect(const vector<double> &listLogs) {
  // simply just direct sum over
  double asum = 0.0;
  for (int i = 0; i < (int)listLogs.size(); ++i) {
    asum += exp(listLogs[i]);
  }
  return log(asum);
}

double GetLogSumOfTwo(double logv1, double logv2) {
  vector<double> vecVals;
  vecVals.push_back(logv1);
  vecVals.push_back(logv2);
  return GetLogSumOfLogs(vecVals);
}

void SumofLogVecs(vector<double> &listLogsAdded,
                  vector<double> &listLogsAdding) {
  YW_ASSERT_INFO(listLogsAdded.size() == listLogsAdding.size(),
                 "Must have the same length");
  for (int i = 0; i < (int)listLogsAdded.size(); ++i) {
    listLogsAdded[i] = GetLogSumOfTwo(listLogsAdded[i], listLogsAdding[i]);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////
// More useful functions

// This is a very useful function, so expose it
int FindMatchedSeqForFounders(const vector<SEQUENCE> &founder,
                              const SEQUENCE &seq, set<int> &endRows,
                              bool fPrefix) {
  // Return the number of crossovers
  // This function computes the minimum recombination weight for the given
  // hapRow when restricted to interval [left, right] in mat
  int res = 0;

  set<int> lastTrackRows; // set of rows that matching the hapRow

  // Ohterwise, we can start from all possible rows
  for (unsigned int i = 0; i < founder.size(); ++i) {
    lastTrackRows.insert(i);
  }

  int curpos = 0;
  int end = seq.size();
  if (fPrefix == false) {
    curpos = seq.size() - 1;
    end = -1;
  }

  while (curpos != end) {
    // Each time, we intersect the set with the sets matching the current bit
    set<int> trackRows;
    for (unsigned int i = 0; i < founder.size(); ++i) {
      if (IsTwoStatesCompatible(founder[i][curpos], seq[curpos]) == true) {
        // Yes, this row matches
        trackRows.insert(i);
      }
    }

    // Now we test if there is intersection, if non-empty, we contiinue
    set<int> sint;
    JoinSets(trackRows, lastTrackRows, sint);
    if (sint.size() == 0) {
      break;
    } else {
      // In this case, we still continue
      lastTrackRows = sint;
    }

    if (fPrefix == true) {
      curpos++;
    } else {
      curpos--;
    }
  }

  endRows = lastTrackRows;

  // what is the length of the prefix/suffix
  if (fPrefix) {
    res = curpos;
  } else {
    res = seq.size() - 1 - curpos;
  }

  return res;
}

int FindNoninformativeRow(const vector<SEQUENCE> &listSeqs, int col) {
  int numZeros = 0, numOnes = 0, numMissing = 0;
  // now we compare these two cols: c1, c2
  // if they match, we put c2 into set
  int res0 = -1, res1 = -1;
  for (unsigned int r = 0; r < listSeqs.size(); ++r) {
    if (listSeqs[r][col] == 0) {
      numZeros++;
      res0 = r;
    } else if (listSeqs[r][col] == 1) {
      numOnes++;
      res1 = r;
    } else if (IsMissingValueBit(listSeqs[r][col]) == true) {
      numMissing++;
    }
    if (numZeros > 1 && numOnes > 1) {
      return -1; // no such row
    }
  }

  // Check to see if this is non-informative
  if (numZeros == 1 && numOnes >= 1) {
    // we find a duplicate
    //			cout << "Site  " << c1+1 << "is  non-informative" <<
    // endl;
    return res0;
  } else if (numOnes == 1 && numZeros >= 1) {
    return res1;
  } else {
    return -1;
  }
}

void ConvVecToArray(const vector<int> &vec, int *arr) {
  // IMPORTANT: ASSUME ARR HAS BEEN ALLOCATED TO PROPER SIZE!!
  for (int i = 0; i < (int)vec.size(); ++i) {
    arr[i] = vec[i];
  }
}
void ConvVecToArray(const vector<double> &vec, double *arr) {
  // IMPORTANT: ASSUME ARR HAS BEEN ALLOCATED TO PROPER SIZE!!
  for (int i = 0; i < (int)vec.size(); ++i) {
    arr[i] = vec[i];
  }
}
void DumpIntArray(int len, int *arr) {
  for (int i = 0; i < len; ++i) {
    cout << arr[i];
    if (i < len - 1) {
      cout << ", ";
    }
  }
  cout << endl;
}
void FlipBinVector(vector<int> &vec) {
  for (int i = 0; i < (int)vec.size(); ++i) {
    if (vec[i] == 0) {
      vec[i] = 1;
    } else {
      vec[i] = 0;
    }
  }
}

void RecoverOrigIndicesAfterDeletion(const vector<int> &removedItems,
                                     const vector<int> &itemsNew,
                                     vector<int> &itemsOrigIndices) {
  // this function is for reconstructing the orignal indices of items
  // after some items were deleted from array, and the passed-in indices are for
  // the NEW positions. We are interested to know the origianl positions
  itemsOrigIndices.clear();

  // first sort the two arrays
  vector<int> removedItemsUse = removedItems;
  vector<int> itemsNewUse = itemsNew;
  SortIntVec(removedItemsUse);
  SortIntVec(itemsNewUse);

  int posNew = 0;
  for (int i = 0; i < (int)removedItemsUse.size(); ++i) {
    //
    // int posDel = removedItemsUse[i];

    // output anything that is smaller or equal to this number
    while (posNew < (int)itemsNewUse.size() &&
           itemsNewUse[posNew] < removedItemsUse[i] - i) {
      // convert it
      itemsOrigIndices.push_back(itemsNewUse[posNew] + i);
      posNew++;
    }

    // stop if nothing left
    if (posNew >= (int)itemsNewUse.size()) {
      break;
    }
  }
  // also output things left over
  for (; posNew < (int)itemsNewUse.size(); ++posNew) {
    //
    itemsOrigIndices.push_back(itemsNewUse[posNew] + removedItemsUse.size());
  }
  // cout << "removedItems = ";
  // DumpIntVec( removedItems );
  // cout << "cur items = ";
  // DumpIntVec(itemsNew);
  // cout << "Converted items = ";
  // DumpIntVec(  itemsOrigIndices );
}

void GetOrigPositionAfterRemoval(int numRemains,
                                 const vector<int> &itemsRemoved,
                                 vector<int> &origPosForRemains) {
  // for now, choose a simple but NOT EFFICIENT way. TBD
  // try to get original positions of the item removal from a list
  // for example, say 3 items remains and itemsREmoved = 1,2 (0-based), then
  // orig pos for remaings = 0, 3,4
  set<int> setItemRemoved;
  PopulateSetByVec(setItemRemoved, itemsRemoved);
  set<int> setItemsOrig;
  PopulateSetWithInterval(setItemsOrig, 0,
                          numRemains + (int)itemsRemoved.size());
  // substract something
  SubtractSets(setItemsOrig, setItemRemoved);
  // now result
  PopulateVecBySet(origPosForRemains, setItemsOrig);

  // for(int i=0; i<numRemains + (int)itemsRemoved.size(); ++i)
  //{
  //		if()
  //	{
  //	}
  //}
}

void ConvOneSideToFullSplit(vector<int> &split, const set<int> &oneside,
                            int numLeaves, int val) {
  split.resize(numLeaves);
  int val0 = 0;
  if (val == 0) {
    val0 = 1;
  }
  for (int i = 0; i < numLeaves; ++i) {
    split[i] = val0;
  }
  for (set<int>::iterator it = oneside.begin(); it != oneside.end(); ++it) {
    split[*it] = val;
  }
}

bool AreTwoMVVecCompat(const vector<int> &vec1, const vector<int> &vec2,
                       int &numTrueMatch) {
  YW_ASSERT_INFO(vec1.size() == vec2.size(), "Fail");
  numTrueMatch = 0;
  int mres = 0;
  for (int i = 0; i < (int)vec1.size(); ++i) {
    if (IsMissingValueBit(vec1[i]) == true ||
        IsMissingValueBit(vec2[i]) == true) {
      // match
      continue;
    } else if (vec1[i] != vec2[i]) {
      return false;
    } else {
      // true match
      mres++;
    }
  }
  numTrueMatch = mres;
  return true;
}

int GetMVNum(const vector<int> &vec) {
  int res = 0;
  for (int i = 0; i < (int)vec.size(); ++i) {
    if (IsMissingValueBit(vec[i]) == true) {
      res++;
    }
  }
  return res;
}

bool AreSeqsOverlap(const vector<int> &vec1, const vector<int> &vec2) {
  for (int i = 0; i < (int)vec1.size(); ++i) {
    if (IsMissingValueBit(vec1[i]) == false &&
        IsMissingValueBit(vec2[i]) == false) {
      return true;
    }
  }
  return false;
}

void InsertOrderedVec(vector<int> &vec, int val) {
#if 0
	// assume vec is already ordered and we will add a new val to keep vec ordered
	// IMPORTANT: remove duplicate copy if any
	vector<int> vecRes;
	int pos = 0;
	for(;  pos  < (int) vec.size() && vec[pos] < val; ++pos)
	{
		vecRes.push_back(  vec[pos] );
	}
	if( pos >= (int) vec.size()  || val != vec[pos]  )
	{
		vecRes.push_back( val );
	}
	// now get the remaining
	for(; pos  < (int) vec.size(); ++pos)
	{
		vecRes.push_back(  vec[pos] );
	}
	vec = vecRes;
#endif

  if (vec.size() == 0) {
    vec.push_back(val);
    return;
  }

  // cout << "In InsertOrderedVec: val = " << val << ", vec = ";
  // DumpIntVec( vec );
  // want to insert the item in space
  // first find the location to add this item
  // doing it in binary search
  int pos = binary_search<int>(vec, 0, vec.size() - 1, val);
  YW_ASSERT_INFO(pos >= 0, "Wrong in binary search");
  // cout << "pos = " << pos << endl;
  if (pos >= (int)vec.size() || val != vec[pos]) {
    // need to add this item in. First shift one item to the right
    vec.push_back(0);
    for (int i = (int)vec.size() - 2; i >= pos; --i) {
      vec[i + 1] = vec[i];
    }
    vec[pos] = val;
  }
}

//! \brief A recursive binary search using STL vectors
//! \param vec The vector whose elements are to be searched
//! \param start The index of the first element in the vector
//! \param end One past the index of the last element in the vector
//! \param key The value being searched for
//! \return The index into the vector where the value is located,
//! or -1 if the value could not be found.
template <typename T>
int binary_search(const std::vector<T> &vec, unsigned start, unsigned end,
                  const T &key) {
  // Termination condition: start index greater than end index
  if (start > end) {
    return start;
  }

  // Find the middle element of the vector and use that for splitting
  // the array into two pieces.
  unsigned middle = (start + ((end - start) / 2));

  if (vec[middle] == key) {
    return middle;
  } else if (vec[middle] > key) {
    return binary_search(vec, start, middle - 1, key);
  }

  return binary_search(vec, middle + 1, end, key);
}

bool ReadIntListFromFile(const char *fname, vector<int> &listInts) {
  // data input
  ifstream inFile(fname);
  if (!inFile) {
    cout << "Can not open " << fname << endl;
    return false;
  }
  listInts.clear();
  while (inFile.eof() == false) {
    const int BUF_SZ = 102400;
    char buffer[BUF_SZ];
    inFile.getline(buffer, BUF_SZ);
    if (strlen(buffer) > 0) {
      // cout << "buffer = " << buffer << endl;
      int val;
      sscanf(buffer, "%d", &val);
      listInts.push_back(val);
    }
  }

  return true;
}

void GetVecPosNotInSet(const vector<int> &vec, const set<int> &s,
                       vector<int> &posDiff) {
  posDiff.clear();
  //
  for (int i = 0; i < (int)vec.size(); ++i) {
    if (s.find(vec[i]) == s.end()) {
      posDiff.push_back(i);
    }
  }
}

// Suppose we have g groups of (indistingishable) items and we want to
// divide each group into numParts colors (distinguishable)
// this support enumerate these choices. For example, we have two segments of 3
// and 4 items each and we have two colors, then the choices will be: [(1,2),
// (2,2)], or [(0.3),(1,3)]
void InitPartitionEnum(const vector<int> &vecSegSizes, int numParts,
                       vector<vector<int> > &parts) {
  // start from each one as the first population type has all the ones in
  // segment
  parts.clear();
  parts.resize(vecSegSizes.size());
  for (int i = 0; i < (int)vecSegSizes.size(); ++i) {
    parts[i].push_back(vecSegSizes[i]);
    for (int j = 1; j < numParts; ++j) {
      parts[i].push_back(0);
    }
    // cout << "InitPartitionEnum: part = ";
    // DumpIntVec(parts[i]);
  }
}

bool GetNextPartitionEnum(const vector<int> &vecSegSizes, int numParts,
                          vector<vector<int> > &parts) {
  // cout << "GetNextPartitionEnum: numParts = " << numParts << ", vecSegSizes =
  // "; DumpIntVec(vecSegSizes);
  // get next partition, return false if done
  // first search for the part where we can change (by moving some item to the
  // front)
  YW_ASSERT_INFO(parts.size() == vecSegSizes.size(),
                 "GetNextPartitionEnum: size mismatch");
  int segChange = -1;
  for (int seg = 0; seg < (int)vecSegSizes.size(); ++seg) {
    YW_ASSERT_INFO((int)parts[seg].size() == numParts,
                   "GetNextPartitionEnum: seg size mismatch");
    // when the part has concerntrated to the last population, this is the sign
    // that this part has changed it partiton to its limit
    if (parts[seg][numParts - 1] != vecSegSizes[seg]) {
      segChange = seg;
      break;
    }
  }
  if (segChange < 0) {
    // done
    return false;
  }
  // cout << "segChange = " << segChange << endl;
  //
  vector<vector<int> > partsNew = parts;
  // the first segments before this seg is re-set
  for (int s = 0; s < segChange; ++s) {
    partsNew[s][0] = vecSegSizes[s];
    for (int j = 1; j < numParts; ++j) {
      partsNew[s][j] = 0;
    }
  }
  // then segChange one gets shift by one
  // this is done by finding the least numbered population and
  // move it out to one larger AND concerntrate all the ones up to this point to
  // the first position
  int pp = -1;
  // int numItemsToi = 0;
  for (int i = 0; i < numParts; ++i) {
    if (parts[segChange][i] > 0) {
      pp = i;
      break;
    }
  }
  // cout << "pp = " << pp << endl;
  YW_ASSERT_INFO(pp >= 0 && pp < numParts - 1, "Can not be true");
  vector<int> segNew = parts[segChange];
  segNew[0] = parts[segChange][pp] - 1;
  if (pp != 0) {
    segNew[pp] = 0;
  }
  segNew[pp + 1]++;

  partsNew[segChange] = segNew;

  // the rest remain the same
  parts = partsNew;
  // cout << "Next parts id = \n";
  // for(int i=0;i<(int)parts.size(); ++i)
  //{
  // DumpIntVec( parts[i] );
  //}
  return true;
}

int GetPartEnumIndex(const vector<int> &vecSegSizes, int numParts,
                     const vector<vector<int> > &parts) {
  // get the index (order in the enumerated list) of the given enumerated
  // partition cout from the right hand side
  YW_ASSERT_INFO(vecSegSizes.size() == parts.size(),
                 "GetPartEnumIndex: size wrong");
  int res = 0;
  for (int i = (int)vecSegSizes.size() - 1; i >= 0; --i) {
    if (i < (int)vecSegSizes.size() - 1) {
      res *= GetPartitionEnumNum(vecSegSizes[i], numParts);
    }
    res += GetPartitionEnumId(vecSegSizes[i], parts[i]);
  }
  return res;
}

// Now allow chaing parts num
void InitPartitionEnumVar(const vector<int> &vecSegSizes,
                          const vector<int> &listNumParts,
                          vector<vector<int> > &parts) {
  // start from each one as the first population type has all the ones in
  // segment
  YW_ASSERT_INFO(vecSegSizes.size() == listNumParts.size(), "Mismatch");
  parts.clear();
  parts.resize(vecSegSizes.size());
  for (int i = 0; i < (int)vecSegSizes.size(); ++i) {
    parts[i].push_back(vecSegSizes[i]);
    for (int j = 1; j < listNumParts[i]; ++j) {
      parts[i].push_back(0);
    }
    // cout << "InitPartitionEnum: part = ";
    // DumpIntVec(parts[i]);
  }
}

bool GetNextPartitionEnumVar(const vector<int> &vecSegSizes,
                             const vector<int> &listNumParts,
                             vector<vector<int> > &parts) {
  // cout << "GetNextPartitionEnumVar: vecSegSizes = ";
  // DumpIntVec(vecSegSizes);
  // cout << "listNumparts: ";
  // DumpIntVec(listNumParts);
  // cout << "parts: ";
  // DumpVecSequences(parts);
  YW_ASSERT_INFO(vecSegSizes.size() == listNumParts.size(), "Mismatch");
  // cout << "GetNextPartitionEnum: numParts = " << numParts << ", vecSegSizes =
  // "; DumpIntVec(vecSegSizes);
  // get next partition, return false if done
  // first search for the part where we can change (by moving some item to the
  // front)
  YW_ASSERT_INFO(parts.size() == vecSegSizes.size(),
                 "GetNextPartitionEnum: size mismatch");
  int segChange = -1;
  for (int seg = 0; seg < (int)vecSegSizes.size(); ++seg) {
    YW_ASSERT_INFO((int)parts[seg].size() == listNumParts[seg],
                   "GetNextPartitionEnum: seg size mismatch");
    // when the part has concerntrated to the last population, this is the sign
    // that this part has changed it partiton to its limit
    if (parts[seg][listNumParts[seg] - 1] != vecSegSizes[seg]) {
      segChange = seg;
      break;
    }
  }
  if (segChange < 0) {
    // done
    // cout << "Done\n";
    return false;
  }
  // cout << "segChange = " << segChange << endl;
  //
  vector<vector<int> > partsNew = parts;
  // the first segments before this seg is re-set
  for (int s = 0; s < segChange; ++s) {
    partsNew[s][0] = vecSegSizes[s];
    for (int j = 1; j < listNumParts[s]; ++j) {
      partsNew[s][j] = 0;
    }
  }
  // then segChange one gets shift by one
  // this is done by finding the least numbered population and
  // move it out to one larger AND concerntrate all the ones up to this point to
  // the first position
  int pp = -1;
  // int numItemsToi = 0;
  for (int i = 0; i < listNumParts[segChange]; ++i) {
    if (parts[segChange][i] > 0) {
      pp = i;
      break;
    }
  }
  // cout << "pp = " << pp << endl;
  YW_ASSERT_INFO(pp >= 0 && pp < listNumParts[segChange] - 1,
                 "Can not be true");
  vector<int> segNew = parts[segChange];
  segNew[0] = parts[segChange][pp] - 1;
  if (pp != 0) {
    segNew[pp] = 0;
  }
  segNew[pp + 1]++;

  partsNew[segChange] = segNew;

  // the rest remain the same
  parts = partsNew;
  // cout << "Next parts id = \n";
  // for(int i=0;i<(int)parts.size(); ++i)
  //{
  // DumpIntVec( parts[i] );
  //}
  return true;
}

int GetPartEnumIndexVar(const vector<int> &vecSegSizes,
                        const vector<int> &listNumParts,
                        const vector<vector<int> > &parts) {
  YW_ASSERT_INFO(vecSegSizes.size() == listNumParts.size(), "Mismatch");
  // get the index (order in the enumerated list) of the given enumerated
  // partition cout from the right hand side
  YW_ASSERT_INFO(vecSegSizes.size() == parts.size(),
                 "GetPartEnumIndex: size wrong");
  int res = 0;
  for (int i = (int)vecSegSizes.size() - 1; i >= 0; --i) {
    if (i < (int)vecSegSizes.size() - 1) {
      res *= GetPartitionEnumNum(vecSegSizes[i], listNumParts[i]);
    }
    res += GetPartitionEnumId(vecSegSizes[i], parts[i]);
  }
  return res;
}

// **************************************************************************************
// code for enumearing partitions (based on Ruhua's code)
// hereis the pre-initied enumeration, format: <num of lins, num of color>,
// enumeration
static map<pair<int, int>, vector_vector_t> mapEnumeratedPartitions;

int GetPartitionEnumNum(int n, int numSPop) {
  // cout << "GetPartitionEnumNum: n = " << n << ", numSPop = " << numSPop;
  if (numSPop == 0) {
    return 0;
  }
  // how many number of partitons of identical balls into p colors
  double resd = 1.0;
  for (int j = 1; j <= numSPop - 1; ++j) {
    resd *= (1.0 * (n + numSPop - j)) / j;
  }
  int res = (int)(resd);
  // cout << ", res = " << res << endl;
  return res;
}

int GetPartitionEnumId(int numItemsTot, const vector<int> &vec) {
  int numColor = vec.size();
  // cout << "numItesmTotl: " << numItemsTot << ", numColor = " << numColor <<
  // ", vec = "; DumpIntVec( vec );
  pair<int, int> pp(numItemsTot, numColor);
  bool fExist =
      mapEnumeratedPartitions.find(pp) != mapEnumeratedPartitions.end();
  if (fExist == false) {
    vector_vector_t tt;
    mapEnumeratedPartitions.insert(
        map<pair<int, int>, vector_vector_t>::value_type(pp, tt));
  }
  int res = -1;
  convert_vector_to_index(fExist, vec, res, mapEnumeratedPartitions[pp]);
  YW_ASSERT_INFO(res >= 0, "Fail in GetPartitioId");
  // cout << "parition id: " << res  << ",numItemsTot: " << numItemsTot << ",
  // vec = "; DumpIntVec( vec );
  return res;

#if 0
	// for this enumerated vector, where does it stand in the enumeration order?
	// find the rightmost non-zero pos
	int posRight = -1;
	for(int p=(int)vec.size()-1; p>=0; --p)
	{
		if( vec[p] > 0 )
		{
			posRight = p;
			break;
		}
	}
	YW_ASSERT_INFO( posRight >= 0, "GetPartitionEnumId: Fail" );
//cout << "posRight = " << posRight << endl;
	// do it recurisvely: breakup into multiple classes: the one with one fewr number
	// and those with this bit removed
	if( posRight == 0 && vec[0] == numItemsTot )
	{
		// this is the first
//cout << "res = 0\n";
		return 0;
	}
	else
	{
		// this shows how many enum possible when the value of the highest order bit is DIFFERENT (ie smaller)
		int numDownShift = 0;
		for(int j=1; j<=vec[posRight]; ++j)
		{
			numDownShift += GetPartitionEnumNum( numItemsTot - vec[posRight] + j, posRight);
		}
		// this consider the situation when the highest bit is the SAME
		vector<int> vecDoneShift;
		for(int j = 0; j<posRight; ++j)
		{
			vecDoneShift.push_back( vec[j] );
		}
		int numSameOrder;
		if( numItemsTot > vec[posRight] )
		{
			numSameOrder = GetPartitionEnumId(numItemsTot-vec[posRight], vecDoneShift);
		}
		else
		{
			numSameOrder = 0;
		}
        int res = numDownShift + numSameOrder;
//cout << "numDownShift : " << numDownShift << ", numSameOrder: " << numSameOrder << ", res = " << res << endl;

		return res;
	}
#endif
}

void GetPartitionEnumPartForId(int numItemsTot, int numParts, int eid,
                               vector<int> &vecres) {
  pair<int, int> pp(numItemsTot, numParts);
  bool fExist =
      mapEnumeratedPartitions.find(pp) != mapEnumeratedPartitions.end();
  if (fExist == false) {
    vector_vector_t tt;
    mapEnumeratedPartitions.insert(
        map<pair<int, int>, vector_vector_t>::value_type(pp, tt));
  }
  convert_index_to_vector(fExist, numParts, numItemsTot, eid, vecres,
                          mapEnumeratedPartitions[pp]);
  // YW_ASSERT_INFO(vecres.size() >= 0, "Fail in GetPartitionEnumPartForId");

  // cout << "ConvPartition to id: parition id: " << eid  << ",numItemsTot: " <<
  // numItemsTot << ", numParts: " << numParts << ", vecres = "; DumpIntVec(
  // vecres );

#if 0

    // YW: TMEP
    //TBD


    //
    //cout << "numItemsTot: " << numItemsTot << ", vec = ";
    //DumpIntVec( vec );
	// for this enumerated vector, where does it stand in the enumeration order?
	// find the rightmost non-zero pos
	int posRight = -1;
	for(int p=(int)vec.size()-1; p>=0; --p)
	{
		if( vec[p] > 0 )
		{
			posRight = p;
			break;
		}
	}
	YW_ASSERT_INFO( posRight >= 0, "GetPartitionEnumId: Fail" );
    //cout << "posRight = " << posRight << endl;
	// do it recurisvely: breakup into multiple classes: the one with one fewr number
	// and those with this bit removed
	if( posRight == 0 )
	{
		// this is the first
		//return 0;
	}
	else
	{
		// this shows how many enum possible when the value of the highest order bit is DIFFERENT (ie smaller)
		int numDownShift = 0;
		for(int j=1; j<=vec[posRight]; ++j)
		{
			numDownShift += GetPartitionEnumNum( numItemsTot - vec[posRight] + j, posRight);
		}
		// this consider the situation when the highest bit is the SAME
		vector<int> vecDoneShift;
		for(int j = 0; j<posRight; ++j)
		{
			vecDoneShift.push_back( vec[j] );
		}
		int numSameOrder;
		if( numItemsTot > vec[posRight] )
		{
			numSameOrder = GetPartitionEnumId(numItemsTot-vec[posRight], vecDoneShift);
		}
		else
		{
			numSameOrder = 0;
		}
        //cout << "numDownShift : " << numDownShift << ", numSameOrder: " << numSameOrder << endl;

		//return numDownShift + numSameOrder;
	}
#endif
}

// **************************************************************************************
void MoveOneItemInPartEnum(const vector<vector<int> > &partsSrc, int part,
                           int psrc, int pdest,
                           vector<vector<int> > &partsDest) {
  YW_ASSERT_INFO(partsSrc.size() > 0, "MoveOneItemInPartEnum: wrong1");
  YW_ASSERT_INFO(part < (int)partsSrc.size(), "MoveOneItemInPartEnum: wrong2");
  YW_ASSERT_INFO(psrc < (int)partsSrc[0].size() &&
                     pdest < (int)partsSrc[0].size(),
                 "MoveOneItemInPartEnum: wrong3");
  partsDest = partsSrc;
  partsDest[part][psrc]--;
  partsDest[part][pdest]++;
}

void ConvIndexToPartEnum(const vector<int> &vecSegSizes, int numParts,
                         int pIndex, vector<vector<int> > &parts) {
  // convert the index of enumeration to a real enumeration
  // parts.clear();

  // it would be nice to implement it, but there is clear use yet. so skip
  // YW_ASSERT_INFO(false, "Not implemented yet. TBD.");
  vector<int> listSizes;
  for (int i = 0; i < (int)vecSegSizes.size(); ++i) {
    listSizes.push_back(numParts);
  }
  ConvIndexToPartEnumVar(vecSegSizes, listSizes, pIndex, parts);
}

void ConvIndexToPartEnumVar(const vector<int> &vecSegSizes,
                            const vector<int> &listNumParts, int pIndex,
                            vector<vector<int> > &parts) {
#if 0
cout << "ConvIndexToPartEnumVar: vecSegSizes: ";
DumpIntVec(vecSegSizes);
cout << "ListNumParts: ";
DumpIntVec(listNumParts);
cout << "pindex: " << pIndex << endl;
#endif
  //
  YW_ASSERT_INFO(vecSegSizes.size() == listNumParts.size(), "Mismatch");
  // get the index (order in the enumerated list) of the given enumerated
  // partition cout from the right hand side
  parts.clear();

  int res = pIndex;
  for (int i = 0; i < (int)vecSegSizes.size(); ++i) {
    int totEnumNumStep = GetPartitionEnumNum(vecSegSizes[i], listNumParts[i]);
    int idStep = (res % totEnumNumStep);

    vector<int> partsStep;
    GetPartitionEnumPartForId(vecSegSizes[i], listNumParts[i], idStep,
                              partsStep);
    parts.push_back(partsStep);

    // reduce res
    res = (res - idStep) / totEnumNumStep;

    // cout << "idStep: " << idStep << ", partsStep: ";
    // DumpIntVec(partsStep);
  }
}

void AddIntVec(vector<int> &vecDest, const vector<int> &vecSrc) {
  YW_ASSERT_INFO(vecDest.size() == vecSrc.size(), "AddIntVec: size mismatch");
  for (int i = 0; i < (int)vecSrc.size(); ++i) {
    vecDest[i] += vecSrc[i];
  }
}

void SubtractIntVec(vector<int> &vecDest, const vector<int> &vecSubtracted) {
  //
  YW_ASSERT_INFO(vecDest.size() == vecSubtracted.size(),
                 "AddIntVec: size mismatch");
  for (int i = 0; i < (int)vecSubtracted.size(); ++i) {
    vecDest[i] -= vecSubtracted[i];
  }
}

void GetItemsInRange(const set<int> &items, int lb, int ub, set<int> &sset) {
  sset.clear();
  //
  for (set<int>::iterator it = items.begin(); it != items.end(); ++it) {
    if (*it >= lb && *it <= ub) {
      sset.insert(*it);
    }
  }
}

void InitRandom(int seed) {
  double randTmp = GetRandFraction();
  cout << "Get one random fraction: " << randTmp
       << ", then initialize random seed to " << seed << endl;
  srand(seed);
}
void PermuatePseudoRandomVec(vector<int> &vecPerm) {
  // take a simple strategy: pick two arbitary positions and exchange them
  int numRounds = vecPerm.size();
  int vecLen = vecPerm.size();
  for (int r = 0; r < numRounds; ++r) {
    int i = (int)((rand() * 1.0 / RAND_MAX) * vecLen);
    int j = (int)((rand() * 1.0 / RAND_MAX) * vecLen);
    // int i = (int) (vecLen * GetRandFraction() );
    // int j = (int) (vecLen * GetRandFraction() );
    int tmp = vecPerm[i];
    vecPerm[i] = vecPerm[j];
    vecPerm[j] = tmp;
  }
}

void UnionMultiset(multiset<int> &setUpdate, const multiset<int> &setAdded) {
  for (multiset<int>::iterator it = setAdded.begin(); it != setAdded.end();
       ++it) {
    setUpdate.insert(*it);
  }
}

void JoinMultiset(const multiset<int> &set1, const multiset<int> &set2,
                  multiset<int> &setInt) {
  for (multiset<int>::iterator it = set1.begin(); it != set1.end(); ++it) {
    if (set2.find(*it) != set2.end()) {
      setInt.insert(*it);
    }
  }
}

void ConvMSetToSet(const multiset<int> &mset, set<int> &ss) {
  ss.clear();
  for (multiset<int>::iterator it = mset.begin(); it != mset.end(); ++it) {
    ss.insert(*it);
  }
}

void DumpMultiset(const multiset<int> &mset) {
  for (multiset<int>::iterator it = mset.begin(); it != mset.end(); ++it) {
    cout << *it << "    ";
  }
  cout << endl;
}

int CalcNumNChooseK(int n, int k) {
  // how many ways to choose k items from n items
  YW_ASSERT_INFO(n >= k, "n must be no smaller than k");
  double res = 1.0;
  int kuse = k;
  if (n - k < kuse) {
    kuse = n - k;
  }
  for (int i = 0; i < kuse; ++i) {
    res *= (1.0 * (n - i)) / (i + 1);
  }
  return (int)res;
}

// the following two functions are used to enumerate all partitions

void InitSubsetPartitionEnum(int numItems, int numParts,
                             vector<vector<int> > &parts) {
  int n = numItems;
  int p = numParts;
  parts.clear();
  parts.push_back(vector<int>());
  for (int i = 0; i <= n - p; i++) {
    parts[0].push_back(i);
  }
  for (int i = n - p + 1; i <= n - 1; i++) {
    parts.push_back(vector<int>());
    parts[parts.size() - 1].push_back(i);
  }
}
bool GetNextSubsetPartitionEnum(int numItems, int numParts,
                                vector<vector<int> > &parts) {
  // assuming all the elements in @parts is distinct and the number of these
  // elements is @numItems
  int n = numItems;
  int p = numParts;
  if (((int)parts.size()) != p)
    return false;
  for (int i = 0; i < (int)parts.size(); i++) {
    if (parts[i].empty())
      return false;
    sort(parts[i].begin(), parts[i].end());
  }
  vector<int> M;
  vector<int> K;
  M.reserve(n);
  K.reserve(n);
  for (int i = 0; i < n; i++) {
    M.push_back(0);
    K.push_back(0);
  }
  int lastmin = -1;
  for (int i = 0; i < (int)parts.size(); ++i) {
    int mmin = n;
    int key = -1;
    for (int j = 0; j < (int)parts.size(); ++j) {
      if (parts[j][0] > lastmin && parts[j][0] < mmin) {
        key = j;
        mmin = parts[j][0];
      }
    }
    lastmin = mmin;
    for (int j = 0; j < (int)parts[key].size(); ++j) {
      K[parts[key][j]] = i;
    }
  }
  M[0] = K[0];
  for (int i = 1; i < n; i++) {
    if (K[i] > M[i - 1])
      M[i] = K[i];
    else
      M[i] = M[i - 1];
  }

  bool success = false;
  for (int i = n - 1; i >= 1; --i) {
    if (K[i] < p - 1 && K[i] <= M[i - 1]) {
      success = true;
      K[i] = K[i] + 1;
      if (K[i] > M[i])
        M[i] = K[i];
      for (int j = i + 1; j <= n - (p - M[i]); ++j) {
        K[j] = 0;
        M[j] = M[i];
      }
      for (int j = n - (p - M[i]) + 1; j <= n - 1; ++j) {
        K[j] = p - (n - j);
        M[j] = p - (n - j);
      }
      break;
    }
  }
  if (!success)
    return false;
  parts.clear();
  for (int i = 0; i < p; i++) {
    parts.push_back(vector<int>());
  }
  for (int i = 0; i < n; i++) {
    parts[K[i]].push_back(i);
  }
  return true;
}

// another enumeration: we have n items, need to consider all possible splits of
// n into k parts where there is a limit of sizes for each of the k parts. E.g.
// n=10, 3 types, bounds=2,4,8 (type 1 has no more than 2, type-2 has no more
// than 4 and type-3 has no more than 8) we assume sum of these bounds >=n.
// Otherwise fatal error. Then we can have [1,3,6],[0,2,8] and so on in the case
// lower bounds are small, we start with the last entry being the highest number
void InitBoundedPartitionEnum(int numItems,
                              const vector<int> &lowerBoundsOnParts,
                              const vector<int> &upperBoundsOnParts,
                              vector<int> &partSizes) {
  YW_ASSERT_INFO(upperBoundsOnParts.size() == lowerBoundsOnParts.size(),
                 "Bound sizes: mismatch");
  YW_ASSERT_INFO(upperBoundsOnParts.size() >= 1,
                 "Must have at least one partition");
  YW_ASSERT_INFO(SumIntVector(upperBoundsOnParts) >= numItems,
                 "InitBoundedPartitionEnum: upper bounds too small");
  int sumLBs = SumIntVector(lowerBoundsOnParts);
  YW_ASSERT_INFO(sumLBs <= numItems,
                 "InitBoundedPartitionEnum: lower bounds too large");
  // now start enumerate
  partSizes = lowerBoundsOnParts;
  partSizes[partSizes.size() - 1] = numItems - sumLBs;
  // cout << "InitBoundedPartitionEnum: partSizes = ";
  // DumpIntVec(partSizes);
}

bool GetNextBoundedPartitionEnum(int numItems,
                                 const vector<int> &lowerBoundsOnParts,
                                 const vector<int> &upperBoundsOnParts,
                                 vector<int> &partSizes) {
#if 0
cout << "numItems = " << numItems << ", LBs = ";
DumpIntVec( lowerBoundsOnParts );
cout << " UBs = ";
DumpIntVec( upperBoundsOnParts );
cout << "Current part sizes = ";
DumpIntVec( partSizes );
#endif
  // in general, try to increase the rightmost (the last part) size unless it is
  // already at the limit that is, search for the second rightmost part (the
  // rightmost one is fixed once the other is fixed)  that is not at its upper
  // bound yet
  int pos = -1;
  int sumRight = 0;
  for (pos = (int)partSizes.size() - 2; pos >= 0; --pos) {
    //
    if (partSizes[pos] < upperBoundsOnParts[pos]) {
      break;
    }
    sumRight += partSizes[pos];
  }
  // cout << "GetNextBoundedPartitionEnum: pos = " << pos << ", sumRight = " <<
  // sumRight << endl;
  // if pos is not found (<0), done
  if (pos < 0) {
    return false;
  }
  // inc the current pos by 1 and reset the positions to its right to lower
  // bound
  partSizes[pos]++;
  sumRight--;
  for (int p = pos + 1; p < (int)partSizes.size() - 1; ++p) {
    partSizes[p] = lowerBoundsOnParts[p];
    sumRight -= lowerBoundsOnParts[p];
  }
  partSizes[(int)partSizes.size() - 1] += sumRight;
  YW_ASSERT_INFO(partSizes[(int)partSizes.size() - 1] <=
                         upperBoundsOnParts[(int)partSizes.size() - 1] &&
                     partSizes[(int)partSizes.size() - 1] >=
                         lowerBoundsOnParts[(int)partSizes.size() - 1],
                 "Part sizes: wrong");
  // cout << "GetNextBoundedPartitionEnum: partSizes = ";
  // DumpIntVec(partSizes);
  return true;
}

void UnionStrings(const set<string> &s1, const set<string> &s2,
                  set<string> &resSet) {
  resSet.clear();
  resSet = s1;
  for (set<string>::iterator it = s2.begin(); it != s2.end(); ++it) {
    resSet.insert(*it);
  }
}
bool AreStringsSubsetOf(const set<string> &s1Contained,
                        const set<string> &s2Container) {
  if (s1Contained.size() > s2Container.size()) {
    return false;
  }
  for (set<string>::iterator it = s1Contained.begin(); it != s1Contained.end();
       ++it) {
    if (s2Container.find(*it) == s2Container.end()) {
      return false;
    }
  }
  return true;
}

int SumIntVector(const vector<int> &vecInts) {
  int res = 0;
  for (int i = 0; i < (int)vecInts.size(); ++i) {
    res += vecInts[i];
  }
  return res;
}

double GetSumOfElements(const vector<double> &listVals) {
  double res = 0.0;
  for (int i = 0; i < (int)listVals.size(); ++i) {
    res += listVals[i];
  }
  return res;
}

void FindAllVectorsKStatesLen(int ks, int lenVec,
                              vector<vector<int> > &listAllVecs,
                              bool fOrderByStates) {
  // find all vectors with certain length and can choose from some states 0 to
  // ks-1 fOrderByStates: means vectors in states must be ordered in their first
  // apearnce that is, 2,3,1,2,3 ==> 1,2,3,1,2
  listAllVecs.clear();
  // recursively: start with a single length
  if (lenVec < 1) {
    // nothing
    return;
  }
  if (lenVec == 1) {
    // have ks states: 0,1,...ks-1
    for (int i = 0; i < ks; ++i) {
      vector<int> vec;
      vec.push_back(i);
      listAllVecs.push_back(vec);
    }
  } else {
    // recurisvely perform it
    vector<vector<int> > listVecsOneLess;
    FindAllVectorsKStatesLen(ks, lenVec - 1, listVecsOneLess);
    for (int jj = 0; jj < (int)listVecsOneLess.size(); ++jj) {
      // for each append one more
      int nsStart = 0;
      if (fOrderByStates == true) {
        // find the largest item so far and start with it
        for (int kk = 0; kk < (int)listVecsOneLess[jj].size(); ++kk) {
          if (listVecsOneLess[jj][kk] > nsStart) {
            nsStart = listVecsOneLess[jj][kk];
          }
        }
      }
      for (int i = nsStart; i < ks; ++i) {
        vector<int> vecnew = listVecsOneLess[jj];
        vecnew.push_back(i);
        listAllVecs.push_back(vecnew);
      }
    }
  }
}

void EraseCommonItemsFrom(vector<int> &listItems1, vector<int> &listItems2) {
  // remove shared common items
  // first sort the list
  SortIntVec(listItems1);
  SortIntVec(listItems2);
  // cout << "Before EraseCommonItemsFrom: \n";
  // DumpIntVec(listItems1);
  // DumpIntVec(listItems2);
  vector<int> listItemNew1, listItemNew2;
  // iterate through the two list concurrently, and avoid one common item when
  // needed
  int pos1 = 0, pos2 = 0;
  while (pos1 < (int)listItems1.size() && pos2 < (int)listItems2.size()) {
    // if one item is bigger than move it
    if (listItems1[pos1] < listItems2[pos2]) {
      // put the item to new list
      listItemNew1.push_back(listItems1[pos1]);
      pos1++;
    } else if (listItems1[pos1] > listItems2[pos2]) {
      listItemNew2.push_back(listItems2[pos2]);
      pos2++;
    } else {
      // move together but skip the common items
      pos1++;
      pos2++;
    }
  }
  // now add whatever left over to the two list
  for (int i = pos1; i < (int)listItems1.size(); ++i) {
    listItemNew1.push_back(listItems1[i]);
  }
  for (int i = pos2; i < (int)listItems2.size(); ++i) {
    listItemNew2.push_back(listItems2[i]);
  }
  listItems1 = listItemNew1;
  listItems2 = listItemNew2;
  // cout << "AFTER EraseCommonItemsFrom: \n";
  // DumpIntVec(listItems1);
  // DumpIntVec(listItems2);
}

void OffsetIntSetBy(set<int> &ss, int offset) {
  //
  set<int> sres;
  for (set<int>::iterator it = ss.begin(); it != ss.end(); ++it) {
    sres.insert((*it) + offset);
  }
  ss = sres;
}

static int QSortComparePairs(const void *arg1, const void *arg2) {
  /* Compare all of both strings: */
  // assume sorting in accending order
  pair<int, void *> p1 = *((pair<int, void *> *)arg1);
  pair<int, void *> p2 = *((pair<int, void *> *)arg2);
  // cout <<"arg1 = " << n1 << ", arg2 = " << n2 << endl;
  if (p1.first > p2.first) {
    return 1;
  } else if (p1.first < p2.first) {
    return -1;
  } else {
    return 0;
  }
}

void SortPairsByNums(vector<pair<int, void *> > &listPairs) {
  //#if 0
  if (listPairs.size() <= 1) {
    // do nothing
    return;
  }
  // cout << "Before sort, double vec = ";
  // DumpDoubleVec( vecVals );
  int sortLen = (int)listPairs.size();

  int start = 0;
  int end = sortLen - 1;
  pair<int, void *> *array = new pair<int, void *>[sortLen];
  for (int i = start; i <= end; ++i) {
    array[i - start] = listPairs[i];
  }
  qsort((void *)array, sortLen, sizeof(pair<int, void *>), QSortComparePairs);
  // Now write back
  for (int i = start; i <= end; ++i) {
    listPairs[i] = array[i - start];
  }

  delete[] array;
  //#endif
  // cout << "After sort, double vec = ";
  // DumpDoubleVec( vecVals );
}

static int QSortComparePairsDouble(const void *arg1, const void *arg2) {
  /* Compare all of both strings: */
  // assume sorting in accending order
  pair<double, void *> p1 = *((pair<double, void *> *)arg1);
  pair<double, void *> p2 = *((pair<double, void *> *)arg2);
  // cout <<"arg1 = " << n1 << ", arg2 = " << n2 << endl;
  if (p1.first > p2.first) {
    return 1;
  } else if (p1.first < p2.first) {
    return -1;
  } else {
    return 0;
  }
}

void SortPairsByNumsDouble(vector<pair<double, void *> > &listPairs) {
  if (listPairs.size() <= 1) {
    // do nothing
    return;
  }
  // cout << "Before sort, double vec = ";
  // DumpDoubleVec( vecVals );
  int sortLen = (int)listPairs.size();

  int start = 0;
  int end = sortLen - 1;
  pair<double, void *> *array = new pair<double, void *>[sortLen];
  for (int i = start; i <= end; ++i) {
    array[i - start] = listPairs[i];
  }
  qsort((void *)array, sortLen, sizeof(pair<double, void *>),
        QSortComparePairsDouble);
  // Now write back
  for (int i = start; i <= end; ++i) {
    listPairs[i] = array[i - start];
  }

  delete[] array;
}

//**************************************************************************************************************
// Ruhua Jiang's code for enumeration

static void convert_index_to_vector_helper(bool store_enum, int query_index,
                                           int color_num, int box_num,
                                           int &count, vector_t &vec,
                                           vector_t &result,
                                           vector_vector_t &enumeration) {
  if (result.size() != 0 && !store_enum)
    return;
  // Base case
  if (color_num == 1) {
    vec.push_back(box_num);
    count++;
    if (store_enum)
      enumeration.push_back(vec);

    if (count - 1 == query_index) {
      // std::cout<<count-1<<"\t";
      for (int k = 0; k < (int)vec.size(); k++) {
        // std::cout<<vec[k]<<" ";
        result.push_back(vec[k]);
      }
    }
    vec.pop_back();
    // std::cout<<endl;
    return;
  }
  // Recursion
  for (int i = box_num; i >= 0; i--) {
    vec.push_back(i);
    convert_index_to_vector_helper(store_enum, query_index, color_num - 1,
                                   box_num - i, count, vec, result,
                                   enumeration);
    vec.pop_back();
  }
}

static void convert_vector_to_int_helper(bool store_enum, vector_t query_vec,
                                         int color_num, int box_num, int &count,
                                         vector_t &vec, bool &find,
                                         vector_vector_t &enumeration) {
  if (find && !store_enum)
    return;
  // Base case
  if (color_num == 1) {
    vec.push_back(box_num);
    count++;
    if (store_enum)
      enumeration.push_back(vec);
    if (vec == query_vec) {
      find = true;
    }
    vec.pop_back();
    // std::cout<<endl;
    return;
  }
  // Recursion
  for (int i = box_num; i >= 0; i--) {
    vec.push_back(i);
    convert_vector_to_int_helper(store_enum, query_vec, color_num - 1,
                                 box_num - i, count, vec, find, enumeration);
    vec.pop_back();
  }
}

// Returns whether enumeration is stored or not. If index is not find,
// result.size() still 0
bool convert_index_to_vector(bool enum_already_set, int color_num, int box_num,
                             int index, vector_t &result,
                             vector_vector_t &enumeration) {
  int count = 0;
  vector_t vec;
  // if enumeration is stored or not, then directly access
  if (enum_already_set) {
    if (index < (int)enumeration.size()) {
      for (int k = 0; k < (int)enumeration[index].size(); k++) {
        result.push_back(enumeration[index][k]);
      }
      // std::cout<<"direct access!";  //uncomments this line if want test
      // whether direct access success or not
    }
    return true;
  } else {

    if (color_num > BOX_NUM_THRESHOLD ||
        box_num > COLOR_NUM_THRESHOLD) // c and n too large, we do not store
                                       // enumeration
    {
      convert_index_to_vector_helper(false, index, color_num, box_num, count,
                                     vec, result, enumeration);
      return false;
    } else {
      convert_index_to_vector_helper(true, index, color_num, box_num, count,
                                     vec, result, enumeration);
      return true;
    }
  }
}

// Returns whether enumeration is stored or not. If query_vec is not find,
// result_index is set to -1
bool convert_vector_to_index(bool enum_already_set, vector_t query_vec,
                             int &result_index, vector_vector_t &enumeration) {
  int color_num = query_vec.size(), box_num = 0, index = 0;
  for (int i = 0; i < (int)query_vec.size(); i++)
    box_num += query_vec[i];
  vector_t vec;
  bool find = false;

  // if enumeration is stored or not, then directly compare
  if (enum_already_set) {
    for (int i = 0; i < (int)enumeration.size(); i++) {
      if (query_vec == enumeration[i]) {
        result_index = i;
        // std::cout<<"direct access!";  //uncomments this line if want test
        // whether direct access success or not
        return enum_already_set;
      }
    }

    result_index = -1;

    // is this correct???
    return false;
  } else {
    if (color_num > BOX_NUM_THRESHOLD || box_num > COLOR_NUM_THRESHOLD) {
      convert_vector_to_int_helper(false, query_vec, color_num, box_num, index,
                                   vec, find, enumeration);
      result_index = index - 1;
      return false;
    } else {
      convert_vector_to_int_helper(true, query_vec, color_num, box_num, index,
                                   vec, find, enumeration);
      result_index = index - 1;
      return false;
    }
  }
}

void ZeroOutVec(vector<int> &vec) {
  for (int i = 0; i < (int)vec.size(); ++i) {
    vec[i] = 0;
  }
}

void GetFourPartsIncompatSplits(const set<int> &setAll, const set<int> &split1,
                                const set<int> &split2, set<int> &part1,
                                set<int> &part2, set<int> &part3,
                                set<int> &part4) {
  //
  set<int> split1b = setAll;
  SubtractSets(split1b, split1);
  set<int> split2b = setAll;
  SubtractSets(split2b, split2);
  JoinSets(split1, split2, part1);
  JoinSets(split1, split2b, part2);
  JoinSets(split1b, split2, part3);
  JoinSets(split1b, split2b, part4);
}

bool IsAllZeroVec(const vector<int> &vec) {
  for (int i = 0; i < (int)vec.size(); ++i) {
    if (vec[i] != 0) {
      return false;
    }
  }
  return true;
}
