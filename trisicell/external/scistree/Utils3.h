#ifndef UTILS3_H
#define UTILS3_H

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace std;

#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "Utils.h"
#include "Utils2.h"

// ***************************************************************************
// HashTable
// ***************************************************************************

// Abstract class for item stored in my hash table
class YWHashItem {
public:
  virtual ~YWHashItem() = 0;
  virtual int Key() = 0;
  virtual bool operator==(const YWHashItem &rhs) = 0;
};

// This is the hash table that mgiht be useful in some applications
// Note that this is a rather static HASH table: you can only add stuff in
// but can not remove. TBD
class YWHashTable {
public:
  YWHashTable(int numBuckets = 100);
  ~YWHashTable(); // NOTE: has to free memory here
  void AddItem(YWHashItem *pItem);
  YWHashItem *GetIdenticalItem(YWHashItem *pItem);
  YWHashItem *GetFirstItem();
  YWHashItem *GetNextItem();
  int GetTotalItemNum() const;
  void Dump() const;

private:
  int numBuckets;
  //    vector< vector<YWHashItem *> > hashTable;

  // Sorry we have not implemented hashing yet
  vector<YWHashItem *> hashTable;

  // TBD. These are for enumeation only. BUT only support single enumeration
  // PLEASE do not use in a double loop
  int curPos;
};

// used to support STL
class SequenceCmp //: public binary_function<SequenceCmp &, const SEQUENCE &,
                  //const SEQUENCE>
{
public:
  bool operator()(const SEQUENCE &seq1, const SEQUENCE &seq2) const;
};

// iterator pattern
class GenericIterator {
public:
  virtual ~GenericIterator() {}
  virtual void First() = 0;
  virtual void Next() = 0;
  virtual bool IsDone() = 0;
  // virtual void *GetCurItem() = 0;
};

// ***************************************************************************
// common utilities
// ***************************************************************************
bool IsIntervalContained(const INTERVAL &iv1, const INTERVAL &iv2);
int GetIntervalLen(const INTERVAL &iv);
int GetRandItemInSet(const set<int> &items);
int GetRandItemInVec(const vector<int> &items);
void GetRandVector(vector<int> &rndVec, int start, int end);
int GetWeightedRandItemInVec(const vector<int> &items,
                             const vector<double> &itemWeights);
int GetWeightedRandItemIndex(const vector<double> &itemWeights);
// this function converts the subset over a vector to the subet over the
// original space
void GetOrigSubset(const vector<int> &origVec, const set<int> &subsetInd,
                   set<int> &subsetOrig);
void MutateSequenceAtSites(SEQUENCE &mutSeq, vector<int> &mutSites);
void DumpDoubleVec(const vector<double> &vecDoubles);
void DumpDoubleVec(const vector<long double> &vecDoubles);
void DumpBoolVec(const vector<bool> &vecBools);
int GetLargestIndiceInDoubleVec(const vector<double> &vecDoubles);
int GetLargestIndiceInDoubleVec(const vector<long double> &vecDoubles);
double FindMedian(const vector<double> &vecVals);
long double FindMedian(const vector<long double> &vecVals);
double FindMaxDouble(const vector<double> &vecVals);
double FindMaxDouble(const vector<long double> &vecVals);
double FindRankedItem(const vector<double> &vecVals, int rank);
void SortDoubleVec(vector<double> &vecVals, int start = 0, int end = -1);
void SortDoubleVec(vector<long double> &vecVals, int start = 0, int end = -1);
void FindUniformColumns(const vector<SEQUENCE> &listSeqs, set<int> &uniSites);
int FindNoninformativeRow(const vector<SEQUENCE> &listSeqs, int col);
void BreakSeqAtBkpt(const SEQUENCE &seq, int bkpt, SEQUENCE &seqLeft,
                    SEQUENCE &seqRight);
bool AreTwoSeqsBroken(const SEQUENCE &seqLeft, const SEQUENCE &seqRight);

// support parition-enumeration
// Suppose we have g groups of (indistingishable) items and we want to
// divide each group into numParts colors (distinguishable)
// this support enumerate these choices. For example, we have two segments of 3
// and 4 items each and we have two colors, then the choices will be: [(1,2),
// (2,2)], or [(0.3),(1,3)]
void InitPartitionEnum(const vector<int> &vecSegSizes, int numParts,
                       vector<vector<int> > &parts);
bool GetNextPartitionEnum(const vector<int> &vecSegSizes, int numParts,
                          vector<vector<int> > &parts);
int GetPartEnumIndex(const vector<int> &vecSegSizes, int numParts,
                     const vector<vector<int> > &parts);
void ConvIndexToPartEnum(const vector<int> &vecSegSizes, int numParts,
                         int pIndex, vector<vector<int> > &parts);
void ConvIndexToPartEnumVar(const vector<int> &vecSegSizes,
                            const vector<int> &numParts, int pIndex,
                            vector<vector<int> > &parts);
void InitPartitionEnumVar(const vector<int> &vecSegSizes,
                          const vector<int> &numParts,
                          vector<vector<int> > &parts);
bool GetNextPartitionEnumVar(const vector<int> &vecSegSizes,
                             const vector<int> &numParts,
                             vector<vector<int> > &parts);
int GetPartEnumIndexVar(const vector<int> &vecSegSizes,
                        const vector<int> &numParts,
                        const vector<vector<int> > &parts);
void MoveOneItemInPartEnum(const vector<vector<int> > &partsSrc, int part,
                           int psrc, int pdest,
                           vector<vector<int> > &partsDest);
int GetPartitionEnumNum(int n, int p);
int GetPartitionEnumId(int numItemsTot, const vector<int> &vec);
void GetPartitionEnumPartForId(int numItemsTot, int numParts, int eid,
                               vector<int> &vec);

// support another version of partiton-enumeration
// Suppose we have n (distinct) items, and we want to partition into k groups
// (each with at least one item) E.g. we have {a,b,c,d} and we want to partition
// into 3 groups. Then choices are: {a,b,cd}, {ab,c,d}, {ac,b,d}, and so on
void InitSubsetPartitionEnum(int numItems, int numParts,
                             vector<vector<int> > &parts);
bool GetNextSubsetPartitionEnum(int numItems, int numParts,
                                vector<vector<int> > &parts);

// another enumeration: we have n items, need to consider all possible splits of
// n into k parts where there is a limit of sizes for each of the k parts. E.g.
// n=10, 3 types, bounds=2,4,8 (type 1 has no more than 2, type-2 has no more
// than 4 and type-3 has no more than 8) we assume sum of these bounds >=n.
// Otherwise fatal error. Then we can have [1,3,6],[0,2,8] and so on
void InitBoundedPartitionEnum(int numItems,
                              const vector<int> &lowerBoundsOnParts,
                              const vector<int> &upperBoundsOnParts,
                              vector<int> &partSizes);
bool GetNextBoundedPartitionEnum(int numItems,
                                 const vector<int> &lowerBoundsOnParts,
                                 const vector<int> &upperBoundsOnParts,
                                 vector<int> &partSizes);

// new things from treeHMM
bool GetFirstMutliChoice(int numStage, int numStageElem,
                         vector<int> &initChoice);
bool GetNextMutliChoice(int numStage, int numStageElem,
                        vector<int> &initChoice);
// void DumpVecSequences( const vector<SEQUENCE> &setSeqs );
void GetVecSequencesIV(const vector<SEQUENCE> &vecSeqs, int left, int right,
                       vector<SEQUENCE> &vecSeqsIV);
int GetNumZerosInSeq(const SEQUENCE &seq);
void GetSeqSplit(const SEQUENCE &seq, set<int> &zeroBits, set<int> &oneBits);
void SortVecIntPairs(vector<pair<int, int> > &listOfPriority);
void ConvVecToArray(const vector<int> &vec, int *arr);
void ConvVecToArray(const vector<double> &vec, double *arr);
void DumpIntArray(int len, int *arr);
void FlipBinVector(vector<int> &vec);
void ConvOneSideToFullSplit(vector<int> &split, const set<int> &oneside,
                            int numLeaves, int val = 1);

// more on missing value
bool AreTwoMVVecCompat(const vector<int> &vec1, const vector<int> &vec2,
                       int &numTrueMatch);
int GetMVNum(const vector<int> &vec);
bool AreSeqsOverlap(const vector<int> &vec1, const vector<int> &vec2);

// ***************************************************************************
// Substring utilities
// ***************************************************************************
// This file contains some extra utilties that are frequently used
typedef pair<INTERVAL, SEQUENCE> INTERVAL_SUBSTRING; // (start,end, sequence)

int GetSubstringLeftPos(const INTERVAL_SUBSTRING &substr);
int GetSubstringRightPos(const INTERVAL_SUBSTRING &substr);
void GetIVSubstringData(const INTERVAL_SUBSTRING &substr, SEQUENCE &seq);
INTERVAL GetSubstringInterval(const INTERVAL_SUBSTRING &substr);
bool GetSubstringSegment(const INTERVAL_SUBSTRING &substr,
                         const INTERVAL &ivToRead, SEQUENCE &segment);
int GetSubstringValAt(const INTERVAL_SUBSTRING &substr, int pos);
bool IsSegmentContained(const INTERVAL_SUBSTRING &seqContained,
                        const INTERVAL_SUBSTRING &seqContainer);
bool AreSegmentsConsistent(const INTERVAL_SUBSTRING &seqContained,
                           const INTERVAL_SUBSTRING &seqContainer);
int GetSegmentsIntersection(const INTERVAL_SUBSTRING &seq1,
                            const INTERVAL_SUBSTRING &seq2, INTERVAL &iv);
bool AreSegmentsNextto(const INTERVAL_SUBSTRING &seq1,
                       const INTERVAL_SUBSTRING &seq2);
void DumpSubstring(const INTERVAL_SUBSTRING &substr);

// ***************************************************************************
// Numerical utilities
// ***************************************************************************
double GetLogSumOfLogs(const vector<double> &listLogs);
double GetLogSumOfLogsDirect(const vector<double> &listLogs);
double GetLogSumOfTwo(double logv1, double logv2);
double GetSumOfElements(const vector<double> &listVals);
void SumofLogVecs(vector<double> &listLogsAdded,
                  vector<double> &listLogsAdding);

// ***************************************************************************
// Other utilities
// ***************************************************************************
int FindMatchedSeqForFounders(const vector<SEQUENCE> &founder,
                              const SEQUENCE &seq, set<int> &endRows,
                              bool fPrefix);
void RecoverOrigIndicesAfterDeletion(const vector<int> &removedItems,
                                     const vector<int> &itemsNew,
                                     vector<int> &itemsOrigIndices);
void GetOrigPositionAfterRemoval(int numRemains,
                                 const vector<int> &itemsRemoved,
                                 vector<int> &origPosForRemains);
void InsertOrderedVec(vector<int> &vec, int val);
template <typename T>
int binary_search(const std::vector<T> &vec, unsigned start, unsigned end,
                  const T &key);
bool ReadIntListFromFile(const char *fname, vector<int> &listInts);
void GetVecPosNotInSet(const vector<int> &vec, const set<int> &s,
                       vector<int> &posDiff);
void AddIntVec(vector<int> &vecDest, const vector<int> &vecSrc);
void SubtractIntVec(vector<int> &vecDest, const vector<int> &vecSubtracted);
void GetItemsInRange(const set<int> &items, int lb, int ub, set<int> &sset);
void InitRandom(int seed);
void PermuatePseudoRandomVec(vector<int> &vecPerm);
void UnionMultiset(multiset<int> &setUpdate, const multiset<int> &setAdded);
void JoinMultiset(const multiset<int> &set1, const multiset<int> &set2,
                  multiset<int> &setInt);
void ConvMSetToSet(const multiset<int> &mset, set<int> &ss);
void DumpMultiset(const multiset<int> &mset);
int CalcNumNChooseK(int n,
                    int k); // how many ways to choose k items from n items
void UnionStrings(const set<string> &s1, const set<string> &s2,
                  set<string> &resSet);
bool AreStringsSubsetOf(const set<string> &s1Contained,
                        const set<string> &s2Container);
int SumIntVector(const vector<int> &vecInts);
void FindAllVectorsKStatesLen(int ks, int ns, vector<vector<int> > &listAllVecs,
                              bool fOrderByStates = false);
void EraseCommonItemsFrom(vector<int> &listItems1, vector<int> &listItems2);
void OffsetIntSetBy(set<int> &ss, int offset);
void SortPairsByNums(vector<pair<int, void *> > &listPairs);
void SortPairsByNumsDouble(vector<pair<double, void *> > &listPairs);
void ZeroOutVec(vector<int> &vec);
void GetFourPartsIncompatSplits(const set<int> &setAll, const set<int> &split1,
                                const set<int> &split2, set<int> &part1,
                                set<int> &part2, set<int> &part3,
                                set<int> &part4);
bool IsAllZeroVec(const vector<int> &vec);

// ***************************************************************************
// template utilties
// ***************************************************************************
template <class TYPE>
void JoinSetsGen(const set<TYPE> &set1, const set<TYPE> &set2,
                 set<TYPE> &sint) {
  //
  sint.clear();
  for (typename set<TYPE>::iterator it = set1.begin(); it != set1.end(); ++it) {
    //
    if (set2.find(*it) != set2.end()) {
      //
      sint.insert(*it);
    }
  }
}

template <class TYPE>
void UnionSetsGen(set<TYPE> &setAdded, const set<TYPE> &setAddin) {
  //
  for (typename set<TYPE>::iterator it = setAddin.begin(); it != setAddin.end();
       ++it) {
    //
    setAdded.insert(*it);
  }
}

template <class TYPE>
void SubtractSetsGen(set<TYPE> &setMain, const set<TYPE> &setSubtracted) {
  //
  for (typename set<TYPE>::iterator it = setSubtracted.begin();
       it != setSubtracted.end(); ++it) {
    //
    setMain.erase(*it);
  }
}

template <class TYPE>
bool AreItemsSimilar(const vector<TYPE> &listItems, const TYPE &tol) {
  // are the number of items within some toleratnce from the average
  TYPE sum = 0;
  for (typename vector<TYPE>::const_iterator it = listItems.begin();
       it != listItems.end(); ++it) {
    //
    sum += *it;
  }
  for (typename vector<TYPE>::const_iterator it = listItems.begin();
       it != listItems.end(); ++it) {
    //
    if ((*it) - sum / listItems.size() > tol * sum / listItems.size() ||
        sum / listItems.size() - (*it) > tol * sum / listItems.size()) {
      return false;
    }
  }
  return true;
}

//**************************************************************************************************************
// Ruhua Jiang's code for enumeration
typedef std::vector<int> vector_t;
typedef std::vector<vector_t> vector_vector_t;
#if 0
typedef struct Enumeration{
	vector_vector_t enumeration;
	bool enumeration_set;
	Enumeration(){enumeration_set = false;};
};
typedef Enumeration Enumeration;
#endif
const int COLOR_NUM_THRESHOLD = 40;
const int BOX_NUM_THRESHOLD = 5;

// Notice that index is zero based number

bool convert_index_to_vector(bool enum_already_set, int color_num, int box_num,
                             int index, vector_t &result,
                             vector_vector_t &enumeration);
bool convert_vector_to_index(bool enum_already_set, vector_t query_vec,
                             int &result_index, vector_vector_t &enumeration);

#endif // UTILS3_H
