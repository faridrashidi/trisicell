#ifndef UTILS2_H
#define UTILS2_H

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>
//#include <limits>
using namespace std;

#include "Utils.h"
#include <ctime>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

// This file contains some extra utilties that are frequently used

// ***************************************************************************
// Common utilities
// ***************************************************************************
long GetCurrentTimeTick();
long GetElapseTime(long lastTime);
void GetCurrentCPUTime(std::clock_t &tmStart);
double GetElapseCPUTime(const std::clock_t &tmStart);
bool IsBoolArrayAllTrue(bool *bArray, int size);
void AppendIntVec(vector<int> &dest, const vector<int> &appending);
bool IsSetContainer(const set<int> &contianer, const set<int> &contained);
bool IsSetContainedInSets(const set<int> &s, const set<set<int> > &sets);
bool IsSetContainingOneOfSets(const set<int> &s, const set<set<int> > &sets);
void ConcatIntVec(vector<int> &vecAdded, const vector<int> &vecToAdd);
int ConvIntSetToPosition(const set<int> &s);
void ConvPositionToIntSet(int val, set<int> &s);
void PopulateSetWithInterval(set<int> &s, int left, int right);
void GetSeqInterval(const SEQUENCE &row, SEQUENCE &rowIV, int left, int right);
bool IsIntervalContained(const set<SEQUENCE> &seqs, int left, int right,
                         const SEQUENCE &seqIV);
void SubtractSequenceSets(set<SEQUENCE> &s1, const set<SEQUENCE> &s2);
void DumpSequence(const SEQUENCE &seq);
void DumpVecSequences(const vector<SEQUENCE> &setSeqs);
void DumpSetSequences(const set<SEQUENCE> &setSeqs);
bool AreTwoInSameSet(int i1, int i2, const set<set<int> > &collections);
int GetItemIndexInVec(const vector<int> &vec, int item);
bool IsIntervalOverlap(const INTERVAL &iv1, const INTERVAL &iv2);
bool GetIntervalOverlap(const INTERVAL &iv1, const INTERVAL &iv2,
                        INTERVAL &ivBoth);
void GenerateRandBinVector(int sz, vector<int> &randVec);
bool IsBinary(int val);
void ReOrderWithRemovedSites(const vector<int> &posAfterRem,
                             const vector<int> &removedPos,
                             vector<int> &posBeforeRemove);
void GetSubsetVec(const vector<int> &vecOriginal, const set<int> &sitesToKeep,
                  vector<int> &vecNew);
void AddMissingVecBits(vector<int> &rowComplete, const set<int> &sitesToAdd,
                       vector<int> &partialRow);

// ***************************************************************************
// Utilies for phasing
// ***************************************************************************
bool IsSequenceHaplotype(const SEQUENCE &seq);
bool IsSequenceGenotype(const SEQUENCE &seq);
bool CanPhaseGenoRow(const SEQUENCE &hap1, const SEQUENCE &hap2,
                     const SEQUENCE &geno);
bool AreHapGenoRowCompatible(const SEQUENCE &hapRow, const SEQUENCE &genoRow,
                             SEQUENCE *pComplement = NULL);
bool AreHapGenoRowsSame(const SEQUENCE &hapRow, const SEQUENCE &genoRow);
bool IsTrivialRow(const SEQUENCE &row, SEQUENCE &resolved1,
                  SEQUENCE &resolved2);
bool IsHapSeqSmaller(const SEQUENCE &hapRow1, const SEQUENCE &hapRow2);
void CreateGenoRowFromHapRows(const SEQUENCE &hapRow1, const SEQUENCE &hapRow2,
                              SEQUENCE &genoRow);

// ***************************************************************************
// Utilies for hypercube related stuff
// ***************************************************************************
// We use the hypercube node index to access an hypercube node sequence
typedef int HCSequence;

void GetHyperCubeSeq(int hcSeq, SEQUENCE &seq, int hcWidth);
int GetHyperCubSeqBit(int hcSeq, int bit, int hcWidth);
int GetSeqIdFromSeq(const SEQUENCE &seq);
void FindNonSegSites(const set<HCSequence> &setSeqs, set<int> &sites,
                     int dataWidth);
void FindNonSegSites(const set<SEQUENCE> &setSeqs, set<int> &sites,
                     int dataLen);
int IsHCSeqsMutPair(HCSequence seq1, HCSequence seq2, int dataWidth);
bool IsHCSeqsMutPairAt(HCSequence seq1, HCSequence seq2, int dataWidth,
                       int pos);
void MutateHCSeqAt(const HCSequence seq, HCSequence &res, int dataWidth,
                   int mutPos);
bool IsHCSeqRecombinnable(HCSequence s1, HCSequence s2, HCSequence st,
                          int dataWidth);
void RecombineHCSeqs(const HCSequence hcSeq1, const HCSequence hcSeq2,
                     HCSequence &res, int dataWidth, int bkpt);

#endif // UTILS2_H
