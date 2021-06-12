#include "Utils2.h"
#include "Utils.h"
#include "cstdio"
#include "cstdlib"
#include "ctime"

//////////////////////////////////////////////////////////////////////////////
// Utility functions
//////////////////////////////////////////////////////////////////////////////

long GetCurrentTimeTick() { return (long)time(NULL); }

long GetElapseTime(long lastTime) {
  return (long)((long)time(NULL) - lastTime);
}

void GetCurrentCPUTime(std::clock_t &tmStart) {
  //
  tmStart = std::clock();
}

double GetElapseCPUTime(const std::clock_t &tmStart) {
  //
  return (std::clock() - tmStart) / (double)CLOCKS_PER_SEC;
}

bool IsBoolArrayAllTrue(bool *bArray, int size) {
  for (int i = 0; i < size; ++i) {
    if (bArray[i] == false) {
      return false;
    }
  }
  return true;
}

void AppendIntVec(vector<int> &dest, const vector<int> &appending) {
  for (unsigned int i = 0; i < appending.size(); ++i) {
    dest.push_back(appending[i]);
  }
}

bool IsSetContainer(const set<int> &container, const set<int> &contained) {
  for (set<int>::iterator it = contained.begin(); it != contained.end(); ++it) {
    if (container.find(*it) == container.end()) {
      return false;
    }
  }
  return true;
}

bool IsSetContainedInSets(const set<int> &s, const set<set<int> > &sets) {
  // This function return true if ONE of the sets in sets contains s, false
  // otherwise
  for (set<set<int> >::iterator it = sets.begin(); it != sets.end(); ++it) {
    if (IsSetContainer(*it, s) == true) {
      return true;
    }
  }
  return false;
}

bool IsSetContainingOneOfSets(const set<int> &s, const set<set<int> > &sets) {
  // This function return true if ONE of the sets in sets contains s, false
  // otherwise
  for (set<set<int> >::iterator it = sets.begin(); it != sets.end(); ++it) {
    if (IsSetContainer(s, *it) == true) {
      return true;
    }
  }
  return false;
}

void ConcatIntVec(vector<int> &vecAdded, const vector<int> &vecToAdd) {
  // Append one vector to another
  for (unsigned int i = 0; i < vecToAdd.size(); i++) {
    vecAdded.push_back(vecToAdd[i]);
  }
}

int ConvIntSetToPosition(const set<int> &s) {
  // cout << "In ConvIntSetToPosition: s = ";
  // DumpIntSet( s );
  // this function convert an integer set to a position index, for example,
  // if range = 8, s={2, 4}, then this converts to 00010100
  int res = 0;
  for (set<int>::iterator it = s.begin(); it != s.end(); ++it) {
    int a = *it;
    int mask = 0x1 << a;
    res = res | mask;
  }
  // cout << "conversion res = " << res << endl;
  return res;
}

void ConvPositionToIntSet(int val, set<int> &s) {
  // inverse to ConvIntSetToPosition: convert an integer back to a set
  s.clear();
  int pos = 0;
  while (val != 0) {
    if ((val & 0x1) != 0) {
      s.insert(pos);
    }
    pos++;
    // left-shift val
    val = (val >> 1);
  }
}

void PopulateSetWithInterval(set<int> &s, int left, int right) {
  s.clear();
  for (int i = left; i <= right; ++i) {
    s.insert(i);
  }
}

void GetSeqInterval(const SEQUENCE &row, SEQUENCE &rowIV, int left, int right) {
  rowIV.clear();
  for (int i = left; i <= right; ++i) {
    rowIV.push_back(row[i]);
  }
}

bool IsIntervalContained(const set<SEQUENCE> &seqs, int left, int right,
                         const SEQUENCE &seqIV) {
  // This function check to see if the seqIV is contained in the seqs
  // when there is missing site, we use COMPABILITY instead of ==
  for (set<SEQUENCE>::iterator it = seqs.begin(); it != seqs.end(); ++it) {
    SEQUENCE substr;
    GetSeqInterval(*it, substr, left, right);
    if (AreTwoSeqsCompatible(substr, seqIV) == true) {
      return true;
    }
  }
  return false;
}

void SubtractSequenceSets(set<SEQUENCE> &s1, const set<SEQUENCE> &s2) {
  if (s2.size() == 0) {
    return;
  }
  set<SEQUENCE> res;
  // this function performs set intersection, i.e. s1=s1 ^s2
  for (set<SEQUENCE>::iterator it = s1.begin(); it != s1.end(); ++it) {
    if (s2.find(*it) == s2.end()) {
      res.insert(*it);
    }
  }
  s1.clear();
  s1 = res;
}

void DumpSequence(const SEQUENCE &seq) {
  for (unsigned int i = 0; i < seq.size(); ++i) {
    if (IsMissingValueBit(seq[i]) == false) {
      cout << seq[i];
    } else {
      cout << "*";
    }
  }
  cout << endl;
}

void DumpVecSequences(const vector<SEQUENCE> &vecSeqs) {
  for (unsigned int i = 0; i < vecSeqs.size(); ++i) {
    DumpSequence(vecSeqs[i]);
  }
}

void DumpSetSequences(const set<SEQUENCE> &setSeqs) {
  for (set<SEQUENCE>::iterator it = setSeqs.begin(); it != setSeqs.end();
       ++it) {
    DumpSequence(*it);
  }
}

bool AreTwoInSameSet(int i1, int i2, const set<set<int> > &collections) {
  // Check to see if i1 and i2 is in same set
  for (set<set<int> >::iterator it = collections.begin();
       it != collections.end(); ++it) {
    bool found1 = false, found2 = false;
    if (it->find(i1) != it->end()) {
      found1 = true;
    }
    if (it->find(i2) != it->end()) {
      found2 = true;
    }
    if (found1 == true && found2 == true) {
      // cout << "i1 = " << i1 << ", i2 = " << i2 << " INDDED in same set.\n";
      return true;
    }
    if (found1 || found2) {
      // cout << "i1 = " << i1 << ", i2 = " << i2 << " not in same set.\n";
      return false;
    }
  }
  // should not need this, in case
  YW_ASSERT_INFO(false, "Bad i1 or i2.");
  return false;
}

int GetItemIndexInVec(const vector<int> &vec, int item) {
  for (unsigned int i = 0; i < vec.size(); ++i) {
    if (vec[i] == item) {
      return i;
    }
  }
  return -1;
}

bool IsIntervalOverlap(const INTERVAL &iv1, const INTERVAL &iv2) {
  if (iv1.second < iv2.first || iv2.second < iv1.first) {
    return false;
  } else {
    return true;
  }
}

bool GetIntervalOverlap(const INTERVAL &iv1, const INTERVAL &iv2,
                        INTERVAL &ivBoth) {
  int left = iv1.first;
  if (left < iv2.first) {
    left = iv2.first;
  }
  int right = iv1.second;
  if (right > iv2.second) {
    right = iv2.second;
  }
  if (left > right) {
    return false;
  } else {
    ivBoth.first = left;
    ivBoth.second = right;
    return true;
  }
}

void GenerateRandBinVector(int sz, vector<int> &randVec) {
  // cout << "GenerateRandBinVector: sz = " << sz << endl;
  // Generate random vector
  randVec.clear();
  for (int i = 0; i < sz; ++i) {
    // cout << " i = " << i << endl;
    double r = GetRandFraction();
    if (r >= 0.5) {
      randVec.push_back(0);
    } else {
      randVec.push_back(1);
    }
  }
}
bool IsBinary(int val) {
  if (val == 0 || val == 1) {
    return true;
  } else {
    return false;
  }
}
void ReOrderWithRemovedSites(const vector<int> &posAfterRem,
                             const vector<int> &removedPos,
                             vector<int> &posBeforeRemove) {
  // THis funciton is often used here
  // For example, we often removed sites from the matrix but then we need to
  // know their original positions this function consider that by adding the
  // removed sites back into order (not directly into posBeforeRem) but rather
  // consider them when adding
  posBeforeRemove.clear();

  unsigned int pos = 0;
  for (unsigned int i = 0; i < posAfterRem.size(); ++i) {
    while (pos < removedPos.size() &&
           posAfterRem[i] + (int)pos >= removedPos[pos]) {
      pos++;
    }
    posBeforeRemove.push_back(posAfterRem[i] + pos);
  }
}

void GetSubsetVec(const vector<int> &vecOriginal, const set<int> &sitesToKeep,
                  vector<int> &vecNew) {
  vecNew.clear();
  for (unsigned int i = 0; i < vecOriginal.size(); ++i) {
    if (sitesToKeep.find(i) != sitesToKeep.end()) {
      vecNew.push_back(vecOriginal[i]);
    }
  }
}

void AddMissingVecBits(vector<int> &rowOrig, const set<int> &sitesToAdd,
                       vector<int> &partialRow) {
  YW_ASSERT_INFO(sitesToAdd.size() == partialRow.size(),
                 "Parameter size mismatch");

  // If there is othing to work, stop
  if (sitesToAdd.size() == 0) {
    return;
  }

  cout << "AddMissingVecBits: rowOrig = ";
  DumpSequence(rowOrig);
  cout << "Append sites ";
  DumpIntSet(sitesToAdd);
  cout << "Missing values = ";
  DumpIntVec(partialRow);
  // Here we try to add back some missing sites
  vector<int> missingSites;
  PopulateVecBySet(missingSites, sitesToAdd);

  vector<int> res;
  int posMiss = 0;
  int posOrig = 0;
  int curpos = 0;

  while (posMiss < (int)partialRow.size() || posOrig < (int)rowOrig.size()) {
    // check to see which bit to use and move
    if (curpos != missingSites[posMiss]) {
      // This bit is original
      YW_ASSERT_INFO(posOrig < (int)rowOrig.size(),
                     "Serious error: not enough bits.");
      res.push_back(rowOrig[posOrig]);
      posOrig++;
    } else {
      // No this is a missing bit
      res.push_back(partialRow[posMiss]);
      posMiss++;
      ;
    }

    // now move on
    curpos++;
  }
  rowOrig = res;
  cout << "AddMissingVecBits: res = ";
  DumpSequence(rowOrig);
}

////////////////////////////////////////////////////////////////////////////////////////
bool IsSequenceHaplotype(const SEQUENCE &seq) {
  // note need to consider missing value!
  for (unsigned int i = 0; i < seq.size(); ++i) {
    if (seq[i] != 0 && seq[i] != 1 && IsMissingValueBit(seq[i]) == false) {
      return false;
    }
  }
  return true;
}

bool IsSequenceGenotype(const SEQUENCE &seq) {
  for (unsigned int i = 0; i < seq.size(); ++i) {
    if (seq[i] != 0 && seq[i] != 1 && seq[i] != 2 &&
        IsMissingValueBit(seq[i]) == false) {
      return false;
    }
  }
  return true;
}

bool CanPhaseGenoRow(const SEQUENCE &hap1, const SEQUENCE &hap2,
                     const SEQUENCE &geno) {
  YW_ASSERT_INFO(IsSequenceHaplotype(hap1), "hap1 is not haplotype row.");
  YW_ASSERT_INFO(IsSequenceHaplotype(hap2), "hap1 is not haplotype row.");
  YW_ASSERT_INFO(IsSequenceGenotype(geno), "hap1 is not haplotype row.");
  YW_ASSERT_INFO(hap1.size() == hap2.size(),
                 "Tow hap rows are not equal length");
  YW_ASSERT_INFO(geno.size() == geno.size(),
                 "Geno row is not the same size as hap row.");
  // for now, do not allow hap1/hap2 contian missing value
  YW_ASSERT_INFO(IsSeqHasMV(hap1) == false && IsSeqHasMV(hap2) == false,
                 "Hap1/Hap2 can not contain missing values");
  // cout << "hap1 = ";
  // DumpIntVec(hap1 );
  // cout << "hap2 = ";
  // DumpIntVec(hap2 );
  // cout << "geno = ";
  // DumpIntVec( geno );
  for (unsigned int i = 0; i < hap1.size(); i++) {
    // a missing vlaue can be phased either way
    if (IsMissingValueBit(geno[i]) == true) {
      continue;
    }

    if (geno[i] == 2) {
      if (hap1[i] + hap2[i] != 1) {
        return false;
      }
    } else {
      if (hap1[i] + hap2[i] != 2 * geno[i]) {
        return false;
      }
    }
  }
  return true;
}

bool AreHapGenoRowCompatible(const SEQUENCE &hapRow, const SEQUENCE &genoRow,
                             SEQUENCE *pComplement) {
  if (pComplement != NULL) {
    pComplement->clear();
  }

  // Check if the haplotype row can be a phasing of the geno row
  YW_ASSERT_INFO(IsSequenceHaplotype(hapRow), "hap is not haplotype row.");
  YW_ASSERT_INFO(IsSequenceGenotype(genoRow), "genorow is not haplotype row.");
  for (unsigned int i = 0; i < hapRow.size(); i++) {
    // if either one is missing value, they match!
    if (IsMissingValueBit(genoRow[i]) == true ||
        IsMissingValueBit(hapRow[i]) == true) {
      continue;
    }

    if (genoRow[i] != 2) {
      if (hapRow[i] != genoRow[i]) {
        return false;
      } else {
        if (pComplement != NULL) {
          pComplement->push_back(genoRow[i]);
        }
      }
    } else {
      if (pComplement != NULL) {
        if (hapRow[i] == 0) {
          pComplement->push_back(1);
        } else {
          pComplement->push_back(0);
        }
      }
    }
  }
  return true;
}

bool AreHapGenoRowsSame(const SEQUENCE &hapRow, const SEQUENCE &genoRow) {
  YW_ASSERT_INFO(IsSequenceHaplotype(hapRow), "hap is not haplotype row.");
  YW_ASSERT_INFO(IsSequenceGenotype(genoRow), "genorow is not haplotype row.");
  return AreTwoSeqsCompatible(hapRow, genoRow);
}

bool IsTrivialRow(const SEQUENCE &row, SEQUENCE &resolved1,
                  SEQUENCE &resolved2) {
  resolved1.clear();
  resolved2.clear();
  // A row is trivial if it contains only a single 2
  YW_ASSERT_INFO(IsSequenceGenotype(row), "hap is not haplotype row.");
  int num2s = 0;
  for (unsigned int i = 0; i < row.size(); ++i) {
    if (row[i] == 2) {
      ++num2s;
      if (num2s > 1) {
        break;
      }
    }
    if (row[i] == 2) {
      resolved1.push_back(0);
      resolved2.push_back(1);
    } else {
      resolved1.push_back(row[i]);
      resolved2.push_back(row[i]);
    }
  }
  if (num2s == 1) {
    // For now, we do not consider a row with no twos as new
    return true;
  } else {
    return false;
  }
}

bool IsHapSeqSmaller(const SEQUENCE &hapRow1, const SEQUENCE &hapRow2) {
  // Decide whether hapRow1 is smaller
  // used in situations when we need to compare two rows
  YW_ASSERT_INFO(IsSequenceHaplotype(hapRow1), "hap1 is not haplotype row.");
  YW_ASSERT_INFO(IsSequenceHaplotype(hapRow2), "hap2 is not haplotype row.");
  YW_ASSERT_INFO(hapRow1.size() == hapRow2.size(),
                 "Tow hap rows are not equal length");
  // do not handle MV in this function
  YW_ASSERT_INFO(IsSeqHasMV(hapRow1) == false && IsSeqHasMV(hapRow2) == false,
                 "Can not handle MV here");

  for (unsigned int i = 0; i < hapRow1.size(); ++i) {
    if (hapRow1[i] < hapRow2[i]) {
      return true;
    }
  }
  return false;
}

void GetHyperCubeSeq(int hcSeq, SEQUENCE &seq, int hcWidth) {
  ConvIntToVecMSB(hcSeq, seq, hcWidth);
}

int GetSeqIdFromSeq(const SEQUENCE &seq) {
  // do not support MV
  YW_ASSERT_INFO(IsSeqHasMV(seq) == false, "Can not support MV");
  return ConvVecToIntMSB(seq);
}

int GetHyperCubSeqBit(int hcSeq, int bit, int hcWidth) {
  // Retrive the bit in the hcSeq at the specified bit
  // But  note that we are assuming bit 0 is on the left (BIG ENDIAN)
  int shiftpos = hcWidth - bit - 1;
  int mask = 0x1 << shiftpos;
  int res = (hcSeq & mask) >> shiftpos;
  YW_ASSERT_INFO(res == 0 || res == 1, "Serious error here.");
  return res;
}

void FindNonSegSites(const set<HCSequence> &setSeqs, set<int> &sites,
                     int dataWidth) {
  // Now we find out whe mutation sites that are not used in the setSeqs
  // Find the set of sites that are not segragating in the set of sequences
  for (int i = 0; i < dataWidth; ++i) {
    bool fZero = false, fOne = false;
    // Check to see if site i is segragating or not
    for (set<int>::iterator it = setSeqs.begin(); it != setSeqs.end(); ++it) {
      int rn = *it;
      // cout <<"In FindNonSegragateSites: rn = " << rn << endl;
      if (GetHyperCubSeqBit(rn, i, dataWidth) == 0) {
        fZero = true;
      } else {
        fOne = true;
      }
      if (fZero && fOne) {
        break;
      }
    }
    if (fZero == false || fOne == false) {
      sites.insert(i);
    }
  }
}

void FindNonSegSites(const set<SEQUENCE> &setSeqs, set<int> &sites,
                     int dataLen) {
  sites.clear();
  if (setSeqs.size() == 0) {
    // Every one is non-segragating
    for (int i = 0; i < dataLen; ++i) {
      sites.insert(i);
    }
    return;
  }

  for (int i = 0; i < dataLen; ++i) {
    bool fZero = false, fOne = false;
    // Check to see if site i is segragating or not
    for (set<SEQUENCE>::iterator it = setSeqs.begin(); it != setSeqs.end();
         ++it) {
      SEQUENCE row = *it;
      YW_ASSERT_INFO(IsSequenceHaplotype(row),
                     "This function only works for haplotype");
      // cout <<"In FindNonSegragateSites: rn = " << rn << endl;
      if (row[i] == 0) {
        fZero = true;
      } else if (row[i] == 1) {
        fOne = true;
      }
      if (fZero && fOne) {
        break;
      }
    }
    if (fZero == false || fOne == false) {
      sites.insert(i);
    }
  }
}

void CreateGenoRowFromHapRows(const SEQUENCE &hapRow1, const SEQUENCE &hapRow2,
                              SEQUENCE &genoRow) {
  // Check if the haplotype row can be a phasing of the geno row
  YW_ASSERT_INFO(IsSequenceHaplotype(hapRow1), "hap1 is not haplotype row.");
  YW_ASSERT_INFO(IsSequenceHaplotype(hapRow2), "hap2 is not haplotype row.");
  // do not allow missing vlaues
  YW_ASSERT_INFO(IsSeqHasMV(hapRow1) == false && IsSeqHasMV(hapRow2) == false,
                 "Can not handle MV");
  genoRow.clear();
  for (unsigned int i = 0; i < hapRow1.size(); i++) {
    if (hapRow1[i] == hapRow2[i]) {
      genoRow.push_back(hapRow1[i]);
    } else {
      genoRow.push_back(2);
    }
  }
}

int IsHCSeqsMutPair(HCSequence seq1, HCSequence seq2, int dataWidth) {
  // This function test if seq1/seq2 is mutation pair, and return the mut site
  // if so
  for (int p = 0; p < dataWidth; p++) {
    int shiftpos = dataWidth - p - 1;
    int mask = 0x1 << shiftpos;

    if ((seq1 | mask) == (seq2 | mask)) {
      return p;
    }
  }

  return -1; // indicate NOT-pair
}

bool IsHCSeqsMutPairAt(HCSequence seq1, HCSequence seq2, int dataWidth,
                       int pos) {
  // Different from the previous function, this check for a specific location
  // instead of trying all possible positions
  int shiftpos = dataWidth - pos - 1;
  int mask = 0x1 << shiftpos;

  if ((seq1 | mask) == (seq2 | mask)) {
    return true;
  }
  return false;
}

void MutateHCSeqAt(const HCSequence seq, HCSequence &res, int dataWidth,
                   int mutPos) {
  int shiftpos = dataWidth - mutPos - 1;
  int mask = 0x1 << shiftpos;

  res = (seq ^ mask);
}

bool IsHCSeqRecombinnable(HCSequence hcSeq1, HCSequence hcSeq2, HCSequence st,
                          int dataWidth) {
  // ASSUME: s1 is LEFT part and s2 is RIGHT part
  // This function test if s1 and s2 can recombine into st
  // Now start recombining
  for (int bkpt = 0; bkpt < dataWidth - 1; ++bkpt) {
    unsigned int maskLower = (0x1 << (bkpt + 1)) - 1;
    unsigned int maskUpper = (0x1 << dataWidth) - 1 - maskLower;

    // Generate s sequence
    int seq1 = ((hcSeq1 & maskLower) | (hcSeq2 & maskUpper));
    if (seq1 == st) {
      return true;
    }
#if 0 // here we
        int seq2 = ((hcSeq1 & maskUpper)  | (hcSeq2 & maskLower) );
        if( seq2 == st)
        {
            return true;
        }
#endif
  }

  return false;
}

void RecombineHCSeqs(const HCSequence hcSeq1, const HCSequence hcSeq2,
                     HCSequence &res, int dataWidth, int bkpt) {
  // ASSUME: s1 is LEFT part and s2 is RIGHT part
  // This function test if s1 and s2 can recombine into st
  // Now start recombining
  unsigned int maskLower = (0x1 << (bkpt + 1)) - 1;
  unsigned int maskUpper = (0x1 << dataWidth) - 1 - maskLower;

  // Generate s sequence
  res = ((hcSeq1 & maskLower) | (hcSeq2 & maskUpper));
}
