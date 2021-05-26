//
//  Utils4.cpp
//
//
//  Created by Yufeng Wu on 9/30/13.
//
//

#include "Utils4.h"
#include <cmath>
#include <fstream>

// template functios goes to header file
#if 0

// utilities
template<class TYPE1, class TYPE2>
void CreateMapFromTwoVec( const vector<TYPE1> &vecKey, const vector<TYPE2> &vecval, map<TYPE1,TYPE2> &mapCreated  )
{
    //
    YW_ASSERT_INFO( vecKey.size() == vecval.size(), "veckey has different size as vecval" );
    mapCreated.clear();
    for( int i=0; i<(int)vecKey.size(); ++i )
    {
        //
        mapCreated.insert(map<TYPE1,TYPE2> :: value_type(vecKey[i], vecval[i] ) );
    }
}

template<class TYPE1,class TYPE2>
void KeepCommonInMaps( map<TYPE1,TYPE2> &mapSubtracted, const map<TYPE1,TYPE2> &mapToSub )
{
    // only keep those that is also in the second map
    map<TYPE1, TYPE2> mapNew;
    for( typename map<TYPE1,TYPE2> :: iterator it = mapSubtracted.begin(); it != mapSubtracted.end(); ++it )
    {
        //
        if( mapToSub.find(it->first) != mapToSub.end() )
        {
            // appear in second map so keep
            mapNew.insert( map<TYPE1,TYPE2> :: value_type(it->first, it->second) );
        }
    }
    mapSubtracted = mapNew;
}

template<class TYPE1,class TYPE2>
void KeepCommonInMapsSet( map<TYPE1,TYPE2> &mapSubtracted, const set<TYPE1> &setKept )
{
    // only keep those that is also in the second map
    map<TYPE1, TYPE2> mapNew;
    for( typename map<TYPE1,TYPE2> :: iterator it = mapSubtracted.begin(); it != mapSubtracted.end(); ++it )
    {
        //
        if( setKept.find(it->first) != setKept.end() )
        {
            // appear in second map so keep
            mapNew.insert( map<TYPE1,TYPE2> :: value_type(it->first, it->second) );
        }
    }
    mapSubtracted = mapNew;

}

template<class TYPE1, class TYPE2>
void CreateTwoVecFromMap(const map<TYPE1,TYPE2> &mapIn, vector<TYPE1> &vecKey, vector<TYPE2> &vecval )
{
    vecKey.clear();
    vecval.clear();
    for( typename map<TYPE1,TYPE2> :: iterator it = mapIn.begin(); it != mapIn.end(); ++it )
    {
        vecKey.push_back( it->first );
        vecval.push_back( it->second );
    }
}

#endif

int GetZeroOneDiff(int x, int y) {
  if (x == y) {
    return 0;
  } else {
    return 1;
  }
}

void GetMatchingPosIntVec(const int val, const vector<int> &listVals,
                          vector<int> &listPos) {
  listPos.clear();
  for (int i = 0; i < (int)listVals.size(); ++i) {
    if (val == listVals[i]) {
      listPos.push_back(i);
    }
  }
}

void FormUnitVector(int numItems, int posUnit, vector<int> &vecUnit) {
  //
  YW_ASSERT_INFO(posUnit < numItems, "Wrong");
  vecUnit.clear();
  for (int i = 0; i < numItems; ++i) {
    //
    vecUnit.push_back(0);
  }
  vecUnit[posUnit] = 1;
}

void FormZeroVector(int numItems, vector<int> &vecZero) {
  //
  vecZero.clear();
  for (int i = 0; i < numItems; ++i) {
    //
    vecZero.push_back(0);
  }
}

bool AreTwoSetsCompatible(const set<int> &set1, const set<int> &set2) {
  // are two sets either disjoint or one contins another
  set<int> sint;
  JoinSets(set1, set2, sint);
  if (sint.size() == 0 || sint.size() == set1.size() ||
      sint.size() == set2.size()) {
    return true;
  }
  return false;
}

bool IsSetCompatibleWithSets(const set<int> &set1,
                             const set<set<int> > &setSets) {
  bool res = true;
  for (set<set<int> >::const_iterator it = setSets.begin(); it != setSets.end();
       ++it) {
    if (AreTwoSetsCompatible(set1, *it) == false) {
      res = false;
      break;
    }
  }
  return res;
}

bool AreTwoSetsCompatible(const set<int> &set1, const set<int> &set2,
                          int numTotalElem) {
  // are two sets either disjoint or one contins another
  set<int> sint;
  JoinSets(set1, set2, sint);
  if (sint.size() == 0 || sint.size() == set1.size() ||
      sint.size() == set2.size()) {
    return true;
  }
  set<int> ssTot = set1;
  UnionSets(ssTot, set2);
  if ((int)ssTot.size() == numTotalElem) {
    return true;
  }
  return false;
}

bool IsSetCompatibleWithSets(const set<int> &set1,
                             const set<set<int> > &setSets, int numTotalElem) {
  bool res = true;
  for (set<set<int> >::const_iterator it = setSets.begin(); it != setSets.end();
       ++it) {
    if (AreTwoSetsCompatible(set1, *it, numTotalElem) == false) {
      res = false;
      break;
    }
  }
  return res;
}

bool IsSignificantFraction(int totNum, int numTypes, int numOneType,
                           double minFrac) {
  // test whether the num of one type occupies a siinficant portion of the
  // totNum (composed of numTypes types)
  if (minFrac >= 0.0) {
    return numOneType >= totNum * minFrac;
  }
  // if not specific fraction is givn, then use the following rule based on the
  // number of types basicallly require appearing at least two times
  return numOneType >= 2;
}

void IncAllNumInSet(set<int> &sint) {
  //
  set<int> res;
  for (set<int>::iterator it = sint.begin(); it != sint.end(); ++it) {
    res.insert(*it + 1);
  }
  sint = res;
}

void DecAllNumInSet(set<int> &sint) {
  //
  set<int> res;
  for (set<int>::iterator it = sint.begin(); it != sint.end(); ++it) {
    res.insert(*it - 1);
  }
  sint = res;
}

void IncAllNumInSets(set<set<int> > &setInts) {
  //
  set<set<int> > res;
  for (set<set<int> >::iterator it = setInts.begin(); it != setInts.end();
       ++it) {
    set<int> sint = *it;
    IncAllNumInSet(sint);
    res.insert(sint);
  }
  setInts = res;
}

void GetNonZeroPosofVec(const vector<int> &vec, set<int> &setpos) {
  //
  setpos.clear();
  for (int i = 0; i < (int)vec.size(); ++i) {
    if (vec[i] != 0) {
      setpos.insert(i);
    }
  }
}

int GetSegIndex(int val, const vector<int> &listSegSizes) {
  //
  int res = -1;
  int szSoFar = 0;
  while (val >= szSoFar && res < (int)listSegSizes.size()) {
    ++res;
    szSoFar += listSegSizes[res];
  }
  return res;
}

// Prob related utilties
double CalcPoisonProb(double rate, int numEvts) {
  //
  double res = exp(-1.0 * rate);
  for (int i = 1; i <= numEvts; ++i) {
    res *= rate / i;
  }
  return res;
}

void GetDiffPosOfTwoVec(const vector<int> &vec1, const vector<int> &vec2,
                        set<int> &setpos) {
  //
  YW_ASSERT_INFO(vec1.size() == vec2.size(), "Size: mismatch");
  setpos.clear();
  for (int i = 0; i < (int)vec1.size(); ++i) {
    if (vec1[i] != vec2[i]) {
      setpos.insert(i);
    }
  }
}

void ComplementBoolVec(vector<bool> &listVals) {
  // T->F and vice versa
  for (int i = 0; i < (int)listVals.size(); ++i) {
    if (listVals[i] == true) {
      listVals[i] = false;
    } else {
      listVals[i] = true;
    }
  }
}

void GetAllGridPoints(int gridLB, int gridUB, int dimGrid,
                      set<vector<int> > &setGridPts) {
  //  get all grid points whose num is within the range [lb,ub]
  YW_ASSERT_INFO(gridLB <= gridUB, "Bounds wrong");
  YW_ASSERT_INFO(dimGrid >= 1, "Dimension must be positive");
  // apply recurrence
  setGridPts.clear();
  if (dimGrid == 1) {
    for (int v = gridLB; v <= gridUB; ++v) {
      //
      vector<int> vec;
      vec.push_back(v);
      setGridPts.insert(vec);
    }
  } else {
    //
    set<vector<int> > setGridPtsSmall;
    GetAllGridPoints(gridLB, gridUB, dimGrid - 1, setGridPtsSmall);
    for (set<vector<int> >::iterator it = setGridPtsSmall.begin();
         it != setGridPtsSmall.end(); ++it) {
      //
      for (int v = gridLB; v <= gridUB; ++v) {
        //
        vector<int> vec = *it;
        vec.push_back(v);
        setGridPts.insert(vec);
      }
    }
  }
}

void MapIntListToAnother(const vector<int> &vec1, const vector<int> &vec2,
                         map<int, int> &mapVec1IndexToVec2) {
  // given two vectors, e.g. vec1 = [2,1,3] and vec2 = [3,2,1]. Create a map
  // from vec1's index to vec2 map = [0,1], [1,2], [2,0] we assume there is no
  // dupllicate for now
  // cout << "MapIntListToAnother: vec1: ";
  // DumpIntVec(vec1);
  // cout << "vec2: ";
  // DumpIntVec(vec2);
  mapVec1IndexToVec2.clear();
  YW_ASSERT_INFO(vec1.size() == vec2.size(), "size: mismatch");
  map<int, int> mapValToIndex1, mapValToIndex2;
  for (int i = 0; i < (int)vec1.size(); ++i) {
    //
    YW_ASSERT_INFO(mapValToIndex1.find(vec1[i]) == mapValToIndex1.end(),
                   "Duplicate found");
    mapValToIndex1.insert(map<int, int>::value_type(vec1[i], i));
    // cout << "mapValToIndex1: " << vec1[i] << ", " << i << endl;
  }
  for (int i = 0; i < (int)vec2.size(); ++i) {
    //
    YW_ASSERT_INFO(mapValToIndex2.find(vec2[i]) == mapValToIndex2.end(),
                   "Duplicate found");
    mapValToIndex2.insert(map<int, int>::value_type(vec2[i], i));
    // cout << "mapValToIndex12 " << vec2[i] << ", " << i << endl;
  }
  for (map<int, int>::iterator it = mapValToIndex1.begin();
       it != mapValToIndex1.end(); ++it) {
    YW_ASSERT_INFO(mapValToIndex2.find(it->first) != mapVec1IndexToVec2.end(),
                   "Two lists: not idential");
    mapVec1IndexToVec2.insert(
        map<int, int>::value_type(it->second, mapVec1IndexToVec2[it->first]));
  }
}

void FindEvenDistriPoints(double valMin, double valMax, double valResolution,
                          int maxNumPoints, vector<double> &listChosenVals) {
  // pick uniformly some number (<= maxNumPoints) of points within [valMin,
  // valMax}, with distance no more than resolution first figure out spacing
  double valSpacing = (valMax - valMin) / maxNumPoints;
  if (valSpacing < valResolution) {
    valSpacing = valResolution;
  }
  for (int i = 0; i < (int)(valMax - valMin) / valSpacing; ++i) {
    //
    double val = (i + 0.5) * valSpacing;
    listChosenVals.push_back(val);
  }
}

// bits operation
bool IsBitSetInt(int val, int posBit) {
  //
  // for an index of AC, which src populaiton is a leave
  int mask = (0x1 << posBit);
  // assume only two populaitons for now
  bool res = false;
  if ((val & mask) != 0) {
    res = true;
  }
  return res;
}

int ToggleBitInt(int val, int posBit) {
  //
  return val ^ (1 << posBit);
}

double StrToDouble(const string &s) {
  double d;
  stringstream ss(s); // turn the string into a stream
  ss >> d;            // convert
  return d;
}

double CalcProductBetween(int lb, int ub) {
  double res = 1.0;
  for (int i = lb; i <= ub; ++i) {
    res *= i;
  }
  return res;
}

void CreateClustersFromMultisets(
    const multiset<multiset<int> > &setMultisets,
    map<multiset<int>, vector<multiset<int> > > &mapMultisetClusters) {
  cout << "CreateClustersFromMultisets: DONOT WORK YET\n";
  // give multisets S1, S2, .... Sn
  // find the ancestral (clustering) relations among them
  // i.e. if S1 contains S2 and S3 (as the smallest enclosing), then we have: S1
  // (S2, S3) and so on
  mapMultisetClusters.clear();

  // this refers to the smallest container set
  map<multiset<int>, multiset<int> > mapSmallestContainer;

  //
  // YW: TBD: issue: there may be duplicate clusters
  // TBDDDDDDDDDDDDDDDDDDD

  for (multiset<multiset<int> >::const_iterator it1 = setMultisets.begin();
       it1 != setMultisets.end(); ++it1) {
    //
    for (multiset<multiset<int> >::const_iterator it2 = setMultisets.begin();
         it2 != setMultisets.end(); ++it2) {
      //
      if (it1 == it2) {
        continue;
      }

      // is s1 contained by s2? and also cannot allow the two becomes the same
      if (IsMultisetContainedIn(*it1, *it2) == true &&
          it1->size() < it2->size()) {
        //
        if (mapSmallestContainer.find(*it1) == mapSmallestContainer.end()) {
          mapSmallestContainer.insert(
              map<multiset<int>, multiset<int> >::value_type(*it1, *it2));
        } else if (mapSmallestContainer[*it1].size() > it2->size()) {
          mapSmallestContainer[*it1] = *it2;
        }
      }
    }
  }
  cout << "here...\n";
  // now from the smallest container, create the clusters
  for (map<multiset<int>, multiset<int> >::iterator it =
           mapSmallestContainer.begin();
       it != mapSmallestContainer.end(); ++it) {
    if (mapMultisetClusters.find(it->second) == mapMultisetClusters.end()) {
      vector<multiset<int> > listMSs;
      mapMultisetClusters.insert(
          map<multiset<int>, vector<multiset<int> > >::value_type(it->second,
                                                                  listMSs));
    }
    mapMultisetClusters[it->second].push_back(it->first);
  }
}

void CountMultiset(const multiset<int> &s1, map<int, int> &msMap) {
  for (multiset<int>::const_iterator it = s1.begin(); it != s1.end(); ++it) {
    if (msMap.find(*it) == msMap.end()) {
      msMap.insert(map<int, int>::value_type(*it, 0));
    }
    ++msMap[*it];
  }
}

bool IsMultisetContainedIn(const multiset<int> &s1, const multiset<int> &s2) {
  map<int, int> msMap1, msMap2;
  CountMultiset(s1, msMap1);
  CountMultiset(s2, msMap2);
  for (map<int, int>::iterator it1 = msMap1.begin(); it1 != msMap1.end();
       ++it1) {
    if (msMap2.find(it1->first) == msMap2.end() ||
        it1->second > msMap2[it1->first]) {
      return false;
    }
  }
  return true;
}

void DumpIntMultiset(const multiset<int> &ms) {
  for (multiset<int>::const_iterator it = ms.begin(); it != ms.end(); ++it) {
    cout << *it << "  ";
  }
  cout << endl;
}

void OutputStringsToFile(const char *filename,
                         const vector<string> &listStrsOut) {
  ofstream outFile(filename);
  if (outFile.is_open() == false) {
    cout << "Fatal error: Can not open output file: " << filename << endl;
    exit(1);
  }

  for (int i = 0; i < (int)listStrsOut.size(); ++i) {
    outFile << listStrsOut[i] << endl;
  }
  outFile.close();
}

unsigned int ConvVecToIntGen(const vector<int> &vec, int base) {
  // assume vec[0] is least siginicant
  unsigned int res = 0;

  for (int i = (int)vec.size() - 1; i >= 0; --i) {
    YW_ASSERT_INFO(vec[i] >= 0 && vec[i] < base,
                   "In ConvVecToIntGen, vector value overflow.");
    // cout << "res = " << res << endl;

    res += vec[i];
    if (i > 0) {
      res = res * base;
    }
  }

  return res;
}

unsigned int ConvVecToIntGenMSB(const vector<int> &vec, int base) {
  vector<int> vecMSB = vec;
  // cout << "vec = ";
  // DumpIntVec( vec );
  ReverseIntVec(vecMSB);
  // cout << "vec = ";
  // DumpIntVec( vec );
  return ConvVecToIntGen(vecMSB, base);
}

int ConvVecToIntGenBounds(const vector<int> &vec, const vector<int> &bounds) {
  // bound[i]: the largest value a digit can reach at position i
  // assume vec[0] is least siginicant
  unsigned int res = 0;

  for (int i = (int)vec.size() - 1; i >= 0; --i) {
    YW_ASSERT_INFO(vec[i] >= 0 && vec[i] <= bounds[i],
                   "In ConvVecToIntGen, vector value overflow.");
    // cout << "res = " << res << endl;

    res += vec[i];
    if (i > 0) {
      res = res * (bounds[i - 1] + 1);
    }
  }

  return res;
}

void ConvIntToVecGen(int val, const vector<int> &bounds, vector<int> &vec) {
  vec.clear();

  int numBits = bounds.size();
  YW_ASSERT_INFO(numBits < 30, "Overflow000");

  // we would store the least significant bit as vec[0]
  for (int i = 0; i < numBits; ++i) {
    int bound0 = bounds[i];
    YW_ASSERT_INFO(bound0 >= 0, "Cannot be too small");
    int val2 = val % (bound0 + 1);
    vec.push_back(val2);
    val = (val - val2) / (bound0 + 1);
  }
}

int ConvRowMajorPosVecToIntGenBounds(const vector<int> &vec,
                                     const vector<int> &bounds) {
  // different from above: bound b means that max value is actaully b-1 (like
  // those) bound[i]: the largest value a digit can reach at position i assume
  // vec[0] is least siginicant
  unsigned int res = 0;

  for (int i = 0; i < (int)vec.size(); ++i) {
    if (i > 0) {
      res = res * (bounds[i]);
    }
    YW_ASSERT_INFO(vec[i] >= 0 && vec[i] <= bounds[i],
                   "In ConvVecToIntGen, vector value overflow.");
    // cout << "res = " << res << endl;
    res += vec[i];
  }

  return res;
}

void ConvRowMajorIntPosToVecGen(int val, const vector<int> &bounds,
                                vector<int> &vec) {
  //
  vec.clear();

  int numBits = bounds.size();
  YW_ASSERT_INFO(numBits < 30, "Overflow000");

  // we would store the least significant bit as vec[0]
  for (int i = numBits - 1; i >= 0; --i) {
    int bound0 = bounds[i];
    YW_ASSERT_INFO(bound0 >= 1, "Cannot be too small");
    int val2 = val % (bound0);
    vec.push_back(val2);
    val = (val - val2) / (bound0);
  }
  ReverseIntVec(vec);
}

// utility
class ClusterPosition {
public:
  ClusterPosition() { pos = 0; }
  ClusterPosition(const ClusterPosition &rhs) : pos(rhs.pos) {}
  ClusterPosition(int posIn) { pos = posIn; }
  int GetPosition() const { return pos; }

private:
  int pos;
};

void ClusterLinearPoints(const vector<double> &listPoints,
                         double ratioMaxInOutCmp, vector<int> &listBkpts) {
  // assume points are sorted!!!
  if (listPoints.size() <= 1) {
    // nothing to cluster
    return;
  }

  // rationInOutCmp: the max ratio btwn inside cluster and outside cluster that
  // we will merge two groups
  map<pair<int, int>, double> mapClusterInfo; // current max distance within
                                              // group
  map<int, pair<int, int> > mapPointMembership;

  // init each point to self
  for (int i = 0; i < (int)listPoints.size(); ++i) {
    pair<int, int> pp(i, i);
    mapClusterInfo.insert(map<pair<int, int>, double>::value_type(pp, 0.0));
    mapPointMembership.insert(map<int, pair<int, int> >::value_type(i, pp));
  }
  // sort the values
  vector<ClusterPosition> vecPosRecords;
  for (int i = 0; i < (int)listPoints.size(); ++i) {
    ClusterPosition cp(i);
    vecPosRecords.push_back(cp);
  }
  vector<pair<double, void *> > listPointsSortedWithPos;
  for (int i = 0; i < (int)listPoints.size() - 1; ++i) {
    pair<double, void *> pp(listPoints[i + 1] - listPoints[i],
                            &vecPosRecords[i]);
    listPointsSortedWithPos.push_back(pp);
  }
  SortPairsByNumsDouble(listPointsSortedWithPos);
  for (int i = 0; i < (int)listPointsSortedWithPos.size(); ++i) {
    //
    double diststep = listPointsSortedWithPos[i].first;
    ClusterPosition *ptr =
        (ClusterPosition *)(listPointsSortedWithPos[i].second);
    int pos = ptr->GetPosition();
    int posNext = pos + 1;
    YW_ASSERT_INFO(mapPointMembership.find(pos) != mapPointMembership.end(),
                   "Fail");
    YW_ASSERT_INFO(mapPointMembership.find(posNext) != mapPointMembership.end(),
                   "Fail");
    pair<int, int> pp1 = mapPointMembership[pos];
    pair<int, int> pp2 = mapPointMembership[posNext];
    // should we merge the two; do so if the current distance
    bool fMerge1 = true;
    if (pp1.second > pp1.first) {
      YW_ASSERT_INFO(mapClusterInfo.find(pp1) != mapClusterInfo.end(),
                     "Fail to find");
      double distCur = mapClusterInfo[pp1];
      if (diststep <= distCur * ratioMaxInOutCmp) {
        fMerge1 = true;
      } else {
        fMerge1 = false;
      }
    }
    bool fMerge2 = true;
    if (pp2.second > pp2.first) {
      YW_ASSERT_INFO(mapClusterInfo.find(pp2) != mapClusterInfo.end(),
                     "Fail to find");
      double distCur = mapClusterInfo[pp2];
      if (diststep <= distCur * ratioMaxInOutCmp) {
        fMerge2 = true;
      } else {
        fMerge2 = false;
      }
    }
    if (fMerge1 && fMerge2) {
      cout << "Merging: (" << pp1.first << ", " << pp1.second << "): and ("
           << pp2.first << "," << pp2.second << ")\n";
      // merge
      pair<int, int> ppnew(pp1.first, pp2.second);
      double distMaxNew = std::max(
          diststep, std::max(mapClusterInfo[pp1], mapClusterInfo[pp2]));
      mapClusterInfo.insert(
          map<pair<int, int>, double>::value_type(ppnew, distMaxNew));
      mapClusterInfo.erase(pp1);
      mapClusterInfo.erase(pp2);
      for (int s = ppnew.first; s <= ppnew.second; ++s) {
        mapPointMembership.erase(s);
        mapPointMembership.insert(
            map<int, pair<int, int> >::value_type(s, ppnew));
      }
    }
  }
  // now insert all segments
  for (map<pair<int, int>, double>::iterator it = mapClusterInfo.begin();
       it != mapClusterInfo.end(); ++it) {
    int bkptRight = it->first.second;
    if (bkptRight < (int)listPoints.size() - 1) {
      listBkpts.push_back(bkptRight);
    }
  }
}

void FindConsecutiveIntervals(const set<int> &setItems,
                              vector<pair<int, int> > &listIVs) {
  listIVs.clear();
  if (setItems.size() == 0) {
    return;
  }
  int itemStart = *setItems.begin();
  int itemPrev = itemStart;
  set<int>::const_iterator it = setItems.begin();
  ++it;
  while (it != setItems.end()) {
    if (*it != itemPrev + 1) {
      // this is an IV
      pair<int, int> pp(itemStart, itemPrev);
      listIVs.push_back(pp);
      itemStart = *it;
    }

    itemPrev = *it;
    ++it;
    if (it == setItems.end()) {
      // ouput the prev
      pair<int, int> pp(itemStart, itemPrev);
      listIVs.push_back(pp);
    }
  }
}

void ComplementIntSet(int numTot, set<int> &setToComp) {
  // YW: assume numbers start from 0 to numTot-1
  set<int> ssTot;
  PopulateSetWithInterval(ssTot, 0, numTot - 1);
  SubtractSets(ssTot, setToComp);
  setToComp = ssTot;
}

void GetCountsItems(int range, const set<int> &listNumbers,
                    vector<int> &listCnts) {
  // count occurance of numbers: listCnts[k] = # of items that is smaller or
  // equal to k in the set
  YW_ASSERT_INFO(range >= 0, "Must be positive");
  listCnts.clear();
  listCnts.resize(range + 1);
  int cntTot = 0;
  int posLast = -1;
  for (set<int>::const_iterator it = listNumbers.begin();
       it != listNumbers.end(); ++it) {
    int val = *it;
    YW_ASSERT_INFO(val <= range, "Wrong");
    for (int i = posLast + 1; i < val; ++i) {
      listCnts[i] = cntTot;
    }
    ++cntTot;
    listCnts[val] = cntTot;
    posLast = val;
  }
}

void FindGapBlocksWithinPosVec(const vector<int> &posvec, int numItemsEnum,
                               int numItemsGap,
                               vector<pair<int, int> > &listSegs) {
  // in a position vector (i.e. subset of positions 0, 1, ..., k; find gaps in
  // between the chosen positions gaps are re-ordered to consecutive from 0, 1,
  // ...
  listSegs.clear();
  vector<int> listGapLens;
  for (int i = 0; i < (int)posvec.size(); ++i) {
    int posLast = -1;
    if (i > 0) {
      posLast = posvec[i - 1];
    }
    int len = posvec[i] - posLast - 1;
    listGapLens.push_back(len);
  }
  // cout << "numItemsEnum: " << numItemsEnum << ", listGapLens: ";
  // DumpIntVec(listGapLens);
  // last segment
  int posFinal = numItemsEnum + numItemsGap - 1;
  int posFirst = 0;
  if (posvec.size() > 0) {
    posFirst = posvec[posvec.size() - 1];
  } else {
    posFinal = numItemsGap;
  }
  int lenFinal = posFinal - posFirst;
  // cout << "posFirst: " << posFirst << ", posFinal: " << posFinal << ",
  // lenFinal: " << lenFinal << endl;
  YW_ASSERT_INFO(lenFinal >= 0, "Cannot be negative");
  listGapLens.push_back(lenFinal);
  int posCur = 0;
  for (int i = 0; i < (int)listGapLens.size(); ++i) {
    pair<int, int> pp;
    pp.first = posCur;
    pp.second = posCur + listGapLens[i];

    if (pp.first > numItemsGap) {
      pp.first = -1;
    }
    if (pp.second > numItemsGap) {
      pp.second = -1;
    }

    listSegs.push_back(pp);

    // note: consecutive IV overlaps
    posCur = pp.second;
  }
}

void GetSetsIntParts(const set<int> &set1, const set<int> &set2,
                     const set<int> &setAll, set<int> &set1Only,
                     set<int> &set2Only, set<int> &set12, set<int> &setNone) {
  //
  set1Only = set1;
  SubtractSets(set1Only, set2);
  set2Only = set2;
  SubtractSets(set2Only, set1);
  set12 = set1;
  UnionSets(set12, set2);
  setNone = setAll;
  SubtractSets(setNone, set12);
}
