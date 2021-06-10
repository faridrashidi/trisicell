//
//  ScistGenotype.cpp
//
//
//  Created by Yufeng Wu on 5/25/18.
//
//

#include "ScistGenotype.hpp"
#include "MarginalTree.h"
#include "PhylogenyTree.h"
#include "RerootTreeUtils.h"
#include "TreeBuilder.h"
#include "Utils3.h"
#include "Utils4.h"
#include <cmath>
#include <iomanip>

// *************************************************************************************
// genotypes: integer matrix

ScistGenGenotypeMat::ScistGenGenotypeMat() : thresSignifcant(0.0) {}

void ScistGenGenotypeMat ::TrimCliquesMaxDiff(
    std::set<std::set<int> > &listCliques, int maxToKeep) const {
  // cout << "Entering trim, number of cliques: " << listCliques.size() << ",
  // maxToKeep: " << maxToKeep << endl;
  // keep only the most different ones
  if ((int)listCliques.size() <= maxToKeep) {
    return;
  }
  // find the distance between two sets
  map<pair<const set<int> *, const set<int> *>, int> mapPairCliqueDiff;
  for (set<set<int> >::iterator it1 = listCliques.begin();
       it1 != listCliques.end(); ++it1) {
    set<set<int> >::iterator it2 = it1;
    ++it2;
    for (; it2 != listCliques.end(); ++it2) {
      //
      set<int> sint;
      JoinSets(*it1, *it2, sint);
      pair<const set<int> *, const set<int> *> pp1(&(*it1), &(*it2)),
          pp2(&(*it2), &(*it1));
      mapPairCliqueDiff[pp1] = it1->size() + it2->size() - 2 * sint.size();
      mapPairCliqueDiff[pp2] = mapPairCliqueDiff[pp1];
    }
  }

  // increamentally add the most different; first add the first clique
  set<const set<int> *> listCliquesNext;
  listCliquesNext.insert(&(*listCliques.begin()));
  while ((int)listCliquesNext.size() < maxToKeep) {
    const set<int> *pcliqueToAdd = NULL;
    int diffMax = 0;
    for (set<set<int> >::iterator it1 = listCliques.begin();
         it1 != listCliques.end(); ++it1) {
      //
      if (listCliquesNext.find(&(*it1)) != listCliquesNext.end()) {
        //
        continue;
      }

      //
      int diffCurr = 0;
      for (set<const set<int> *>::iterator it2 = listCliquesNext.begin();
           it2 != listCliquesNext.end(); ++it2) {
        pair<const set<int> *, const set<int> *> pp(*it2, &(*it1));
        YW_ASSERT_INFO(mapPairCliqueDiff.find(pp) != mapPairCliqueDiff.end(),
                       "Fail to find");
        diffCurr += mapPairCliqueDiff[pp];
      }
      if (diffCurr > diffMax) {
        diffMax = diffCurr;
        pcliqueToAdd = &(*it1);
      }
    }
    YW_ASSERT_INFO(pcliqueToAdd != NULL, "Cannot be null");
    listCliquesNext.insert(pcliqueToAdd);
    // cout << "In TrimCliquesMaxDiff: adding clique: ";
    // DumpIntSet(*pcliqueToAdd);
  }
  set<set<int> > listCliquesNextUse;
  for (set<const set<int> *>::iterator it = listCliquesNext.begin();
       it != listCliquesNext.end(); ++it) {
    listCliquesNextUse.insert(*(*it));
  }
  listCliques = listCliquesNextUse;
}

ScistGenGenotypeMat *
ScistGenGenotypeMat ::SubMatrix(const std::set<int> &setRows,
                                const std::set<int> &setSites) const {
  ScistGenGenotypeMat *pMatNew = CreateNewMat();
  pMatNew->SetSize(setRows.size(), setSites.size());
  // set row name
  int rowCurr = 0;
  for (set<int>::iterator it = setRows.begin(); it != setRows.end(); ++it) {
    int siteCurr = 0;
    for (set<int>::iterator it2 = setSites.begin(); it2 != setSites.end();
         ++it2) {
      pMatNew->SetGenotypeAt(rowCurr, siteCurr, GetGenotypeAt(*it, *it2));
      pMatNew->SetGenotypeProbAt(rowCurr, siteCurr,
                                 GetGenotypeProbAllele0At(*it, *it2));
      ++siteCurr;
    }

    pMatNew->SetGenotypeName(rowCurr, GetGenotypeName(*it));
    ++rowCurr;
  }
  return pMatNew;
}

std::string ScistGenGenotypeMat ::ConsNJTree() const {
  //
  PhyloDistance dist;
  // setup pairwise hamming distance
  for (int i = 0; i < GetNumHaps(); ++i) {
    for (int j = i + 1; j < GetNumHaps(); ++j) {
      //
      double d = CalcHammingDistBetwHaps(i, j);
      dist.SetDistance(i, j, d);
      // cout << "Distance between (" << i << "," << j << "): " << d << endl;
    }
  }
  DistanceTreeBuilder dtb(dist);
  for (int i = 0; i < GetNumHaps(); ++i) {
    int indexUse = i + 1;
    string strIndexToUse = std::to_string(indexUse);
    dtb.SetTaxonName(i, strIndexToUse);
  }
  return dtb.NJ();
}

std::string ScistGenGenotypeMat ::ConsNJTreeZeroRoot() const {
  //
  PhyloDistance dist;
  // setup pairwise hamming distance
  for (int i = 0; i < GetNumHaps(); ++i) {
    for (int j = i + 1; j < GetNumHaps(); ++j) {
      //
      double d = CalcHammingDistBetwHaps(i, j);
      dist.SetDistance(i, j, d);
      // cout << "Distance between (" << i << "," << j << "): " << d << endl;
    }
  }
  // add one more hap: all-0
  for (int i = 0; i < GetNumHaps(); ++i) {
    //
    double d = 0.0;
    for (int s = 0; s < GetNumSites(); ++s) {
      if (GetGenotypeAt(i, s) != 0) {
        d += 1.0;
      }
    }
    d = d / GetNumSites();
    dist.SetDistance(i, GetNumHaps(), d);
  }

  DistanceTreeBuilder dtb(dist);
  for (int i = 0; i <= GetNumHaps(); ++i) {
    int indexUse = i + 1;
    string strIndexToUse = std::to_string(indexUse);
    dtb.SetTaxonName(i, strIndexToUse);
  }
  string strNJWithRoot = dtb.NJ();
  // cout << "strNJWithRoot: " << strNJWithRoot << endl;
  // reroot
  string strIdRoot = std::to_string(GetNumHaps() + 1);
  char strNJWithRootBuf[102400];
  strcpy(strNJWithRootBuf, strNJWithRoot.c_str());
  char strIdRootBuf[102400];
  strcpy(strIdRootBuf, strIdRoot.c_str());
  string strNJWithRootReroot = ReRootTreeNewick(strNJWithRootBuf, strIdRootBuf);
  // cout << "strNJWithRootReroot: " << strNJWithRootReroot << endl;
  // remove the root
  MarginalTree mtree;
  ReadinMarginalTreesNewickWLenString(strNJWithRootReroot,
                                      this->GetNumHaps() + 1, mtree);
  mtree.BuildDescendantInfo();
  int posRootLeaf = mtree.GetPosForLabel(this->GetNumHaps() + 1);
  YW_ASSERT_INFO(posRootLeaf >= 0, "Fail to find the root");
  mtree.RemoveLeafNodeFromBinaryTree(posRootLeaf);
  mtree.BuildDescendantInfo();
  // cout << "Aftre removing reoot: " << mtree.GetNewickSorted(false) << endl;
  return mtree.GetNewickSorted(false);
}

std::string ScistGenGenotypeMat ::ConsNJTreeNoInc() const {
  PhyloDistance dist;
  // setup pairwise hamming distance
  for (int i = 0; i < GetNumHaps(); ++i) {
    for (int j = i + 1; j < GetNumHaps(); ++j) {
      //
      double d = CalcHammingDistBetwHaps(i, j);
      dist.SetDistance(i, j, d);
      // cout << "Distance between (" << i << "," << j << "): " << d << endl;
    }
  }
  DistanceTreeBuilder dtb(dist);
  return dtb.NJ();
}

double ScistGenGenotypeMat ::CalcHammingDistBetwHaps(int h1, int h2) const {
  int numDiffs = 0;
  for (int c = 0; c < GetNumSites(); ++c) {
    if (GetGenotypeAt(h1, c) != GetGenotypeAt(h2, c) &&
        IsProbAtCellPosSignificant(h1, c, GetSignificanceThres()) &&
        IsProbAtCellPosSignificant(h2, c, GetSignificanceThres())) {
      ++numDiffs;
    }
  }
  return (1.0 * numDiffs) / GetNumSites();
}

void ScistGenGenotypeMat ::ConsCompatMap(
    std::set<std::pair<int, int> > &setCompatPairs) const {
  //
  setCompatPairs.clear();
  for (int s1 = 0; s1 < GetNumSites(); ++s1) {
    for (int s2 = s1 + 1; s2 < GetNumSites(); ++s2) {
      if (IsCompatible(s1, s2)) {
        pair<int, int> pp(s1, s2);
        setCompatPairs.insert(pp);
      }
    }
  }
}

bool ScistGenGenotypeMat ::AreSitesCompatInMap(
    const std::set<std::pair<int, int> > &setCompatPairs, int s1, int s2) {
  //
  pair<int, int> pp(s1, s2);
  OrderInt(pp.first, pp.second);
  return setCompatPairs.find(pp) != setCompatPairs.end();
}

int ScistGenGenotypeMat ::GetGenotypeNumOf(int geno) const {
  int res = 0;
  for (int i = 0; i < GetNumHaps(); ++i) {
    for (int j = 0; j < GetNumSites(); ++j) {
      if (GetGenotypeAt(i, j) == geno) {
        ++res;
      }
    }
  }
  return res;
}

int ScistGenGenotypeMat ::FindCellByName(const std::string &strName) const {
  //
  for (int i = 0; i < GetNumHaps(); ++i) {
    if (GetGenotypeName(i) == strName) {
      return i;
    }
  }
  return -1;
}

void ScistGenGenotypeMat ::Dump() const {
  cout << "Genotype names: ";
  for (int i = 0; i < GetNumHaps(); ++i) {
    cout << GetGenotypeName(i) << "  ";
  }
  cout << endl;
}

void ScistGenGenotypeMat ::ChangeGenosAtPositions(
    const std::set<std::pair<std::pair<int, int>, int> > &listChangedPlaces) {
  //
  for (std::set<std::pair<std::pair<int, int>, int> >::const_iterator it =
           listChangedPlaces.begin();
       it != listChangedPlaces.end(); ++it) {
    //
    SetGenotypeAt(it->first.first, it->first.second, it->second);
  }
}

// *************************************************************************************
// genotypes: binary matrix

ScistHaplotypeMat ::ScistHaplotypeMat() {}

ScistGenGenotypeMat *ScistHaplotypeMat ::Copy() const {
  //
  ScistHaplotypeMat *pMatCopy = new ScistHaplotypeMat();

  for (int i = 0; i < GetNumNames(); ++i) {
    pMatCopy->AddGenotypeName(GetGenotypeName(i));
  }

  pMatCopy->SetSize(GetNumHaps(), GetNumSites());

  //
  for (int i = 0; i < GetNumHaps(); ++i) {
    for (int j = 0; j < GetNumSites(); ++j) {
      pMatCopy->SetGenotypeAt(i, j, GetGenotypeAt(i, j));
      pMatCopy->SetGenotypeProbAt(i, j, GetGenotypeProbAllele0At(i, j));
    }
  }

  return pMatCopy;
}

bool ScistHaplotypeMat ::ReadFromFile(std::ifstream &infile, int numSites,
                                      int numSCs, bool fSiteName) {
  // cout << "ScistHaplotypeMat :: ReadFromFile: numSites: " << numSites << ",
  // numSCs: " << numSCs << endl;
  //
  // assume each site is independent
  SetSize(numSCs, numSites);
  for (int i = 0; i < numSites; ++i) {
    string strName;
    if (fSiteName) {
      infile >> strName;
    } else {
      strName = std::to_string(i + 1);
    }
    AddSiteName(strName);

    // cout << "Read in site: " << i << endl;
    for (int j = 0; j < numSCs; ++j) {
      double prob0 = 0.0;
      bool res = ReadFromFileHapProb(infile, prob0);
      if (res == false) {
        return false;
      }
      // choose the allele w/ higher prob
      int allele = 0;
      if (prob0 < 0.5) {
        allele = 1;
      }
      SetGenotypeAt(j, i, allele);

      matHaplotypesProb0[j][i] = prob0;
    }
  }

  // cout << "Input matrix: ";
  // this->matHaplotypes.Dump();

  return true;
}

bool ScistHaplotypeMat ::ReadFromFileHapProb(std::ifstream &infile,
                                             double &prob0) {
  // read in the prob of haploid allele: 0.6 means prob of 0 is 0.6
  // assume prob of 0 + prob of 1 = 1
  infile >> prob0;
  return true;
}

void ScistHaplotypeMat ::SetSize(int numHaps, int numSites) {
  matHaplotypes.SetSize(numHaps, numSites);

  matHaplotypesProb0.clear();
  matHaplotypesProb0.resize(numHaps);

  bool fNameInit = GetNumNames() > 0;

  for (int i = 0; i < numHaps; ++i) {
    matHaplotypesProb0[i].resize(numSites);

    // by default, use the numericals, starting from one
    if (fNameInit == false) {
      string str = std::to_string(i + 1);
      AddGenotypeName(str);

      // cout << "Init name: " << str << endl;
    }
  }
}

void ScistHaplotypeMat ::SetGenotypeAt(int sc, int site, int geno) {
  matHaplotypes(sc, site) = geno;
}

void ScistHaplotypeMat ::AddGenotypeAt(int sc, int site, int geno) {
  // append the genotype into it
  int genoThis = GetGenotypeAt(sc, site);
  if (genoThis == 0 && geno == 1) {
    SetGenotypeAt(sc, site, 1);
  }
}

int ScistHaplotypeMat ::GetAltGenotypeAt(int sc, int site) const {
  int genoThis = GetGenotypeAt(sc, site);
  if (genoThis == 0) {
    return 1;
  } else {
    return 0;
  }
}

double ScistHaplotypeMat ::GetGenotypeProbAllele0At(int sc, int site) const {
  // return proble of allele 0
  return this->matHaplotypesProb0[sc][site];
}

void ScistHaplotypeMat ::SetGenotypeProbAt(int sc, int site, double prob) {
  this->matHaplotypesProb0[sc][site] = prob;
}

void ScistHaplotypeMat ::SetGenotypeProbOfGenoAt(int sc, int site, int geno,
                                                 double prob) {
  if (geno == 0) {
    SetGenotypeProbAt(sc, site, prob);
  } else {
    SetGenotypeProbAt(sc, site, 1.0 - prob);
  }
}

int ScistHaplotypeMat ::GetGenotypeAt(int sc, int site) const {
  return matHaplotypes(sc, site);
}

void ScistHaplotypeMat ::FindMaximalCompatSites(
    const std::vector<double> &wtSites,
    std::vector<std::map<int, std::set<int> > > &listSetSitesCompat,
    int maxNumSets,
    const std::set<std::pair<int, int> > *pSetCompatPairs) const {
  //#if 0
  // const double DEF_MIN_FRAC = 0.5;

  // we find the maximum weightd clique of compatible pairs
  // construct compat pairs if not done yet
  set<std::pair<int, int> > *pSetCompatPairsUse =
      const_cast<set<std::pair<int, int> > *>(pSetCompatPairs);
  set<pair<int, int> > setCompatPairsAlt;
  if (pSetCompatPairsUse == NULL) {
    ConsCompatMap(setCompatPairsAlt);
    pSetCompatPairsUse = &setCompatPairsAlt;
  }

  // implement the simple heuristics by Johnson 1974
  // BinaryMatrix &matHaplotypesUse = const_cast<BinaryMatrix &>(
  // this->matHaplotypes );

  //
  listSetSitesCompat.clear();
  // vector<vector<bool> > vecHapsFullyCompat( GetNumSites() );
  // for(int i=0; i<GetNumSites(); ++i)
  //{
  //    vecHapsFullyCompat[i].resize( GetNumSites() );
  //}

  // for(int s1 = 0; s1<GetNumSites(); ++s1)
  //{
  //    vecHapsFullyCompat[s1][s1] = true;
  //    for(int s2=s1+1; s2<GetNumSites(); ++s2)
  //    {
  //        // root allele: 0
  //        bool fCompat = matHaplotypesUse.IsCompatibleRooted(s1, s2, 0, 0);
  //        vecHapsFullyCompat[s1][s2] = fCompat;
  //        vecHapsFullyCompat[s2][s1] = fCompat;
  // cout << "Sites " << s1 << "," << s2 << ": ";
  // if(fCompat)
  //{
  // cout << " compatible\n";
  //}
  // else
  //{
  // cout << " not compatible\n";
  //}
  //    }
  //}

  //
  set<pair<set<int>, set<int> > > listSetMaxCompatChosen;
  // init
  set<int> ss;
  set<int> setSitesRemainInit;
  PopulateSetWithInterval(setSitesRemainInit, 0, GetNumSites() - 1);
  pair<set<int>, set<int> > pp(ss, setSitesRemainInit);
  listSetMaxCompatChosen.insert(pp);

  while (true) {
    //
    set<pair<set<int>, set<int> > > listSetMaxCompatChosenNext;

    for (set<pair<set<int>, set<int> > >::iterator it =
             listSetMaxCompatChosen.begin();
         it != listSetMaxCompatChosen.end(); ++it) {
      set<int> setSitesRemain = it->second;

      if (setSitesRemain.size() == 0) {
        continue;
      }

      set<int> setMaxCompatChosen = it->first;

      // find the one that is the most compatible with remaining sites
      // int maxNumCompat = -1;
      double wtSiteMax = -1.0 * HAP_MAX_INT;

      vector<int> listSitesNext;
      for (set<int>::iterator it = setSitesRemain.begin();
           it != setSitesRemain.end(); ++it) {
        //    int numCompat = 0;
        //    for(set<int> :: iterator it2 = setSitesRemain.begin(); it2 !=
        //    setSitesRemain.end(); ++it2)
        //    {
        //        if( AreSitesCompatInMap(*pSetCompatPairsUse, *it,*it2) )
        //        {
        //            ++numCompat;
        //        }
        //    }
        double wtCur = wtSites[*it];
        //    if( numCompat > maxNumCompat )
        if (wtCur > wtSiteMax) {
          listSitesNext.clear();
          listSitesNext.push_back(*it);
          wtSiteMax = wtCur;
          // maxNumCompat = numCompat;
        }
        //    else if( numCompat == maxNumCompat )
        if (wtCur == wtSiteMax) {
          listSitesNext.push_back(*it);
        }
      }

      // if weight is too small now, stop if we have already get enough
      // if( wtSiteMax < 1.0)
      //{
      //    if( ((int)DEF_MIN_FRAC*GetNumSites()) <=
      //    (int)setMaxCompatChosen.size() )
      //    {
      //        break;
      //    }
      //}

      for (int jj = 0; jj < (int)listSitesNext.size(); ++jj) {
        // don't continue adding if we are at the limit
        if ((int)listSetMaxCompatChosenNext.size() > maxNumSets) {
          continue;
        }

        int sChose = listSitesNext[jj];
        set<int> setMaxCompatChosenNew = setMaxCompatChosen;
        setMaxCompatChosenNew.insert(sChose);

        // remove any sites that are incompatible with the chosen sites
        set<int> setSitesRemainNew;
        for (set<int>::iterator it = setSitesRemain.begin();
             it != setSitesRemain.end(); ++it) {
          if (AreSitesCompatInMap(*pSetCompatPairsUse, sChose, *it) == true) {
            setSitesRemainNew.insert(*it);
          }
        }
        setSitesRemainNew.erase(sChose);

        pair<set<int>, set<int> > pp(setMaxCompatChosenNew, setSitesRemainNew);

        listSetMaxCompatChosenNext.insert(pp);
      }
    }

    //
    if (listSetMaxCompatChosenNext.size() == 0) {
      //
      break;
    } else {
      listSetMaxCompatChosen = listSetMaxCompatChosenNext;
    }
  }

  YW_ASSERT_INFO(listSetMaxCompatChosen.size() > 0, "Cannot be empty");
  for (set<pair<set<int>, set<int> > >::iterator it =
           listSetMaxCompatChosen.begin();
       it != listSetMaxCompatChosen.end(); ++it) {
    // cout << "Maximum clique found by the heuristic: ";
    // DumpIntSet( it->first );
    map<int, set<int> > mm;
    for (set<int>::iterator it2 = it->first.begin(); it2 != it->first.end();
         ++it2) {
      set<int> ss;
      GetMutRowsHapAtSite(*it2, ss);
      mm[*it2] = ss;
    }
    listSetSitesCompat.push_back(mm);
  }

  //#endif

#if 0
    BinaryMatrix &matHaplotypesUse = const_cast<BinaryMatrix &>( this->matHaplotypes );

    //
    listSetSitesCompat.clear();
    vector<vector<bool> > vecHapsFullyCompat( GetNumSites() );
    for(int i=0; i<GetNumSites(); ++i)
    {
        vecHapsFullyCompat[i].resize( GetNumSites() );
    }

    for(int s1 = 0; s1<GetNumSites(); ++s1)
    {
        for(int s2=s1+1; s2<GetNumSites(); ++s2)
        {
            // root allele: 0
            bool fCompat = matHaplotypesUse.IsCompatibleRooted(s1, s2, 0, 0);
            vecHapsFullyCompat[s1][s2] = fCompat;
            vecHapsFullyCompat[s2][s1] = fCompat;
//cout << "Sites " << s1 << "," << s2 << ": ";
//if(fCompat)
//{
//cout << " compatible\n";
//}
//else
//{
//cout << " not compatible\n";
//}
        }
    }
    // find maximal compatible components
    set< set<int> > setMaximalComps;
    // start by putting all compatible pairs
    for(int s1 = 0; s1<GetNumSites(); ++s1)
    {
        set<int> ss;
        ss.insert(s1);
        setMaximalComps.insert(ss);
    }
    // find larger
    while(true)
    {
        // every time, make sure size is not too large
        TrimCliquesMaxDiff( setMaximalComps, maxNumSets );
//cout << "Size of current cliques to grow: " << setMaximalComps.size() << endl;
//for( set<set<int> > :: iterator it = setMaximalComps.begin(); it != setMaximalComps.end(); ++it)
//{
//DumpIntSet(*it);
//}

        set< set<int> > setMaximalCompsNext;
        // try to grow by adding one more
        for( set<set<int> > :: iterator it = setMaximalComps.begin(); it != setMaximalComps.end(); ++it )
        {
            for(int s=0; s<GetNumSites(); ++s)
            {
                if(  it->find(s) == it->end() )
                {
                    bool fCompat = true;
                    for(set<int> :: iterator it2 = it->begin(); it2 != it->end(); ++it2 )
                    {
                        if( vecHapsFullyCompat[ s ][ *it2 ]  == false )
                        {
                            fCompat = false;
                            break;
                        }
                    }
                    if( fCompat )
                    {
                        set<int> ss = *it;
                        ss.insert( s );
                        setMaximalCompsNext.insert(ss);
//cout << "Growing a subset: ";
//DumpIntSet(ss);
                    }
                }
            }
        }
        if( setMaximalCompsNext.size() == 0 )
        {
            break;
        }
        else
        {
            setMaximalComps = setMaximalCompsNext;
        }
    }
    //
    //TrimCliquesMaxDiff( setMaximalComps, maxNumSets );

    YW_ASSERT_INFO( setMaximalComps.size() > 0, "Cannot be empty" );
    for( set<set<int> > :: iterator it = setMaximalComps.begin(); it != setMaximalComps.end(); ++it )
    {
cout << "Clique found: ";
DumpIntSet(*it);
        map<int, std::set<int> >  setSitesCompat;

        set<int> ssChosen = *it;
        for(set<int> :: iterator it = ssChosen.begin(); it != ssChosen.end(); ++it)
        {
            set<int> ss;
            GetMutRowsHapAtSite(*it, ss);
            setSitesCompat[*it] = ss;
        }
        listSetSitesCompat.push_back(setSitesCompat);
    }
#endif
}

int ScistHaplotypeMat ::GetNumSites() const {
  return matHaplotypes.GetColNum();
}

int ScistHaplotypeMat ::GetNumHaps() const { return matHaplotypes.GetRowNum(); }

void ScistHaplotypeMat ::GetMutRowsHapAtSite(int site,
                                             std::set<int> &setRows) const {
  // any allele w/ non-zero is mutant
  setRows.clear();
  for (int r = 0; r < matHaplotypes.GetRowNum(); ++r) {
    if (matHaplotypes(r, site) == 1) {
      setRows.insert(r);
    }
  }
}

void ScistHaplotypeMat ::GetRowsWithGenoAtSite(int site, int geno,
                                               std::set<int> &setRows) const {
  setRows.clear();
  if (geno == 1) {
    GetMutRowsHapAtSite(site, setRows);
  } else if (geno == 0) {
    // get the complement
    setRows.clear();
    PopulateSetWithInterval(setRows, 0, GetNumHaps() - 1);
    set<int> setRows1;
    GetMutRowsHapAtSite(site, setRows1);
    SubtractSets(setRows, setRows1);
  }
}

double ScistHaplotypeMat ::GetScoreForGeno(int scIndex, int site,
                                           int genotype) const {
  int allele = this->matHaplotypes(scIndex, site);
  if (allele == genotype) {
    // when greeing, score is 0
    return 0.0;
  }

  // for now, only use default scoring
  double res = 0.0;
  double prob0 = this->matHaplotypesProb0[scIndex][site];
  double prob1 = 1.0 - prob0;
  if (genotype == 1) {
    // change from 0 to 1
    if (prob1 <= 0.0) {
      res = HAP_MAX_INT * 1.0;
    } else {
      res = log(prob0 / prob1);
    }
  } else {
    if (prob0 <= 0.0) {
      res = HAP_MAX_INT * 1.0;
    } else {
      res = log(prob1 / prob0);
    }
  }
  if (res < 0.0) {
    this->Dump();
    cout << "cell: " << scIndex << ", site: " << site
         << ", genotype: " << genotype << ", prob0: " << prob0 << endl;
  }
  YW_ASSERT_INFO(res >= 0.0, "Prob: wrong");
  return res;
}

bool ScistHaplotypeMat ::IsNoninformative(int site) const {
  //
  BinaryMatrix &matHaplotypesUse =
      const_cast<BinaryMatrix &>(this->matHaplotypes);
  return matHaplotypesUse.IsColNonInformative(site);
}

bool ScistHaplotypeMat ::IsCompatible(int s1, int s2) const {
  //
  BinaryMatrix &matHaplotypesUse =
      const_cast<BinaryMatrix &>(this->matHaplotypes);
  return matHaplotypesUse.IsCompatible(s1, s2);
}

std::string ScistHaplotypeMat ::ConsTree() const {
  //
  // construct phylogeny
  vector<int> rootZero;
  for (int i = 0; i < GetNumSites(); ++i) {
    rootZero.push_back(0);
  }
  PhylogenyTree phTree;
  phTree.SetRoot(rootZero);
  phTree.ConsOnBinMatrix(this->matHaplotypes);
  phTree.RemoveDegreeTwoNodes();

  // now assign leaf labels
  map<string, string> mapIdToLabels;
  for (int i = 0; i < GetNumHaps(); ++i) {
    // cout << "i: " << i << ", name: " << this->genosInput.GetGenotypeName(i)
    // << endl;
    string str = "(" + std::to_string(i) + ")";
    mapIdToLabels[str] = GetGenotypeName(i);
  }
  phTree.ReassignLeafLabels(mapIdToLabels);

  string res;
  phTree.ConsNewickSorted(res);
  // phTree.ConsNewick(res, false, 0.0, true);
  return res;
}

double ScistHaplotypeMat ::SumLogProbs() const {
  //
  double res = 0.0;
  for (int i = 0; i < (int)matHaplotypesProb0.size(); ++i) {
    res += GetSumOfVecElements(matHaplotypesProb0[i]);
  }
  return res;
}

void ScistHaplotypeMat ::Dump() const {
  ScistGenGenotypeMat ::Dump();

  //
  cout << "Matrix: [" << GetNumHaps() << "," << GetNumSites() << "]" << endl;
  this->matHaplotypes.Dump();
#if 0
    cout << "Clusters\n";
    for(int c=0; c<GetNumSites(); ++c)
    {
        cout << "Site " << c+1 << ": ";
        set<int> rowsMut;
        this->matHaplotypes.GetRowsWithAllele(c, 1, rowsMut);
        DumpIntSet(rowsMut);
    }
#endif
  cout << "Probabilities: \n";
  for (int i = 0; i < (int)matHaplotypesProb0.size(); ++i) {
    DumpDoubleVec(matHaplotypesProb0[i]);
  }
}

void ScistHaplotypeMat ::OutputImput(const string *pStrDesc) const {
  //
  cout << "Lineages: ";
  for (int i = 0; i < GetNumNames(); ++i) {
    cout << GetGenotypeName(i) << "  ";
  }
  cout << endl;
  if (pStrDesc != NULL) {
    cout << *pStrDesc << endl;
  } else {
    cout << "Imputed genotypes: \n";
  }
  for (int s = 0; s < GetNumSites(); ++s) {
    cout << "Site " << setw(6) << s + 1 << ":\t";

    for (int i = 0; i < GetNumHaps(); ++i) {
      cout << GetGenotypeAt(i, s) << " ";
    }
    cout << endl;
  }
}

bool ScistHaplotypeMat ::IsProbSignificant(double prob, double thresVal) const {
  //
  const double probConst = 0.5;
  if (prob < probConst && prob > (probConst - thresVal / 2)) {
    return false;
  }
  if (prob > probConst && prob < (probConst + thresVal / 2)) {
    return false;
  }
  return true;
}

// *************************************************************************************
// genotypes: ternary matrix

ScistTernaryMat ::ScistTernaryMat() {}

ScistGenGenotypeMat *ScistTernaryMat ::Copy() const {
  //
  ScistTernaryMat *pMatCopy = new ScistTernaryMat();

  for (int i = 0; i < GetNumNames(); ++i) {
    pMatCopy->AddGenotypeName(GetGenotypeName(i));
  }

  pMatCopy->SetSize(GetNumHaps(), GetNumSites());

  //
  for (int i = 0; i < GetNumHaps(); ++i) {
    for (int j = 0; j < GetNumSites(); ++j) {
      pMatCopy->SetGenotypeAt(i, j, GetGenotypeAt(i, j));
      pMatCopy->SetGenotypeProbOfGenoAt(i, j, 0, GetGenotypeProbAt(i, j, 0));
      pMatCopy->SetGenotypeProbOfGenoAt(i, j, 1, GetGenotypeProbAt(i, j, 1));
    }
  }

  return pMatCopy;
}

bool ScistTernaryMat ::ReadFromFile(std::ifstream &infile, int numSites,
                                    int numSCs, bool fSiteName) {
  //
  // assume each site is independent
  SetSize(numSCs, numSites);
  for (int i = 0; i < numSites; ++i) {
    string strName;
    if (fSiteName) {
      infile >> strName;
    } else {
      strName = std::to_string(i + 1);
    }
    AddSiteName(strName);

    // cout << "Read in site: " << i << endl;
    for (int j = 0; j < numSCs; ++j) {
      double prob0 = 0.0, prob1 = 0.0;
      bool res = ReadFromFileTernaryProb(infile, prob0, prob1);
      if (res == false) {
        return false;
      }

      SetGenotypeProbOfGenoAt(j, i, 0, prob0);
      SetGenotypeProbOfGenoAt(j, i, 1, prob1);

      // choose the allele w/ higher prob
      int allele = 0;
      double probMax = GetGenotypeProbAt(j, i, 0);
      if (probMax < GetGenotypeProbAt(j, i, 1)) {
        probMax = GetGenotypeProbAt(j, i, 1);
        allele = 1;
      }
      if (probMax < GetGenotypeProbAt(j, i, 2)) {
        probMax = GetGenotypeProbAt(j, i, 2);
        allele = 2;
      }
      SetGenotypeAt(j, i, allele);
    }
  }

  cout << "Input matrix: ";
  this->matTernary.Dump();

  return true;
}

bool ScistTernaryMat ::ReadFromFileTernaryProb(std::ifstream &infile,
                                               double &prob0, double &prob1) {
  // read in the prob of allele: (0.6,0.1) 0.6 means prob of 0 is 0.6 and prob
  // of 1 is 0.1 assume prob of 0 + 1 + 2 = 1
  infile >> prob0 >> prob1;
  return true;
}

void ScistTernaryMat ::SetSize(int numSCs, int numSites) {
  matTernary.SetSize(numSCs, numSites);

  matTernaryProbs.clear();
  matTernaryProbs.resize(numSCs);

  bool fNameInit = GetNumNames() > 0;

  for (int i = 0; i < numSCs; ++i) {
    matTernaryProbs[i].resize(numSites);
    for (int s = 0; s < numSites; ++s) {
      SetGenotypeProbOfGenoAt(i, s, 0, 1.0);
      SetGenotypeProbOfGenoAt(i, s, 1, 0.0);
    }

    // by default, use the numericals, starting from one
    if (fNameInit == false) {
      string str = std::to_string(i + 1);
      AddGenotypeName(str);
      // cout << "Init name: " << str << endl;
    }
  }
}

int ScistTernaryMat ::GetGenotypeAt(int sc, int site) const {
  return matTernary(sc, site);
}

int ScistTernaryMat ::GetAltGenotypeAt(int sc, int site) const {
  YW_ASSERT_INFO(false, "Not supported1");
  return 1;
}

void ScistTernaryMat ::SetGenotypeAt(int sc, int site, int geno) {
  matTernary(sc, site) = geno;
}

void ScistTernaryMat ::AddGenotypeAt(int sc, int site, int geno) {
  // append the genotype into it
  int genoThis = GetGenotypeAt(sc, site);
  if (genoThis != geno) {
    SetGenotypeAt(sc, site, geno);
  }
}

double ScistTernaryMat ::GetGenotypeProbAllele0At(int sc, int site) const {
  return GetGenotypeProbAt(sc, site, 0);
}

double ScistTernaryMat ::GetGenotypeProbAt(int sc, int site, int geno) const {
  if (geno == 0) {
    return this->matTernaryProbs[sc][site].first;
  } else if (geno == 1) {
    return this->matTernaryProbs[sc][site].second;
  } else {
    return 1.0 - GetGenotypeProbAt(sc, site, 0) -
           GetGenotypeProbAt(sc, site, 1);
  }
}

void ScistTernaryMat ::SetGenotypeProbAt(int sc, int site, double prob) {
  YW_ASSERT_INFO(false, "Not impelemented");
}

void ScistTernaryMat ::SetGenotypeProbOfGenoAt(int sc, int site, int geno,
                                               double prob) {
  if (geno == 0) {
    matTernaryProbs[sc][site].first = prob;
  } else if (geno == 1) {
    matTernaryProbs[sc][site].second = prob;
  } else {
    YW_ASSERT_INFO(false, "Cannot only set the homozygous mutant probility");
  }
}

void ScistTernaryMat ::FindMaximalCompatSites(
    const std::vector<double> &wtSites,
    std::vector<std::map<int, std::set<int> > > &listSetSitesCompat,
    int maxNumSets,
    const std::set<std::pair<int, int> > *pSetCompatPairs) const {
  YW_ASSERT_INFO(false, "Not implemented");
}

int ScistTernaryMat ::GetNumSites() const { return matTernary.GetColNum(); }

int ScistTernaryMat ::GetNumHaps() const { return matTernary.GetRowNum(); }

void ScistTernaryMat ::GetMutRowsHapAtSite(int site,
                                           std::set<int> &setRows) const {
  // YW_ASSERT_INFO(false, "Not supported2");
  // for now, use both 1/2 rows
  GetRowsWithGenoAtSite(site, 1, setRows);
  set<int> setRows2;
  GetRowsWithGenoAtSite(site, 2, setRows2);
  UnionSets(setRows, setRows2);
}

void ScistTernaryMat ::GetRowsWithGenoAtSite(int site, int geno,
                                             std::set<int> &setRows) const {
  setRows.clear();
  for (int h = 0; h < GetNumHaps(); ++h) {
    if (GetGenotypeAt(h, site) == geno) {
      setRows.insert(h);
    }
  }
}

double ScistTernaryMat ::GetScoreForGeno(int scIndex, int site,
                                         int genotype) const {
  YW_ASSERT_INFO(false, "Not supported3");
  return 0.0;
}

bool ScistTernaryMat ::IsNoninformative(int site) const {
  YW_ASSERT_INFO(false, "Not supported4");
  return false;
}

bool ScistTernaryMat ::IsCompatible(int s1, int s2) const {
  YW_ASSERT_INFO(false, "Not supported5");
  return false;
}

std::string ScistTernaryMat ::ConsTree() const {
  // construct phylogeny
  vector<int> rootZero;
  for (int i = 0; i < GetNumSites(); ++i) {
    rootZero.push_back(0);
  }

  // construct binary matrix for distance computation
  BinaryMatrix binMat;
  ConsHapMatForDistCalc(binMat);

  PhylogenyTree phTree;
  phTree.SetRoot(rootZero);
  phTree.ConsOnBinMatrix(binMat);
  phTree.RemoveDegreeTwoNodes();

  // now assign leaf labels
  map<string, string> mapIdToLabels;
  for (int i = 0; i < GetNumHaps(); ++i) {
    // cout << "i: " << i << ", name: " << this->genosInput.GetGenotypeName(i)
    // << endl;
    string str = "(" + std::to_string(i) + ")";
    mapIdToLabels[str] = GetGenotypeName(i);
  }
  phTree.ReassignLeafLabels(mapIdToLabels);

  string res;
  phTree.ConsNewickSorted(res);
  // phTree.ConsNewick(res, false, 0.0, true);
  return res;
}

double ScistTernaryMat ::SumLogProbs() const {
  YW_ASSERT_INFO(false, "Not impelemtned");
  return 0.0;
}

void ScistTernaryMat ::Dump() const {
  ScistGenGenotypeMat::Dump();

  //
  cout << "Matrix: [" << GetNumHaps() << "," << GetNumSites() << "]" << endl;
  this->matTernary.Dump();

  cout << "Probabilities: \n";
  for (int i = 0; i < (int)matTernaryProbs.size(); ++i) {
    for (int j = 0; j < (int)matTernaryProbs[i].size(); ++j) {
      cout << "(" << matTernaryProbs[i][j].first << ","
           << matTernaryProbs[i][j].second << ")  ";
    }
    cout << endl;
  }
}

void ScistTernaryMat ::OutputImput(const string *pStrDesc) const {
  //
  cout << "Lineages: ";
  for (int i = 0; i < GetNumNames(); ++i) {
    cout << GetGenotypeName(i) << "  ";
  }
  cout << endl;
  if (pStrDesc != NULL) {
    cout << *pStrDesc << endl;
  } else {
    cout << "Imputed genotypes: \n";
  }
  for (int s = 0; s < GetNumSites(); ++s) {
    cout << "Site " << setw(6) << s + 1 << ":\t";

    for (int i = 0; i < GetNumHaps(); ++i) {
      cout << GetGenotypeAt(i, s) << " ";
    }
    cout << endl;
  }
}

void ScistTernaryMat ::ConsHapMatForDistCalc(
    BinaryMatrix &matHaplotypes) const {
  matHaplotypes.SetSize(GetNumHaps(), 2 * GetNumSites());
  for (int r = 0; r < GetNumHaps(); ++r) {
    for (int s = 0; s < GetNumSites(); ++s) {
      int geno = GetGenotypeAt(r, s);
      int allele0 = 0, allele1 = 0;
      if (geno != 0) {
        allele0 = 1;
      }
      if (geno == 2) {
        allele1 = 1;
      }
      matHaplotypes(r, 2 * s) = allele0;
      matHaplotypes(r, 2 * s + 1) = allele1;
    }
  }
}

bool ScistTernaryMat ::IsProbSignificant(double prob, double thresVal) const {
  //
  const double probConst = 0.3333333;
  if (prob < probConst && prob > (probConst - thresVal / 2)) {
    return false;
  }
  if (prob > probConst && prob < (probConst + thresVal / 2)) {
    return false;
  }
  return true;
}
