//
//  ScistDoublet.cpp
//
//
//  Created by Yufeng Wu on 6/2/18.
//
//

#include "ScistDoublet.hpp"
#include "PhylogenyTree.h"
#include "PhylogenyTreeBasic.h"
#include "ScistGenotype.hpp"
#include "ScistPerfPhyImp.hpp"
#include "Utils3.h"
#include <iomanip>

// *************************************************************************************
// DP backtrace info

ScistDoubletDPTraceback ::ScistDoubletDPTraceback()
    : indexChild1(-1), phaseChild1(-1), indexChild2(-1), phaseChild2(-1) {}

ScistDoubletDPTraceback ::ScistDoubletDPTraceback(
    const ScistDoubletDPTraceback &rhs)
    : indexChild1(rhs.indexChild1), phaseChild1(rhs.phaseChild1),
      indexChild2(rhs.indexChild2), phaseChild2(rhs.phaseChild2) {}

ScistDoubletDPTraceback &
ScistDoubletDPTraceback ::operator=(const ScistDoubletDPTraceback &rhs) {
  indexChild1 = rhs.indexChild1;
  phaseChild1 = rhs.phaseChild1;
  indexChild2 = rhs.indexChild2;
  phaseChild2 = rhs.phaseChild2;
  return *this;
}

// *************************************************************************************
// Deal with doublet

ScistDoublet ::ScistDoublet(const ScistGenGenotypeMat &genosInputIn)
    : genosInput(genosInputIn) {}

double ScistDoublet ::EvalGenoDoublet(const set<int> &setTemplateRows,
                                      int genoDoublet,
                                      vector<int> &genoDoublePhase1,
                                      vector<int> &genoDoublePhase2) const {
  // construct cluster trees
  map<int, ScistPerfPhyCluster> setTemplateSites;
  std::map<const ScistPerfPhyCluster *, int> mapClusToSiteIndex;
  ConsClustersForTemplates(setTemplateRows, setTemplateSites,
                           mapClusToSiteIndex);

  ScistPerfPhyClusTreeNode *pClusTreeRoot =
      ScistPerfPhyClusTreeNode::ConsClusterTree(setTemplateSites);

  // construct solution based on this
  std::map<ScistPerfPhyClusTreeNode *,
           std::vector<std::pair<double, ScistDoubletDPTraceback> > >
      mapNodeVals;
  ConsDPTblDoubletNodes(setTemplateSites, mapClusToSiteIndex, genoDoublet,
                        pClusTreeRoot, mapNodeVals);

  //
  double minCost = mapNodeVals[pClusTreeRoot][3].first;
  // cout << "The min-cost phasing has optimal cost: " << minCost << endl;

  vector<int> vecPhasing;
  ConsPhasing(mapClusToSiteIndex, genoDoublet, pClusTreeRoot, mapNodeVals,
              vecPhasing);
  // cout << "Phasing vector: ";
  // DumpIntVec( vecPhasing);

  // now construct phasing
  ConsPhasingVec(vecPhasing, genoDoublePhase1, genoDoublePhase2);

  delete pClusTreeRoot;
  return minCost;
}

void ScistDoublet ::ConsClustersForTemplates(
    const set<int> &setTemplateRows,
    std::map<int, ScistPerfPhyCluster> &setTemplateSites,
    std::map<const ScistPerfPhyCluster *, int> &mapClusToSiteIndex) const {
  // only use those rows
  setTemplateSites.clear();

  for (int s = 0; s < genosInput.GetNumSites(); ++s) {
    set<int> rowsMut;
    genosInput.GetMutRowsHapAtSite(s, rowsMut);
    set<int> rowsMutInTemp;
    JoinSets(rowsMut, setTemplateRows, rowsMutInTemp);

    // ignore any singleton
    if (rowsMutInTemp.size() == 0) {
      continue;
    }

    ScistPerfPhyCluster clus(rowsMutInTemp);
    setTemplateSites[s] = rowsMutInTemp;
  }

  // construct reverse mapping
  for (map<int, ScistPerfPhyCluster>::iterator it = setTemplateSites.begin();
       it != setTemplateSites.end(); ++it) {
    mapClusToSiteIndex[&(it->second)] = it->first;
  }

  // cout << "ConsClustersForTemplates: template rows\n";
  // for( map<int, ScistPerfPhyCluster > :: iterator it =
  // setTemplateSites.begin(); it != setTemplateSites.end(); ++it )
  //{
  // cout << "Site " << it->first << ": mut rows within template: ";
  // it->second.Dump();
  //}
}

void ScistDoublet ::ConsDPTblDoubletNodes(
    const std::map<int, ScistPerfPhyCluster> &setTemplateSites,
    const std::map<const ScistPerfPhyCluster *, int> &mapClusToSiteIndex,
    int genoDoublet, ScistPerfPhyClusTreeNode *pNodeCurr,
    std::map<ScistPerfPhyClusTreeNode *,
             vector<pair<double, ScistDoubletDPTraceback> > > &mapNodeVals)
    const {
  // cons DP table for doublet recursively from bottom up
  //
  const ScistPerfPhyCluster *pClus = pNodeCurr->GetClus();

#if 0
    // work with all sites
    vector<pair<double, ScistDoubletDPTraceback> > vecThis(4);
    vecThis[0].first = 0.0;
    vecThis[1].first = 0.0;
    vecThis[2].first = 0.0;
    vecThis[3].first = 0.0;
    if( pNodeCurr->IsLeaf() )
    {
        YW_ASSERT_INFO( pClus != NULL, "Leaf: cluster cannot be null" );
        for(int s=0; s<this->genosInput.GetNumSites(); ++s)
        {
            double prob0Orig = this->genosInput.GetGenotypeProbAllele0At(genoDoublet, s);
            double prob0 = -1.0*log(prob0Orig);
            double prob1 = -1.0*log(1.0-prob0Orig);

            vecThis[0].first += prob0;
            vecThis[1].first += prob1;
            vecThis[2].first += prob1;
            vecThis[3].first += prob1;
        }
    }
    mapNodeVals[ pNodeCurr ] = vecThis;
#endif
  //#if 0
  if (pClus != NULL) {
    map<const ScistPerfPhyCluster *, int>::const_iterator it =
        mapClusToSiteIndex.find(pClus);
    YW_ASSERT_INFO(it != mapClusToSiteIndex.end(), "Fail to find the cluster2");
    int site = it->second;

    //
    // double prob0 = this->genosInput.GetScoreForGeno( genoDoublet, site, 0 );
    // double prob1 = this->genosInput.GetScoreForGeno( genoDoublet, site, 1 );
    double prob0Orig =
        this->genosInput.GetGenotypeProbAllele0At(genoDoublet, site);
    double prob0 = -1.0 * log(prob0Orig);
    double prob1 = -1.0 * log(1.0 - prob0Orig);
    vector<pair<double, ScistDoubletDPTraceback> > vecThis(4);
    vecThis[0].first = prob0;
    vecThis[1].first = prob1;
    vecThis[2].first = prob1;
    vecThis[3].first = prob1;
    mapNodeVals[pNodeCurr] = vecThis;
  } else {
    // otherwise everything is zero
    vector<pair<double, ScistDoubletDPTraceback> > vecThis(4);
    vecThis[0].first = 0.0;
    vecThis[1].first = 0.0;
    vecThis[2].first = 0.0;
    vecThis[3].first = 0.0;
    mapNodeVals[pNodeCurr] = vecThis;
  }
  //#endif

  if (pNodeCurr->IsLeaf()) {
    return;
  }

  // internal node: first construct all the descendents
  for (int i = 0; i < pNodeCurr->GetNumChildren(); ++i) {
    ScistPerfPhyClusTreeNode *pChild = pNodeCurr->GetChild(i);
    ConsDPTblDoubletNodes(setTemplateSites, mapClusToSiteIndex, genoDoublet,
                          pChild, mapNodeVals);
  }

  // now setup the values for the current node
  vector<pair<double, ScistDoubletDPTraceback> > vec;

  // phasing 00
  pair<double, ScistDoubletDPTraceback> mv00;
  mv00.first = 0.0;
  for (int i = 0; i < pNodeCurr->GetNumChildren(); ++i) {
    ScistPerfPhyClusTreeNode *pChild = pNodeCurr->GetChild(i);
    mv00.first += mapNodeVals[pChild][0].first;
  }
  // use default traceback
  vec.push_back(mv00);

  // phasing 01
  pair<double, ScistDoubletDPTraceback> mv01;
  mv01.first = mv00.first;
  for (int i = 0; i < pNodeCurr->GetNumChildren(); ++i) {
    ScistPerfPhyClusTreeNode *pChild = pNodeCurr->GetChild(i);

    //
    double mv01i = mv00.first - mapNodeVals[pChild][0].first +
                   mapNodeVals[pChild][1].first;
    if (mv01i < mv01.first) {
      mv01.first = mv01i;
      mv01.second.SetChild1(i);
      mv01.second.SetPhase1(1);
    }
  }
  vec.push_back(mv01);

  // phasing 10
  pair<double, ScistDoubletDPTraceback> mv10;
  mv10.first = mv00.first;
  for (int i = 0; i < pNodeCurr->GetNumChildren(); ++i) {
    ScistPerfPhyClusTreeNode *pChild = pNodeCurr->GetChild(i);

    //
    double mv10i = mv00.first - mapNodeVals[pChild][0].first +
                   mapNodeVals[pChild][2].first;
    if (mv10i < mv10.first) {
      mv10.first = mv10i;
      mv10.second.SetChild1(i);
      mv10.second.SetPhase1(2);
    }
  }
  vec.push_back(mv10);

  // phasing 11
  pair<double, ScistDoubletDPTraceback> mv11;
  mv11.first = std::min(mv01.first, mv10.first);
  // setup trace back
  if (mv11.first == mv01.first) {
    mv11.second = mv01.second;
  } else {
    mv11.second = mv10.second;
  }

  // consider exatly one is 11
  for (int i = 0; i < pNodeCurr->GetNumChildren(); ++i) {
    ScistPerfPhyClusTreeNode *pChild = pNodeCurr->GetChild(i);

    //
    double mv11i = mv00.first - mapNodeVals[pChild][0].first +
                   mapNodeVals[pChild][3].first;
    if (mv11i < mv11.first) {
      mv11.first = mv11i;
      mv11.second.SetChild1(i);
      mv11.second.SetPhase1(3);
    }
  }
  // consider a pair of i and j
  for (int i = 0; i < pNodeCurr->GetNumChildren(); ++i) {
    ScistPerfPhyClusTreeNode *pChildi = pNodeCurr->GetChild(i);

    for (int j = 0; j < pNodeCurr->GetNumChildren(); ++j) {
      if (i == j) {
        continue;
      }

      ScistPerfPhyClusTreeNode *pChildj = pNodeCurr->GetChild(j);

      //
      double mv11i = mv00.first - mapNodeVals[pChildi][0].first -
                     mapNodeVals[pChildj][0].first +
                     mapNodeVals[pChildi][1].first +
                     mapNodeVals[pChildj][2].first;
      if (mv11i < mv11.first) {
        mv11.first = mv11i;
        mv11.second.SetChild1(i);
        mv11.second.SetPhase1(1);
        mv11.second.SetChild2(j);
        mv11.second.SetPhase2(2);
      }
    }
  }
  vec.push_back(mv11);

  // add the current cost
  for (int i = 0; i < (int)vec.size(); ++i) {
    vec[i].first += mapNodeVals[pNodeCurr][i].first;
  }
  mapNodeVals[pNodeCurr] = vec;
}

void ScistDoublet ::ConsPhasing(
    const std::map<const ScistPerfPhyCluster *, int> &mapClusToSiteIndex,
    int genoDoublet, ScistPerfPhyClusTreeNode *pNodeRoot,
    const std::map<ScistPerfPhyClusTreeNode *,
                   std::vector<std::pair<double, ScistDoubletDPTraceback> > >
        &mapNodeVals,
    vector<int> &vecPhasing) const {
  //
  vecPhasing.resize(this->genosInput.GetNumSites());

  // init all phasing to be 00 for genotype 0 and 01 for genotype 1
  for (int i = 0; i < this->genosInput.GetNumSites(); ++i) {
    int geno = this->genosInput.GetGenotypeAt(genoDoublet, i);
    if (geno == 0) {
      vecPhasing[i] = 0;
    } else {
      vecPhasing[i] = 1;
    }
  }
  const int ROOT_PHASING = 3;
  TracePhasingAtNode(mapClusToSiteIndex, genoDoublet, pNodeRoot, ROOT_PHASING,
                     mapNodeVals, vecPhasing);
}

void ScistDoublet ::TracePhasingAtNode(
    const std::map<const ScistPerfPhyCluster *, int> &mapClusToSiteIndex,
    int genoDoublet, ScistPerfPhyClusTreeNode *pNodeCurr, int phasingCurr,
    const std::map<ScistPerfPhyClusTreeNode *,
                   std::vector<std::pair<double, ScistDoubletDPTraceback> > >
        &mapNodeVals,
    vector<int> &vecPhasing) const {
  //
  const ScistPerfPhyCluster *pClus = pNodeCurr->GetClus();
  if (pClus != NULL) {
    map<const ScistPerfPhyCluster *, int>::const_iterator it =
        mapClusToSiteIndex.find(pClus);
    YW_ASSERT_INFO(it != mapClusToSiteIndex.end(), "Fail to find the cluster2");
    int site = it->second;

    // record this phasing
    vecPhasing[site] = phasingCurr;
  }

  // consider all children
  std::map<ScistPerfPhyClusTreeNode *,
           std::vector<std::pair<double, ScistDoubletDPTraceback> > >::
      const_iterator it = mapNodeVals.find(pNodeCurr);
  YW_ASSERT_INFO(it != mapNodeVals.end(), "Fail to find");
  for (int i = 0; i < pNodeCurr->GetNumChildren(); ++i) {
    ScistPerfPhyClusTreeNode *pChild = pNodeCurr->GetChild(i);
    int phasingChild = 0;
    if (it->second[phasingCurr].second.GetChild1() == i) {
      phasingChild = it->second[phasingCurr].second.GetPhase1();
    } else if (it->second[phasingCurr].second.GetChild2() == i) {
      phasingChild = it->second[phasingCurr].second.GetPhase2();
    }
    TracePhasingAtNode(mapClusToSiteIndex, genoDoublet, pChild, phasingChild,
                       mapNodeVals, vecPhasing);
  }
}

void ScistDoublet ::ConsPhasingVec(const std::vector<int> &vecPhasing,
                                   std::vector<int> &genoDoublePhase1,
                                   std::vector<int> &genoDoublePhase2) const {
  //
  genoDoublePhase1.clear();
  genoDoublePhase2.clear();
  for (int i = 0; i < (int)vecPhasing.size(); ++i) {
    int p = vecPhasing[i];
    int a1, a2;
    if (p == 0) {
      a1 = 0;
      a2 = 0;
    } else if (p == 1) {
      a1 = 0;
      a2 = 1;
    } else if (p == 2) {
      a1 = 1;
      a2 = 0;
    } else {
      a1 = 1;
      a2 = 1;
    }
    genoDoublePhase1.push_back(a1);
    genoDoublePhase2.push_back(a2);
  }
}

// *************************************************************************************
// Deal with doublet (search)

const double DEF_DOUBLET_COST = 0.0;

ScistDoubletSearch ::ScistDoubletSearch(const ScistGenGenotypeMat &genosInputIn,
                                        int maxDoubletSubsetSzIn)
    : genosInput(genosInputIn), maxDoubletSubsetSz(maxDoubletSubsetSzIn),
      costDoublet(DEF_DOUBLET_COST), fVerbose(false),
      fOutputPPWithEdgeLabels(false) {}

void ScistDoubletSearch ::Search() {
  // cout << "Matrix: ";
  // this->genosInput.Dump();
  set<int> setCandidates;
  FindDoubletCandidates(setCandidates);
  // cout << "Candidates: ";
  // DumpIntSet(setCandidates);
  int szDoublets = this->maxDoubletSubsetSz;
  if (szDoublets > (int)setCandidates.size()) {
    szDoublets = (int)setCandidates.size();
  }
  YW_ASSERT_INFO(szDoublets > 0, "Wrong: no doublets to work with. Consider "
                                 "run without specifying doublets");

  // try all subset up to a level
  double opt = HAP_MAX_INT * 1.0;
  ScistGenGenotypeMat *pMatRes = NULL;
  for (int szDoubletsStep = 0; szDoubletsStep <= szDoublets; ++szDoubletsStep) {
    vector<int> posvec;
    GetFirstCombo(szDoubletsStep, (int)setCandidates.size(), posvec);
    while (true) {
      // now work with the chosen subset
      set<int> rowsDoubles;
      PopulateSetByVec(rowsDoubles, posvec);
      // cout << "Processing doublets: ";
      // DumpIntSet(rowsDoubles);
      //
      double optStep = 0.0;
      ScistGenGenotypeMat *pMatStep =
          EvalGenoDoubletSet(this->genosInput, rowsDoubles, optStep);
      YW_ASSERT_INFO(pMatStep != NULL, "Canot be null");
      // cout << "optStep: " << optStep << endl;
      if (optStep < opt) {
        // cout << "BETTER\n";
        opt = optStep;
        if (pMatRes != NULL) {
          delete pMatRes;
        }
        pMatRes = pMatStep;
      } else {
        delete pMatStep;
      }

      if (GetNextCombo(szDoubletsStep, (int)setCandidates.size(), posvec) ==
          false) {
        break;
      }
    }
  }
  YW_ASSERT_INFO(pMatRes != NULL, "Resulting matrix: not found");
  cout << "**** Optimal cost for doublet resoultion: " << opt << endl;
  if (fVerbose) {
    pMatRes->OutputImput();
  }
  string strTree = pMatRes->ConsTree();
  cout << "Constructed single cell phylogeny: " << strTree << endl;

  if (this->fVerbose) {
    // keep track of imputation results
    ScistGenGenotypeMat *pMatImpute = genosInput.Copy();
    std::map<int, std::set<int> > mapDoublets;
    FindOrigImputedGeno(*pMatRes, *pMatImpute, mapDoublets);

    //
    cout << "Doublet genotypes (1-based)): <original dobule genotype> : <list "
            "of expanded doublet rows in imputed matrix>\n";
    for (map<int, set<int> >::iterator it = mapDoublets.begin();
         it != mapDoublets.end(); ++it) {
      cout << it->first << " : ";
      for (set<int>::const_iterator it2 = it->second.begin();
           it2 != it->second.end(); ++it2) {
        cout << *it2 + 1 << " ";
      }
      cout << endl;
    }

    // also output the imputaton results
    cout << "Imputed genotypes: \n";
    pMatImpute->OutputImput();

    set<pair<pair<int, int>, int> > listChangedPlaces;
    for (int i = 0; i < genosInput.GetNumHaps(); ++i) {
      for (int j = 0; j < genosInput.GetNumSites(); ++j) {
        if (genosInput.GetGenotypeAt(i, j) != pMatImpute->GetGenotypeAt(i, j)) {
          pair<int, int> pp(i, j);
          pair<pair<int, int>, int> pp1(pp, pMatImpute->GetGenotypeAt(i, j));
          listChangedPlaces.insert(pp1);
        }
      }
    }
    cout << "List of corrected genotypes (site, cell, new genotype) in base-1: "
            "\n";
    for (set<pair<pair<int, int>, int> >::iterator it =
             listChangedPlaces.begin();
         it != listChangedPlaces.end(); ++it) {
      cout << "[ " << setw(6) << it->first.second + 1 << " " << setw(6)
           << it->first.first + 1 << " ]: " << it->second << endl;
    }

    delete pMatImpute;
  }

  delete pMatRes;
}

void ScistDoubletSearch ::SearchInc() {
  // search incrementally for doublets
  ScistGenGenotypeMat *pMatRes = this->genosInput.Copy();
  double optFinal = 1.0 * HAP_MAX_INT;
  bool fInit = false;

  int numDoublesUsed = 0;
  while (numDoublesUsed < this->maxDoubletSubsetSz) {
    double opt = HAP_MAX_INT * 1.0;
    set<int> rowsDoublesEmpty;
    ScistGenGenotypeMat *pMatInitDump =
        EvalGenoDoubletSet(*pMatRes, rowsDoublesEmpty, opt);
    YW_ASSERT_INFO(pMatInitDump != NULL, "Cannot be null");
    // cout << "pMatInitDump: ";
    // pMatInitDump->Dump();
    // ScistHaplotypeMat *pMatResHap0 = dynamic_cast<ScistHaplotypeMat
    // *>(pMatInitDump); string strTreeEdgeLabel0 =
    // ConsRootedPerfectPhylogenyFromMat(pMatResHap0->GetHapMat(), true, true);
    // cout << "Stepwise tree: " << strTreeEdgeLabel0 << endl;
    delete pMatInitDump;

    if (fInit == false) {
      fInit = true;
      optFinal = opt;
    }

    // cout << "Finding doublet: opt=" << opt << ", num of doublet so far: " <<
    // numDoublesUsed+1 << ", current matrix: "; pMatRes->Dump();
    // try to find the best single doublet row to expand
    double optLoop = HAP_MAX_INT * 1.0;
    ScistGenGenotypeMat *pMatLoop = NULL;
    int indexDouble = -1;
    for (int i = 0; i < pMatRes->GetNumHaps(); ++i) {
      // cout << "i = " << i << endl;
      // now work with the chosen subset
      set<int> rowsDoubles;
      rowsDoubles.insert(i);
      //
      double optStep = 0.0;
      ScistGenGenotypeMat *pMatStep =
          EvalGenoDoubletSet(*pMatRes, rowsDoubles, optStep);
      if (pMatStep != NULL) {
        // cout << "Stepwise matrix: ";
        // pMatStep->Dump();
        // ScistHaplotypeMat *pMatResHap = dynamic_cast<ScistHaplotypeMat
        // *>(pMatStep); string strTreeEdgeLabel1 =
        // ConsRootedPerfectPhylogenyFromMat(pMatResHap->GetHapMat(), true,
        // true); cout << "Stepwise tree: " << strTreeEdgeLabel1 << endl; cout
        // << "for genotype: " << i << ", optStep: " << optStep << endl;
        if (optStep < optLoop) {
          // cout << "BETTER\n";
          optLoop = optStep;
          if (pMatLoop != NULL) {
            delete pMatLoop;
          }
          pMatLoop = pMatStep;
          indexDouble = i;
        } else {
          delete pMatStep;
        }
      }
    }
    if (indexDouble < 0) {
      break;
    }
    if (optLoop >= opt) {
      // YW: 08/22/18, now force to have the same number of doublets
      // break;
    }
    if (pMatLoop == NULL) {
      break;
    }

    opt = optLoop;
    optFinal = optLoop;
    YW_ASSERT_INFO(pMatLoop != NULL, "Cannot be null");
    YW_ASSERT_INFO(indexDouble >= 0, "Wrong");
    // cout << "pMatLoop: ";
    // pMatLoop->Dump();
    ScistGenGenotypeMat *pMatLoopConv =
        CreateGnoesWithDouble(*pMatRes, indexDouble, *pMatLoop);
    // cout << "Converted matrix: ";
    // pMatLoopConv->Dump();

    delete pMatLoop;

    if (IsOverImpute(*pMatLoopConv) == true) {
      delete pMatLoopConv;
      break;
    }

    if (pMatRes != NULL) {
      delete pMatRes;
    }
    pMatRes = pMatLoopConv;

    ++numDoublesUsed;
  }

  YW_ASSERT_INFO(pMatRes != NULL, "Resulting matrix: not found");
  cout << "**** Optimal cost for doublet resoultion: " << optFinal << endl;
  if (fVerbose) {
    pMatRes->OutputImput();

    // analyze doublets
    int numDoublets = 0;
    for (int h = 0; h < pMatRes->GetNumHaps(); ++h) {
      string strName = pMatRes->GetGenotypeName(h);
      string strLastChar = strName.substr(strName.length() - 1, 1);
      if (strLastChar == "'") {
        //
        string strNameOrig = GetGenoDoubleRowName(strName);
        cout << "Doublet: imputed haplotype " << h + 1
             << " (with assigned name " << strName
             << ") is a doublet from cell " << strNameOrig << endl;
        ++numDoublets;
      }
    }
    cout << "Number of found doublets: " << numDoublets << endl;
  }
  if (fOutputPPWithEdgeLabels) {
    // cout << "Imputed genotypes: ";
    // pMatRes->Dump();
    OutputMutTree(*pMatRes);

#if 0
        ScistHaplotypeMat *pMatResHap = dynamic_cast<ScistHaplotypeMat *>(pMatRes);
        if(pMatResHap == NULL)
        {
            cout << "** Right now, only output perfect phylogeny for binary genotypes\n";
        }
        else
        {
            string strTreeEdgeLabel = ConsRootedPerfectPhylogenyFromMat(pMatResHap->GetHapMat(), true, true);
            //cout << "** Perfect phylogeny (with sites labeled on edges) from the imputed genotypes: " << strTreeEdgeLabel << endl;
            string strMutTree = ConsEdgeLabeTree(strTreeEdgeLabel);
            string strMutTreeConv = ConvMutTreeStr(strMutTree);
            cout << "^^ Mutation tree: " << strMutTreeConv << endl;
            // output mutation tree file
            OutputMutationTree( this->strMutTreeFileName.c_str(), strMutTreeConv, true );
        }
#endif
  }

  // string strTree = pMatRes->ConsTree();
  // cout << "Constructed single cell phylogeny: " << strTree << endl;
  string strNW;
  double likeliOpt = ConsTree(*pMatRes, strNW);
  // cout << "Optimal log-likelihood is " << likeliOpt << endl;
  cout << "**** Maximum log-likelihood: " << likeliOpt << endl;
  cout << "Constructed single cell phylogeny: " << strNW << endl;

#if 0
    if( this->fVerbose )
    {
        // keep track of imputation results
        ScistGenGenotypeMat *pMatImpute = genosInput.Copy();
        std::map<int,std::set<int> > mapDoublets;
        FindOrigImputedGeno( *pMatRes, *pMatImpute, mapDoublets );

        //
        cout << "Doublet genotypes (1-based)): <original dobule genotype> : <list of expanded doublet rows in imputed matrix>\n";
        for(map<int,set<int> > :: iterator it = mapDoublets.begin(); it != mapDoublets.end(); ++it)
        {
            cout << it->first << " : ";
            for(set<int> :: const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
            {
                cout << *it2 + 1 << " ";
            }
            cout << endl;
        }

        // also output the imputaton results
        cout << "Imputed genotypes: \n";
        pMatImpute->OutputImput();

        set< pair<pair<int,int>, int> > listChangedPlaces;
        for(int i=0; i<genosInput.GetNumHaps(); ++i)
        {
            for(int j=0; j<genosInput.GetNumSites(); ++j)
            {
                if( genosInput.GetGenotypeAt(i,j) != pMatImpute->GetGenotypeAt(i,j) )
                {
                    pair<int,int> pp(i,j);
                    pair<pair<int,int>,int> pp1(pp, pMatImpute->GetGenotypeAt(i,j));
                    listChangedPlaces.insert(pp1);
                }
            }
        }
        cout << "List of corrected genotypes (site, cell, new genotype) in base-1: \n";
        for(set<pair<pair<int,int>,int> > :: iterator it = listChangedPlaces.begin(); it != listChangedPlaces.end(); ++it )
        {
            cout << "[ " << setw(6) << it->first.second+1 << " " << setw(6) << it->first.first+1 << " ]: " << it->second << endl;
        }

        delete pMatImpute;
    }
#endif

  delete pMatRes;
}

double ScistDoubletSearch ::ConsTree(ScistGenGenotypeMat &genosNoDoublets,
                                     std::string &strNW) const {
  //
  ScistPerfPhyMLE sciInf1(genosNoDoublets);
  sciInf1.SetOutput(false);
  sciInf1.SetVerbose(false);
  std::set<std::pair<std::pair<int, int>, int> > listChangedPlaces;
  std::string strTreeNW;
  double opt = sciInf1.Infer(&listChangedPlaces, &strTreeNW);
  // cout << "Before mapping: inferred tree is " << strTreeNW << endl;
  // now remap
  map<string, string> mapIdToOrig;
  for (int h = 0; h < genosNoDoublets.GetNumHaps(); ++h) {
    string idCur = std::to_string(h + 1);
    string idMapped = genosNoDoublets.GetGenotypeName(h);
    mapIdToOrig[idCur] = idMapped;
    // cout << idCur << " mapped to " << idMapped << endl;
  }
  strNW = strTreeNW;
  NewickUtils::UpdateLabells(strNW, mapIdToOrig);
  // cout << "After mapping, inferred tree is: " << strNW << endl;
  return opt;
}

static string GetNonDoubleName(const string &strTaxon) {
  int posLast = (int)strTaxon.length() - 1;
  while (posLast >= 0) {
    string str = strTaxon.substr(posLast, 1);
    if (str == "'") {
      break;
    }
    --posLast;
  }
  //
  YW_ASSERT_INFO(posLast >= 0, "Fail111");
  return strTaxon.substr(0, posLast + 1);
}

bool ScistDoubletSearch ::IsOverImpute(
    const ScistGenGenotypeMat &genosDbl) const {
  // simple rule: if it use the same row again, then it overimputes
  for (int h = 0; h < genosDbl.GetNumHaps(); ++h) {
    string strName = genosDbl.GetGenotypeName(h);
    string strLastChar = strName.substr(strName.length() - 1, 1);
    string str2ndLastChar;
    if (strName.length() >= 2) {
      str2ndLastChar = strName.substr(strName.length() - 2, 1);
    }
    if (strLastChar == "'" && str2ndLastChar == "'") {
      //
      return true;
    }
  }
  return false;
}

void ScistDoubletSearch ::FindDoubletHapsInMat(
    const ScistGenGenotypeMat &genosDbl, std::set<int> &setHapsDoubles) const {
  //
  setHapsDoubles.clear();
  set<string> setDoubles;
  for (int h = 0; h < genosDbl.GetNumHaps(); ++h) {
    string strName = genosDbl.GetGenotypeName(h);
    string strLastChar = strName.substr(strName.length() - 1, 1);
    if (strLastChar == "'") {
      //
      string strNameOrig = GetGenoDoubleRowName(strName);
      setDoubles.insert(strNameOrig);
      setHapsDoubles.insert(h);
    }
  }
  for (int h = 0; h < genosDbl.GetNumHaps(); ++h) {
    string strName = genosDbl.GetGenotypeName(h);
    if (setDoubles.find(strName) != setDoubles.end()) {
      //
      setHapsDoubles.insert(h);
    }
  }
}

void ScistDoubletSearch ::OutputMutTree(
    ScistGenGenotypeMat &genosNoDoublets) const {
  // output the matrix
  ScistGenGenotypeMat *pMatRes = genosNoDoublets.Copy();

  // YW: 05/16/19 try to make tree inference with doublet more accurate
  // set all doublets haplotypes to be uncertain
  // analyze doublets
  set<int> setHapsDoubles;
  FindDoubletHapsInMat(*pMatRes, setHapsDoubles);
  // cout << "Set of doublet haplotypes: ";
  // DumpIntSet(setHapsDoubles);

  // now set uncertain haps to those positions
  for (set<int>::iterator it = setHapsDoubles.begin();
       it != setHapsDoubles.end(); ++it) {
    for (int s = 0; s < pMatRes->GetNumSites(); ++s) {
      double probOld = pMatRes->GetGenotypeProbAllele0At(*it, s);
      if (probOld < 0.5) {
        pMatRes->SetGenotypeProbAt(*it, s, probOld / 2 + 0.25);
      }

      // pMatRes->SetGenotypeProbAt(*it, s, 0.5);
      //}
      // else
      //{
      //    pMatRes->SetGenotypeProbAt(*it, s, 0.7);
      //}
    }
  }
  // cout << "After revision, genotype matrix: ";
  // pMatRes->Dump();

  //
  ScistPerfPhyMLE sciInf1(*pMatRes);
  sciInf1.SetOutput(false);
  sciInf1.SetVerbose(false);
  std::set<std::pair<std::pair<int, int>, int> > listChangedPlaces;
  std::string strTreeNW;
  // double opt =
  sciInf1.Infer(&listChangedPlaces, &strTreeNW);
  // cout << "Before mapping: inferred tree is " << strTreeNW << endl;

  pMatRes->ChangeGenosAtPositions(listChangedPlaces);
  // if( fVerbose )
  //{
  //    cout << "Called genotypes\n";
  //    pMatRes->OutputImput();
  //}
  ScistHaplotypeMat *pMatResHap = dynamic_cast<ScistHaplotypeMat *>(pMatRes);
  if (pMatResHap == NULL) {
    cout
        << "** Right now, only output perfect phylogeny for binary genotypes\n";
  } else {
    string strTreeEdgeLabel =
        ConsRootedPerfectPhylogenyFromMat(pMatResHap->GetHapMat(), true, true);
    // cout << "** Perfect phylogeny (with sites labeled on edges) from the
    // imputed genotypes: " << strTreeEdgeLabel << endl;

    string strMutTree = ConsEdgeLabeTree(strTreeEdgeLabel);
    string strMutTreeConv = ConvMutTreeStr(strMutTree);
    cout << "^^ Mutation tree: " << strMutTreeConv << endl;

    // output mutation tree file
    OutputMutationTree(this->strMutTreeFileName.c_str(), strMutTreeConv, true);
  }

  delete pMatRes;
}

ScistGenGenotypeMat *ScistDoubletSearch ::CreateGnoesWithDouble(
    const ScistGenGenotypeMat &genosOrig, int indexDouble,
    const ScistGenGenotypeMat &genosDoubleInfer) const {
  // cout << "CreateGnoesWithDouble: genosOrig: ";
  // genosOrig.Dump();
  // cout << "indexDouble: " << indexDouble << endl;
  // cout << "genosDoubleInfer: ";
  // genosDoubleInfer.Dump();

  // create a new genotype matrix w/ doublets
  ScistGenGenotypeMat *pResMat = genosOrig.CreateNewMat();
  pResMat->SetSize(genosOrig.GetNumHaps() + 1, genosOrig.GetNumSites());

  // fill in old values
  for (int i = 0; i < genosOrig.GetNumHaps(); ++i) {
    pResMat->SetGenotypeName(i, genosOrig.GetGenotypeName(i));
    for (int j = 0; j < genosOrig.GetNumSites(); ++j) {
      pResMat->SetGenotypeAt(i, j, genosOrig.GetGenotypeAt(i, j));
      pResMat->SetGenotypeProbAt(i, j,
                                 genosOrig.GetGenotypeProbAllele0At(i, j));
    }
  }
#if 0
    // fill in imputed values values
    for(int i=0; i<genosDoubleInfer.GetNumHaps()-1; ++i)
    {
        int hapUse=i;
        if(i >= indexDouble)
        {
            hapUse=i+1;
        }

        for(int j=0; j<genosDoubleInfer.GetNumSites(); ++j)
        {
            pResMat->SetGenotypeAt( hapUse, j, genosDoubleInfer.GetGenotypeAt(i,j) );
            pResMat->SetGenotypeProbAt( hapUse, j, genosDoubleInfer.GetGenotypeProbAllele0At(i,j) );
        }
    }
#endif

  // fill in imputed dobulet genos (two last rows)
  pResMat->SetGenotypeName(genosOrig.GetNumHaps(),
                           GetNewGenoDoubleRowName(genosOrig, indexDouble));
  for (int s = 0; s < genosOrig.GetNumSites(); ++s) {
    double p0 = genosOrig.GetGenotypeProbAllele0At(indexDouble, s);
    int g1 =
        genosDoubleInfer.GetGenotypeAt(genosDoubleInfer.GetNumHaps() - 2, s);
    pResMat->SetGenotypeAt(indexDouble, s, g1);
    double p0Use1 = p0;
    if ((g1 == 0 && p0 < 0.5) || (g1 == 1 && p0 > 0.5)) {
      p0Use1 = 1.0 - p0;
    }
    pResMat->SetGenotypeProbAt(indexDouble, s, p0Use1);
    int g2 =
        genosDoubleInfer.GetGenotypeAt(genosDoubleInfer.GetNumHaps() - 1, s);
    pResMat->SetGenotypeAt(genosOrig.GetNumHaps(), s, g2);
    double p0Use2 = p0;
    if ((g2 == 0 && p0 < 0.5) || (g2 == 1 && p0 > 0.5)) {
      p0Use2 = 1.0 - p0;
    }
    pResMat->SetGenotypeProbAt(genosOrig.GetNumHaps(), s, p0Use2);
  }

  return pResMat;
}

// construct matrix that is constructed from doublet result
void ScistDoubletSearch ::FindOrigImputedGeno(
    const ScistGenGenotypeMat &genosDoubletRes,
    ScistGenGenotypeMat &genosImpute,
    std::map<int, std::set<int> > &mapDoublets) const {
  // cout << "FindOrigImputedGeno: genosDoubletRes: ";
  // genosDoubletRes.Dump();
  mapDoublets.clear();
  // match any row
  map<string, int> mapNameToRowIndexDouble;
  for (int i = 0; i < genosImpute.GetNumHaps(); ++i) {
    //
    mapNameToRowIndexDouble[genosImpute.GetGenotypeName(i)] = i;
  }

  // first copy any row that is not double
  set<int> rowsDouble;
  for (int i = 0; i < genosDoubletRes.GetNumHaps(); ++i) {
    if (mapNameToRowIndexDouble.find(genosDoubletRes.GetGenotypeName(i)) !=
        mapNameToRowIndexDouble.end()) {
      // copy
      int index = mapNameToRowIndexDouble[genosDoubletRes.GetGenotypeName(i)];
      for (int j = 0; j < genosDoubletRes.GetNumSites(); ++j) {
        genosImpute.SetGenotypeAt(index, j,
                                  genosDoubletRes.GetGenotypeAt(i, j));
      }
    } else {
      rowsDouble.insert(i);
    }
  }
  // cout << "RowsDouble: ";
  // DumpIntSet(rowsDouble);
  // now add those
  for (set<int>::iterator it = rowsDouble.begin(); it != rowsDouble.end();
       ++it) {
    int i = *it;
    string strName = GetGenoDoubleRowName(genosDoubletRes.GetGenotypeName(i));
    YW_ASSERT_INFO(mapNameToRowIndexDouble.find(strName) !=
                       mapNameToRowIndexDouble.end(),
                   "Fail to find the row");

    // copy
    int index = mapNameToRowIndexDouble[strName];
    for (int j = 0; j < genosDoubletRes.GetNumSites(); ++j) {
      genosImpute.AddGenotypeAt(index, j, genosDoubletRes.GetGenotypeAt(i, j));
    }

    // record it
    int strNameInt = std::stoi(strName);
    mapDoublets[strNameInt].insert(i);
    mapDoublets[strNameInt].insert(genosDoubletRes.FindCellByName(strName));
  }
}

string ScistDoubletSearch ::GetGenoDoubleRowName(const string &strName) const {
  // if last character is '
  if (strName.length() > 0 && strName.substr(strName.length() - 1, 1) == "'") {
    // return the portion that doesn't have trailing '
    int posNone = strName.find_last_not_of("'");
    return strName.substr(0, posNone + 1);
  }
  YW_ASSERT_INFO(false, "The row is doublet");
  string strDummy;
  return strDummy;
}

ScistGenGenotypeMat *
ScistDoubletSearch ::EvalGenoDoubletSet(const ScistGenGenotypeMat &matToSearch,
                                        const set<int> &setDoubletRows,
                                        double &resOpt) {
  //
  resOpt = setDoubletRows.size() * this->costDoublet;
  set<int> setDoubleRowsConv;
  double costInit = 0.0;
  ScistGenGenotypeMat *pMatDouble = InitSearchGenotypes(
      matToSearch, setDoubletRows, setDoubleRowsConv, costInit);
  resOpt += costInit;
  // cout << "costInit: " << costInit << ", matrixDouble: ";
  // pMatDouble->Dump();

  if (setDoubletRows.size() == 0) {
    return pMatDouble;
  }

  // now score doublet
  set<int> rowsTemplate;
  PopulateSetWithInterval(rowsTemplate, 0, pMatDouble->GetNumHaps() - 1);
  SubtractSets(rowsTemplate, setDoubleRowsConv);

  // each time pick the lowest cost change to resolve doublets
  while (setDoubleRowsConv.size() > 0) {
    // cout << "setDoubleRowsConv: ";
    // DumpIntSet(setDoubleRowsConv);
    // cout << "rowsTemplate: ";
    // DumpIntSet(rowsTemplate);
    // evaluate each
    set<int> rowsDone;

    double optBest = HAP_MAX_INT * 1.0;
    vector<int> vecHap1, vecHap2;
    int rowBest = -1;
    for (set<int>::iterator it = setDoubleRowsConv.begin();
         it != setDoubleRowsConv.end(); ++it) {
      if (rowsDone.find(*it) != rowsDone.end()) {
        continue;
      }

      vector<int> vecHap1Step, vecHap2Step;
      double optStep = ScoreDoubletRow(pMatDouble, rowsTemplate, *it,
                                       vecHap1Step, vecHap2Step);

      // cout << "ScoreDoubleRow for row " << *it << ", two resolved haplotypes:
      // "; DumpIntVec(vecHap1Step); DumpIntVec(vecHap2Step);

      // if there is no change of doublets, stop
      if (IsAllZeroVec(vecHap1Step) || IsAllZeroVec(vecHap2Step) ||
          vecHap1Step == vecHap2Step) {
        // this is trivial doublet, stop
        break;
      }

      if (optStep < optBest) {
        optBest = optStep;
        vecHap1 = vecHap1Step;
        vecHap2 = vecHap2Step;
        rowBest = *it;
        // cout << "better....\n";
      }

      rowsDone.insert(*it);
      rowsDone.insert(*it + 1);
    }

    if (rowBest < 0) {
      delete pMatDouble;
      pMatDouble = NULL;
      break;
    }

    // take the best one
    YW_ASSERT_INFO(rowBest >= 0, "Wrong");
    resOpt += optBest;
    // cout << "**Resolve double: optBest: " << optBest << ", rowBest: " <<
    // rowBest << endl; cout << "vecHap1: "; DumpIntVec(vecHap1); cout <<
    // "vecHap2: "; DumpIntVec(vecHap2);
    UpdateSearchGenotypes(pMatDouble, rowBest, vecHap1, vecHap2);

    // cout << "Evl step matrix: ";
    // pMatDouble->Dump();

    // ScistHaplotypeMat *pMatResHap0 = dynamic_cast<ScistHaplotypeMat
    // *>(pMatDouble); string strTreeEdgeLabel0 =
    // ConsRootedPerfectPhylogenyFromMat(pMatResHap0->GetHapMat(), true, true);
    // cout << "EvalGenoDoubletSet tree (step): " << strTreeEdgeLabel0 << endl;

    setDoubleRowsConv.erase(rowBest);
    setDoubleRowsConv.erase(rowBest + 1);
    rowsTemplate.insert(rowBest);
    rowsTemplate.insert(rowBest + 1);
  }

  return pMatDouble;
}

void ScistDoubletSearch ::FindDoubletCandidates(set<int> &candidatesDoublet) {
  // for now, each row can be a doublet
  candidatesDoublet.clear();
  PopulateSetWithInterval(candidatesDoublet, 0,
                          this->genosInput.GetNumHaps() - 1);
}

ScistGenGenotypeMat *
ScistDoubletSearch ::InitSearchGenotypes(const ScistGenGenotypeMat &matToSearch,
                                         const set<int> &candidatesDoubletCurr,
                                         set<int> &setDoubletRows,
                                         double &costInit) {
  // cout << "candidatesDoubletCurr: ";
  // DumpIntSet(candidatesDoubletCurr);
  // cout << "matToSearch: ";
  // matToSearch.Dump();
  // in the new matrix to work with, put the single genotype together, and then
  // put the doublets later
  ScistGenGenotypeMat *pMatToProc = new ScistHaplotypeMat();
  int numHapsNew = matToSearch.GetNumHaps() + (int)candidatesDoubletCurr.size();
  pMatToProc->SetSize(numHapsNew, matToSearch.GetNumSites());

  // fill single rows
  set<int> setTemplateRows;
  int hapCur = 0;
  for (int i = 0; i < matToSearch.GetNumHaps(); ++i) {
    if (candidatesDoubletCurr.find(i) != candidatesDoubletCurr.end()) {
      continue;
    }
    // copy it
    for (int s = 0; s < matToSearch.GetNumSites(); ++s) {
      pMatToProc->SetGenotypeAt(hapCur, s, matToSearch.GetGenotypeAt(i, s));
      pMatToProc->SetGenotypeProbAt(hapCur, s,
                                    matToSearch.GetGenotypeProbAllele0At(i, s));
    }
    // set name
    pMatToProc->SetGenotypeName(hapCur, matToSearch.GetGenotypeName(i));
    setTemplateRows.insert(hapCur);

    ++hapCur;
  }
  // cout << "After filling single rows: pMatToProc: ";
  // pMatToProc->Dump();
  // now  copy the doublet rows
  for (int i = 0; i < matToSearch.GetNumHaps(); ++i) {
    if (candidatesDoubletCurr.find(i) == candidatesDoubletCurr.end()) {
      continue;
    }
    // copy it
    for (int s = 0; s < matToSearch.GetNumSites(); ++s) {
      pMatToProc->SetGenotypeAt(hapCur, s, matToSearch.GetGenotypeAt(i, s));
      pMatToProc->SetGenotypeProbAt(hapCur, s,
                                    matToSearch.GetGenotypeProbAllele0At(i, s));
      pMatToProc->SetGenotypeAt(hapCur + 1, s, matToSearch.GetGenotypeAt(i, s));
      pMatToProc->SetGenotypeProbAt(hapCur + 1, s,
                                    matToSearch.GetGenotypeProbAllele0At(i, s));
    }
    // set name
    pMatToProc->SetGenotypeName(hapCur, matToSearch.GetGenotypeName(i));
    string strName1 = GetNewGenoDoubleRowName(matToSearch, i);
    pMatToProc->SetGenotypeName(hapCur + 1, strName1);

    setDoubletRows.insert(hapCur);
    setDoubletRows.insert(hapCur + 1);

    hapCur += 2;
  }
  // cout << "After filling double rows: ";
  // pMatToProc->Dump();

  // now fit perfect phylogeny
  costInit = FitPerfPhyFor(pMatToProc, setTemplateRows);

  // cout << "Inflated genotype matrix: ";
  // pMatToProc->Dump();

  return pMatToProc;
}

std::string ScistDoubletSearch ::GetNewGenoDoubleRowName(
    const ScistGenGenotypeMat &matToSearch, int index) const {
  // find a new name for the doublet s.t. it is new
  string strName1 = matToSearch.GetGenotypeName(index) + "'";
  while (matToSearch.FindCellByName(strName1) >= 0) {
    strName1 = strName1 + "'";
  }
  return strName1;
}

void ScistDoubletSearch ::UpdateSearchGenotypes(
    ScistGenGenotypeMat *pMatCurr, int genoDoublet,
    const vector<int> &genoDoublePhase1, const vector<int> &genoDoublePhase2) {
  // fill in the new values the two rows genoDouble and the next row
  YW_ASSERT_INFO(pMatCurr->GetNumSites() == (int)genoDoublePhase1.size(),
                 "Wrong size");
  for (int s = 0; s < pMatCurr->GetNumSites(); ++s) {
    pMatCurr->SetGenotypeAt(genoDoublet, s, genoDoublePhase1[s]);
    pMatCurr->SetGenotypeAt(genoDoublet + 1, s, genoDoublePhase2[s]);
  }
}

double ScistDoubletSearch ::ScoreDoubletRow(ScistGenGenotypeMat *pMatCurr,
                                            const set<int> &rowsTemplate,
                                            int rowDouble,
                                            vector<int> &genoDoublePhase1,
                                            vector<int> &genoDoublePhase2) {
  // cout << "ScistDoubletSearch :: ScoreDoubletRow: curr mat: ";
  // pMatCurr->Dump();
  //
  ScistDoublet sciDouble(*pMatCurr);
  return sciDouble.EvalGenoDoublet(rowsTemplate, rowDouble, genoDoublePhase1,
                                   genoDoublePhase2);
}

double
ScistDoubletSearch ::FitPerfPhyFor(ScistGenGenotypeMat *pMatCurr,
                                   const std::set<int> &setTemplateRows) {
  // cout << "template rows: ";
  // DumpIntSet(setTemplateRows);
  // cout << "Current matrix: ";
  // pMatCurr->Dump();
  // Make the chosen rows to be perfect phylogeny
  set<int> sitesUse;
  PopulateSetWithInterval(sitesUse, 0, this->genosInput.GetNumSites() - 1);

  // create a submatrix to fit perfect phylogeny
  ScistGenGenotypeMat *pMatSub = pMatCurr->SubMatrix(setTemplateRows, sitesUse);
  // cout << "Submatrix: ";
  // pMatSub->Dump();
  ScistPerfPhyMLE sciInf1(*pMatSub);
  sciInf1.SetOutput(false);
  sciInf1.SetVerbose(false);
  double opt = -1.0 * sciInf1.Infer();

  // update genotype
  // cout << "After perfect phylogeny fitting: genotypes are: opt = " << opt <<
  // ": "; pMatSub->Dump();
  int rowCur = 0;
  for (set<int>::const_iterator it = setTemplateRows.begin();
       it != setTemplateRows.end(); ++it) {
    for (int s = 0; s < pMatCurr->GetNumSites(); ++s) {
      pMatCurr->SetGenotypeAt(*it, s, pMatSub->GetGenotypeAt(rowCur, s));
    }
    ++rowCur;
  }
  // cout << "After perfect phylogeny fitting, current matrix: ";
  // pMatCurr->Dump();

  delete pMatSub;
  return opt;
}

std::string
ScistDoubletSearch ::ConvMutTreeStr(const std::string &strTree) const {
  //
  if (this->listSiteNames.size() == 0) {
    // no conversion if no cell names specified
    return strTree;
  }

  TaxaMapper taxaMapper;
  for (int i = 0; i < (int)listSiteNames.size(); ++i) {
    taxaMapper.AddTaxaStringWithId(i + 1, listSiteNames[i]);
  }
  //
  return taxaMapper.ConvIdStringWithOrigTaxa(strTree);
}

// *************************************************************************************

void ScistDoubletTest() {
  ScistHaplotypeMat genoMat;
  const int numSCs = 5, numSites = 3;
  genoMat.SetSize(numSCs, numSites);
  genoMat.SetGenotypeAt(0, 0, 0);
  genoMat.SetGenotypeAt(0, 1, 0);
  genoMat.SetGenotypeAt(0, 2, 1);
  genoMat.SetGenotypeAt(1, 0, 0);
  genoMat.SetGenotypeAt(1, 1, 1);
  genoMat.SetGenotypeAt(1, 2, 0);
  genoMat.SetGenotypeAt(2, 0, 1);
  genoMat.SetGenotypeAt(2, 1, 1);
  genoMat.SetGenotypeAt(2, 2, 0);
  genoMat.SetGenotypeAt(3, 0, 1);
  genoMat.SetGenotypeAt(3, 1, 1);
  genoMat.SetGenotypeAt(3, 2, 0);
  genoMat.SetGenotypeAt(4, 0, 1);
  genoMat.SetGenotypeAt(4, 1, 0);
  genoMat.SetGenotypeAt(4, 2, 1);

  // genoMat.SetSize(numSCs, numSites);
  genoMat.SetGenotypeProbAt(0, 0, 0.8);
  genoMat.SetGenotypeProbAt(0, 1, 0.8);
  genoMat.SetGenotypeProbAt(0, 2, 0.1);
  genoMat.SetGenotypeProbAt(1, 0, 0.8);
  genoMat.SetGenotypeProbAt(1, 1, 0.1);
  genoMat.SetGenotypeProbAt(1, 2, 0.8);
  genoMat.SetGenotypeProbAt(2, 0, 0.1);
  genoMat.SetGenotypeProbAt(2, 1, 0.1);
  genoMat.SetGenotypeProbAt(2, 2, 0.8);
  genoMat.SetGenotypeProbAt(3, 0, 0.1);
  genoMat.SetGenotypeProbAt(3, 1, 0.1);
  genoMat.SetGenotypeProbAt(3, 2, 0.8);
  genoMat.SetGenotypeProbAt(4, 0, 0.3);
  genoMat.SetGenotypeProbAt(4, 1, 0.8);
  genoMat.SetGenotypeProbAt(4, 2, 0.1);

  const int SZ_DOUBLETS = 2;
  ScistDoubletSearch sds(genoMat, SZ_DOUBLETS);
  sds.Search();

#if 0
    set<int> setTemplateRows;
    setTemplateRows.insert(0);
    setTemplateRows.insert(1);
    setTemplateRows.insert(2);
    setTemplateRows.insert(3);
    int genoDoublet = 4;

    ScistDoublet ppInf( genoMat );
    vector<int> genoDouble1, genoDouble2;
    double minCost = ppInf.EvalGenoDoublet( setTemplateRows, genoDoublet, genoDouble1, genoDouble2 );
    cout << "Min-cost phasing: min-cost is " << minCost << endl;
    cout << "One haplotype: ";
    DumpIntVec(genoDouble1);
    cout << "Second haplotype: ";
    DumpIntVec(genoDouble2);
#endif
}
