//
//  ScistPerfPhyImp.cpp
//
//
//  Created by Yufeng Wu on 7/27/18.
//
//

#include "ScistPerfPhyImp.hpp"
#include "MarginalTree.h"
#include "PhylogenyTree.h"
#include "RBT.h"
#include "ScistGenotype.hpp"
#include "TreeBuilder.h"
#include "Utils3.h"
#include "Utils4.h"
#include "UtilsNumerical.h"
#include <cmath>
#include <iomanip>

const int MAX_SPR_OP = 1;

// *************************************************************************************
// Utiltiies

void OutputMutationTree(const char *filenameMT, const string &strMutTree,
                        bool fLabel) {
  PhylogenyTreeBasic treeMut;
  treeMut.ConsOnNewickEdgeLabelTree(strMutTree);
  if (fLabel) {
    treeMut.OutputGML(filenameMT);
  } else {
    treeMut.OutputGMLNoLabel(filenameMT);
  }
}

// *************************************************************************************
// Build phylogeny by tree search with branch length

ScistFullPerfPhyMLE ::ScistFullPerfPhyMLE(ScistGenGenotypeMat &genos)
    : genosInput(genos), fVerbose(false), pMargTreeOptBrLen(NULL),
      brOptIndex(-1) {
  Init();
}

void ScistFullPerfPhyMLE ::Infer() {
  set<ScistPerfPhyCluster> setClusAllGuide;
  this->treeGuide.GetAllClusters(setClusAllGuide);
  string strTreeOpt = ConsTreeFromSetClusters(setClusAllGuide);

  MarginalTree treeOpt;
  ReadinMarginalTreesNewickWLenString(strTreeOpt, this->genosInput.GetNumHaps(),
                                      treeOpt);
  treeOpt.InitUnitEdgelen();

  // double loglikeliOptInit = CalcLikelihoodOf(treeOpt);

  // optimize branch length
  double loglikeliOptBr = OptBranchLens(treeOpt);
  strTreeOpt = treeOpt.GetNewickSorted(true);
  // cout << "Initial tree: "  << treeOpt.GetNewick() << ", log-likelihood: " <<
  // loglikeliOptBr << endl;

  set<string> setTreeSearchedBefore;
  setTreeSearchedBefore.insert(strTreeOpt);

  // now search for neighborhood of the current tree to optimize the tree
  while (true) {
    set<string> setNgbrTrees;
    // GetNgbrTreesFromSPR( this->genosInput.GetNumHaps(), strTreeOpt,
    // setNgbrTrees );
    ScistPerfPhyMLE ::GetNgbrTreesFrom(this->genosInput.GetNumHaps(),
                                       strTreeOpt, setNgbrTrees);
    if (fVerbose) {
      cout << "Current best likelihood: " << loglikeliOptBr
           << ", current tree: " << treeOpt.GetNewickSorted(true)
           << ", tree neighborhood size: " << setNgbrTrees.size() << endl;
    }
    bool fCont = false;
    for (set<string>::iterator it = setNgbrTrees.begin();
         it != setNgbrTrees.end(); ++it) {
      if (setTreeSearchedBefore.find(*it) != setTreeSearchedBefore.end()) {
        continue;
      }
      setTreeSearchedBefore.insert(*it);

      // cout << "Neighbor tree: " << *it << endl;
      MarginalTree treeStep;
      ReadinMarginalTreesNewickWLenString(*it, this->genosInput.GetNumHaps(),
                                          treeStep);
      // treeStep.InitUnitEdgelen();
      // cout << "treeStep: " << treeStep.GetNewick() << endl;
      double loglikeliStep = OptBranchLens(treeStep);
      // double loglikeliStep = CalcLikelihoodOf( treeStep );
      // cout << ", loglikeliStep (w/ branch length optimization): " <<
      // loglikeliStep << endl;
      if (loglikeliStep > loglikeliOptBr) {
        // cout << "BETTER.\n";
        loglikeliOptBr = loglikeliStep;
        strTreeOpt = *it;
        treeOpt = treeStep;
        fCont = true;
      }
    }
    if (fCont == false) {
      break;
    }
  }

  cout << "**** Optimal cost: " << loglikeliOptBr << endl;
  cout << "Constructed single cell phylogeny: "
       << treeOpt.GetNewickSorted(false) << endl;
  cout << "With branch length: " << treeOpt.GetNewickSorted(true) << endl;
}

void ScistFullPerfPhyMLE ::Init() {
  //
  cacheProbMutClades.resize(genosInput.GetNumSites());
  // get all clusters
  // listClusMutsInput.clear();
  // for(int s=0; s<genosInput.GetNumSites(); ++s)
  //{
  //    set<int> muts;
  //    genosInput.GetMutRowsHapAtSite(s, muts);
  //    ScistPerfPhyCluster clus(muts);
  //    listClusMutsInput.push_back(clus);
  //}
  listClusMutsInputHetero.clear();
  listClusMutsInputHomo.clear();
  for (int s = 0; s < genosInput.GetNumSites(); ++s) {
    set<int> muts;
    genosInput.GetRowsWithGenoAtSite(s, 1, muts);
    ScistPerfPhyCluster clus(muts);
    listClusMutsInputHetero.push_back(clus);

    set<int> muts2;
    genosInput.GetRowsWithGenoAtSite(s, 2, muts2);
    ScistPerfPhyCluster clus2(muts2);
    listClusMutsInputHomo.push_back(clus2);
  }

  this->genosInput.GetColMultiplicityMap(listInputColMulti);

  // construct NJ tree as the initial tree
  string strNJ = this->genosInput.ConsNJTreeZeroRoot();
  this->treeGuide.Init(strNJ);
}

double ScistFullPerfPhyMLE ::OptBranchLens(MarginalTree &tree) {
  //
  this->pMargTreeOptBrLen = &tree;

  const double MIN_BR_LEN = 0.01;
  const double MAX_BR_LEN = 10.0;
  const double TOLNUM = 0.2;

  double loglikeliRes = -1.0 * HAP_MAX_INT;

  // optimize branch of each once and only once
  for (int br = 0; br < tree.GetTotNodesNum(); ++br) {
    if (br == tree.GetRoot()) {
      continue;
    }
    this->brOptIndex = br;
    double brLen = tree.GetEdgeLen(br);
    double brNew = brLen;
    double likeliMax =
        -1.0 * Func1DMinBrent(MIN_BR_LEN, brLen, MAX_BR_LEN, TOLNUM, &brNew);
    if (likeliMax > loglikeliRes) {
      loglikeliRes = likeliMax;
      tree.SetBranchLen(br, brNew);
    } else {
      tree.SetBranchLen(br, brLen);
    }
  }
  return loglikeliRes;
}

double ScistFullPerfPhyMLE ::EvaluateAt(double pt, void *pParam) {
  //
  YW_ASSERT_INFO(pMargTreeOptBrLen != NULL, "Tree to opt branch: null");
  YW_ASSERT_INFO(brOptIndex >= 0, "Branch opt not set");
  pMargTreeOptBrLen->SetBranchLen(brOptIndex, pt);
  return -1.0 * CalcLikelihoodOf(*pMargTreeOptBrLen);
}

double ScistFullPerfPhyMLE ::CalcLikelihoodOf(MarginalTree &tree) const {
  set<pair<ScistPerfPhyCluster, ScistPerfPhyCluster> > setClusDone;
  double res = 0.0;

  vector<set<int> > listClades;
  tree.ConsDecedentLeavesInfoLabels(listClades);
  for (int i = 0; i < (int)listClades.size(); ++i) {
    DecAllNumInSet(listClades[i]);
    // cout << "Tree clade: ";
    // DumpIntSet(listClades[i]);
  }
  double totEdgeLen = tree.GetTotEdgeLen();
  ScistPerfPhyProbOnTree sppp(this->genosInput, tree);

  for (int site = 0; site < genosInput.GetNumSites(); ++site) {
    pair<ScistPerfPhyCluster, ScistPerfPhyCluster> pp(
        listClusMutsInputHetero[site], listClusMutsInputHomo[site]);
    if (setClusDone.find(pp) != setClusDone.end()) {
      continue;
    }
    int multi = this->listInputColMulti[site];
    double loglikeliSite =
        CalcLikelihoodOf(sppp, site, tree, totEdgeLen, listClades);
    res += loglikeliSite * multi;
    setClusDone.insert(pp);
  }

  return res;
}

double ScistFullPerfPhyMLE ::CalcLikelihoodOf(
    ScistPerfPhyProbOnTree &sppp, int site, MarginalTree &tree,
    double totEdgeLen, const vector<set<int> > &listClades) const {
  return sppp.CalcProbForSite(site, totEdgeLen, listClades);
}

std::string ScistFullPerfPhyMLE ::ConsTreeFromSetClusters(
    const std::set<ScistPerfPhyCluster> &setClusters) const {
  //
  // now construct tree
  ScistInfPerfPhyUtils treeBuild;
  map<int, ScistPerfPhyCluster> mapPickedClus;
  int s = 0;
  for (set<ScistPerfPhyCluster>::iterator it = setClusters.begin();
       it != setClusters.end(); ++it) {
    mapPickedClus[s] = *it;
    ++s;
  }
  string strTree =
      treeBuild.ConsTreeWCombDistClus(this->genosInput, mapPickedClus);
  return strTree;
}

// *************************************************************************************
// Build phylogeny by tree search

ScistPerfPhyMLE ::ScistPerfPhyMLE(ScistGenGenotypeMat &genos)
    : genosInput(genos), fVerbose(false), fOptBrLen(false), fOutput(true),
      fOutputPPWithEdgeLabels(false), fOutputLabel(true), fSPR(false),
      maxSPRNum(MAX_SPR_OP) {
  Init();
}

double ScistPerfPhyMLE ::Infer(
    std::set<std::pair<std::pair<int, int>, int> > *plistChangedPlaces,
    std::string *pstrTreeNW) {
  // cout << "ScistPerfPhyMLE :: Infer\n";
  //
  set<ScistPerfPhyCluster> setClusAllGuide;
  this->treeGuide.GetAllClusters(setClusAllGuide);
  // cout << "Number of clusters: " << setClusAllGuide.size() << endl;
  string strTreeOpt = ConsTreeFromSetClusters(setClusAllGuide);
  // cout << "strTreeOpt: " << strTreeOpt << endl;
  // set<ScistPerfPhyCluster> setClusAllGuideUse;
  // GetClustersFromTree(strTreeOpt, setClusAllGuideUse);
  std::vector<pair<ScistPerfPhyCluster, ScistPerfPhyCluster> >
      listChangedClustersOpt;
  // double loglikeliBest = ScoreSetClusters( setClusAllGuideUse,
  // listChangedClustersOpt );
  double loglikeliBest = ScoreTree(strTreeOpt, listChangedClustersOpt);
  // cout << "Init likelihood: " << loglikeliBest << endl;
  set<string> setTreeSearchedBefore;
  setTreeSearchedBefore.insert(strTreeOpt);

  // now search for neighborhood of the current tree to optimize the tree
  int numSPRPerformed = 0;
  bool fNNI = true;
  while (true) {
    // if(fNNI)
    //{
    // cout << "NNI mode\n";
    //}
    // else
    //{
    // cout << "SPR mode\n";
    //}

    set<string> setNgbrTrees;
    if (fNNI == false && fSPR && numSPRPerformed <= maxSPRNum) {
      GetNgbrTreesFromSPR(this->genosInput.GetNumHaps(), strTreeOpt,
                          setNgbrTrees);
      ++numSPRPerformed;
      // fNNI = true;
    } else if (fNNI == true) {
      GetNgbrTreesFrom(this->genosInput.GetNumHaps(), strTreeOpt, setNgbrTrees);
    } else {
      break;
    }
    if (fVerbose) {
      cout << "Current best likelihood: " << loglikeliBest
           << ", current cost: " << CalcMaxProbUpperBound() - loglikeliBest
           << ", opt tree: " << strTreeOpt
           << ", tree neighborhood size: " << setNgbrTrees.size() << endl;
    }
    // cout << "Current opt tree: " << strTreeOpt << endl;
    bool fCont = false;
    for (set<string>::iterator it = setNgbrTrees.begin();
         it != setNgbrTrees.end(); ++it) {
      if (setTreeSearchedBefore.find(*it) != setTreeSearchedBefore.end()) {
        continue;
      }
      setTreeSearchedBefore.insert(*it);

      // cout << "Neighbor tree: " << *it << endl;
      // set<ScistPerfPhyCluster> setClus;
      // GetClustersFromTree(*it, setClus);
      vector<pair<ScistPerfPhyCluster, ScistPerfPhyCluster> >
          listChangedClustersStep;
      // double loglikeliStep = ScoreSetClusters( setClus,
      // listChangedClustersStep);
      double loglikeliStep = ScoreTree(*it, listChangedClustersStep);
      // cout << ", loglikeliStep: " << loglikeliStep << ", cost: " <<
      // CalcMaxProbUpperBound()- loglikeliStep << endl; if( loglikeliStep <
      // loglikeliBest )
      if (loglikeliStep > loglikeliBest) {
        // cout << "BETTER.\n";
        loglikeliBest = loglikeliStep;
        strTreeOpt = *it;
        listChangedClustersOpt = listChangedClustersStep;
        fCont = true;
      }
    }
    if (fCont == false) {
      if (fNNI == false) {
        break;
      }

      fNNI = false;
      // break;
    } else {
      fNNI = true;
    }
  }
  // output the final tree
  std::set<std::pair<std::pair<int, int>, int> > listChangedPlaces;
  for (int site = 0; site < this->genosInput.GetNumSites(); ++site) {
    FindChangedGenos(site, listChangedClustersOpt[site], listChangedPlaces);
  }
  if (plistChangedPlaces != NULL) {
    *plistChangedPlaces = listChangedPlaces;
  }
  if (pstrTreeNW != NULL) {
    *pstrTreeNW = strTreeOpt;
  }

  if (fVerbose) {
    if (fOutput) {
      cout << "Genotypes called by maximal single position probability\n";
      const string strDesc = "Single-site maximal probability genotypes";
      this->genosInput.OutputImput(&strDesc);
    }

    cout << "List of corrected genotypes (site, cell, new genotype) in base-1: "
            "\n";
    for (set<pair<pair<int, int>, int> >::iterator it =
             listChangedPlaces.begin();
         it != listChangedPlaces.end(); ++it) {
      cout << "[ " << setw(6) << it->first.second + 1 << " " << setw(6)
           << it->first.first + 1 << " ]: " << it->second << endl;
    }
  }

  if (fOutput) {
    // output the matrix
    ScistGenGenotypeMat *pMatRes = this->genosInput.Copy();
    pMatRes->ChangeGenosAtPositions(listChangedPlaces);
    if (fVerbose) {
      cout << "Called genotypes\n";
      pMatRes->OutputImput();
    }
    if (fOutputPPWithEdgeLabels) {
      ScistHaplotypeMat *pMatResHap =
          dynamic_cast<ScistHaplotypeMat *>(pMatRes);
      if (pMatResHap == NULL) {
        cout << "** Right now, only output perfect phylogeny for binary "
                "genotypes\n";
      } else {
        string strTreeEdgeLabel = ConsRootedPerfectPhylogenyFromMat(
            pMatResHap->GetHapMat(), true, true);
        // cout << "** Perfect phylogeny (with sites labeled on edges) from the
        // imputed genotypes: " << strTreeEdgeLabel << endl;

        string strMutTree = ConsEdgeLabeTree(strTreeEdgeLabel);
        string strMutTreeConv = ConvMutTreeStr(strMutTree);
        cout << "^^ Mutation tree: " << strMutTreeConv << endl;

        // output mutation tree file
        OutputMutationTree(this->strMutTreeFileName.c_str(), strMutTreeConv,
                           this->fOutputLabel);
      }
    }

    delete pMatRes;
  }

  // change genotype
  for (set<pair<pair<int, int>, int> >::iterator it = listChangedPlaces.begin();
       it != listChangedPlaces.end(); ++it) {
    this->genosInput.SetGenotypeAt(it->first.first, it->first.second,
                                   it->second);
  }

  double res = loglikeliBest;

  if (fOutput) {
    cout << "**** Maximum log-likelihood: " << loglikeliBest
         << ", number of changed genotypes: " << listChangedPlaces.size()
         << endl;
    cout << "Computed log-lielihood from changed genotypes: "
         << CalcChangedGenosProb(listChangedPlaces) << endl;
    // cout << "Minimum cost: " << CalcMaxProbUpperBound() - loglikeliBest <<
    // endl;

    string strTreeOptOut = ConvCellTreeStr(strTreeOpt);
    cout << "Constructed single cell phylogeny: " << strTreeOptOut << endl;
  }
  if (fOptBrLen) {
    string strTreeBrOpt;
    double loglikeliBestBr = OptBranchLens(strTreeOpt, strTreeBrOpt);
    res = loglikeliBestBr;
    if (fOutput) {
      cout << "**** Maximum log-likelihood (with branch length optimization): "
           << loglikeliBestBr << endl;
      string strTreeBrOptOut = ConvCellTreeStr(strTreeBrOpt);
      cout << "Single cell phylogeny with branch length: " << strTreeBrOptOut
           << endl;
    }
  }
  return res;
}

double ScistPerfPhyMLE ::OptBranchLens(const std::string &strTree,
                                       std::string &strTreeBrOpt) {
  //
  MarginalTree treeBrOpt;
  ReadinMarginalTreesNewickWLenString(strTree, this->genosInput.GetNumHaps(),
                                      treeBrOpt);
  ScistFullPerfPhyMLE sfpp(this->genosInput);
  double res = sfpp.OptBranchLens(treeBrOpt);
  strTreeBrOpt = treeBrOpt.GetNewickSorted(true);
  return res;
}

void ScistPerfPhyMLE ::Init() {
  //
  // get all clusters
  listClusMutsInputHetero.clear();
  listClusMutsInputHomo.clear();
  for (int s = 0; s < genosInput.GetNumSites(); ++s) {
    set<int> muts;
    genosInput.GetRowsWithGenoAtSite(s, 1, muts);
    ScistPerfPhyCluster clus(muts);
    listClusMutsInputHetero.push_back(clus);

    set<int> muts2;
    genosInput.GetRowsWithGenoAtSite(s, 2, muts2);
    ScistPerfPhyCluster clus2(muts2);
    listClusMutsInputHomo.push_back(clus2);
  }

  this->genosInput.GetColMultiplicityMap(listInputColMulti);

  // construct NJ tree as the initial tree
  string strNJ = this->genosInput.ConsNJTreeZeroRoot();
  // cout << "Guide tree: " << strNJ << endl;
  // string strNJ = this->genosInput.ConsNJTree();
  // cout << "Zero-rooted initial tree: " << strNJ << endl;
  // cout << "Genotype input: \n";
  // this->genosInput.Dump();
  //
  this->treeGuide.Init(strNJ);

  // set the prior score to be zero
  listSitePriorScore.clear();
  for (int i = 0; i < this->genosInput.GetNumSites(); ++i) {
    double logprobInit = 0.0;
    for (int h = 0; h < this->genosInput.GetNumHaps(); ++h) {
      double p = this->genosInput.GetGenotypeProbAllele0At(h, i);
      logprobInit += log(p);
    }
    listSitePriorScore.push_back(logprobInit);
  }
}

std::string ScistPerfPhyMLE ::ConsTreeFromSetClusters(
    const std::set<ScistPerfPhyCluster> &setClusters) const {
  // cout << "All the clusters: \n";
  // for(set<ScistPerfPhyCluster> :: const_iterator it = setClusters.begin(); it
  // != setClusters.end(); ++it)
  //{
  // it->Dump();
  //}
  //
  // now construct tree
  ScistInfPerfPhyUtils treeBuild;
  map<int, ScistPerfPhyCluster> mapPickedClus;
  int s = 0;
  for (set<ScistPerfPhyCluster>::iterator it = setClusters.begin();
       it != setClusters.end(); ++it) {
    mapPickedClus[s] = *it;
    ++s;
  }
  string strTree =
      treeBuild.ConsTreeWCombDistClus(this->genosInput, mapPickedClus, false);
  return strTree;
}

void ScistPerfPhyMLE ::GetNgbrTreesFrom(int numHaps, const std::string &strTree,
                                        std::set<std::string> &setNgbrTrees) {
  // cout << "GetNgbrTreesFrom: numHaps: " << numHaps << ", tree: " << strTree
  // << endl;
  //
  setNgbrTrees.clear();
  MarginalTree treeCurr;
  ReadinMarginalTreesNewickWLenString(strTree, numHaps, treeCurr);
  vector<MarginalTree> listNgbrTrees;
  FindOneNNIMTreesFrom(treeCurr, listNgbrTrees);
  for (int i = 0; i < (int)listNgbrTrees.size(); ++i) {
    string strTree = listNgbrTrees[i].GetNewickSorted(false);
    setNgbrTrees.insert(strTree);
  }
}

void ScistPerfPhyMLE ::GetNgbrTreesFromSPR(
    int numHaps, const std::string &strTree,
    std::set<std::string> &setNgbrTrees) {
  //
  setNgbrTrees.clear();
  MarginalTree treeCurr;
  ReadinMarginalTreesNewickWLenString(strTree, numHaps, treeCurr);
  string strSelf = treeCurr.GetNewickSorted(false);
  // cout << "strTree: " << strTree << ", strSelf: " << strSelf << endl;

  // map to consecutive order as required by RBT
  vector<int> listLeafLblsOld;
  // treeCurr.MapLeafLblConsecutiveOrder( listLeafLblsOld );
  treeCurr.GetLabelList(listLeafLblsOld);
  // cout << "Mapped leaves: ";
  // DumpIntVec(listLeafLblsOld);
  // cout << "Changed tree: " << treeCurr.GetNewick() << endl;

  // use RBT utility
  vector<int> listLbls;
  treeCurr.GetLabelList(listLbls);
  // cout << "listLbss: ";
  // DumpIntVec(listLbls);
  vector<int> parPosList;
  treeCurr.GetParPosInfo(parPosList);
  // cout << "parPosList: ";
  // DumpIntVec(parPosList);
  vector<double> listEdgeDistOut;
  treeCurr.GetTreeEdgeLen(listEdgeDistOut);
  RBT treeCurrRBT(numHaps, listLbls, parPosList, listEdgeDistOut);
  vector<RBT *> ngbrTrees;
  treeCurrRBT.FindSPRDistOneNgbrs(ngbrTrees);

  // cout << "GetNgbrTreesFromSPR: init tree: " << strTree << endl;
  for (int i = 0; i < (int)ngbrTrees.size(); ++i) {
    string strNW = ngbrTrees[i]->GetNewick();
    string strNWBack = RemapLeafLbls(numHaps, strNW, listLeafLblsOld);
    setNgbrTrees.insert(strNWBack);
    // cout << "strNW: " << strNW  << ", SPR tree: " << strNWBack << endl;
  }
  // remove self
  setNgbrTrees.erase(strSelf);

  for (int i = 0; i < (int)ngbrTrees.size(); ++i) {
    delete ngbrTrees[i];
  }
}

std::string ScistPerfPhyMLE ::RemapLeafLbls(int numHaps,
                                            const std::string &strTree0Based,
                                            const vector<int> &listLblsOld) {
  //
  MarginalTree treeCurr;
  ReadinMarginalTreesNewickWLenString(strTree0Based, numHaps, treeCurr);
  map<int, int> mapLblsBack;
  for (int i = 0; i < (int)listLblsOld.size(); ++i) {
    mapLblsBack[i] = listLblsOld[i];
  }
  treeCurr.RemapLeafLabels(mapLblsBack);
  return treeCurr.GetNewickSorted(false);
}

std::string
ScistPerfPhyMLE ::RemapLeafLbls(int numHaps, const std::string &strTree,
                                const std::map<int, int> &mapLabels) {
  //
  MarginalTree treeCurr;
  ReadinMarginalTreesNewickWLenString(strTree, numHaps, treeCurr);
  treeCurr.RemapLeafLabels(mapLabels);
  return treeCurr.GetNewickSorted(false);
}

std::string
ScistPerfPhyMLE ::ConvCellTreeStr(const std::string &strTree) const {
  //
  if (this->listCellNames.size() == 0) {
    // no conversion if no cell names specified
    return strTree;
  }

  TaxaMapper taxaMapper;
  for (int i = 0; i < (int)listCellNames.size(); ++i) {
    taxaMapper.AddTaxaStringWithId(i + 1, listCellNames[i]);
  }
  //
  return taxaMapper.ConvIdStringWithOrigTaxa(strTree);
}

std::string ScistPerfPhyMLE ::ConvMutTreeStr(const std::string &strTree) const {
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

void ScistPerfPhyMLE ::FindChangedGenos(
    int siteToAdd,
    const pair<ScistPerfPhyCluster, ScistPerfPhyCluster> &clusToAdd,
    set<pair<pair<int, int>, int> > &listChangedPlaces) const {
  // find list of positions where the genos are changed
  ScistPerfPhyCluster clusInt, clusThisOnly, clusRHSOnly;
  clusToAdd.first.IntersectWith(listClusMutsInputHetero[siteToAdd], clusInt,
                                clusThisOnly, clusRHSOnly);
  ScistPerfPhyCluster clusInt2, clusThisOnly2, clusRHSOnly2;
  clusToAdd.second.IntersectWith(listClusMutsInputHomo[siteToAdd], clusInt2,
                                 clusThisOnly2, clusRHSOnly2);
  // get changed 0
  set<int> setss;
  PopulateSetWithInterval(setss, 0, this->genosInput.GetNumHaps() - 1);
  set<int> rows0Orig;
  this->genosInput.GetRowsWithGenoAtSite(siteToAdd, 0, rows0Orig);
  SubtractSets(setss, rows0Orig);
  ScistPerfPhyCluster clus0(setss);
  clus0.SubtractFrom(clusToAdd.first);
  clus0.SubtractFrom(clusToAdd.second);

  ScistPerfPhyClusterItor itor0(clus0);
  itor0.First();
  while (itor0.IsDone() == false) {
    int sc = itor0.GetCurrentSC();
    pair<int, int> pp(sc, siteToAdd);
    pair<pair<int, int>, int> pp0(pp, 0);
    listChangedPlaces.insert(pp0);
    itor0.Next();
  }

  // This only: new mutants
  ScistPerfPhyClusterItor itor1(clusThisOnly);
  itor1.First();
  while (itor1.IsDone() == false) {
    int sc = itor1.GetCurrentSC();
    pair<int, int> pp(sc, siteToAdd);
    pair<pair<int, int>, int> pp1(pp, 1);
    listChangedPlaces.insert(pp1);
    itor1.Next();
  }
  // RHS only: new wildtype
  ScistPerfPhyClusterItor itor2(clusThisOnly2);
  itor2.First();
  while (itor2.IsDone() == false) {
    int sc = itor2.GetCurrentSC();
    pair<int, int> pp(sc, siteToAdd);
    pair<pair<int, int>, int> pp2(pp, 2);
    listChangedPlaces.insert(pp2);
    itor2.Next();
  }
}

double ScistPerfPhyMLE ::ScoreTree(
    const string &strTree,
    std::vector<std::pair<ScistPerfPhyCluster, ScistPerfPhyCluster> >
        &listChangedCluster) const {
  // cout << "ScoreTree: tree: " << strTree << endl;
  // score the current tree
  MarginalTree treeToScore;
  ReadinMarginalTreesNewickWLenString(strTree, this->genosInput.GetNumHaps(),
                                      treeToScore);
  // cout << "Score tree: " << treeToScore.GetNewick() << endl;
  set<pair<ScistPerfPhyCluster, ScistPerfPhyCluster> > setClusDone;
  map<pair<ScistPerfPhyCluster, ScistPerfPhyCluster>,
      pair<ScistPerfPhyCluster, ScistPerfPhyCluster> >
      mapChangedClus;
  double res = 0.0;
  ScistPerfPhyProbOnTree probTree(this->genosInput, treeToScore);

  for (int site = 0; site < genosInput.GetNumSites(); ++site) {
    // cout << "ScoreTree: site " << site << " multi:" <<
    // this->listInputColMulti[site] << endl; cout << "Heterozygote clus: ";
    // listClusMutsInputHetero[site].Dump();
    // cout << "Homozygous clus: ";
    // listClusMutsInputHomo[site].Dump();
    pair<ScistPerfPhyCluster, ScistPerfPhyCluster> pp0(
        listClusMutsInputHetero[site], listClusMutsInputHomo[site]);
    if (setClusDone.find(pp0) != setClusDone.end()) {
      listChangedCluster.push_back(mapChangedClus[pp0]);
      continue;
    }
    int multi = this->listInputColMulti[site];
    pair<ScistPerfPhyCluster, ScistPerfPhyCluster> clusChanged;
    double loglikeliSite = ScoreTreeWithSite(
        probTree, treeToScore, site, clusChanged.first, clusChanged.second);
    mapChangedClus[pp0] = clusChanged;
    listChangedCluster.push_back(clusChanged);
    res += loglikeliSite * multi;
    setClusDone.insert(pp0);
    // cout << "site prob: " << loglikeliSite << ": clusChanged: ";
    // clusChanged.first.Dump();
    // cout << "  and ";
    // clusChanged.second.Dump();
  }

  return res;
}

double
ScistPerfPhyMLE ::ScoreTreeWithSite(ScistPerfPhyProbOnTree &probTree,
                                    MarginalTree &tree, int site,
                                    ScistPerfPhyCluster &clusChanged1,
                                    ScistPerfPhyCluster &clusChanged2) const {
  // cout << "site: " << site << ", tree: " << tree.GetNewickSorted(false) <<
  // endl;
  return probTree.CalcProbMaxForSite(site, clusChanged1, clusChanged2);
}

double ScistPerfPhyMLE ::CalcMaxProbUpperBound() const {
  //
  double res = 0.0;
  for (int s = 0; s < this->genosInput.GetNumSites(); ++s) {
    for (int h = 0; h < this->genosInput.GetNumHaps(); ++h) {
      double p0 = this->genosInput.GetGenotypeProbAllele0At(h, s);
      double p1 = 1 - p0;
      if (p0 >= p1) {
        res += log(p0);
      } else {
        res += log(p1);
      }
    }
  }
  return res;
}

double ScistPerfPhyMLE ::CalcChangedGenosProb(
    const std::set<std::pair<std::pair<int, int>, int> > &listChangedPlaces)
    const {
  //
  double res = 0.0;
  map<pair<int, int>, int> mapChangedPlaces;
  for (std::set<std::pair<std::pair<int, int>, int> >::const_iterator it =
           listChangedPlaces.begin();
       it != listChangedPlaces.end(); ++it) {
    mapChangedPlaces[it->first] = it->second;
  }

  for (int s = 0; s < this->genosInput.GetNumSites(); ++s) {
    for (int h = 0; h < this->genosInput.GetNumHaps(); ++h) {
      pair<int, int> pp(h, s);
      int allele = this->genosInput.GetGenotypeAt(h, s);
      std::map<std::pair<int, int>, int>::const_iterator it =
          mapChangedPlaces.find(pp);
      if (it != mapChangedPlaces.end()) {
        int alleleAlt = it->second;
        YW_ASSERT_INFO(allele == alleleAlt, "Wrong");
        allele = alleleAlt;
      }

      double p0 = this->genosInput.GetGenotypeProbAllele0At(h, s);
      double p1 = 1 - p0;
      if (allele == 0) {
        res += log(p0);
      } else {
        res += log(p1);
      }
    }
  }

  return res;
}

// *************************************************************************************
// Tree probability

ScistPerfPhyProbOnTree ::ScistPerfPhyProbOnTree(ScistGenGenotypeMat &genos,
                                                MarginalTree &mtreeIn)
    : genosInput(genos), mtree(mtreeIn) {
  // set the prior score to be zero
  listSitePriorScore.clear();
  for (int i = 0; i < this->genosInput.GetNumSites(); ++i) {
    double logprobInit = 0.0;
    for (int h = 0; h < this->genosInput.GetNumHaps(); ++h) {
      double p = this->genosInput.GetGenotypeProbAt(h, i, 0);
      logprobInit += log(p);
    }
    listSitePriorScore.push_back(logprobInit);
  }
  Init();
}

void ScistPerfPhyProbOnTree ::Init() {
  //
  ScistTernaryMat *pGenoMat =
      dynamic_cast<ScistTernaryMat *>(&this->genosInput);
  if (pGenoMat == NULL) {
    return; // only work with genotype data
  }
  this->genosInputHap.SetSize(this->genosInput.GetNumHaps(),
                              this->genosInput.GetNumSites() * 2);
  for (int h = 0; h < this->genosInput.GetNumHaps(); ++h) {
    for (int s = 0; s < this->genosInput.GetNumSites(); ++s) {
      double p0 = pGenoMat->GetGenotypeProbAt(h, s, 0);
      double p1 = pGenoMat->GetGenotypeProbAt(h, s, 1);
      double p2 = pGenoMat->GetGenotypeProbAt(h, s, 2);
      double p12 = p1 + p2;
      double p01 = p0 + p1;
      int allele0 = 0;
      if (p0 < p12) {
        allele0 = 1;
      }
      this->genosInputHap.SetGenotypeAt(h, 2 * s, allele0);
      this->genosInputHap.SetGenotypeProbAt(h, 2 * s, p0);
      int allele1 = 0;
      if (p01 < p2) {
        allele1 = 1;
      }
      this->genosInputHap.SetGenotypeAt(h, 2 * s + 1, allele1);
      this->genosInputHap.SetGenotypeProbAt(h, 2 * s + 1, p01);
    }
  }
}

double ScistPerfPhyProbOnTree ::CalcProbMaxForSite(
    int site, ScistPerfPhyCluster &clusChangedMut,
    ScistPerfPhyCluster &clusChangedHomoMut) const {
  ScistHaplotypeMat *pHapMat =
      dynamic_cast<ScistHaplotypeMat *>(&this->genosInput);

  if (pHapMat != NULL) {
    clusChangedHomoMut.Clear();
    return CalcProbMaxForSiteHap(site, clusChangedMut);
  } else // right now, must be of genotype
  {
    return CalcProbMaxForSiteGeno(site, clusChangedMut, clusChangedHomoMut);
  }
}

double ScistPerfPhyProbOnTree ::CalcProbMaxForSiteHap(
    int site, ScistPerfPhyCluster &clusChanged) const {
  // cout << "ScoreTreeWithSite: tree: " << tree.GetNewick() << ", site: " <<
  // site << endl;
  // score the site wrt the tree (i.e. find the best split of the tree for this
  // site)
  double res = -1.0 * HAP_MAX_INT;
  // do a bottom up
  vector<double> listNodeSplitProb;
  // init to be bad
  for (int node = 0; node < mtree.GetTotNodesNum(); ++node) {
    listNodeSplitProb.push_back(-1.0 * HAP_MAX_INT);
  }

  // cout << "CalcProbMaxForSiteHap: mtree: " << mtree.GetNewickSorted(false) <<
  // endl; mtree.Dump();

  int nodeOpt = -1;
  for (int node = 0; node < mtree.GetTotNodesNum(); ++node) {
    // cout << "node " << node << endl;
    if (node == mtree.GetRoot()) {
      // continue;
    }
    double logpStep;
    if (mtree.IsLeaf(node)) {
      // a single leaf in the split
      int lvlbl = mtree.GetLabel(node) - 1;
      // cout << "Leaf: " << lvlbl << endl;
      double p0 = this->genosInput.GetGenotypeProbAllele0At(lvlbl, site);
      if (p0 < YW_VERY_SMALL_FRACTION) {
        p0 = YW_VERY_SMALL_FRACTION;
      } else if (p0 > 1.0 - YW_VERY_SMALL_FRACTION) {
        p0 = 1.0 - YW_VERY_SMALL_FRACTION;
      }
      logpStep = log((1 - p0) / p0);
      // cout << "Set leaf " << node << " log prob to: " << logpStep << ", p0="
      // << p0 << endl;
    } else {
      // get the two children and add them up
      int childLeft = mtree.GetLeftDescendant(node);
      int childRight = mtree.GetRightDescendant(node);
      // cout << "node: " << node << ", childLeft: " << childLeft << ",
      // childRight: " << childRight << endl; cout << "childLeft: " << childLeft
      // << ", right: " << childRight << endl;

      YW_ASSERT_INFO(listNodeSplitProb[childLeft] > -1.0 * HAP_MAX_INT,
                     "Bad left");
      YW_ASSERT_INFO(listNodeSplitProb[childRight] > -1.0 * HAP_MAX_INT,
                     "Bad right1");
      logpStep = listNodeSplitProb[childLeft] + listNodeSplitProb[childRight];
    }
    // cout << "log prob: " << logpStep << " for node: " << node << endl;
    listNodeSplitProb[node] = logpStep;
    if (logpStep > res) {
      // cout << "Better at node: " << node << endl;
      res = logpStep;
      nodeOpt = node;
    }
  }

  set<int> nodeOptSplitLbls;

  // if nothing is good, just take all-0
  if (res < 0.0) {
    //
    res = 0;
    nodeOpt = -1;
  } else {
    YW_ASSERT_INFO(nodeOpt >= 0, "Node not found");
    set<int> nodeOptSplit;
    mtree.GetLeavesUnder(nodeOpt, nodeOptSplit);
    mtree.GetlabelsFor(nodeOptSplit, nodeOptSplitLbls);
    DecAllNumInSet(nodeOptSplitLbls);
  }
  ScistPerfPhyCluster clus(nodeOptSplitLbls);
  clusChanged = clus;
  // cout << "Max prob at this site: " << res + this->listSitePriorScore[site]
  // << " at site " << nodeOpt << endl; cout << "clust changed: ";
  // clusChanged.Dump();
  return res + this->listSitePriorScore[site];
}

double ScistPerfPhyProbOnTree ::CalcProbMaxForSiteGeno(
    int site, ScistPerfPhyCluster &clusChangedHetero,
    ScistPerfPhyCluster &clusChangedHomo) const {
  //
  set<int> setSC0, setSC1, setSC2;
  this->genosInput.GetRowsWithGenoAtSite(site, 0, setSC0);
  this->genosInput.GetRowsWithGenoAtSite(site, 1, setSC1);
  this->genosInput.GetRowsWithGenoAtSite(site, 2, setSC2);

  // first accumulate for each node, the sum of diff p1/p0
  vector<double> vecSumDiffP10, vecSumDiffP21;
  vector<double> vecMaxSumDiff21;
  vector<int> vecMaxSumDiff21Node;
  for (int node = 0; node < mtree.GetTotNodesNum(); ++node) {
    double logpStep, logpStep2;
    if (mtree.IsLeaf(node)) {
      // a single leaf in the split
      int lvlbl = mtree.GetLabel(node) - 1;
      // cout << "Leaf: " << lvlbl << endl;
      double p0 = this->genosInput.GetGenotypeProbAt(lvlbl, site, 0);
      double p1 = this->genosInput.GetGenotypeProbAt(lvlbl, site, 1);
      double p2 = this->genosInput.GetGenotypeProbAt(lvlbl, site, 2);
      logpStep = log(p1 / p0);
      logpStep2 = log(p2 / p1);
      vecMaxSumDiff21.push_back(logpStep2);
      vecMaxSumDiff21Node.push_back(node);
    } else {
      // get the two children and add them up
      int childLeft = mtree.GetLeftDescendant(node);
      int childRight = mtree.GetRightDescendant(node);
      // cout << "childLeft: " << childLeft << ", right: " << childRight <<
      // endl;

      YW_ASSERT_INFO(vecSumDiffP10[childLeft] > -1.0 * HAP_MAX_INT,
                     "Bad left (geno)");
      YW_ASSERT_INFO(vecSumDiffP10[childRight] > -1.0 * HAP_MAX_INT,
                     "Bad right2");
      logpStep = vecSumDiffP10[childLeft] + vecSumDiffP10[childRight];
      logpStep2 = vecSumDiffP21[childLeft] + vecSumDiffP21[childRight];

      double maxSumLogp21 = logpStep2;
      int nodeMax = node;
      if (vecSumDiffP21[childLeft] > maxSumLogp21) {
        maxSumLogp21 = vecSumDiffP21[childLeft];
        nodeMax = vecMaxSumDiff21Node[childLeft];
      }
      if (vecSumDiffP21[childRight] > maxSumLogp21) {
        maxSumLogp21 = vecSumDiffP21[childRight];
        nodeMax = vecMaxSumDiff21Node[childRight];
      }
      vecMaxSumDiff21.push_back(maxSumLogp21);
      vecMaxSumDiff21Node.push_back(nodeMax);
    }
    // cout << "log prob: " << logpStep << endl;
    vecSumDiffP10.push_back(logpStep);
    vecSumDiffP21.push_back(logpStep2);
  }

  // do another scan to find the best
  double res = -1.0 * HAP_MAX_INT;
  int node1 = -1, node2 = -1;
  for (int node = 0; node < mtree.GetTotNodesNum(); ++node) {
    double p2Part = 0.0;
    double node2MaxUse = -1;
    if (vecMaxSumDiff21[node] > 0.0) {
      p2Part = vecMaxSumDiff21[node];
      node2MaxUse = vecMaxSumDiff21Node[node];
    }
    if (vecSumDiffP10[node] + p2Part > res) {
      res = vecSumDiffP10[node] + p2Part;

      node1 = node;
      node2 = node2MaxUse;
    }
  }

  // figure out the genos
  set<int> dummy;
  ScistPerfPhyCluster clusDummy(dummy);
  if (res < 0.0) {
    //
    clusChangedHetero = clusDummy;
    clusChangedHomo = clusDummy;
  } else {
    YW_ASSERT_INFO(node1 >= 0, "Wrong");
    set<int> nodeOptSplit, nodeOptSplitLbls;
    mtree.GetLeavesUnder(node1, nodeOptSplit);
    mtree.GetlabelsFor(nodeOptSplit, nodeOptSplitLbls);
    DecAllNumInSet(nodeOptSplitLbls);
    set<int> nodeOptSplitLbls2;
    if (node2 >= 0) {
      set<int> nodeOptSplit2;
      mtree.GetLeavesUnder(node2, nodeOptSplit2);
      mtree.GetlabelsFor(nodeOptSplit2, nodeOptSplitLbls2);
      DecAllNumInSet(nodeOptSplitLbls2);
    }
    SubtractSets(nodeOptSplitLbls, nodeOptSplitLbls2);

    ScistPerfPhyCluster clus1(nodeOptSplitLbls);
    clusChangedHetero = clus1;
    ScistPerfPhyCluster clus2(nodeOptSplitLbls2);
    clusChangedHomo = clus2;
  }

  return res + this->listSitePriorScore[site];
}

double ScistPerfPhyProbOnTree ::CalcProbForSite(
    int site, double totEdgeLen, const vector<set<int> > &listClades) const {
  ScistHaplotypeMat *pHapMat =
      dynamic_cast<ScistHaplotypeMat *>(&this->genosInput);

  if (pHapMat != NULL) {
    return CalcProbForSiteHap(site, totEdgeLen, listClades);
  } else // right now, must be of genotype
  {
    return CalcProbForSiteGeno(site, totEdgeLen, listClades);
  }
}

double ScistPerfPhyProbOnTree ::CalcProbForSiteHap(
    int site, double totEdgeLen, const vector<set<int> > &listClades) const {
  vector<double> listCladeProb;
  for (int i = 0; i < mtree.GetTotNodesNum(); ++i) {
    listCladeProb.push_back(-1.0 * HAP_MAX_INT);
  }

  // get the sum of prob0
  double sumProb0 = 0.0;
  for (int h = 0; h < this->genosInput.GetNumHaps(); ++h) {
    sumProb0 += log(this->genosInput.GetGenotypeProbAllele0At(h, site));
  }

  double loglikeTot = -1.0 * HAP_MAX_INT;
  for (int i = 0; i < mtree.GetTotNodesNum(); ++i) {
    if (i == mtree.GetRoot()) {
      continue;
    }
    double brLen = mtree.GetEdgeLen(i);
    double probPrior = brLen / totEdgeLen;
    double probCladeOnly = 0.0;
    if (mtree.IsLeaf(i)) {
      int lbl = *listClades[i].begin();
      double p0 = this->genosInput.GetGenotypeProbAllele0At(lbl, site);
      double p1 = 1 - p0;
      double pr = log(p1 / p0);
      probCladeOnly = pr;
    } else {
      int childLeft = mtree.GetLeftDescendant(i);
      int childRight = mtree.GetRightDescendant(i);
      probCladeOnly = listCladeProb[childLeft] + listCladeProb[childRight];
    }
    // cout << "probPrior: " << probPrior << endl;
    listCladeProb[i] = probCladeOnly + log(probPrior);
    // double probClade = CalcProbMutClade(site, listClades[i] );
    // YW: need to check this
    // loglikeTot = GetLogSumOfTwo(loglikeTot, log(probPrior) + probCladeOnly);
    if (loglikeTot < listCladeProb[i]) {
      loglikeTot = listCladeProb[i];
    }
  }
  // return loglikeTot + sumProb0;
  double res = loglikeTot + sumProb0;
  // cout << "log prob at site: " << res << endl;
  return res;
}

double ScistPerfPhyProbOnTree ::CalcProbForSiteGeno(
    int site, double totEdgeLen, const vector<set<int> > &listClades) const {
  ScistPerfPhyProbOnTree *pthis = const_cast<ScistPerfPhyProbOnTree *>(this);
  ScistPerfPhyProbOnTree spppt(pthis->genosInputHap, this->mtree);
  return spppt.CalcProbForSite(2 * site, totEdgeLen, listClades) +
         spppt.CalcProbForSite(2 * site + 1, totEdgeLen, listClades);
}
