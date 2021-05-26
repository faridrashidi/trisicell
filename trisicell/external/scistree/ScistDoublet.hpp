//
//  ScistDoublet.hpp
//
//
//  Created by Yufeng Wu on 6/2/18.
//
//

#ifndef ScistDoublet_hpp
#define ScistDoublet_hpp

#include <map>
#include <set>
#include <string>
#include <vector>

class ScistGenGenotypeMat;
class ScistPerfPhyCluster;
class ScistPerfPhyClusTreeNode;

// *************************************************************************************
// DP backtrace info

class ScistDoubletDPTraceback {
public:
  ScistDoubletDPTraceback();
  ScistDoubletDPTraceback(const ScistDoubletDPTraceback &rhs);
  ScistDoubletDPTraceback &operator=(const ScistDoubletDPTraceback &rhs);

  void AddTraceback(int indexChild, int phase);
  int GetChild1() const { return indexChild1; }
  int GetPhase1() const { return phaseChild1; }
  int GetChild2() const { return indexChild2; }
  int GetPhase2() const { return phaseChild2; }
  void SetChild1(int c) { indexChild1 = c; }
  void SetPhase1(int p) { phaseChild1 = p; }
  void SetChild2(int c) { indexChild2 = c; }
  void SetPhase2(int p) { phaseChild2 = p; }

private:
  int indexChild1;
  int phaseChild1;
  int indexChild2;
  int phaseChild2;
};

// *************************************************************************************
// Deal with doublet (single genotype row)

class ScistDoublet {
public:
  ScistDoublet(const ScistGenGenotypeMat &genosInputIn);
  double EvalGenoDoublet(const std::set<int> &setTemplateRows, int genoDoublet,
                         std::vector<int> &genoDoublePhase1,
                         std::vector<int> &genoDoublePhase2) const;

private:
  void ConsDPTblDoubletNodes(
      const std::map<int, ScistPerfPhyCluster> &setTemplateSites,
      const std::map<const ScistPerfPhyCluster *, int> &mapClusToSiteIndex,
      int genoDoublet, ScistPerfPhyClusTreeNode *pNodeCurr,
      std::map<ScistPerfPhyClusTreeNode *,
               std::vector<std::pair<double, ScistDoubletDPTraceback> > >
          &mapNodeVals) const;
  void ConsClustersForTemplates(
      const std::set<int> &setTemplateRows,
      std::map<int, ScistPerfPhyCluster> &setTemplateSites,
      std::map<const ScistPerfPhyCluster *, int> &mapClusToSiteIndex) const;
  void ConsPhasing(
      const std::map<const ScistPerfPhyCluster *, int> &mapClusToSiteIndex,
      int genoDoublet, ScistPerfPhyClusTreeNode *pNodeRoot,
      const std::map<ScistPerfPhyClusTreeNode *,
                     std::vector<std::pair<double, ScistDoubletDPTraceback> > >
          &mapNodeVals,
      std::vector<int> &vecPhasing) const;
  void TracePhasingAtNode(
      const std::map<const ScistPerfPhyCluster *, int> &mapClusToSiteIndex,
      int genoDoublet, ScistPerfPhyClusTreeNode *pNodeCurr, int phasingCurr,
      const std::map<ScistPerfPhyClusTreeNode *,
                     std::vector<std::pair<double, ScistDoubletDPTraceback> > >
          &mapNodeVals,
      std::vector<int> &vecPhasing) const;
  void ConsPhasingVec(const std::vector<int> &vecPhasing,
                      std::vector<int> &genoDoublePhase1,
                      std::vector<int> &genoDoublePhase2) const;

  const ScistGenGenotypeMat &genosInput;
};

// *************************************************************************************
// Deal with doublet (search)

class ScistDoubletSearch {
public:
  ScistDoubletSearch(const ScistGenGenotypeMat &genosInputIn,
                     int maxDoubletSubsetSzIn);
  void Search();
  void SearchInc();
  void SetDouletCost(double c) { costDoublet = c; }
  void SetVerbose(bool f) { fVerbose = f; }
  void SetMutTreeOut(bool f) { fOutputPPWithEdgeLabels = f; }
  void SetCellNames(const std::vector<std::string> &listCellNamesIn) {
    listCellNames = listCellNamesIn;
  }
  void SetSiteNames(const std::vector<std::string> &listSiteNamesIn) {
    listSiteNames = listSiteNamesIn;
  }
  void SetMutTreeFileName(const std::string &strMutTreeFileNameIn) {
    this->strMutTreeFileName = strMutTreeFileNameIn;
  }
  static void GetNgbrTreesFrom(int numHaps, const std::string &strTree,
                               std::set<std::string> &setNgbrTrees);

private:
  ScistGenGenotypeMat *
  EvalGenoDoubletSet(const ScistGenGenotypeMat &matToSearch,
                     const std::set<int> &setDoubletRows, double &optCost);
  double FitPerfPhyFor(ScistGenGenotypeMat *pMatCurr,
                       const std::set<int> &setTemplateRows);
  double ScoreDoubletRow(ScistGenGenotypeMat *pMatCurr,
                         const std::set<int> &rowsTemplate, int rowDouble,
                         std::vector<int> &genoDoublePhase1,
                         std::vector<int> &genoDoublePhase2);
  void FindDoubletCandidates(std::set<int> &candidatesDoublet);
  ScistGenGenotypeMat *
  InitSearchGenotypes(const ScistGenGenotypeMat &matToSearch,
                      const std::set<int> &candidatesDoublet,
                      std::set<int> &setDoubletRows, double &costInit);
  void UpdateSearchGenotypes(ScistGenGenotypeMat *pMatCurr, int genoDoublet,
                             const std::vector<int> &genoDoublePhase1,
                             const std::vector<int> &genoDoublePhase2);
  void FindOrigImputedGeno(const ScistGenGenotypeMat &genosPhasingRes,
                           ScistGenGenotypeMat &genosImpute,
                           std::map<int, std::set<int> > &mapDoublets) const;
  std::string GetGenoDoubleRowName(const std::string &strName) const;
  std::string GetNewGenoDoubleRowName(const ScistGenGenotypeMat &matToSearch,
                                      int index) const;
  ScistGenGenotypeMat *
  CreateGnoesWithDouble(const ScistGenGenotypeMat &genosOrig, int indexDobule,
                        const ScistGenGenotypeMat &genosDoubleInfer) const;
  double ConsTree(ScistGenGenotypeMat &genosNoDoublets,
                  std::string &strTreeNW) const;
  void OutputMutTree(ScistGenGenotypeMat &genosNoDoublets) const;
  std::string ConvMutTreeStr(const std::string &strTree) const;
  void FindDoubletHapsInMat(const ScistGenGenotypeMat &genosDbl,
                            std::set<int> &setHapsDoubles) const;
  bool IsOverImpute(const ScistGenGenotypeMat &genosDbl) const;

  const ScistGenGenotypeMat &genosInput;
  int maxDoubletSubsetSz;
  double costDoublet;
  bool fVerbose;
  bool fOutputPPWithEdgeLabels;
  std::vector<std::string> listCellNames;
  std::vector<std::string> listSiteNames;
  std::string strMutTreeFileName;
};

// *************************************************************************************
void ScistDoubletTest();

#endif /* ScistDoublet_hpp */
