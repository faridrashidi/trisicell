#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <vector>

using namespace std;

#include "ScistDoublet.hpp"
#include "ScistErrRateInf.hpp"
#include "ScistGenotype.hpp"
#include "ScistPerfPhyImp.hpp"
#include "ScistPerfPhyUtils.hpp"
#include "Utils2.h"
#include "Utils3.h"

//*****************************************************************************
// Main driving functions
// ***************************************************************************

// ***************************************************************************
// Main for computing lower bound
// ***************************************************************************
static void Usage() {
  cout << "Usage: ./scistree <options> <input file> " << endl;
  cout << "Options:\n";
  // cout << "\t -d <dn>           dn: number of allowed doublet genotypes; dc:
  // cost of having a doublet\n";
  cout << "\t -d <dn>           dn: number of allowed doublet genotypes\n";
  cout << "\t -v                Turn on verbose mode  \n";
  // cout << "\t -p                Find optimal false positive rate and false
  // negative rate\n"; cout << "\t -l                Find cell tree with branch
  // length (by default, constructed cell trees don't have branch length\n";
  cout << "\t -n                Only build simple neighbor joining tree (may "
          "be useful for very large data)\n";
  cout << "\t -e                Output mutation tree (may not be binary tree) "
          "from called genotypes branch labels.\n";
  cout << "\t -e0               Output mutation tree but don't output labels "
          "(for visualizing large trees).\n";
  // cout << "\t -s <level>        Use SPR tree search (this will be slower);
  // level: # of SPRs to allow (default is 1)\n";
  cout << "\t -o <output-file>  Set output file (used for mutation tree output "
          "(in GML) format; should have suffix .gml (default: "
          "mutation-tree.gml)\n";
  cout << "\t -t <threshold>    Discard somewhat ambigous genotyeps when "
          "constructing intial trees: \n\t\t\t genotypes discarded if the "
          "prob. of alternative genotypes is less than <threshold> "
          "\n\t\t\t(default is 0, i.e. use all genotypes)\n";
  exit(1);
}

// settings
static int fileInArgIndex = 1;
static int numDoublets = 0;
static double costDoublet = 0.0;
static bool fVerbose = false;
static bool fOptParam = false;
static bool fOptBrLen = false;
static bool fNJOnly = false;
static bool fSPR = false;
static int numSPR = 1;
static double thresProbSignificance = 0.0;
static vector<string> listCellNames;
static vector<string> listSiteNames;
static int numSites = 0;
static int numSCs = 0;
static string strMutTreeOutFile = "mutation-tree.gml";
static bool fOutPPEdgeLabel = false;
static bool fOutputLabel = true;
// GLobal variables

// Local functions
static bool CheckArguments(int argc, char **argv) {
  if (argc <= 1) {
    return false;
  }

  // Check argument one by one
  // int argpos = 1;
  for (int i = 1; i < argc; ++i) {
    if (argv[i][0] == '-' && argv[i][1] == 'l') {
      YW_ASSERT_INFO(i < argc - 1, "Check input");
      fOptBrLen = true;
      cout << "Turn on branch optimization. " << endl;
    } else if (argv[i][0] == '-' && argv[i][1] == 'd') {
      YW_ASSERT_INFO(i < argc - 1, "Check input");
      ++i;
      sscanf(argv[i], "%d", &numDoublets);
      // YW_ASSERT_INFO( i <argc-1, "Check input" );
      //++i;
      // float costDoubletThis;
      // sscanf(argv[i], "%f", &costDoubletThis);
      // costDoublet = costDoubletThis;
      // cout << "Setting doublet number to " << numDoublets << ", and doublet
      // cost to " << costDoublet << endl;
      cout << "Setting doublet number to " << numDoublets << endl;
    } else if (argv[i][0] == '-' && argv[i][1] == 'v') {
      YW_ASSERT_INFO(i < argc - 1, "Check input");
      fVerbose = true;
      cout << "Turn on verbose mode" << endl;
    } else if (argv[i][0] == '-' && argv[i][1] == 'n') {
      YW_ASSERT_INFO(i < argc - 1, "Check input");
      fNJOnly = true;
      cout << "Only build neighbor joining tree." << endl;
    } else if (argv[i][0] == '-' && argv[i][1] == 'p') {
      YW_ASSERT_INFO(i < argc - 1, "Check input");
      fOptParam = true;
      cout << "Search for optimal genotype error rates" << endl;
    } else if (argv[i][0] == '-' && argv[i][1] == 'e') {
      YW_ASSERT_INFO(i < argc - 1, "Check input");
      fOutPPEdgeLabel = true;
      cout << "Output perfect phylogeny with edge labels" << endl;

      string strOpt = argv[i];
      if (strOpt.length() >= 3 && strOpt[2] == '0') {
        cout << "  -- no labels in mutation tree\n";
        fOutputLabel = false;
      }
    } else if (argv[i][0] == '-' && argv[i][1] == 's') {
      YW_ASSERT_INFO(i < argc - 1, "Check input");
      fSPR = true;
      ++i;
      sscanf(argv[i], "%d", &numSPR);
      cout << "Use SPR tree search: level set to " << numSPR << endl;
    } else if (argv[i][0] == '-' && argv[i][1] == 't') {
      YW_ASSERT_INFO(i < argc - 1, "Check input");
      ++i;
      float thresUse = 0.0;
      sscanf(argv[i], "%f", &thresUse);
      thresProbSignificance = thresUse;
      cout << "Threshold for probability significance: set to "
           << thresProbSignificance << endl;
    } else if (argv[i][0] == '-' && argv[i][1] == 'o') {
      YW_ASSERT_INFO(i < argc - 1, "Check input");
      ++i;
      strMutTreeOutFile = argv[i];
      cout << "Use mutation tree file name to " << strMutTreeOutFile << endl;
    }

    else if (argv[i][0] != '-') {
      // not an option one. Right now the only one is file
      fileInArgIndex = i;
      // filenameGMLPrefix = argv[i];
    } else {
      return false;
    }
  }

  return true;
}

// input handling
static ScistGenGenotypeMat *ReadsInput(const char *filename) {
  //
  ifstream inFile(filename);
  if (!inFile) {
    cout << "Can not open " << filename << endl;
    YW_ASSERT_INFO(false, "Stop");
  }
  ScistGenGenotypeMat *pMatIn = NULL;
  while (inFile.eof() == false) {
    const int BUF_SZ = 102400;
    char buffer[BUF_SZ];
    inFile.getline(buffer, BUF_SZ);
    if (strlen(buffer) > 0) {
      // cout << "read one line: " << buffer << endl;
      // now try to read alleles
      std::istringstream is(buffer);

      // looking for keyword
      string strKey;
      is >> strKey;
      if (strKey == "HAPLOTYPES" || strKey == "HAPLOID") {
        is >> numSites >> numSCs;
        // cout << "numSites: " << numSites << ", numSCs: " << numSCs << endl;
        YW_ASSERT_INFO(numSites > 0 && numSCs > 0,
                       "Site and single cells numbers: Cannot be zeros");

        // read in names if specified
        while (is.eof() == false) {
          string strName;
          is >> strName;
          if (strName.length() > 0) {
            listCellNames.push_back(strName);
            // cout << "One lineage name: " << strName << endl;
          }
          if ((int)listCellNames.size() > numSCs) {
            break;
          }
        }
        if (listCellNames.size() > 0 && (int)listCellNames.size() != numSCs) {
          YW_ASSERT_INFO(
              false, "Fatal error: you must provide names for each lineage");
        }
        bool fSiteName = false;
        if (listCellNames.size() > 0) {
          fSiteName = true;
          ;
        }

        pMatIn = new ScistHaplotypeMat;
        for (int i = 0; i < (int)listCellNames.size(); ++i) {
          pMatIn->AddGenotypeName(listCellNames[i]);
        }

        pMatIn->ReadFromFile(inFile, numSites, numSCs, fSiteName);

#if 0
if( fSiteName )
{
cout << "List of site names: ";
for(int i=0; i<numSites; ++i)
{
cout << pMatIn->GetSiteName(i) << " ";
}
cout << endl;
}
#endif

        break;
      } else if (strKey == "TERNARY") {
        is >> numSites >> numSCs;
        // cout << "numSites: " << numSites << ", numSCs: " << numSCs << endl;
        YW_ASSERT_INFO(numSites > 0 && numSCs > 0,
                       "Site and single cells numbers: Cannot be zeros");

        // read in names if specified
        while (is.eof() == false) {
          string strName;
          is >> strName;
          if (strName.length() > 0) {
            listCellNames.push_back(strName);
            // cout << "One lineage name: " << strName << endl;
          }
          if ((int)listCellNames.size() > numSCs) {
            break;
          }
        }
        if (listCellNames.size() > 0 && (int)listCellNames.size() != numSCs) {
          YW_ASSERT_INFO(
              false, "Fatal error: you must provide names for each lineage");
        }
        bool fSiteName = false;
        if (listCellNames.size() > 0) {
          fSiteName = true;
          ;
        }

        pMatIn = new ScistTernaryMat;
        for (int i = 0; i < (int)listCellNames.size(); ++i) {
          pMatIn->AddGenotypeName(listCellNames[i]);
        }

        pMatIn->ReadFromFile(inFile, numSites, numSCs, fSiteName);

        break;
      }
    }
  }
  pMatIn->SetSignificantThres(thresProbSignificance);

  // initialize cell names to plain 1, 2, ... if not specified
  if (listCellNames.size() == 0) {
    YW_ASSERT_INFO(numSCs > 0, "Number of SCs: not intiialized");
    for (int c = 1; c <= numSCs; ++c) {
      string str = std::to_string(c);
      listCellNames.push_back(str);
    }
  }
  pMatIn->GetSiteNamesAll(listSiteNames);

  return pMatIn;
}

// test code
static void TestCode(const char *filename) {
  //

  ScistGenGenotypeMat *pMatInput = ReadsInput(filename);
  string filenameUse = filename;
  pMatInput->SetFileName(filenameUse);

  // cout << "Input genotype matrix:\n";
  // pMatInput->Dump();
  // string strNJ2 = pMatInput->ConsNJTree();
  // cout << "NJ tree: " << strNJ2 << endl;
  // delete pMatInput;
  // exit(1);

  if (fOptParam) {
    cout << "Now searching for optimal genotype error rates...\n";
    ScistErrRateInf serInf(*pMatInput);
    serInf.SetVerbose(fVerbose);
    serInf.Infer();
  } else {
    string treeNJ = pMatInput->ConsNJTree();
    if (fVerbose) {
      cout << "Neighbor joining tree from noisy genotypes: " << treeNJ << endl;
    }
    if (fNJOnly) {
      delete pMatInput;
      return;
    }

    // ScistInfPerfPhyTest();
    // plain mode if no double is allowed
    if (numDoublets == 0) {
#if 0
            ScistFullPerfPhyMLE ppInfHeu(*pMatInput);
            ppInfHeu.SetVerbose(fVerbose);
            ppInfHeu.Infer();
#endif

      ScistPerfPhyMLE ppInfHeu(*pMatInput);
      ppInfHeu.SetBrOpt(fOptBrLen);
      ppInfHeu.SetVerbose(fVerbose);
      ppInfHeu.SetPPOut(fOutPPEdgeLabel);
      ppInfHeu.SetPPOutLabel(fOutputLabel);
      ppInfHeu.SetSPR(fSPR);
      ppInfHeu.SetSPRNum(numSPR);
      ppInfHeu.SetCellNames(listCellNames);
      ppInfHeu.SetSiteNames(listSiteNames);
      ppInfHeu.SetMutTreeFileName(strMutTreeOutFile);
      ppInfHeu.Infer();
    } else {
      // right now only work with haplotype matrix
      ScistHaplotypeMat *pMatInputUse =
          dynamic_cast<ScistHaplotypeMat *>(pMatInput);
      YW_ASSERT_INFO(
          pMatInputUse != NULL,
          "At present, doublet feature only works for binary genotype matrix.");

      cout << "SEARCHING FOR DOUBLETS...\n";
      ScistDoubletSearch sds(*pMatInput, numDoublets);
      sds.SetVerbose(fVerbose);
      sds.SetDouletCost(costDoublet);
      sds.SetMutTreeOut(fOutPPEdgeLabel);
      sds.SetCellNames(listCellNames);
      sds.SetSiteNames(listSiteNames);
      sds.SetMutTreeFileName(strMutTreeOutFile);
      sds.SearchInc();
    }
  }

  delete pMatInput;
}

////////////////////////////////////////////////////////////////////////////////////////

const char *CODE_VER_INFO = "*** SCISTREE ver. 1.2.0.6, May 19, 2019 ***";

//******************************************************************
int main_in_c(int argc, char **argv) {
  //    int seq = 0x001;
  //    int seqMut;
  //    MutateHCSeqAt(seq, seqMut, 4, 2);
  // cout << "mutated seq = " << seqMut << endl;

  string outputfile = argv[argc - 1];
  string str2 = "scistree.input";
  string str3 = "scistree.output";
  outputfile.replace(outputfile.find(str2), str2.length(), str3);

  std::ofstream out(outputfile);
  std::streambuf *coutbuf = std::cout.rdbuf(); // save old buf
  std::cout.rdbuf(out.rdbuf());                // redirect std::cout to out.txt!

  cout << CODE_VER_INFO << endl << endl;

  // first verify usage
  if (CheckArguments(argc, argv) == false) {
    Usage();
  }

  // cout << "here0\n";
  long tstart1 = GetCurrentTimeTick();

  TestCode(argv[fileInArgIndex]);

  cout << "Elapsed time = " << GetElapseTime(tstart1) << " seconds." << endl;

  // dump out stats
  // ApproxGTPStats::Instance().DumpStats();

  std::cout.rdbuf(coutbuf); // reset to standard output again

  return 0;
}
