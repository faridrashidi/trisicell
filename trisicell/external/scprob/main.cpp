// Not exactly a script, but I have to do it this way since AWK does not work
// May 12, 2019: change genotype probability format to Scistree v.1.2.0
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <cmath>
#include <map>
#include <set>
#include <queue>
#include <random>
#include <string.h>
using namespace std;


// parameters
static double rateDropout = 0.2;
static double rateDropoutstd = 0.0;
static double rateErr = 0.002;
static double fracMissing = 0.0;
static bool seedRndUserDef = true;
static double minAlleFreq = 0.0;
const double SMALL_NUM = 0.00001;
const double SMALL_PROB = 0.000000001;
const double LARGE_PROB = 0.999999999;
bool fBinary;
// the following refers to the single-cell sequencing realted aspects
static double aveReadDepth = 4.0;
static double stdReadDepth = 2.0;
static double aveStrandBias = 0.0;
static double stdStrandBias = 0.0;
static vector<double> listCellDropoutRates;
static bool fNoReadsAsMiss = true;
static double thresMinFracDiffNoMiss = 0.8;

//*******************************************************************************************************
// Utilies

static int QSortCompareInt( const void *arg1, const void *arg2 )
{
    /* Compare all of both strings: */
    // assume sorting in accending order
    int n1 = *((int *) arg1);
    int n2 = *((int *) arg2);
    //cout <<"arg1 = " << n1 << ", arg2 = " << n2 << endl;
    if( n1 > n2)
    {
        return 1;
    }
    else if( n1 < n2)
    {
        return -1;
    }
    else
    {
        return 0;
    }
}
void SortIntVec( vector<int> &vecVals, int start, int end )
{
    //#if 0
    if( vecVals.size() <= 1)
    {
        // do nothing
        return;
    }
    if (end < 0 )
    {
        end = vecVals.size() - 1;
    }
    int sortLen = end - start +1;
    int *array = new int[sortLen];
    for(int i=start; i<= end; ++i)
    {
        array[i-start] = vecVals[i];
    }
    qsort( (void *)array, sortLen, sizeof( int ), QSortCompareInt );
    // Now write back
    for(int i=start; i<=end; ++i)
    {
        vecVals[i] = array[i-start];
    }
    delete [] array;
}
static int FindTwoClusterBoundary( const vector<int> &listVals )
{
    // goal: find the cutoff point s.t. the two parts are as similar as possible
    // return the largest value in part1 (i.e. anything <= this value is the second part)
    vector<int> listValsUse = listVals;
    SortIntVec(listValsUse, 0, listValsUse.size()-1);
    //
    int res = -1;
    double minCenterSum = 1000000000000000000.0;
    for(int i=1; i<(int)listValsUse.size()-1; ++i)
    {
        // compute the partition
        int szPart1 = i, szPart2 = (int)listValsUse.size()-szPart1;
        if(szPart2 == 0)
        {
            szPart2 = 1;
        }
        int sumPart1 = 0, sumPart2 = 0;
        for(int j=0; j<i; ++j)
        {
            sumPart1 += listValsUse[j];
        }
        for(int j=i; j<(int)listValsUse.size(); ++j)
        {
            sumPart2 += listValsUse[j];
        }
        double partCenter1 = ((double)sumPart1/szPart1);
        double partCenter2 = ((double)sumPart2/szPart2);
        double sumSqDiff1 = 0.0, sumSqDiff2 = 0.0;
        for(int j=0; j<i; ++j)
        {
            double diff = listValsUse[j] - partCenter1;
            sumSqDiff1 += diff*diff;
        }
        for(int j=i; j<(int)listValsUse.size(); ++j)
        {
            double diff = listValsUse[j] - partCenter2;
            sumSqDiff2 += diff*diff;
        }
        double meanSqDiff1 = sumSqDiff1/szPart1;
        double meanSqDiff2 = sumSqDiff2/szPart2;
        double meanSqDiffSum = meanSqDiff1 + meanSqDiff2;
        if( meanSqDiffSum < minCenterSum )
        {
            minCenterSum = meanSqDiffSum;
            res = listValsUse[i-1];
        }
    }
    return res;
}


//******************************************************************************************************
// Prob computation

static void MakeProbInRange(double &prob)
{
    if(prob < SMALL_PROB)
    {
        prob = SMALL_PROB;
    }
    if(prob> LARGE_PROB)
    {
        prob = LARGE_PROB;
    }
}

static double CalcProb0(int cell)
{
    // allele 0: may be a result of dropout or genotype error
    double rateDropoutCell = listCellDropoutRates[cell];
    double res = 1.0-rateDropoutCell-rateErr;
    MakeProbInRange(res);
    return res;
}
static double CalcProb1()
{
    // allele 1: error is only cause
    return 1.0-rateErr;
}


//******************************************************************************************************
// sequence reads sampling from genotypes


static double CalcReadAlleleForGeno(int allele, int g)
{
    // allele = 0 or 1, g=0,1,2
    // use errRate as the reads error rate, also drop rate matters too
    if( g == 0 )
    {
        if( allele == 0 )
        {
            return 1.0 - rateErr;
        }
        else
        {
            return rateErr;
        }
    }
    else if(g == 2)
    {
        if( allele == 1 )
        {
            return 1.0 - rateErr;
        }
        else
        {
            return rateErr;
        }
    }
    else
    {
        return 0.5;
    }
}

static double CalcReadAlleleForGenoWithOneDrop0(int n0, int n1)
{
    double p0 = 1.0 - rateErr;
    double p1 = rateErr;
    return pow(p0, n0) * pow(p1, n1);
}
static double CalcReadAlleleForGenoWithOneDrop2(int n0, int n1)
{
    double p1 = 1.0 - rateErr;
    double p0 = rateErr;
    return pow(p0, n0) * pow(p1, n1);
}

static double CalcReadAlleleForGenoWithOneDrop(int g, int n0, int n1)
{
    // droput can occur at either allele
    if( g== 0 )
    {
        return CalcReadAlleleForGenoWithOneDrop0(n0, n1);
    }
    else if(g == 2)
    {
        return CalcReadAlleleForGenoWithOneDrop2(n0, n1);
    }
    else
    {
        // assume dropout occurs w/ 50/50 chance at one allele
        double p00 = CalcReadAlleleForGenoWithOneDrop0(n0, n1);
        double p11 = CalcReadAlleleForGenoWithOneDrop2(n0, n1);
        return (p00+p11)*0.5;
    }
}

static void CalcGenotypeProbOfReads(const pair<int,int> &readsCount01, double allele0Freq, vector<double> &listGenoProbs)
{
//#if 0
    //double p0 = rateDropout + (1.0-rateDropout) * allele0Freq*allele0Freq, p1 = (1.0-rateDropout) * 2*allele0Freq*(1.0-allele0Freq), p2 = (1.0-rateDropout) * (1.0-allele0Freq)*(1.0-allele0Freq);
    double p00 = allele0Freq*allele0Freq, p10 = 2*allele0Freq*(1.0-allele0Freq), p20 = (1.0-allele0Freq)*(1.0-allele0Freq);
    double p0 = 1.0, p1 = 1.0, p2 = 1.0;
    double pa0g0 = CalcReadAlleleForGeno(0, 0);
    double pa0g1 = CalcReadAlleleForGeno(0, 1);
    double pa0g2 = CalcReadAlleleForGeno(0, 2);
    double pa1g0 = CalcReadAlleleForGeno(1, 0);
    double pa1g1 = CalcReadAlleleForGeno(1, 1);
    double pa1g2 = CalcReadAlleleForGeno(1, 2);
//cout << "p00:" << p00 << ", p10:" << p10 << ",p20:" << p20 << ",pa0g0: " << pa0g0 << ", pa0g1: " << pa0g1 << ", pa0g2: " << pa0g2 << ", pa1g0:" << pa1g0 << ", pa1g1:" << pa1g1 << ", pa1g2:" << pa1g2 << endl;
    p0 = pow(pa0g0, readsCount01.first) * pow(pa1g0, readsCount01.second) * p00;
    p1 = pow(pa0g1, readsCount01.first) * pow(pa1g1, readsCount01.second) * p10;
    p2 = pow(pa0g2, readsCount01.first) * pow(pa1g2, readsCount01.second) * p20;
//cout << "Before normalize: p0: " << p0 << ", p1: " << p1 << ", p2: " << p2 << endl;
    double sum = p0+p1+p2;
    p0 = p0/sum;
    p1 = p1/sum;
    p2 = p2/sum;

#if 0
    // add droput
    double p0all = rateDropout + (1.0-rateDropout)*p0;
    double p1all = (1.0-rateDropout)*p1;
    double p2all = (1.0-rateDropout)*p2;
#endif

    double p0d = CalcReadAlleleForGenoWithOneDrop(0, readsCount01.first, readsCount01.second);
    double p1d = CalcReadAlleleForGenoWithOneDrop(1, readsCount01.first, readsCount01.second);
    double p2d = CalcReadAlleleForGenoWithOneDrop(2, readsCount01.first, readsCount01.second);
    double p0all = rateDropout*p0d + (1.0-rateDropout)*p0;
    double p1all = rateDropout*p1d + (1.0-rateDropout)*p1;
    double p2all = rateDropout*p2d + (1.0-rateDropout)*p2;

    listGenoProbs.clear();
    listGenoProbs.push_back(p0all);
    listGenoProbs.push_back(p1all);
    listGenoProbs.push_back(p2all);

    // make sure the prob is not too small
    MakeProbInRange( listGenoProbs[0] );
    MakeProbInRange( listGenoProbs[1] );
    MakeProbInRange( listGenoProbs[2] );
//#endif

//cout <<"List of simulated genotype probability (after normalize): " << listGenoProbs[0] << "," << listGenoProbs[1] << "," << listGenoProbs[2] << endl;
}


static double CalcProbReadDepthNormDist(double numReads, double aveReadDepthSite, double stdReadDepthSite )
{
    // assume normal distribution
    double diff = numReads - aveReadDepthSite;
    return exp( -1.0*diff*diff/(2.0*stdReadDepthSite*stdReadDepthSite) );
}
static double CalcGenotypeProbForGenoSingleStrand( const pair<int,int> &readsCount01,  double aveReadDepthSite, double stdReadDepthSite, int geno )
{

    // single strand: no droput,
    // if the geno is different from the read, must be errors
    double probErr = 1.0;
    if( geno == 0 )
    {
        probErr = pow(rateErr, readsCount01.second);
    }
    else
    {
        probErr = pow(rateErr, readsCount01.first);
    }
    int readNum = readsCount01.first+readsCount01.second;
    double probDepth = CalcProbReadDepthNormDist(readNum, aveReadDepthSite, stdReadDepthSite);
    double res = probDepth * probErr;
//cout << "readNum: " << readNum << ", probDepth: " << probDepth << ", probErr: " << probErr << endl;
//cout << "CalcGenotypeProbForGenoSingleStrand: read: " << readsCount01.first << "," << readsCount01.second << ", geno: " << geno << ", aveRead: " << aveReadDepthSite << ", stdRead: " << stdReadDepthSite << ", prob=" << res << endl;
    return res;
}
static double CalcGenotypeProbForGenoDrop( const pair<int,int> &readsCount01,  double aveReadDepthSite, double stdReadDepthSite, int geno0, int geno1, bool fdrop0, bool fdrop1 )
{
    if( fdrop0 == true && fdrop1 == true )
    {
        // all errors
        double res = pow(rateErr, readsCount01.first+readsCount01.second);
        return res;
    }
    else if( fdrop0 == true && fdrop1 == false )
    {
        return CalcGenotypeProbForGenoSingleStrand( readsCount01, aveReadDepthSite, stdReadDepthSite, geno1 );
    }
    else if( fdrop0 == false && fdrop1 == true )
    {
        return CalcGenotypeProbForGenoSingleStrand( readsCount01, aveReadDepthSite, stdReadDepthSite, geno0 );;
    }
    else
    {
        // consider cases: if homozygote, then take half
        if( geno0 == geno1)
        {
            pair<int,int> pp1(readsCount01.first/2, readsCount01.second/2);
            pair<int,int> pp2(readsCount01.first-pp1.first, readsCount01.second -pp1.second);
            double p1=CalcGenotypeProbForGenoSingleStrand( pp1, aveReadDepthSite, stdReadDepthSite, geno0 );
            double p2=CalcGenotypeProbForGenoSingleStrand( pp2, aveReadDepthSite, stdReadDepthSite, geno1 );
            return p1*p2;
        }
        else
        {
            // heterozygote: just use the read count
            pair<int,int> pp1(readsCount01.first, 0), pp2(0, readsCount01.second);
            double p1=CalcGenotypeProbForGenoSingleStrand( pp1, aveReadDepthSite, stdReadDepthSite, 0 );
            double p2=CalcGenotypeProbForGenoSingleStrand( pp2, aveReadDepthSite, stdReadDepthSite, 1 );
            return p1*p2;
        }
    }
}

static void CalcGenotypeProbOfReadsNew( int cell, double aveReadDepthSite, double stdReadDepthSite, const pair<int,int> &readsCount01, double allele0Freq, int numInfDrops, vector<double> &listGenoProbs)
{
//cout << "CalcGenotypeProbOfReadsNew: cell: " << cell << ", aveReadDepthSite: " << aveReadDepthSite << ", stdReadDepthSite: " << stdReadDepthSite << ", read count: " << readsCount01.first << ", " << readsCount01.second  << ", allele0Freq: " << allele0Freq << ", numInfDrops: " << numInfDrops << endl;

//#if 0
    // assume normal distribution of read counts; calculate prob of 0
    double p00d00 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 0, false, false );
    double p00d01 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 0, false, true );
    double p00d10 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 0, true, false );
    double p00d11 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 0, true, true );
    double p01d00 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 1, false, false );
    double p01d01 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 1, false, true );
    double p01d10 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 1, true, false );
    double p01d11 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 0, 1, true, true );
    double p11d00 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 1, 1, false, false );
    double p11d01 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 1, 1, false, true );
    double p11d10 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 1, 1, true, false );
    double p11d11 = CalcGenotypeProbForGenoDrop( readsCount01, aveReadDepthSite, stdReadDepthSite, 1, 1, true, true );
//cout << "p00d00:" << p00d00 << ", p00d01:" << p00d01 << ", p00d10:" << p00d10 << ", p00d11: " << p00d11 << endl;
//cout << "p01d00:" << p01d00 << ", p01d01:" << p01d01 << ", p01d10:" << p01d10 << ", p01d11: " << p01d11 << endl;
//cout << "p11d00:" << p11d00 << ", p11d01:" << p11d01 << ", p11d10:" << p11d10 << ", p11d11: " << p11d11 << endl;
    double rateDropoutCell = listCellDropoutRates[cell];
//cout << "rateDropoutCell: " << rateDropoutCell << endl;
    double p00d00prior = allele0Freq*allele0Freq*(1.0-rateDropoutCell)*(1.0-rateDropoutCell);
    double p00d01prior = allele0Freq*allele0Freq*rateDropoutCell*(1.0-rateDropoutCell);
    double p00d11prior = allele0Freq*allele0Freq*rateDropoutCell*rateDropoutCell;
    double p01d00prior = 2.0*(1.0-allele0Freq)*allele0Freq*(1.0-rateDropoutCell)*(1.0-rateDropoutCell);
    double p01d01prior = 2.0*(1.0-allele0Freq)*allele0Freq*rateDropoutCell*(1.0-rateDropoutCell);
    double p01d11prior = 2.0*(1.0-allele0Freq)*allele0Freq*rateDropoutCell*rateDropoutCell;
    double p11d00prior = (1.0-allele0Freq)*(1.0-allele0Freq)*(1.0-rateDropoutCell)*(1.0-rateDropoutCell);
    double p11d01prior = (1.0-allele0Freq)*(1.0-allele0Freq)*rateDropoutCell*(1.0-rateDropoutCell);
    double p11d11prior = (1.0-allele0Freq)*(1.0-allele0Freq)*rateDropoutCell*rateDropoutCell;
//cout << "p00d00prior:" << p00d00prior << ", p00d01prior:" << p00d01prior << ", p00d11prior:" << p00d11prior << endl;
//cout << "p01d00prior:" << p01d00prior << ", p01d01prior:" << p00d01prior << ", p01d11prior:" << p01d11prior << endl;
//cout << "p11d00prior:" << p11d00prior << ", p11d01prior:" << p11d01prior << ", p11d11prior:" << p11d11prior << endl;
    double facNorm = p00d00*p00d00prior + p00d01*p00d01prior + p00d10*p00d01prior + p00d11*p00d11prior + p01d00*p01d00prior + p01d01*p01d01prior + p01d10*p01d01prior + p01d11*p01d11prior + p11d00*p11d00prior + p11d01*p11d01prior + p11d10*p11d01prior + p11d11*p11d11prior;
//cout << "facNorm: " << facNorm << endl;
    listGenoProbs.clear();
    double p00 = (p00d00*p00d00prior + p00d01*p00d01prior + p00d10*p00d01prior + p00d11*p00d11prior )/ facNorm;
    double p01 = ( p01d00*p01d00prior + p01d01*p01d01prior + p01d10*p01d01prior + p01d11*p01d11prior )/facNorm;
    double p11 = ( p11d00*p11d00prior + p11d01*p11d01prior + p11d10*p11d01prior + p11d11*p11d11prior )/facNorm;
    listGenoProbs.push_back(p00);
    listGenoProbs.push_back(p01);
    listGenoProbs.push_back(p11);
    MakeProbInRange( listGenoProbs[0] );
    MakeProbInRange( listGenoProbs[1] );
    MakeProbInRange( listGenoProbs[2] );
//cout << "probability: " << listGenoProbs[0] << ", " << listGenoProbs[1] << ", " << listGenoProbs[2] << endl;
//#endif
}

static void ConvGenoProbToOututProb(const vector<double> &listGenoProbs, bool fBinary, vector<double> &listOutProbs)
{
    listOutProbs.clear();
    if( fBinary == false )
    {
        for(int i=0; i<(int)listGenoProbs.size()-1; ++i)
        {
            listOutProbs.push_back(listGenoProbs[i]);
        }
    }
    else
    {
        // only output the zero prob
        listOutProbs.push_back(listGenoProbs[0]);
    }
}

static int CalcAllele0FreqForGenos(const vector<int> &listGenos)
{
    //
    int res = (int)(listGenos.size())*2;
    //if( fBinary == false )
    //{
    //    res *= 2;
    //}
    for( int i=0; i<(int)listGenos.size(); ++i )
    {
        res -= listGenos[i];
    }
    return res;
}
static int CountTotAlleleForGenos(const map<int,pair<int,int> > &mapGenosOfCells)
{
    //
    int res = 0;

    //if( fBinary )
    //{
    //    res = mapGenosOfCells.size();
    //}
    //else
    //{
    for( map<int,pair<int,int> > :: const_iterator it = mapGenosOfCells.begin(); it != mapGenosOfCells.end(); ++it )
    {
        res += it->second.first;
        res += it->second.second;
    }
    // }
    return res;
}
static int CallGenoWithProb(const vector<double> &listProbs)
{
    int res = -1;
    double probMax = -1.0;
    double probSum = 0.0;
    for(int i=0; i<(int)listProbs.size(); ++i)
    {
        if( probMax < listProbs[i])
        {
            res = i;
            probMax = listProbs[i];
        }
        probSum += listProbs[i];
    }
    double probLast = 1.0-probSum;
    if( probLast < 0.0 )
    {
        probLast = 0.0;
//cout << "probSum: " << probSum << endl;
//        cout << "FATAL ERROR: genotype probability should be between 0 and 1, and all probability should sum to 1.0 for each genotype.\n";
//        exit(1);
    }
    if( probMax < probLast )
    {
        res = listProbs.size();
    }
    return res;
}
static void ConvReadBasedProbToGeno( const vector<vector<vector<double> > > &listReadBasedProb, vector<vector<int> > &listGenos)
{
    listGenos.resize( listReadBasedProb.size() );
    for(int i=0; i<(int)listReadBasedProb.size(); ++i)
    {
        listGenos[i].resize( listReadBasedProb[i].size() );
        for(int j=0; j<(int)listReadBasedProb[i].size(); ++j)
        {
            int geno = CallGenoWithProb( listReadBasedProb[i][j] );
            listGenos[i][j] = geno;
        }
    }
}

static void AdjustReadCountForDrop(int cell, int site, const map<pair<int,int>, pair<int,int> > &setDropPos, pair<int,int> &readCount)
{
cout << "At cell " << cell << ", site: " << site << ", read count: " << readCount.first << " " << readCount.second << endl;
    pair<int,int> pp(cell, site);
    if( setDropPos.find(pp) != setDropPos.end() )
    {
        map<pair<int,int>,pair<int,int> > :: const_iterator it = setDropPos.find(pp);
        // make read count of 1 to be 0 if dropout
//cout << "Setting read count to be 0 at site " << site << ", cell " << cell << endl;
        if( it->second.first == 1 )
        {
cout << "Setting read count of 0-allele to be 0 at site " << site << ", cell " << cell << endl;
            readCount.first = 0;
        }
        if( it->second.second == 1 )
        {
cout << "Setting read count of 1-allele to be 0 at site " << site << ", cell " << cell << endl;
            readCount.second = 0;
        }
    }
}

// analyze reads counts
static void AnalyzeReadsAtSites( const vector<pair<int,int> > &listReadCnts, double &fracAllele0, double &aveReadDepthSite, double &stdReadDepthSite)
{
//cout << "AnalyzeReadsAtSites: listReadCnts size: " << listReadCnts.size() << endl;
    //
    int numNonDoubleDrops = 0;
    int numZeros = 0, numAlleleTot = 0;
    double sumReadDepthSingleStrand = 0.0;
    vector<double> listReadSingle;
    for(int i=0; i<(int)listReadCnts.size(); ++i)
    {
        if( listReadCnts[i].first+listReadCnts[i].second>0)
        {
            int numReads = listReadCnts[i].first+listReadCnts[i].second;
            double rateDropoutCell = listCellDropoutRates[i];
            double readSingleExp = (1.0-rateDropoutCell)*numReads/2 + rateDropoutCell*numReads;
            sumReadDepthSingleStrand += readSingleExp;
            listReadSingle.push_back(readSingleExp);
            ++numNonDoubleDrops;
        }
        numZeros += listReadCnts[i].first;
        numAlleleTot += listReadCnts[i].first + listReadCnts[i].second;
    }
    fracAllele0 = ((double)numZeros+1)/((double)numAlleleTot+1);
    aveReadDepthSite = (sumReadDepthSingleStrand+1.0)/(numNonDoubleDrops+1.0);
//cout << "numZeros: " << numZeros << ", numOnes: " << numAlleleTot << ", numReadDepth; " << sumReadDepthSingleStrand << ", numNonDoubleDrop: "  << numNonDoubleDrops << endl;
    // std
    double sum = 0.0;
    for(int i=0; i<(int)listReadSingle.size(); ++i)
    {
        double diff = listReadSingle[i] - aveReadDepthSite;
        sum += diff*diff;
    }
    stdReadDepthSite = sqrt((sum+1.0)/(listReadSingle.size()+1.0));
}

static double EstNoneZeroFracGenoFromReadsAtSite( const vector<pair<int,int> > &listReadCnts )
{
    //
    int numNoneZeroGenos = 0, numZeroGenos = 0;
    for(int i=0; i<(int)listReadCnts.size(); ++i)
    {
        if( listReadCnts[i].first > 0 || listReadCnts[i].second > 0 )
        {
            if( listReadCnts[i].second > 0)
            {
                ++numNoneZeroGenos;
            }
            else
            {
                ++numZeroGenos;
            }
        }
    }
    return ( numNoneZeroGenos+1.0 )/( numNoneZeroGenos + numZeroGenos + 1.0 );
}

static void InfDroputsByReadsAtSite( const vector<pair<int,int> > &listReadCnts, vector<int> &listNumDropouts)
{
    // based on the read counts on the genotypes of this site, infer the number of dropouts
    // approach: if read counts are zero, then double droput; otherwise, find the bipartiton
    vector<int> listTotNonZeroReadCounts;
    for(int i=0; i<(int)listReadCnts.size(); ++i)
    {
        int rc = listReadCnts[i].first + listReadCnts[i].second;
        if( rc > 0)
        {
            listTotNonZeroReadCounts.push_back(rc);
        }
    }
    int bound0and1 = FindTwoClusterBoundary(listTotNonZeroReadCounts);
    listNumDropouts.clear();
    for(int i=0; i<(int)listReadCnts.size(); ++i)
    {
        int rc = listReadCnts[i].first + listReadCnts[i].second;
        int ndrops = 0;
        if( rc == 0)
        {
            ndrops = 2;
        }
        else if( rc <= bound0and1)
        {
            ndrops = 1;
        }
        listNumDropouts.push_back(ndrops);
    }
//cout << "bound0and1: " << bound0and1 << ", and number of inferred dropouts: " << endl;
//cout << "InfDroputsByReadsAtSite: listReadCounts: ";
//for(int i=0; i<(int)listReadCnts.size(); ++i)
//{
//cout << "(" << listReadCnts[i].first << "," << listReadCnts[i].second << ") d: " << listNumDropouts[i] << " ";
//}
//cout << endl;
}

static void SkipAmbiguousSites( vector<vector<vector<double> > > &listGenoProbs )
{
    // noise cancelation: don't use ambiguous sites
    // do nothing if not binary
    if( fBinary == false )
    {
        return;
    }
    for(int i=0; i<(int)listGenoProbs.size(); ++i)
    {
        for(int j=0; j<(int)listGenoProbs[i].size(); ++j)
        {
            if( listGenoProbs[i][j].size() != 1 )
            {
                cout << "Wrong size\n";
                exit(1);
            }
            double genoProb = listGenoProbs[i][j][0];
            if( genoProb < 0.5 && genoProb > 0.5-thresMinFracDiffNoMiss/2 )
            {
                genoProb = 0.499999;
            }
            else if( genoProb > 0.5 && genoProb < 0.5+thresMinFracDiffNoMiss/2 )
            {
                genoProb = 0.500001;
            }
            listGenoProbs[i][j][0] = genoProb;
        }
    }
}


//******************************************************************************************************
// Given: list of read counts at a site

static void CalcGenoProbForReadCountsAtSite( const vector<pair<int,int> > &listReadCntsAtSite, vector< vector<double> > &listGenoProbsSite )
{
//cout << "List of alleles: ";
//for(int i=0; i<(int)listReadCntsAtSite.size(); ++i)
//{
//cout << listReadCntsAtSite[i].first << "," << listReadCntsAtSite[i].second << "  ";
//}
//cout << endl;

    // now calculate the probability
    double allele0Freq = 0.0, aveReadDepthSite = 0.0, stdReadDepthSite = 0.0;
    AnalyzeReadsAtSites( listReadCntsAtSite, allele0Freq, aveReadDepthSite, stdReadDepthSite  );
//cout << "Allelefreq: " << allele0Freq << ", aveDepth: " << aveReadDepthSite << ", std: " << stdReadDepthSite << endl;

    double fracNoneZeroGenos = EstNoneZeroFracGenoFromReadsAtSite(listReadCntsAtSite);

    // also impute dropouts
    vector<int> listInfDropouts;
    InfDroputsByReadsAtSite(listReadCntsAtSite, listInfDropouts);

    //
    listGenoProbsSite.clear();
    for(int j=0; j<(int)listReadCntsAtSite.size(); ++j)
    {
        // simulate prob
        //cout << "Prob at site: " << i << " cell " << j << endl;
        pair<int,int> readCnt = listReadCntsAtSite[j];
        vector<double> probGenoReads;

        if( readCnt.first ==0 && readCnt.second == 0 && fBinary )
        {
            //
            if( fNoReadsAsMiss == false )
            {
                probGenoReads.push_back(1.0-fracNoneZeroGenos);
                probGenoReads.push_back(fracNoneZeroGenos);
            }
            else
            {
                // treat as missing values
                double probMiss0 = 0.5, probMiss1 = 0.5;
                if( fracNoneZeroGenos < 0.5 )
                {
                    probMiss0 = 0.500001;
                    probMiss1 = 0.499999;
                }
                else
                {
                    probMiss1 = 0.500001;
                    probMiss0 = 0.499999;
                }
                probGenoReads.push_back(probMiss0);
                probGenoReads.push_back(probMiss1);
            }
        }
        else
        {
            //cout << "Site: " << i << ", cell: " << j << endl;
            CalcGenotypeProbOfReadsNew( j, aveReadDepthSite, stdReadDepthSite, readCnt, allele0Freq, listInfDropouts[j], probGenoReads );
        }

        //CalcGenotypeProbOfReads( readCnt, allele0Freq, probGenoReads );
        vector<double> probGenoReadsOut;
        ConvGenoProbToOututProb(probGenoReads, fBinary, probGenoReadsOut);

        listGenoProbsSite.push_back(probGenoReadsOut);
        //cout << "[" << i << "," << j << "]: prob=" << probGenoReadsOut[0] << ", readCnt: " << readCnt.first << "," << readCnt.second << endl;
    }
}


//******************************************************************************************************
// Usage: ./xxxx allele-counts-file
// allele-count-file: each line is for a SNV site;

const char *VERSION = "scprob ver 1.2.0, May 12, 2019";

//
int main_in_c(int argc, char **argv)
{
    if( argc >= 4 || argc <= 1  )
    {
        cout << VERSION << endl;
        cout << "Usage: ./scprob-v1.2.0  <allele-count-file>\n";
        exit(1);
    }
    if( strcmp(argv[2],"0")==0 )
    {
        fBinary = true;
    } else {
        fBinary = false;
    }

    string outputfile = argv[argc - 2];
    string str2 = "rscistree.counts";
    string str3 = "rscistree.input";
    outputfile.replace(outputfile.find(str2), str2.length(), str3);

    ofstream out(outputfile);
    streambuf *coutbuf = cout.rdbuf(); // save old buf
    cout.rdbuf(out.rdbuf());          // redirect cout to out.txt!

    // read the allele count file
    vector<vector<pair<int,int> > > listSCAlleles;
    ifstream inFile(argv[1]);
    bool fHeader = true;
    int nc = 0, ns = 0;
    while( !inFile.eof() )
    {
        char buf[102400];
        inFile.getline(buf, 102400);
        //inFile.close();
        string strBuf = buf;

        if( strBuf.length() == 0 )
        {
            break;
        }

        stringstream ssinput(strBuf);
//cout << "Line: " << strBuf << endl;

        if( fHeader )
        {
            // read in the num of sites and cells
            // for now simply skip header
            fHeader = false;

            // read in the number of cells
            string strTemp;
            ssinput >> strTemp;
            ssinput >> ns;
            ssinput >> nc;
//cout << "Number of sites: " << ns << ", number of cells: " << nc << endl;
        }
        else
        {
//cout << "* A new site\n";
            // read in the counts
            vector<pair<int,int> > listAllels;
            for(int i=0; i<nc; ++i)
            {
                int na0=0, na1=0;
                ssinput >> na0 >> na1;
                pair<int,int> pp(na0, na1);
                listAllels.push_back(pp);
//cout << "Alleles: " << na0 << "," << na1 << endl;
            }
            listSCAlleles.push_back(listAllels);
        }

    }

    inFile.close();

    // initialize dropout rates
    listCellDropoutRates.clear();
    for(int c=0; c<nc; ++c)
    {
        listCellDropoutRates.push_back( rateDropout );
    }

    // calculate geno prob
    vector< vector< vector<double> > > listReadBasedProb;
    for(int s=0; s<(int)listSCAlleles.size(); ++s)
    {
//cout << "Processing site: " << s << endl;
        vector<vector<double> > listGenoProbsSite;
        CalcGenoProbForReadCountsAtSite( listSCAlleles[s], listGenoProbsSite );
        listReadBasedProb.push_back(listGenoProbsSite);
    }

    // apply filter here
    SkipAmbiguousSites(listReadBasedProb);

    //
    int numSC = 0;
    if( fBinary )
    {
        numSC = (int)listReadBasedProb[0].size();
        cout << "HAPLOID ";
    }
    else
    {
        numSC = (int)(listReadBasedProb[0].size());
        cout << "TERNARY ";
    }
    cout << listReadBasedProb.size() << " " << numSC << " ";
    for(int i=0; i<numSC; ++i)
    {
      cout << "c" << i+1 << " ";
    }
    cout << endl;
    cout.precision(10);
    for(int i=0; i<(int)listReadBasedProb.size(); ++i)
    {
        cout << "s" << i+1 << " ";
        for( int j=0; j<(int)listReadBasedProb[i].size(); ++j)
        {
            for(int k=0; k<(int)listReadBasedProb[i][j].size(); ++k)
            {
                cout << listReadBasedProb[i][j][k] << " " ;
            }
            cout << " ";
        }
        cout << endl;
    }

    cout.rdbuf(coutbuf); // reset to standard output again
    out.close();

    return 0;
}
