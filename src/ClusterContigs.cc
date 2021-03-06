#include <string>
#include <omp.h>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/Logger.h"
#include "ContigClusterUnit.h"


int main(int argc,char** argv)
{

    commandArg<string> aCmmd("-i","FASTA file containing input reads");
    commandArg<string> bCmmd("-o","File to Output overlaps", "contigClusters.out");
    commandArg<string> b0Cmmd("-t","Temporary file for pair/size info to be used by layout search", "pairSz.tmp");
    commandArg<string> b1Cmmd("-C","File to Output consensus reads", "contigConsesReads.out");
    commandArg<bool>   cCmmd("-s","Single-strand data", false);
    commandArg<int>    dCmmd("-b","Subread block step", 50);
    commandArg<int>    eCmmd("-S","Seed size for choosing candidates", 30);
    commandArg<double> fCmmd("-I","Minimum acceptable identity for overlapping sequences", 0.95);
    commandArg<double> gCmmd("-c","Minimum Coverage of the read ends for an acceptable extension", 0.96);
    commandArg<int>    hCmmd("-O","Minimum overlap of an acceptable extension", 300);
    commandArg<int>    iCmmd("-B","Bandwidth for local alignments", 4);
    commandArg<int>    jCmmd("-l","Minimum length of assembled scaffold", 400);
    commandArg<string> kCmmd("-a","Auxillary information output file", "aux.out");
    commandArg<string> lCmmd("-L","Application logging file","application.log");
    commandArg<int>    threadCmmd("-T","Number of Cores to run with", 2);
    commandArg<string> readGroupFileCmmd("-g","read groupin information file if available","");

    commandLineParser P(argc,argv);
    P.SetDescription("Cluster the Assembled contigs");
    
    P.registerArg(aCmmd);
    P.registerArg(bCmmd);
    P.registerArg(b0Cmmd);
    P.registerArg(b1Cmmd);
    P.registerArg(cCmmd);
    P.registerArg(dCmmd);
    P.registerArg(eCmmd);
    P.registerArg(fCmmd);
    P.registerArg(gCmmd);
    P.registerArg(hCmmd);
    P.registerArg(iCmmd);
    P.registerArg(jCmmd);
    P.registerArg(kCmmd);
    P.registerArg(lCmmd);
    P.registerArg(threadCmmd);
    P.registerArg(readGroupFileCmmd);
    P.parse();

    string inputFile       = P.GetStringValueFor(aCmmd);
    string overlapFile     = P.GetStringValueFor(bCmmd);
    string pairSzFile      = P.GetStringValueFor(b0Cmmd);
    string consensFile     = P.GetStringValueFor(b1Cmmd);
    bool   singleStrand    = P.GetBoolValueFor(cCmmd);
    int    readBlockSize   = P.GetIntValueFor(dCmmd);
    int    seedSize        = P.GetIntValueFor(eCmmd);
    double minIdent        = P.GetDoubleValueFor(fCmmd);
    double minCoverage     = P.GetDoubleValueFor(gCmmd);
    int    minOverlap      = P.GetIntValueFor(hCmmd);
    int    alignBand       = P.GetIntValueFor(iCmmd);
    int    minBasePerScaf  = P.GetIntValueFor(jCmmd);
    string auxFile         = P.GetStringValueFor(kCmmd);
    string applicationFile = P.GetStringValueFor(lCmmd);
    int    numOfCores      = P.GetIntValueFor(threadCmmd);
    int    numOfThreads    = P.GetIntValueFor(threadCmmd);
    string readGroupFile   = P.GetStringValueFor(readGroupFileCmmd);
    
    FILE* pFile = fopen(applicationFile.c_str(), "w");
    Output2FILE::Stream()     = pFile;
    FILELog::ReportingLevel() = logINFO; 
    
#if defined(OPEN_MP)
    omp_set_num_threads(numOfThreads); //The sort functions still use OMP
#endif
    
    AssemblyParams params(singleStrand, readBlockSize, seedSize,
                          minIdent, minCoverage, minOverlap,
                          alignBand, minBasePerScaf);
    ContigClusterUnit CLUnit(params, inputFile);
    CLUnit.clusterContigs(numOfThreads);
    CLUnit.writeContigClusters(overlapFile);
    CLUnit.writeContigPairs("pairs"+overlapFile);
    return 0;
}



