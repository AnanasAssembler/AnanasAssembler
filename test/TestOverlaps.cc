#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <string>
#include <omp.h>
#include "base/CommandLineParser.h"
#include "extern/logger/log.h"
#include "src/Ananas/ConsensOverlapUnit.h"
#include "test/ReadSimulator.h"


int main(int argc,char** argv)
{

    commandArg<string> aCmmd("-i","FASTA file containing input ref sequences to simulate reads from");
    commandArg<string> bCmmd("-o","File to Output overlaps", "allOverlaps.out");
    commandArg<string> b0Cmmd("-t","Temporary file for pair/size info to be used by layout search", "pairSz.tmp");
    commandArg<string> b1Cmmd("-C","File to Output consensus reads", "consensusReads.out");
    commandArg<bool>   cCmmd("-s","Single-strand data", true);
    commandArg<int>    dCmmd("-b","Subread block step", 10);
    commandArg<int>    eCmmd("-S","Seed size for choosing candidates", 15);
    commandArg<double> fCmmd("-I","Minimum acceptable identity for overlapping sequences", 0.99);
    commandArg<double> gCmmd("-c","Minimum Coverage of the read ends for an acceptable extension", 0.98);
    commandArg<int>    hCmmd("-O","Minimum overlap of an acceptable extension", 40);
    commandArg<int>    iCmmd("-B","Bandwidth for local alignments", 3);
    commandArg<int>    jCmmd("-l","Minimum length of assembled scaffold", 400);
    commandArg<string> kCmmd("-a","Auxillary information output file", "aux.out");
    commandArg<string> lCmmd("-L","Application logging file","application.log");
    commandArg<int>    threadCmmd("-T","Number of Cores to run with", 2);
    commandArg<double> readGroupIdentThreshCmmd("-d","read grouping threshold for identity",0.99);
    commandArg<string> readGroupFileCmmd("-g","read grouping information file if available","");
    commandArg<string> nOutCmmd("-prefix","simulation outputName prefix", "SIM");
    commandArg<int> intervalCmmd("-step","step/interval size", 35);
    commandArg<double> indelCmmd("-indel","indel error rate", 0.);
    commandArg<double> subCmmd("-sub","substitution error rate", 0.);

    commandLineParser P(argc,argv);
    P.SetDescription("Testing Overlap finder for assembly");
    
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
    P.registerArg(readGroupIdentThreshCmmd);
    P.registerArg(readGroupFileCmmd);

    P.registerArg(nOutCmmd);
    P.registerArg(intervalCmmd);
    P.registerArg(indelCmmd);
    P.registerArg(subCmmd);
 
    P.parse();

    string inputFile          = P.GetStringValueFor(aCmmd);
    string overlapFile        = P.GetStringValueFor(bCmmd);
    string pairSzFile         = P.GetStringValueFor(b0Cmmd);
    string consensFile        = P.GetStringValueFor(b1Cmmd);
    bool   singleStrand       = P.GetBoolValueFor(cCmmd);
    int    readBlockSize      = P.GetIntValueFor(dCmmd);
    int    seedSize           = P.GetIntValueFor(eCmmd);
    double minIdent           = P.GetDoubleValueFor(fCmmd);
    double minCoverage        = P.GetDoubleValueFor(gCmmd);
    int    minOverlap         = P.GetIntValueFor(hCmmd);
    int    alignBand          = P.GetIntValueFor(iCmmd);
    int    minBasePerScaf     = P.GetIntValueFor(jCmmd);
    string auxFile            = P.GetStringValueFor(kCmmd);
    string applicationFile    = P.GetStringValueFor(lCmmd);
    int    numOfCores         = P.GetIntValueFor(threadCmmd);
    int    numOfThreads       = P.GetIntValueFor(threadCmmd);
    double readGroupThresh    = P.GetDoubleValueFor(readGroupIdentThreshCmmd);
    string readGroupFile      = P.GetStringValueFor(readGroupFileCmmd);

    string outPrefix          = P.GetStringValueFor(nOutCmmd);
    int interval              = P.GetIntValueFor(intervalCmmd);
    double indel              = P.GetDoubleValueFor(indelCmmd);
    double sub                = P.GetDoubleValueFor(subCmmd);
    
    FILE* pFile               = fopen(applicationFile.c_str(), "w");
    Output2FILE::Stream()     = pFile;
    FILELog::ReportingLevel() = logINFO; 

#if defined(OPEN_MP)
    omp_set_num_threads(numOfThreads); //The sort functions still use OMP
#endif

    vecDNAVector inputSeqs;
    inputSeqs.Read(inputFile);

    cout<<"Simulating reads..."<<endl;
    ReadSimulator  readSim;
    readSim.generatePEReads(inputSeqs, interval, indel, sub);  

    ofstream fout1;
    string simReadOutputFile = outPrefix + "_reads.fa";
    fout1.open(simReadOutputFile.c_str());
    readSim.writeReads(fout1);
    fout1.close();

    cout<<"Finding Simulated read overlaps..."<<endl;
    AllReadOverlaps allSimOverlaps(readSim.getReadCnt());
    readSim.findAllOverlaps(allSimOverlaps, minOverlap);
    string simOverlapOutputFile = outPrefix + "_readOverlaps.out";
    allSimOverlaps.write(simOverlapOutputFile, 1);

    cout<<"Finding Overlaps ext..."<<endl;
    AssemblyParams params(singleStrand, readBlockSize, seedSize,
                          minIdent, minCoverage, minOverlap,
                          alignBand, minBasePerScaf);
    ConsensOverlapUnit COUnit(params, simReadOutputFile);
    COUnit.findOverlaps(numOfThreads, 0, 1, readGroupThresh, readGroupFile);
    COUnit.writePairSzInfo(pairSzFile);
    COUnit.writeOverlaps(overlapFile, 1);

    AllReadOverlaps allOverlaps = COUnit.getOverlaps();
    for(int i=0; i<allSimOverlaps.getSize(); i++) {
       ReadInfo rInfo = allSimOverlaps[i];
       for(int j=0; j<rInfo.getNumLaps(); j++) {
           if(!allOverlaps.hasOverlap(i, rInfo.getLap(j).getOverlapIndex())) {
               cout << i << "\t" <<allSimOverlaps[i].getLap(j).toString()<<endl;
           }
       }
    }

    return 0;
}
