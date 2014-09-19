#include <string>
#include <omp.h>
#include "base/CommandLineParser.h"
#include "extern/logger/log.h"
#include "src/Ananas/ConsensOverlapUnit.h"


int main(int argc,char** argv)
{

  commandArg<string> aCmmd("-i","FASTA file containing input reads");
  commandArg<string> a1Cmmd("-p","FASTA file containing paired input reads", "");
  commandArg<string> bCmmd("-o","File to Output overlaps", "allOverlaps.out");
  commandArg<string> b0Cmmd("-t","Temporary file for pair/size info to be used by layout search", "pairSz.tmp");
  commandArg<string> b1Cmmd("-C","File to Output consensus reads", "consensusReads.out");
  commandArg<bool>   cCmmd("-s","Single-strand data", true);
  commandArg<int>    dCmmd("-b","Subread block step", 10);
  commandArg<int>    eCmmd("-S","Seed size for choosing candidates", 15);
  commandArg<double> fCmmd("-I","Minimum acceptable identity for overlapping sequences", 0.99);
  commandArg<double> gCmmd("-c","Minimum Coverage of the read ends for an acceptable extensiont", 0.98);
  commandArg<int>    hCmmd("-O","Minimum overlap of an acceptable extension", 40);
  commandArg<int>    iCmmd("-B","Bandwidth for local alignments", 3);
  commandArg<int>    jCmmd("-l","Minimum length of assembled scaffold", 400);
  commandArg<string> kCmmd("-a","Auxillary information output file", "aux.out");
  commandArg<string> lCmmd("-L","Application logging file","application.log");
  commandArg<int>    threadCmmd("-T","Number of Cores to run with", 2);

  commandLineParser P(argc,argv);
  P.SetDescription("Overlap finder for assembly");

  P.registerArg(aCmmd);
  P.registerArg(a1Cmmd);
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

  P.parse();

  string inputFile       = P.GetStringValueFor(aCmmd);
  string pairFile        = P.GetStringValueFor(a1Cmmd);
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

  FILE* pFile = fopen(applicationFile.c_str(), "w");
  Output2FILE::Stream()     = pFile;
  FILELog::ReportingLevel() = logINFO; 

  #if defined(OPEN_MP)
    omp_set_num_threads(numOfThreads); //The sort functions still use MPI
  #endif

  AssemblyParams params(singleStrand, readBlockSize, seedSize,
                        minIdent, minCoverage, minOverlap,
                        alignBand, minBasePerScaf);
  ConsensOverlapUnit COUnit(params, inputFile);
  COUnit.findOverlaps(numOfThreads);

  COUnit.writePairSzInfo(pairSzFile);
  COUnit.writeOverlaps(overlapFile, 1);
  COUnit.writeConsensInfo(consensFile, 1);
  //COUnit.writeConsensReads("temp_consensRead.fa");

  return 0;
}
  


