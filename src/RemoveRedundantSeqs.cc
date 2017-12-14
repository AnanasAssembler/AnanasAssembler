#ifndef FORCE_DEBUG
#define NDEBUG
#endif


#include "ryggrad/src/base/CommandLineParser.h"
#include "RedundanceRemoval.h"

int main(int argc,char** argv)
{
  commandArg<string> fastaFileCmmd("-if","input fasta file");
  commandArg<string> layoutFileCmmd("-il","input layout file");
  commandArg<double> minIdentCmmd("-ident","Minium acceptable identity", 0.99);
  commandArg<int>    bandedCmmd("-b", "The bandwidth for banded mode", 3);
  commandArg<int>  threadCmmd("-n","Number of Cores to run with", 1);
  commandArg<string> appLogCmmd("-L","Application logging file","application.log");

  commandLineParser P(argc,argv);
  P.SetDescription("Removes redundant transcripts if too similar.");
  P.registerArg(fastaFileCmmd);
  P.registerArg(layoutFileCmmd);
  P.registerArg(minIdentCmmd);
  P.registerArg(bandedCmmd);
  P.registerArg(threadCmmd);
  P.registerArg(appLogCmmd);

  P.parse();

  string fastaFileName         = P.GetStringValueFor(fastaFileCmmd);
  string layoutFileName        = P.GetStringValueFor(layoutFileCmmd);
  double      minIdent         = P.GetDoubleValueFor(minIdentCmmd);
  int         banded           = P.GetIntValueFor(bandedCmmd);
  int         numThreads       = P.GetIntValueFor(threadCmmd);
  string      applicationFile  = P.GetStringValueFor(appLogCmmd);

  FILE* pFile               = fopen(applicationFile.c_str(), "w");
  Output2FILE::Stream()     = pFile;
  FILELog::ReportingLevel() = logINFO; 
 
  AlignmentParams params(1, 20, minIdent, banded, 0.75);
  RedundanceRemoval rr;
  rr.removeRedundants(params, fastaFileName, layoutFileName, numThreads); 
  return 0;
}

