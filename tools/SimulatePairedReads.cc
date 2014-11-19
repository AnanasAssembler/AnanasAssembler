#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "base/CommandLineParser.h"
#include "test/ReadSimulator.h"


int main( int argc, char** argv )
{
  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> nOutCmmd("-o","outputName", "SIM");
  commandArg<int> iCmmd("-s","step/interval size", 35);
  commandArg<double> indelCmmd("-indel","indel error rate", 0.);
  commandArg<double> subCmmd("-S","substitution error rate", 0.);
  commandArg<int> olCmmd("-O","Accepted Overlap", 15);
  commandArg<string> appCmmd("-L","Application logging file","application.log");
  commandLineParser P(argc,argv);
  P.SetDescription("Simulating paired-end reads for assembly");
  P.registerArg(fileCmmd);
  P.registerArg(nOutCmmd);
  P.registerArg(iCmmd);
  P.registerArg(indelCmmd);
  P.registerArg(subCmmd);
  P.registerArg(olCmmd);
  P.registerArg(appCmmd);
  
  P.parse();
  
  string fileName        = P.GetStringValueFor(fileCmmd);
  string outName         = P.GetStringValueFor(nOutCmmd);
  int interval           = P.GetIntValueFor(iCmmd);
  double indel           = P.GetDoubleValueFor(indelCmmd);
  double sub             = P.GetDoubleValueFor(subCmmd);
  double minOverlap      = P.GetIntValueFor(olCmmd);
  string applicationFile = P.GetStringValueFor(appCmmd);
  vecDNAVector inputSeqs;
  inputSeqs.Read(fileName);

  FILE* pFile = fopen(applicationFile.c_str(), "w");
  Output2FILE::Stream()     = pFile;
  FILELog::ReportingLevel() = logINFO; 

  cout<<"Simulating reads..."<<endl;
  ReadSimulator  readSim;
  readSim.generatePEReads(inputSeqs, interval, indel, sub);  

  ofstream fout1;
  fout1.open((outName+"_reads.fa").c_str());
  readSim.writeReads(fout1);
  fout1.close();

  cout<<"Getting read overlaps..."<<endl;
  AllReadOverlaps allOverlaps(readSim.getReadCnt());
  readSim.findAllOverlaps(allOverlaps, minOverlap);
  string outputFile = outName+"_readOverlaps.out";
  allOverlaps.write(outputFile, 1);

  return 0;
}
