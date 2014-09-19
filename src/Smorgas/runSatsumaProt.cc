#ifndef FORCE_DEBUG
#define NDEBUG
#endif


#include "base/CommandLineParser.h"
#include "src/Smorgas/SatsumaProt.h"
#include <sstream>

int main( int argc, char** argv )
{

  commandArg<string> aStringCmmd("-t","target protein fasta", "");
  commandArg<string> dStringCmmd("-db","protein database in Spines binary format", "");
  commandArg<string> bStringCmmd("-q","query protein fasta");
  commandArg<int>    kCmmd("-k","k-mer size", 4);
  commandArg<int>    kStepCmmd("-kmerStep","Step size to be used in generating filteration kmers", 1);
  commandArg<int>    filterCmmd("-filter","Type of prefilter of hits to use- 1:fixed distance k-mer based 2:max k-mer based", 1);
  commandArg<int>    failCntCmmd("-allowFails","Number of failures to allow before stopping the search", 50);
  commandArg<int>    blockCmmd("-block","search only this subset (requires -e)", 0);
  commandArg<int>    nBlocksCmmd("-n_blocks","number of blocks (w/ -e)", 0);
  commandArg<int>    manyCmmd("-m","number of results to display", 50);
  commandArg<bool>   exCmmd("-e","exhaustive search (no filtering)", false);
  commandArg<bool>   selfCmmd("-self","self-alignments", false);
  commandArg<bool>   sameCmmd("-same","same alignments", false);
  commandArg<bool>   rnaCmmd("-rna","Do RNA alignments", false);
  commandArg<bool>   timeCmmd("-timestamps","Print time stamps", false);
  commandArg<double> cutoffCmmd("-cutoff","show only alignments at this (ingapped) identity or higher", 0.);
  commandArg<int>    kmerWindowSlideCmmd("-wSlide","Filter Window slide, if set to 1 the window sliding will cover all kmers.", 2);
  commandArg<double> eThreshCmmd("-E-value","show only alignments with better E-value", 0.001);
  commandArg<string> appLogCmmd("-l", "Application logging file","application.log");

  commandLineParser P(argc,argv);
  P.SetDescription("Satsuma-based (cross-correlation) protein alignment tool.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(dStringCmmd);
  P.registerArg(exCmmd);
  P.registerArg(kCmmd);
  P.registerArg(kStepCmmd);
  P.registerArg(filterCmmd);
  P.registerArg(failCntCmmd);
  P.registerArg(eThreshCmmd);
  P.registerArg(selfCmmd);
  P.registerArg(sameCmmd);
  P.registerArg(blockCmmd);
  P.registerArg(nBlocksCmmd);
  P.registerArg(rnaCmmd);
  P.registerArg(timeCmmd);
  P.registerArg(manyCmmd);
  P.registerArg(cutoffCmmd);
  P.registerArg(kmerWindowSlideCmmd);
  P.registerArg(appLogCmmd);


  P.parse();

  string REF      = P.GetStringValueFor(aStringCmmd);
  string DB       = P.GetStringValueFor(dStringCmmd);
  string QUERY    = P.GetStringValueFor(bStringCmmd);
  bool bExhaust   = P.GetBoolValueFor(exCmmd);
  bool bSelf      = P.GetBoolValueFor(selfCmmd);
  bool bSame      = P.GetBoolValueFor(sameCmmd);
  int protein_K   = P.GetIntValueFor(kCmmd);
  int kmerStep    = P.GetIntValueFor(kStepCmmd);
  int filterType  = P.GetIntValueFor(filterCmmd);
  int failAllowed = P.GetIntValueFor(failCntCmmd);
  int block       = P.GetIntValueFor(blockCmmd);
  int n_blocks    = P.GetIntValueFor(nBlocksCmmd);
  int limit       = P.GetIntValueFor(manyCmmd);
  bool bRNA       = P.GetBoolValueFor(rnaCmmd);
  double cutoff   = P.GetDoubleValueFor(cutoffCmmd);
  int kmerWSlide  = P.GetIntValueFor(kmerWindowSlideCmmd);
  string logFile  = P.GetStringValueFor(appLogCmmd);
  bool bQuiet     = true;
  double eThresh  = P.GetDoubleValueFor(eThreshCmmd);

  FILE* pFile = fopen(logFile.c_str(), "w");
  Output2FILE::Stream()     = pFile;
  FILELog::ReportingLevel() = logINFO; 
  FILE_LOG(logINFO) <<"Running Smörgås protein aligner ";

  if (P.GetBoolValueFor(timeCmmd))
    bQuiet = false;
  if (bSame)
    bExhaust = true;

  SatsumaProtParams params(bExhaust, bSelf, bSame, protein_K, kmerStep,
                           filterType, failAllowed, block, n_blocks, limit,
                           bRNA, bQuiet, cutoff, kmerWSlide, eThresh);
  SatsumaProt protAligner(REF, DB, QUERY, params);
  protAligner.alignAll();
}
