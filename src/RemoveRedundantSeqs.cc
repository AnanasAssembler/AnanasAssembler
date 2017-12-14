#include <string>
#include <set>
#include "cola/src/fastAlign/FastAlignUnit.h"
#include "ryggrad/src/base/CommandLineParser.h"
#include "cola/src/cola/Cola.h"
#include "ContScaff.h"
#include "ContScaffIO.h"


void handleAlignments(svec<AlignmentInfo> currAlignments, set<string>& contigsToRemove, float identRatio) {
  for(AlignmentInfo& algn:currAlignments) {
    if(algn.getTargetBaseAligned() >= identRatio*algn.getTargetLength() &&
      algn.getIdentity()>=identRatio) {
      contigsToRemove.insert(algn.getTargetName());
    }
    if(algn.getQueryBaseAligned() >= identRatio*algn.getQueryLength() &&
      algn.getIdentity()>=identRatio) {
      contigsToRemove.insert(algn.getQueryName());
      return; //The current sequence being queried has been found as redundant
    }
  }
}

void removeRedundants(const string& fastaFileName, const string& layoutFileName,
                      set<string>& contigsToRemove) {
  ofstream fout;
  fout.open(fastaFileName + ".rr");
  vecDNAVector dna;
  dna.Read(fastaFileName);
  for (int i=0; i<dna.isize(); i++) {
    string name = ">";
    name += dna.NameClean(i);
    const DNAVector & d = dna[i];
    if (contigsToRemove.count(name)==0) { //Keep
      fout << name << endl;
      for (int j=0; j<d.isize(); j++)
        fout << d[j];
      fout << endl;
    }
  }
  fout.close();

  Assembled assembly;
  ContigScaffoldIO io;
  io.Read(assembly, layoutFileName);
  // Main loop over scaffolds/contigs to remove redundant contigs
  for (int scaffCnt=0; scaffCnt<assembly.isize(); scaffCnt++) {
    Scaffold & currScaff  = assembly[scaffCnt];
    for (int contigCnt=0; contigCnt<assembly[scaffCnt].isize(); contigCnt++) {
      string name = currScaff[contigCnt].Name();
      if (contigsToRemove.count(name)!=0) { // Do not keep
        currScaff[contigCnt].SetDiscard(true);
      } 
    }
  }
  io.Write(assembly, layoutFileName+".rr");
}
 

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
  FastAlignTargetUnit qUnit(fastaFileName, 1);

  set<string> contigsToRemove; 
  FastAlignUnit FAUnit(fastaFileName, qUnit, params, numThreads);
  FastAlignUnit FAUnitRev(fastaFileName, qUnit, params, numThreads, true);
  int progCount = 0;
  int totSize   = FAUnit.getNumQuerySeqs();
  cout<<"Total number of query sequences: " << totSize <<endl;
  int inc       = totSize/10000;
  if (inc < 1) inc = 1;
  for(int qIdx=0; qIdx<FAUnit.getNumQuerySeqs(); qIdx++) {
    if(!contigsToRemove.count(FAUnit.getQuerySeqName(qIdx))) {
      svec<AlignmentInfo> currAlignments;
      FAUnit.alignSequence(qIdx, currAlignments);
      handleAlignments(currAlignments, contigsToRemove, minIdent);
      svec<AlignmentInfo> currAlignments_rev;
      FAUnitRev.alignSequence(qIdx, currAlignments_rev);
      handleAlignments(currAlignments_rev, contigsToRemove, minIdent);
    }
    progCount++;
   if (progCount % inc == 0) 
     cout << "\r===================== " << 100.0*progCount/totSize 
          << "%  " << flush; 
  }
  removeRedundants(fastaFileName, layoutFileName, contigsToRemove);
  return 0;

}

