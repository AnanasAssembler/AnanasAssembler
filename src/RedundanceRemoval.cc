#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <set>
#include "RedundanceRemoval.h"
#include "cola/src/fastAlign/FastAlignUnit.h"
#include "ContScaff.h"
#include "ContScaffIO.h"

void RedundanceRemoval::removeRedundants(const AlignmentParams& aParams, const string& fastaFileName,
                                         const string& layoutFileName, int numThreads) {
  FastAlignTargetUnit qUnit(fastaFileName, 1);
  set<string> contigsToRemove; 
  FastAlignUnit FAUnit(fastaFileName, qUnit, aParams, numThreads);
  FastAlignUnit FAUnitRev(fastaFileName, qUnit, aParams, numThreads, true);
  int progCount = 0;
  int totSize   = FAUnit.getNumQuerySeqs();
  cout<<"Total number of query sequences: " << totSize <<endl;
  int inc       = totSize/10000;
  if (inc < 1) inc = 1;
  for(int qIdx=0; qIdx<FAUnit.getNumQuerySeqs(); qIdx++) {
    if(!contigsToRemove.count(FAUnit.getQuerySeqName(qIdx))) {
      svec<AlignmentInfo> currAlignments;
      FAUnit.alignSequence(qIdx, currAlignments);
      handleAlignments(currAlignments, contigsToRemove, aParams.getMinIdentity());
      svec<AlignmentInfo> currAlignments_rev;
      FAUnitRev.alignSequence(qIdx, currAlignments_rev);
      handleAlignments(currAlignments_rev, contigsToRemove, aParams.getMinIdentity());
    }
    progCount++;
   if (progCount % inc == 0) 
     cout << "\r===================== " << 100.0*progCount/totSize 
          << "%  " << flush; 
  }
  cout << endl;
  removeItems(fastaFileName, layoutFileName, contigsToRemove);
}

void RedundanceRemoval::handleAlignments(svec<AlignmentInfo> currAlignments, set<string>& contigsToRemove, float identRatio) {
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

void RedundanceRemoval::removeItems(const string& fastaFileName, const string& layoutFileName,
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
 
