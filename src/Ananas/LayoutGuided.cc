#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "base/CommandLineParser.h"
#include "src/Ananas/SearchOverlaps.h"
#include "src/Ananas/ContScaff.h"
#include "src/Ananas/ContScaffIO.h"
#include "src/Ananas/ConsensOverlapUnit.h"



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input read pair/size info file");
  commandArg<string> lapCmmd("-l","input read overlap file");
  commandArg<string> consCmmd("-g","input read consensus group file");
  commandArg<string> scaffCmmd("-s","scaffolds file");
  commandArg<string> layoutCmmd("-o","output layout file", "contigs_guided.layout");
  commandArg<string> dirCmmd("-dir","direction of pairs (fr or ff)");
  commandArg<double> minCmmd("-m","minimum overlap identity", 0.985);
  commandArg<int> sizeCmmd("-size","minimum length of a single-contig scaffold to report", 300);
  commandArg<bool> exCmmd("-e","DO NOT DO exhaustive search (report top-n)", false);
  commandArg<int> numCmmd("-num","Process # (for parallel runs)", 0);
  commandArg<int> ofCmmd("-of","Out of # processes (for parallel runs)", 1);
  commandLineParser P(argc,argv);
  P.SetDescription("Assembles COUnit from overlaps.");
  P.registerArg(fileCmmd);
  P.registerArg(lapCmmd);
  P.registerArg(consCmmd);
  P.registerArg(scaffCmmd);
  P.registerArg(layoutCmmd);
  P.registerArg(minCmmd);
  P.registerArg(exCmmd);
  P.registerArg(dirCmmd);
  P.registerArg(sizeCmmd);
  P.registerArg(numCmmd);
  P.registerArg(ofCmmd);
 
  P.parse();
  
  string pairSzFileName = P.GetStringValueFor(fileCmmd);
  string lapName = P.GetStringValueFor(lapCmmd);
  string layoutName = P.GetStringValueFor(layoutCmmd);
  string scaffName = P.GetStringValueFor(scaffCmmd);
  double mI = P.GetDoubleValueFor(minCmmd);
  bool bEx2 = P.GetBoolValueFor(exCmmd);
  int minSize = P.GetIntValueFor(sizeCmmd);
  int num = P.GetIntValueFor(numCmmd);
  int of = P.GetIntValueFor(ofCmmd);
  string consName = P.GetStringValueFor(consCmmd);
 
  string dir = P.GetStringValueFor(dirCmmd);
  bool bEx = true;
  if (bEx2)
    bEx = false;

  ConsensOverlapUnit COUnit(pairSzFileName, consName, "");

  int i, j, k, l;
 
  Assembled assembly;
  ContigScaffoldIO io;

  io.Read(assembly, scaffName);
  
  svec<int> ids;
  ids.resize(COUnit.GetNumReads(), 0);
  //=========================================================
  // Collect all COUnit (WARNING: Duplicated code!!)
  for (l=num; l<assembly.isize(); l += of) {
    const Scaffold & s = assembly[l];
    if (s.isize() == 1) {
      if (s[0].Highest() < minSize) {
	continue;
      }
    }
    for (i=0; i<s.isize(); i++) {
      const Contig & c = s[i];
      for (j=0; j<c.isize(); j++) {
        const ReadPlacement & r = c[j]; 
        int id = r.Read();
	ids[id] = 1;
      }
    }
  }
  //=========================================================
  // Free up some unused space
  // COUnit.Prune(ids);

  COUnit.ReadOverlaps(lapName, ids);

  Search search;
  search.SetExhaustive(bEx);
  search.SetDir(dir);

  if (of > 1) {
    char tmp[256];
    sprintf(tmp, ".%d", num);
    layoutName += tmp;
  }

  search.SetOutput(layoutName);
  search.SetIndex(num);
  search.SetMinAltKeep(minSize);

  for (l=num; l<assembly.isize(); l += of) {
    const Scaffold & s = assembly[l];
    if (s.isize() == 1) {
      if (s[0].Highest() < minSize) {
	cout << "Skipping scaffold " << l << endl;
	continue;
      }
    }

    cout << "Processing scaffold " << l << endl;
    search.SetUsedAll(COUnit);
    for (i=0; i<s.isize(); i++) {
      const Contig & c = s[i];
      for (j=0; j<c.isize(); j++) {
        const ReadPlacement & r = c[j]; 
        int id = r.Read();
	search.SetUsed(id, false);
      }
    }
    cout << "Searching... " << l << endl;
    search.SetExhaustive(true);
    if (s.isize() == 1)
      search.SetExhaustive(false);
      

    search.DoSearchAll(COUnit);
  }
 
  cout << "All done!!!" << endl;
  return 0;
}
