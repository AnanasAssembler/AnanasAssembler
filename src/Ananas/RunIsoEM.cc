#define FORCE_DEBUG

#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "base/CommandLineParser.h"
#include "src/Ananas/ContScaff.h"
#include "src/Ananas/ContScaffIO.h"
#include "src/Ananas/ConsensOverlapUnit.h"


bool DoRemove(const Contig & major, const Contig & minor)
{
  int i;
  svec<ReadPlacement> all;
  svec<ReadPlacement> unique;
  int hi = 1;
  
  //  cout << "Major: " << major.isize() << " minor: " << minor.isize() << endl;

  for (i=0; i<major.isize(); i++) {
    all.push_back(major[i]);
    if (major[i].Stop() > hi)
      hi = major[i].Stop();
  }

  Sort(all);
  int n = 0;
  for (i=0; i<minor.isize(); i++) {
    const ReadPlacement & r = minor[i];
    int index = BinSearch(all, r);
    if (index < 0) 
      unique.push_back(r);
  }
  Sort(unique);
  for (i=1; i<unique.isize(); i++) {
    if (unique[i-1].Stop() <= unique[i].Start())
      continue;
    n++;
  }
  double cov = (double)major.isize()/(double)hi;
  double diff = (double)n/(double)minor.isize();
  
  //cout << "Ratio: " << diff/cov << " " << diff << " " << cov << endl;
  if (diff/cov < 0.1) {
    //cout << "Remove " << endl;
    return true;
  } else {
    //cout << "Keep" << endl;
    return false;
  }

}

/*void Clean(svec<int> & scaff) 
{
  int i, j;
  j = 1;
  for (i=0; i<scaff.isize(); i++) {
    if (j >= scaff.isize())
      break;
    int n = scaff[i];
    if (n == 0) {
      while (j < scaff.isize()) {
	if (scaff[j] == 0)
	  j++;
	else
	  break;
      }
      if (j >= scaff.isize())
	break;
      scaff[i] = scaff[j];
      scaff[j] = 0;
    }
  }
  scaff.resize(i);
  }*/

void CleanScaffold(Scaffold & scaff) 
{
  int i, j;
  j = 1;
  for (i=0; i<scaff.isize(); i++) {
    const Contig & c = scaff[i];
    //cout << "Contig " << i << " size " << c.isize() << endl;
    if (c.isize() == 0) {
      j = i+1;
      while (j < scaff.isize()) {
	if (scaff[j].isize() == 0)
	  j++;
	else
	  break;
      }
      if (j < scaff.isize()) {       
	scaff[i] = scaff[j];
	scaff[j].clear();
      }
    }
  }
  scaff.resize(i);
}

int main( int argc, char** argv )
{

    commandArg<string> fileCmmd("-i","input scaffold/contig/layout file");
    //commandArg<string> lapCmmd("-l","input read overlap file");
    //commandArg<string> fastaCmmd("-f","read fasta file"); // Do we really need this???
    //commandArg<string> consCmmd("-g","input read consensus group file");
    commandArg<string> outCmmd("-o","output file");
    //commandArg<string> pairCmmd("-dir","pairing (fr, ff, or na)");

  
    commandLineParser P(argc,argv);
    P.SetDescription("Joins contigs into scaffolds.");
    P.registerArg(fileCmmd);
    //P.registerArg(lapCmmd);
    //P.registerArg(consCmmd);
    //P.registerArg(fastaCmmd);
    P.registerArg(outCmmd);
    //P.registerArg(pairCmmd);
  
    P.parse();
  
    string fileName = P.GetStringValueFor(fileCmmd);
    //string lapName = P.GetStringValueFor(lapCmmd);
    //string fastaName = P.GetStringValueFor(fastaCmmd);
    string outName = P.GetStringValueFor(outCmmd);
    //string pairing = P.GetStringValueFor(pairCmmd);
    //string consName = P.GetStringValueFor(consCmmd); 
 
    Assembled assembly;

    ContigScaffoldIO io;

    io.Read(assembly, fileName);

    //ConsensOverlapUnit laps(fastaName, consName, "");


    
    int i, j, k;
    
    for (i=0; i<assembly.isize(); i++) {
      Scaffold & scaff = assembly[i];
      cout << "Scaffold " << i << endl;
      for (j=scaff.isize()-1; j>0; j--) {     
	Contig & minor = scaff[j];
	//cout << "j " << j << endl;
	if (minor.isize() == 0)
	  continue;
	for (k=j-1; k>=0; k--) {
	  //cout << "k " << k << endl;
	  const Contig & major = scaff[k];
	  if (major.isize() == 0)
	    continue;
	  if (DoRemove(major, minor)) {
	    minor.clear();
	    break;
	  }
	}
	
      }
      CleanScaffold(scaff);
    }

    io.Write(assembly, outName);

    return 0;
}

