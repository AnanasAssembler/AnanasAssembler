#include <string>
#include "base/CommandLineParser.h"
#include "src/DNAVector.h"
#include "base/RandomStuff.h"


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<int> nCmmd("-n", "size");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);
  P.registerArg(nCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  int n = P.GetIntValueFor(nCmmd);
 
  vecDNAVector dna;
  dna.Read(fileName, false, false, false);


  int howmany = 10000;

  DNAVector & d = dna[0];
  int i;

  vecDNAVector a, b;

  a.resize(howmany);
  b.resize(howmany);

  while (howmany > 0) {
    int s = RandomInt(d.isize()-n);
    bool bGood = true;
    for (i=s; i<s+n; i++) {
      if (d[i] == 'N' || d[i] >= 'a') {
	bGood = false;
      }
    }
    if (!bGood)
      continue;
    howmany--;
    a[howmany].SetToSubOf(d, s, n);
    b[howmany] = a[howmany];
    b[howmany].ReverseComplement();

  }

  char out[256];
  sprintf(out, "FDRTestFW%d.fasta", n);
  a.Write(out);
  sprintf(out, "FDRTestRC%d.fasta", n);
  b.Write(out);
  
  return 0;
}
