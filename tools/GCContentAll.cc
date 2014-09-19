#include <string>
#include "base/CommandLineParser.h"
#include "src/DNAVector.h"

//#include "analysis/SlideSearch.h"

int main( int argc, char** argv )
{

  commandArg<string> leftCmmd("-i","input fasta file");
  //commandArg<string> rightCmmd("-r","right input fasta file", "");
  //commandArg<string> primerCmmd("-p","primer input fasta file");
  //commandArg<int> misCmmd("-mis","maximum number of mismatches", 3);
  commandLineParser P(argc,argv);
  P.SetDescription("Prints out all reads and mates that hit the primers.");
  P.registerArg(leftCmmd);
  //P.registerArg(rightCmmd);
  //P.registerArg(primerCmmd);
  //P.registerArg(misCmmd);
  
  P.parse();
 
  vecDNAVector left, right;
 
  //primer.Read(P.GetStringValueFor(primerCmmd));
  left.Read(P.GetStringValueFor(leftCmmd));
  //right.Read(P.GetStringValueFor(rightCmmd));
  
  int i, j;
  double n = 0.;
  double cg = 0.;
  for (i=0; i<left.isize(); i++) {
    const DNAVector & a = left[i];
 
    double div = a.isize();
    for (j=0; j<a.isize(); j++) {
      if (a[j] == 'C' || a[j] == 'G') {
	cg += 1.;
	n += 1.;
      }
     if (a[j] == 'A' || a[j] == 'T') {
	n += 1.;
      }
    }
    //cout << left.Name(i) << "\t" << cg / div << endl;
  
  }
  cout << "GC content: " << 100 * cg/n << " %" << endl;
  return 0;
}
