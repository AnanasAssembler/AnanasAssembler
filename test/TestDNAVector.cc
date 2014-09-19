#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
 
  vecDNAVector dna;
  dna.Read(fileName);

  int i;
  for (i=0; i<dna.isize(); i++) {
    cout << dna.Name(i) << endl;
    const DNAVector & d = dna[i];
    for (int j=0; j<d.isize(); j++)
      cout << d[j];
    cout << endl;
  }


  return 0;
}
