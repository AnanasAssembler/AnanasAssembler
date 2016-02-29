#include <string>
#include "base/CommandLineParser.h"
#include "src/DNAVector.h"



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> labelCmmd("-l","label fasta header", "");
  commandLineParser P(argc,argv);
  P.SetDescription("Keeps only the most likely isoform.");
  P.registerArg(fileCmmd); 
  P.registerArg(labelCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string label = P.GetStringValueFor(labelCmmd);
 
  vecDNAVector dna;
  dna.Read(fileName);

  int i, j;

  for (i=0; i<dna.isize(); i++) {
    const string name = dna.Name(i);
    if (name[name.length()-1] != '0')
      continue;
    if (name[name.length()-2] != '0')
      continue;
    if (name[name.length()-3] != '0')
      continue;

    
    if (label == "") {
      cout << name << endl;
    } else {
      StringParser p;
      p.SetLine(name, "_");
      for (j=0; j<p.GetItemCount(); j++) {
	if (j > 0)
	  cout << "_";
	if (j == 1) 
	  cout << label;
	else 
	  cout << p.AsString(j);
      }
      cout << endl;
    }

    for (j=0; j<dna[i].isize(); j++) {
      if (j > 0 && j % 80 == 0)
	cout << endl;
      cout << (dna[i])[j];
    }
    cout << endl;
  }

  return 0;
}
