#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"
void Fill(char * p) {
  int l = strlen(p);
  for (int i=0; i<l; i++) {
    if (p[i] == ' ')
      p[i] = '0';
  }
}


int main( int argc, char** argv )
{

  commandArg<string> fastaCmmd("-f","input fasta file");
  commandArg<string> fileCmmd("-i","input GTF file");
  commandLineParser P(argc,argv);
  P.SetDescription("Prints a non-redundant fasta file from a GTF.");
  P.registerArg(fileCmmd);
  P.registerArg(fastaCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string fastaName = P.GetStringValueFor(fastaCmmd);
  
  vecDNAVector dna;
  dna.Read(fastaName);

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  int i, j;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    const string & chr = parser.AsString(0);
    int start = parser.AsInt(3);
    int stop  = parser.AsInt(4);

    DNAVector & d = dna(chr);
    for (i=start; i<=stop; i++) {
      if (i<1)
	continue;
      if (i>=d.isize()-1)
	continue;
      //cout << "Inc qual " << i << endl;
      if (d[i] != 'N')
	d.SetQual(i, 1);
    }
  }
  int x;
  int count = 0;
  for (i=0; i<dna.isize(); i++) {
    const DNAVector & d = dna[i];
    for (j=1; j<d.isize(); j++) {
      if (d.Qual(j-1) == 0 && d.Qual(j) > 0) {
	for (x = j; x<d.isize()-1; x++) {
	  if (d.Qual(x) > 0 && d.Qual(x+1) == 0) {
	    break;
	  }
	}
	char name[512];
	sprintf(name, "RUID%8d", count);
	Fill(name);
	count++;
	cout << dna.Name(i) << "\t" << j << "\t" << x << "\t" << name << endl;
	/*cout << dna.Name(i) << ":" << j << "-" << x << endl;
	for (int y=j; y<x; y++) {
	  cout << d[y];
	}
	cout << endl;*/
      }

    }
  }

  return 0;
}
