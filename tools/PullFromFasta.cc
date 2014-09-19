#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"


int main( int argc, char** argv )
{

  commandArg<string> fastaCmmd("-f","fasta file");
  commandArg<string> fileCmmd("-i","list of names");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
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

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    int index = dna.NameIndex(parser.AsString(0));
    if (index < 0)
      continue;
    const DNAVector & d = dna[index];
    cout << ">" << parser.AsString(0) << endl;
    for (int i=0; i<d.isize(); i++)
      cout << d[i];
    cout << endl;
  }
  return 0;
}
