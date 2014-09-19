#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> fastaCmmd("-f","fasta file");
  commandArg<string> headCmmd ("-header","fasta header", "");
  commandLineParser P(argc,argv);
  P.SetDescription("Prints out top hits.");
  P.registerArg(fileCmmd);
  P.registerArg(fastaCmmd);
  P.registerArg(headCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string fastaName = P.GetStringValueFor(fastaCmmd);
  string head = P.GetStringValueFor(headCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  vecDNAVector dna;
  dna.Read(fastaName);

  string last;
  int k = 0;
  svec<string> found;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;

    char tmp[256];
    strcpy(tmp, parser.AsString(0).c_str());
    tmp[strlen(tmp)-4] = 0;
    if (last == tmp) {
      continue;
    }
    last = tmp;

    if (parser.AsInt(2) - parser.AsInt(1) < 150)
      continue;
    if (parser.AsFloat(6) < 0.54)
      continue;
    DNAVector d = dna(parser.AsString(0));
    if (d.isize() < 400)
      continue;
    
    bool b = false;
    for (int x=0; x<found.isize(); x++) {
      if (found[x] == parser.AsString(0)) {
	b = true;
	break;
      }
    }
    if (b)
      continue;

    found.push_back(parser.AsString(0));

    if (parser.AsString(7) == "-")
      d.ReverseComplement();
    if (head == "") {
      cout << ">" << parser.AsString(0) << endl;
    } else {
      cout << ">" << head << k << endl;
    }
    k++;
    for (int i=0; i<d.isize(); i++)
      cout << d[i];
    cout << endl;
      
  }
  return 0;
}
