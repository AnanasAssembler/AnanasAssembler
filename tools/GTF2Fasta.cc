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

string Clean(const string & s) {
  //string tt = s;
  //return tt;
  char tmp[256];
  strcpy(tmp, &s.c_str()[1]);
  tmp[strlen(tmp)-2] = 0;
  string out = tmp;
  return out;
}

int main( int argc, char** argv )
{

  commandArg<string> fastaCmmd("-f","input fasta file");
  commandArg<string> outCmmd("-o","output fasta file");
  commandArg<string> fileCmmd("-i","input GTF file");
  commandLineParser P(argc,argv);
  P.SetDescription("Prints a fasta file from a GTF.");
  P.registerArg(fileCmmd);
  P.registerArg(outCmmd);
  P.registerArg(fastaCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string fastaName = P.GetStringValueFor(fastaCmmd);
  string outName = P.GetStringValueFor(outCmmd);
  
  vecDNAVector dna;
  dna.Read(fastaName);
  vecDNAVector out;

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  int i, j;

  
  string last;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;

    for (i=0; i<parser.GetItemCount(); i++) {
      if (parser.AsString(i) == "transcript_name")
	break;
    }
    string name = ">";
    name += Clean(parser.AsString(i+1));
    const string & chr = parser.AsString(0);


    int start = parser.AsInt(3) - 1;
    int stop  = parser.AsInt(4);

    const DNAVector & d = dna(chr);
    DNAVector tmp;
    tmp.SetToSubOf(d, start, stop-start);
    if (parser.AsString(6) == "-")
      tmp.ReverseComplement();

    out.push_back(tmp, name);

    //for (i=start; i<=stop; i++) {
    //  cout << d[i];
    //}
    //cout << endl;
  }

  out.Write(outName);
  return 0;
}
