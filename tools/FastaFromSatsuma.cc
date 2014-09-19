#define FORCE_DEBUG

#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","satsuma summary input file");
  commandArg<string> fastaCmmd("-f","fasta file");
  commandArg<string> outCmmd("-o","output file");
  commandArg<int> minCmmd("-min","minimum alignment length", 300);
  commandLineParser P(argc,argv);
  P.SetDescription("Prints aligining sequences from fasta.");
  P.registerArg(fileCmmd);
  P.registerArg(fastaCmmd);
  P.registerArg(outCmmd);
  P.registerArg(minCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string fastaName = P.GetStringValueFor(fastaCmmd);
  string outName = P.GetStringValueFor(outCmmd);
  int min = P.GetIntValueFor(minCmmd);
  
  
  vecDNAVector dna;
  dna.Read(fastaName);

  //comment. ???
  FlatFileParser parser;
  
  //int min = 250;

  parser.Open(fileName);
  svec<int> fw, rc;
  fw.resize(dna.isize(), 0);
  rc.resize(dna.isize(), 0);
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;

    const string & n = parser.AsString(0);
    int index = dna.NameIndex(n);
    int len = parser.AsInt(2) - parser.AsInt(1);
    if (parser.AsString(7) == "-") {
      if (len > rc[index])
	rc[index] = len;      
    } else {
     if (len > fw[index])
	fw[index] = len;      
    }  
    
  }

  vecDNAVector out;
  for (int i=0; i<dna.isize(); i++) {
    if (fw[i] < min && rc[i] < min)
      continue;
    DNAVector tmp = dna[i];
    if (rc[i] > fw[i])
      tmp.ReverseComplement();

    out.push_back(tmp, dna.Name(i));
  }
  out.Write(outName);

  return 0;
}
