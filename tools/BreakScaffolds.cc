#define FORCE_DEBUG

#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"
#include "base/StringUtil.h"

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> fastaCmmd("-f","fasta input file");
  commandArg<string> outCmmd("-o","fasta output file");
  commandLineParser P(argc,argv);
  P.SetDescription("Breaks scaffolds.");
  P.registerArg(fileCmmd);
  P.registerArg(fastaCmmd);
  P.registerArg(outCmmd);
  
  P.parse();
  



  string fileName = P.GetStringValueFor(fileCmmd);
  string fastaName = P.GetStringValueFor(fastaCmmd);
  string outName = P.GetStringValueFor(outCmmd);
  
  vecDNAVector dna, out;
  dna.Read(fastaName);

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  int i;
  svec<int> offset;
  offset.resize(dna.isize(), 0);
  int k = 0;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    
    int c1 = parser.AsInt(7);
    int c2 = parser.AsInt(9);

    if (c1 < 15 || c2 < 15)
      continue;

    int index = dna.NameIndex(parser.AsString(0));
    DNAVector & tmp = dna[index];
    int off = offset[index];

    int from = parser.AsInt(2) - off;
    int to = parser.AsInt(4) - off;
    
    bool bYes = false;
    cout << parser.Line() << endl;
    if (parser.AsInt(5) < 0) {
      bYes = true;
      cout << "Overlap is negative!" << endl;
    }

    if (tmp.isize() == 0)
      continue;
    int n = 0;
    for (i=from; i<to; i++) {
      if (tmp[i] == 'N')
	n++;
    }
    cout << "N's: " << n << endl;
    if (n > 1)
      bYes = true;

    if (bYes) {
      cout << "Breaking." << endl;
      int m = (from + to)/2;
      DNAVector a, b;
      cout << "First part... " << m << endl;
      a.SetToSubOf(tmp, 0, m);
      cout << "Second part... " << tmp.isize() << " " << endl;
      b.SetToSubOf(tmp, m, tmp.isize()-m);
      string nA = ">" + parser.AsString(0) + "_part_" + Stringify(k);
      k++;
      //string nB = ">" + parser.AsString(0) + "_b";
      cout << "Push and done." << endl;
      out.push_back(a, nA);
      //out.push_back(b, nB);
      tmp = b;
      offset[index] += a.isize();
    }
  }

  for (i=0; i<dna.isize(); i++) {
    const DNAVector & d = dna[i];
    if (d.isize() == 0)
      continue;
    out.push_back(d, dna.Name(i));
  }
  
  out.Write(outName);

  return 0;
}
