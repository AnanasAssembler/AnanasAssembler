#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"

#include "src/DNAVector.h"
#include "util/mutil.h"


void Coords(int & start, int & stop, const string & s)
{
  CMTokenizer t;
  t.AddDelimiter("[");
  t.AddDelimiter("-");
  t.AddDelimiter("]");
  const char * p = s.c_str();
  CMPtrStringList out;
  t.Tokenize(out, p);

  start = atol(*out(0));
  stop = atol(*out(1));
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","transcript fasta file");
  commandArg<string> refCmmd("-r","reference fasta file");
  commandArg<string> pCmmd("-p","papaya file");
  commandLineParser P(argc,argv);
  P.SetDescription("Prints out poly-A tails.");
  P.registerArg(fileCmmd);
  P.registerArg(pCmmd);
  P.registerArg(refCmmd);
  
  P.parse();
  
  string papaya = P.GetStringValueFor(pCmmd);
  string fileName = P.GetStringValueFor(fileCmmd);
  string refName = P.GetStringValueFor(refCmmd);
  
  vecDNAVector dna;
  
  dna.Read(fileName);

  vecDNAVector ref;
  
  ref.Read(refName);

  FlatFileParser parser;
  
  parser.Open(papaya);

  svec<string> best;
  best.resize(ref.isize());
  svec<int> count;
  count.resize(ref.isize(), 0);
  svec<int> refBest;
  refBest.resize(ref.isize(), 0);
  svec<int> refBestStart;
  refBestStart.resize(ref.isize(), 0);
  svec<int> txBest;
  txBest.resize(ref.isize(), 0);

  int i;
  parser.ParseLine();
  parser.ParseLine();
  parser.ParseLine();
  parser.ParseLine();
  parser.ParseLine();

  while (parser.ParseLine()) {
    int refStart = parser.AsInt(15);
    int refStop = parser.AsInt(16);
    const string & refName = parser.AsString(13);
    //Coords(refStart, refStop, parser.AsString(1));


    int txStart = parser.AsInt(11);
    int txStop = parser.AsInt(12);
    const string & txName = parser.AsString(9);
    //Coords(txStart, txStop, parser.AsString(4));

    const string & ori = parser.AsString(8);

    if (ori == "-")
      continue;

    int refIdx = ref.NameIndex(refName);
    //const DNAVector & t = dna(txName);
    //cout << refName << " " << refIdx << endl;
    if (refIdx == -1)
      break;


    if (refStop - refStart > count[refIdx]) {
      if (refStop - refStart > 200 
	  && refStop - refStart > count[refIdx]
	  && refStop > ref[refIdx].isize() - 1000
	  && refStart < 1200) {
	count[refIdx] = refStop - refStart;
	best[refIdx] = txName;
	refBest[refIdx] = refStop;
	refBestStart[refIdx] = refStart;
	txBest[refIdx] = txStop;
      }
    }

  }
  
  int k = 0;


  for (i=0; i<ref.isize(); i++) {
    const DNAVector & d = ref[i];
    if (best[i] == "")
      continue;


    //if (d.isize() - txBest[i] <= 1000) {
    const DNAVector & t = dna(best[i]);
    if (t.isize() - txBest[i] > 3 && t.isize() - txBest[i] < 50) {
      cout << ">tail_" << k << "_ext_" << 1000 - (d.isize() - refBest[i]);
      cout << "_" << best[i] << "_" << ref.NameClean(i) << endl;
      k++;
      for (int x = txBest[i]-40; x<t.isize(); x++)
	cout << t[x];
      cout << endl;
    
    }
  }

  /*
  k = 0;
  for (i=0; i<ref.isize(); i++) {
    const DNAVector & d = ref[i];
    if (best[i] == "")
      continue;
    
    const DNAVector & t = dna(best[i]);
    if (refBestStart[i] < 1000) {
      int len = 1000 - refBestStart[i];
      cout << ">head_" << k << "_ext_" << len;
      cout << "_" << best[i] << "_" << ref.NameClean(i) << endl;
      k++;
      for (int x = 0; x<len+40; x++)
	cout << t[x];
      cout << endl;
          
    }
    }*/

  return 0;
}
