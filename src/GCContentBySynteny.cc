#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","satsuma file");
  commandArg<string> tCmmd("-t","target fasta");
  commandArg<string> qCmmd("-q","query fasta");
  commandLineParser P(argc,argv);
  P.SetDescription("Computes GC content via Satsuma alignments.");
  P.registerArg(fileCmmd);
  P.registerArg(tCmmd);
  P.registerArg(qCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  

  vecDNAVector target, query;
  target.Read(P.GetStringValueFor(tCmmd));
  query.Read(P.GetStringValueFor(qCmmd));

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  int i;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    const DNAVector & t = target(parser.AsString(0));
    const DNAVector & q = query(parser.AsString(3));
    
    int tStart = parser.AsInt(1);
    int tStop = parser.AsInt(2);
    int qStart = parser.AsInt(4);
    int qStop = parser.AsInt(5);

    int gc = 0;
    int at = 0;
    int cpg = 0;
    for (i=tStart+1; i<tStop; i++) {
      if (t[i] == 'C' || t[i] == 'G')
	gc++;
      if (t[i] == 'A' || t[i] == 'T')
	at++;
      if (t[i-1] == 'C' && t[i] == 'G')
	cpg++;
    }

    cout << (double)gc/(double)(gc + at) << "\t";
    cout << (double)cpg/(double)(gc + at) << "\t";

    gc = 0;
    at = 0;
    cpg = 0;
    for (i=qStart+1; i<qStop; i++) {
      if (q[i] == 'C' || q[i] == 'G')
	gc++;
      if (q[i] == 'A' || q[i] == 'T')
	at++;
      if (q[i-1] == 'C' && q[i] == 'G')
	cpg++;
    }

    cout << (double)gc/(double)(gc + at) << "\t";
    cout << (double)cpg/(double)(gc + at) << endl;
   

  }
  return 0;
}
