#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"


int main( int argc, char** argv )
{

  commandArg<string> inCmmd("-i","input fasta file");
  commandArg<string> outCmmd("-o","output fasta file");
  commandArg<string> satCmmd("-s","SatsumaProt file");
  commandLineParser P(argc,argv);
  P.SetDescription("Re-assigns protein names based on SatsumaProt alignments.");
  P.registerArg(inCmmd);
  P.registerArg(outCmmd);
  P.registerArg(satCmmd);
   
  P.parse();
  
  string fileName = P.GetStringValueFor(satCmmd);
  string in = P.GetStringValueFor(inCmmd);
  string out = P.GetStringValueFor(outCmmd);
  

  vecDNAVector dna;
  dna.Read(in);

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  svec<string> names;
  svec<double> scores;

  names.resize(dna.isize());
  scores.resize(dna.isize(), 10.);

  double lastIdent = 0;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 4 && parser.AsString(0) == "Identity:") {
      lastIdent = parser.AsFloat(1);
      continue;
    }

    if (parser.GetItemCount() < 6)
      continue;
    if (parser.AsString(0) != "Summary")
      continue;
    int off = 0;
    if (parser.GetItemCount() == 6)
      off = 1;
    const string & seq = parser.AsString(1);
    string frame = "frame:(n/a)";

    if (off == 0)
      frame = parser.AsString(2);
    
    const string & name = parser.AsString(4-off);
    double score = parser.AsFloat(6-off);
    
    cout << "Score: " << score << endl;

    //if (score < 150.)
    //continue;

    int index = dna.NameIndex(seq);
    //if (lastIdent  scores[index]) {
    if (score < scores[index]) {
      scores[index] = score;
      names[index] = seq + " " + name + " " + frame;
    }
      
    lastIdent = 10.;
  }

  for (int i=0; i<dna.isize(); i++) {
    if (names[i] != "") {
      dna.SetName(i, ">" + names[i]);
    }
  }

  dna.Write(out);

  return 0;
}
