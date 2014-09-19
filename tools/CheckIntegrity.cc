#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  string lastT, lastQ;
  int last;
  int lastEnd = 0;
  int same = 0;
  int next = 0;
  bool bBreak = false;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    const string & q = parser.AsString(0);
    const string & t = parser.AsString(3);
    int pos = parser.AsInt(1);
    if (q != lastQ) {
      lastQ = q;
      lastT = t;
      last = pos;
      lastEnd = parser.AsInt(2);
      if (bBreak)
	cout << "\tafter: " << same << endl; 
      same = 0;
      bBreak = false;
      continue;
    }
    if (t != lastT) {
      if (bBreak)
	cout << "\tafter: " << same << endl; 	
      //cout << "\tbefore: " << same << "\t"; 
      cout << q << "\t" << lastT << "\t" << last << "\t" << t << "\t" << pos;
      cout << "\t" << pos - lastEnd;
      cout << "\tbefore: " << same; 
      bBreak = true;
      same = 0;;
    } else {
      same++;
    }
    lastQ = q;
    lastT = t;
    last = pos;
    lastEnd = parser.AsInt(2);
  }
  cout << "\tafter: " << same << endl; 
  return 0;
}
