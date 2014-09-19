#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<bool> leftCmmd("-l","left", false);
  commandArg<bool> rightCmmd("-r","right", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Makes a left & right file.");
  P.registerArg(fileCmmd);
  P.registerArg(leftCmmd);
  P.registerArg(rightCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  bool l = P.GetBoolValueFor(leftCmmd);
  bool r = P.GetBoolValueFor(rightCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (l) {
      cout << parser.Line() << endl;
      parser.ParseLine();
      cout << parser.Line() << endl;
      parser.ParseLine();
      parser.ParseLine();
    }
    if (r) {
      parser.ParseLine();
      parser.ParseLine();
      cout << parser.Line() << endl;
      parser.ParseLine();
      cout << parser.Line() << endl;      
    }
 }
  return 0;
}
