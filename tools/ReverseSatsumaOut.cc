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

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    cout << parser.AsString(3) << "\t";
    cout << parser.AsInt(4) << "\t";
    cout << parser.AsInt(5) << "\t";
    cout << parser.AsString(0) << "\t";
    cout << parser.AsInt(1) << "\t";
    cout << parser.AsInt(2) << "\t";
    cout << parser.AsString(6) << "\t";
    cout << parser.AsString(7) << endl;
  }
  
  return 0;
}
