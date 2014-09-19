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
    for (int i=0; i<parser.GetItemCount(); i++) {
      if (strstr(parser.AsString(i).c_str(), "ENS") != NULL)
	cout << parser.AsString(i) << endl;
    }
  }
  return 0;
}
