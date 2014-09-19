#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Replaces blanks in a file with tabs.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  while (parser.ParseLine()) {
    for (int i=0; i<parser.GetItemCount(); i++) {
      if (i>0)
	cout << "\t";
      cout << parser.AsString(i);
    }
    cout << endl;

  }
  return 0;
}
