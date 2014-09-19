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
    cout << "Items: " << parser.GetItemCount() << " ";
    for (int i=0; i<parser.GetItemCount(); i++) {
      if (parser.IsInt(i)) {
	cout << " int=";
      } else {
	if (parser.IsFloat(i)) {
	  cout << " float=";
	} else {
	  cout << " string=";
	}
      }
    
      cout << parser.AsString(i);
    }
  
    cout << endl;
  }
  return 0;
}
