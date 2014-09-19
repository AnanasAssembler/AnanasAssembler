#include <string>
#include "base/CommandLineParser.h"
#include "base/StreamParser.h"


int main( int argc, char** argv )
{

  stringstream ss;
  ss << "This is a number " << 6 << endl;
  ss << "This is a second number " << 12 << endl;

 
  /*  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the stream parser.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  */

  
  //comment. ???
  StreamParser parser;
  parser.Set(ss);


  while (parser.ParseLine()) {
    cout << "Line: " << parser.Line() << endl;
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
