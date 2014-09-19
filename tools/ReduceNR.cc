#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  // commandArg<int> bCmmd("-from","from column");
  //commandArg<int> eCmmd("-to","to column");
  //commandArg<bool> nCmmd("-newline","add newline", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Parses the fasta headers of the nr database");
  P.registerArg(fileCmmd);
  //P.registerArg(bCmmd);
  //P.registerArg(eCmmd);
  //P.registerArg(nCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  //int from = P.GetIntValueFor(bCmmd);
  //int to = P.GetIntValueFor(eCmmd);
  //bool bN = P.GetBoolValueFor(nCmmd);
 

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  char tmp[64];
  tmp[0] = 1;
  tmp[1] = 0;
  string del = tmp;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.AsString(0)[0] == '>') {
      //StringParser p;
      //p.SetLine(parser.Line(), del);
      //char name[4096];
      
      //char dd[4096];
      //strcpy(dd, parser.Line().c_str());
      for (int i=0; i<parser.Line().size(); i++) {
	if (parser.Line()[i] == 1)
	  break;
	cout <<  parser.Line()[i];
      }
      cout << endl;
      
      //cout << name << endl;
    } else {
      cout << parser.Line() << endl;
    }
  }
  return 0;
}
