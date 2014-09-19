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
    string cons, ens;
    for (int i=0; i<parser.GetItemCount(); i++) {
      if (parser.AsString(i) == "transcript_id")
	cons = parser.AsString(i+1);
     if (parser.AsString(i) == "gene_id")
	ens = parser.AsString(i+1);

    }
     //if (ens != "" && strstr(ens.c_str(), "ENS") != NULL) 
    cout << cons << "\t" << ens << endl;
    

  }
  return 0;
}
