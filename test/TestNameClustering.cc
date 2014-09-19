#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"

#include "src/Smorgas/ProtNameCluster.h"


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the protein name clustering.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  ProteinNameClusterer cl;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0) {
      continue;
    }
    if (parser.AsString(0) != "Summary")
      continue;
    
    const string & s = parser.AsString(4);
    cl.Add(s);
  }

  cl.Cluster();
  cout << "Printing clusters:" << endl;
  for (int i=0; i<cl.Num(); i++) {
    cout << cl.Class(i) << "\t" << cl.Count(i) << endl;
  }

  return 0;
}
