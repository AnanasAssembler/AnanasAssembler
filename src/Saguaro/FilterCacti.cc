#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/Saguaro/Cactus.h"

#define FORCE_DEBUG

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","LocalTree file");
  commandArg<string> cCmmd("-c","cacti file");
  commandLineParser P(argc,argv);
  P.SetDescription("Eliminates cacti that are not used.");
  P.registerArg(fileCmmd);
  P.registerArg(cCmmd);
  
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string cName = P.GetStringValueFor(cCmmd);
  
  int i, j;
  svec<Cactus> cacti;
  LoadCacti(cacti, cName);

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  bool b = false;
  svec<string> used;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    
    if (parser.AsString(0) == "REPORTING")
      b = true;

    if (!b)
      continue;
    if (parser.GetItemCount() > 6 && parser.AsString(5) == "length:") {
      used.push_back(parser.AsString(0));
    }    
  }
  //cout << "Added" << endl;
  UniqueSort(used);

  for (i=0; i<cacti.isize(); i++) {
    const string & name = cacti[i].Name();
    for (j=0; j<used.isize(); j++) {
      if (used[j] == name) {
	cacti[i].Print();
	break;
      }
    }
  }
  

  return 0;
}
