#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> sCmmd("-s","species file");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);
  P.registerArg(sCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string sName = P.GetStringValueFor(sCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(sName);

  svec<string> species;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    species.push_back(parser.AsString(0));
  }

  int i, j;
  for (i=0; i<species.isize(); i++) {
    for (j=0; j<species.isize(); j++) {
      if (i == j)
	continue;
      string cmd = "grep ^" + species[i];
      cmd += " " + fileName;
      cmd += " | grep " + species[j] + " > tmp";
      cout << cmd << endl;
      int r = system(cmd.c_str());
      cmd = "/home/mangr224/test/ryggrad/ERVMatrix -o tmp.ps -i tmp ";
      cmd += " |sort > tmp2";
      cout << cmd << endl;
      r = system(cmd.c_str());
      cmd = "/home/mangr224/test/ryggrad/ERVMatrix -o tmp.ps -i tmp2 ";
      cmd += " |sort > " + species[i] + "_" + species[j] + ".final";
      cout << cmd << endl;
      r = system(cmd.c_str());
    }
  }


  return 0;
}
