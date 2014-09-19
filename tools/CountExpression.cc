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

  svec<string> labels;
  parser.ParseLine();

  int i;
  for (i=0; i<parser.GetItemCount(); i++) {
    labels.push_back(parser.AsString(i));
  }

  int p = 0;
  int d = 0;
  int b = 0;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    bool polyA = false;
    bool dsn = false;
    for (i=0; i<parser.GetItemCount(); i++) {
      if (strstr(labels[i].c_str(), "dsn_FPKM") != NULL) {
	if (parser.AsFloat(i) > 0.1)
	  dsn = true;
      }
      if (strstr(labels[i].c_str(), "polyA_FPKM") != NULL) {
	if (parser.AsFloat(i) > 0.1)
	  polyA = true;
      }
    }
    // if (polyA && dsn)
    // cout << parser.AsString(4) << endl;

    if (polyA && dsn) {
      b++;
    } else {
      if (polyA)
	p++;
      if (dsn)
	d++;
    }
  }

  cout << "polyA: " << p << " dsn: " << d << " both: " << b << endl;

  return 0;
}
