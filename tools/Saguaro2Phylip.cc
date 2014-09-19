#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","cactus input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Makes phyloP compatible distance matrices.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  
  FILE * p = NULL;
  int i;
  int k = 0;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.GetItemCount() == 1) {
      if (p != NULL)
	fclose(p);
      p = fopen(parser.AsString(0).c_str(), "w");
      parser.ParseLine();
      fprintf(p, "%d\n", parser.GetItemCount());
      k = 0;
      continue;
    }
    fprintf(p, "%s", parser.AsString(0).c_str());
    for (i= strlen(parser.AsString(0).c_str()); i<10; i++)
      fprintf(p, " ");
    for (i=1; i<parser.GetItemCount(); i++) {
      if (i-1 == k)
	fprintf(p, "0.000  ");
      else
	fprintf(p, "%1.3f ", parser.AsFloat(i));
    }
    fprintf(p, "\n");
    k++;
  }
  fclose(p);
  return 0;
}
