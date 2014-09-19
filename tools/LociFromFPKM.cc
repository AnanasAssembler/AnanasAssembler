#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "math/Spearman.h"

void Read(svec<string> & names, const string & fileName)
{
  FlatFileParser parser;
  
  parser.Open(fileName);

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    names.push_back(parser.AsString(0));
  }
  Sort(names);
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","FPKM file");
  commandArg<string> listCmmd("-l","list file");
  commandLineParser P(argc,argv);
  P.SetDescription("Prints FPKM stats.");
  P.registerArg(fileCmmd);
  P.registerArg(listCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string listName = P.GetStringValueFor(listCmmd);
  
  svec<string> names;
  Read(names, listName);

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  svec<double> polyA_all, dsn_all;
  int i;
  int s = 4;
  
  double polyA_spec = 0.;
  double dsn_spec = 0.;
  double nn = 0;

  parser.ParseLine();

  //cout << "Read " << names.isize() << " entries." << endl;

  svec<string> label;
  for (i=0; i<parser.GetItemCount(); i++) {
    label.push_back(parser.AsString(i));
  }

  svec<double> a_brainA, a_brainA2, a_brainD, a_kidneyA, a_kidneyA2, a_kidneyD;
  
  double ex_brain = 0.;
  double ex_brain_div = 0.;

  double polyA_count = 0.;
  double dsn_count = 0.;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    const string & n = parser.AsString(0);
    int index = BinSearch(names, n);
    if (index < 0)
      continue;

    cout << parser.AsString(3) << endl;
  } 
  return 0;
}
