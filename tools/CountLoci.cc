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

void Strip(string &s) {

  char tmp[512];
  const char * p = s.c_str();
  strcpy(tmp, &p[1]);
  if (tmp[strlen(tmp)-1] == ';')
    tmp[strlen(tmp)-1] = 0;
  if (tmp[strlen(tmp)-1] == '"')
    tmp[strlen(tmp)-1] = 0;
  s = tmp;
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","GTF file");
  commandArg<string> listCmmd("-l","list file");
  commandLineParser P(argc,argv);
  P.SetDescription("Prints loci.");
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

 
  parser.ParseLine();

  //cout << "Read " << names.isize() << " entries." << endl;
  int i;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;

    string trans, gene;

    for (int i=0; i<parser.GetItemCount(); i++) {
      if (parser.AsString(i) == "transcript_id")
	trans = parser.AsString(i+1);
    }
    for (int i=0; i<parser.GetItemCount(); i++) {
      if (parser.AsString(i) == "gene_id")
	gene = parser.AsString(i+1);
    }
    Strip(trans);
    Strip(gene);
    //cout << trans << " " << gene << endl;
    int index = BinSearch(names, trans);
    if (index < 0)
      continue;

    cout << gene << endl;
  } 
  return 0;
}
