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

  commandArg<string> fileCmmd("-i","GTF file");
  commandArg<string> listCmmd("-l","list file");
  commandArg<bool> keepCmmd("-k","keep, not dump", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Prints FPKM stats.");
  P.registerArg(fileCmmd);
  P.registerArg(listCmmd);
  P.registerArg(keepCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string listName = P.GetStringValueFor(listCmmd);
  bool bKeep = P.GetBoolValueFor(keepCmmd);

  svec<string> names;
  Read(names, listName);

  //comment. ???
  FlatFileParser parser;

  int i;
  
  parser.Open(fileName);

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;

    bool b = false;

    for (i=0; i<parser.GetItemCount(); i++) {
       if (parser.AsString(i) == "transcript_id") {
	
	char tmp[1024];
	const char * p = parser.AsString(i+1).c_str();
	strcpy(tmp, &p[1]);
	if (tmp[strlen(tmp)-1] == ';') 
	  tmp[strlen(tmp)-1] = 0;
	if (tmp[strlen(tmp)-1] == '"') 
	  tmp[strlen(tmp)-1] = 0;
	string n = tmp;
	int index = BinSearch(names, n);
	if (index >= 0) {
	  b = true;
	}
	break;
      }
    }	
    if (bKeep) {
      if (b)
	cout << parser.Line() << endl;
    } else {
      if (!b)
	cout << parser.Line() << endl;
    }
  }
}
