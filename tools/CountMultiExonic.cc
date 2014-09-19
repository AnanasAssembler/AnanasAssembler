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


  cout << "Read " << names.isize() << " entries." << endl;
  string last;
  int spliced = 0;
  int total = 0;
  int count = 0;
  int max = 0;
  string ID;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    string n = parser.AsString(11);
    Strip(n);
    int index = BinSearch(names, n);
    if (index < 0)
      continue;
    if (n == last) {
      count++;
    } else {
      if (count > 1)
	spliced++;
      if (count > max) {
	max = count;
	ID = n;
      }
      count = 1;
      total++;
    }
    last = n;
  }
  cout << "Total: " << total << " Spliced: " << spliced << " Max: " << max << " ID: " << ID << endl;


  return 0;
}
