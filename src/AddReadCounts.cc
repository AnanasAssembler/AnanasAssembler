#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/FileParser.h"



string StripNumber(const string & s) {
  int i;
  //char tmp[512];
  //strcpy(tmp, s.c_str());
  int k = 0;
  string ret = s;
  for (i=s.length()-1; i>=0; i--) {
    if (s[i] == '_')
      k++;
    if (k == 2) {
      ret = &s[i+1];
    }
  } 
  return ret;
}


class Info
{
public:
  string id;
  int r;
  int p;

  bool operator < (const Info & f) const {
    return id < f.id;
  }

};


void ReadInfo(svec<Info> & f, const string & fileName)
{    
  FlatFileParser parser;
  
  parser.Open(fileName);

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    Info info;
    info.id = StripNumber(parser.AsString(0));
    info.r = parser.AsInt(1);
    info.p = parser.AsInt(2);
    f.push_back(info);
  }
  Sort(f);
}


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> countCmmd("-s","count summary file");
  commandArg<int> bCmmd("-c","column containing the contig ids");
  commandArg<int> eCmmd("-o","which column to write to");
  commandArg<bool> rCmmd("-r","use read counts, not pair counts", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Adds read or pair counts.");
  P.registerArg(fileCmmd);
  P.registerArg(countCmmd);
  P.registerArg(bCmmd);
  P.registerArg(eCmmd);
  P.registerArg(rCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string countName = P.GetStringValueFor(countCmmd);
  int ids = P.GetIntValueFor(bCmmd);
  int to = P.GetIntValueFor(eCmmd);
  bool bReads = P.GetBoolValueFor(rCmmd);
  
  svec<Info> info;
  ReadInfo(info, countName);

  //comment. ???
  FlatFileParser parser1;
  
  parser1.Open(fileName);

  while (parser1.ParseLine()) {
    if (parser1.GetItemCount() == 0)
      continue;

    StringParser parser;
    parser.SetLine(parser1.Line(), "\t");
    
    string key = StripNumber(parser.AsString(ids));
    Info tmp;
    tmp.id = key;
    int index = BinSearch(info, tmp);
    if (index < 0) {
      cout << "ERROR: Not found: " << parser.AsString(ids) << " " << key << endl;
      exit(-1);
    }

    int n = info[index].p;
    if (bReads)
      n = info[index].r;

    string sep = "";
    for (int i=0; i<parser.GetItemCount(); i++) {
      if (i == to) {
	cout << sep << n;
	sep = "\t";
      }
      cout << sep << parser.AsString(i);
      sep = "\t";
    }
    cout << endl;
  }
  return 0;
}
