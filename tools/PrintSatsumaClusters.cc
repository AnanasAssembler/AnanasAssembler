#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"

bool Spring(const string & s) {
  char tmp[1024];
  strcpy(tmp, s.c_str());
  if (strstr(tmp, "sample03") != NULL)
    return true;
  if (strstr(tmp, "sample04") != NULL)
    return true;
  if (strstr(tmp, "sample05") != NULL)
    return true;
  if (strstr(tmp, "sample09") != NULL)
    return true;
  if (strstr(tmp, "sample12") != NULL)
  return true;
  if (strstr(tmp, "sample15") != NULL)
    return true;

  return false;
}

bool Fall(const string & s) {
  char tmp[1024];
  strcpy(tmp, s.c_str());
  if (strstr(tmp, "sample20") != NULL)
    return true;
  if (strstr(tmp, "sample21") != NULL)
    return true;
  if (strstr(tmp, "sample22") != NULL)
    return true;
  if (strstr(tmp, "sample23") != NULL)
    return true;
  if (strstr(tmp, "sample09") != NULL)
    return true;
  if (strstr(tmp, "sample12") != NULL)
    return true;
  if (strstr(tmp, "sample15") != NULL)
    return true;
  //if (strstr(tmp, "sample12") != NULL)
  //return true;

  return false;
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","satsuma input file");
  commandArg<string> listCmmd("-f","fasta file");
  commandLineParser P(argc,argv);
  P.SetDescription("Clusters sequences by sequence similarity.");
  P.registerArg(fileCmmd);
  P.registerArg(listCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string listName = P.GetStringValueFor(listCmmd);
  
  vecDNAVector dna;
  dna.Read(listName);

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  int k = 0;
  int i;
  vecDNAVector out;
  int spring = 0;
  int fall = 0;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.AsString(0) == "Cluster:") {
      continue;
    }
    if (parser.AsString(0) == "Members:") {
      char tmp[1024];
      char tmpTop[1024];
      sprintf(tmp, "cluster_%d_members_%d.fasta", k, out.isize());
      sprintf(tmpTop, "top_%d_members_%d.fasta", k, out.isize());
            //cout << out.isize() << " " << spring << " " << fall << endl;
      if (out.isize() > 3 /*&& out.isize() == spring*/) {
	int x = 0;
	int l = 0;
	for (int y=0; y<out.isize(); y++) {
	  if (out[y].isize() > l) {
	    l = out[y].isize();
	    x = y;
	  }
	}
	out.Write(tmp);
	if (x > 0) {
	  out[0] = out[x];
	  out.SetName(0, out.Name(x));
	}
	out.resize(1);
	out.Write(tmpTop);
      }
      k++;
      spring = fall = 0;
      out.resize(0);
      continue;
    }
    const DNAVector & d = dna(parser.AsString(0));
    int index = dna.NameIndex(parser.AsString(0));
    out.push_back(d, dna.Name(index));
    if (Spring(parser.AsString(0)))
      spring++;
    if (Fall(parser.AsString(0)))
      fall++;
  }

  return 0;
}
