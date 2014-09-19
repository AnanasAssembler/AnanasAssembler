#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "base/SVector.h"

class Info
{
public:
  Info(const string & a, const string & b) {
    m_a = a;
    m_b = b;
  }

  const string & A() const {return m_a;}
  const string & B() const {return m_b;}


  bool operator < (const Info & f) const {
    return m_a < f.m_a;
  }

private:
  string m_a;
  string m_b;
};

void Read(svec<Info> & out, const string & file) {
  FlatFileParser parser;
  
  parser.Open(file);

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.GetItemCount() == 2) {
      out.push_back(Info(parser.AsString(0), parser.AsString(1)));
    } else
      out.push_back(Info(parser.AsString(0), ""));
  }

  Sort(out);
}

string Parse(const string & s) {
  char tmp[1024];
  strcpy(tmp, s.c_str());
  tmp[0] = ' ';
  tmp[strlen(tmp)-1] = ' ';
  tmp[strlen(tmp)-2] = ' ';
  string o = tmp;
  return o;
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i1","species 1");
  commandArg<string> specCmmd("-i2","species 2");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);
  P.registerArg(specCmmd);
  
  P.parse();
  
  string specName1 = P.GetStringValueFor(fileCmmd);
  string specName2 = P.GetStringValueFor(specCmmd);
 
  svec<Info> proteins;
  svec<Info> mapping;
  svec<Info> loc1;
  svec<Info> loc2;

  string fileName = "map_" + specName1 + "_" + specName2 + ".full";

  cout << "Reading map." << endl;
  Read(mapping, specName1 + ".map");
  cout << "Reading proteins." << endl;
  Read(proteins, specName1 + "_protein_coding");
  //Read(proteins, specName1 + "_singletons.list");
  //Read(proteins, specName1 + "_families.list");
  //Read(proteins, "gene_families/homo_apo.list");

  cout << "Reading loci." << endl;
  Read(loc1, specName1 + "_loci.map");
  cout << "Reading loci." << endl;
  Read(loc2, specName2 + "_loci.map");
  cout << "Parsing." << endl;

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  int n = 0;
  int u = 0;
  int p = 0;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (strstr(parser.AsString(1).c_str(), "INTRON") != NULL ||
	strstr(parser.AsString(1).c_str(), "OTHER") != NULL)
      continue;
    n++;

    string ss = "\"";
    ss += parser.AsString(0);
    ss += "\";";
    int index = BinSearch(mapping, Info(ss, ""));
    if (index < 0) {
      u++;
      continue;
    }
    const string & b = mapping[index].B();
    index = BinSearch(proteins, Info(b, ""));
    if (index >= 0)
      p++;
    if (index >= 0) {
      index = BinSearch(loc1, Info(ss, ""));
      if (index < 0)
	cout << "ERR1" << endl;
      ss = "\"";
      ss += parser.AsString(2);
      ss += "\";";
      cout << "grep " << Parse(loc1[index].B()) << specName1 << "_loci.fpkm ";
      
      index = BinSearch(loc2, Info(ss, ""));
      if (index < 0)
	cout << "ERR2" << endl;
      cout << "grep " << Parse(loc2[index].B()) << specName2 << "_loci.fpkm " << endl;
    }

  }

  cout << "Total: " << n << " Ens: " << n-u << " Prot: " << p << endl;

  return 0;



}
