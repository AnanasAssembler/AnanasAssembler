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
 
  svec<Info> loc1;
  svec<Info> loc2;

  Read(loc1, specName1);
  Read(loc2, specName2);
  
  for (int i=0; i<loc1.isize(); i++) {
    int index = BinSearch(loc2, Info(loc1[i].A(), ""));
    if (index < 0)
      continue;
    cout << loc1[i].A() << endl;
    

  }

 
  return 0;



}
