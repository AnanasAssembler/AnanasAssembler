#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"

void Parse(string & spec, const string & in) 
{
  StringParser a, b;
  a.SetLine(in, "_[");
  //a.SetLine(in, "Tax=");
  if (a.GetItemCount() < 2)
    return;
  //cout << a.AsString(0) << endl;
  string c = a.AsString(1);
  b.SetLine(c, "_");
  spec = b.AsString(0);
}

class Species
{
public:
  Species() {}

  void Read(FlatFileParser & parser) {
    //cout << parser.Line() << endl;
    m_name = parser.AsString(1);
    int i = 3;
    while (strstr(parser.AsString(i).c_str(), ";") == NULL && i < parser.GetItemCount()) {
      //cout << parser.AsString(i) << endl;
      i++;
    }
    if (i >= parser.GetItemCount()-1) {
      m_name = "";
      return;
    }
    m_tax = parser.AsString(i);
    m_tax += " ";
    m_tax += parser.AsString(i+1);
    if (parser.GetItemCount() > 2) {
      m_tax += " ";
      m_tax += parser.AsString(i+2);    
    }
    //cout << m_tax << endl;
    
  }

  const string & Name() {return m_name;}
  const string & Tax()  {return m_tax;}

private:
  string m_name;
  string m_tax;
};

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> outCmmd("-o","output file");
  commandLineParser P(argc,argv);
  P.SetDescription("Adds taxonomic info to a fasta file.");
  P.registerArg(fileCmmd);
  P.registerArg(outCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string outName = P.GetStringValueFor(outCmmd);
  
  vecDNAVector dna;
  dna.Read(fileName);

  svec<Species> spec;
  FlatFileParser parser1;
  parser1.Open("/references/taxonomy/taxonomy.txt");

  while (parser1.ParseLine()) {
    if (parser1.GetItemCount() == 0)
      continue;
    Species s;
    if (parser1.GetItemCount() > 9) {
      s.Read(parser1);
      spec.push_back(s);
    }
  }
  int i;

  for (int j=0; j<dna.isize(); j++) {
    string s;
    Parse(s, dna.Name(j));

    string tax;
    for (int i=0; i<spec.isize(); i++) {
      if (spec[i].Name() == s) {
        tax = spec[i].Tax();
        break;
      }
    }
    string nn = dna.Name(j);
    nn += " ";
    nn += tax;
    dna.SetName(j, nn);
  }
  dna.Write(outName);
  return 0;
}
