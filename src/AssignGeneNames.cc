#include <string>

#include "base/CommandLineParser.h"
#include "src/DNAVector.h"
#include "base/FileParser.h"



class Names
{
public:
  Names() {}

  void Clear() {
    m_name.clear();
    m_len.clear();
  }

  void Push(const string & name, int len) {
    for (int i=0; i<m_name.isize(); i++) {
      if (m_name[i] == name) {
	m_len[i] += len;
	return;
      }
    }
    m_name.push_back(name);
    m_len.push_back(len);
  }

  int Best(string & name) {
    int max = 0;
    int best = -1;
    for (int i=0; i<m_name.isize(); i++) {
      if (m_len[i] > max) {
	max = m_len[i];
	best = i;
      }
    }
    name = m_name[best];
    return max;
  }


private:
  svec<int> m_len;
  svec<string> m_name;
};


string Strip(const string & s) 
{
  char o[1024];
  strcpy(o, s.c_str());
  int n = strlen(s.c_str());
  if (s[n-3] == 'O' && s[n-2] == 'R' && s[n-1] == 'F') {
    o[n-4] = 0;
  }
  return o;
}

int main(int argc,char** argv)
{

  
  commandArg<string> aStringCmmd("-f","fasta file");
  commandArg<string> papayaCmmd("-i","Papaya alignment file");
  commandArg<string> bStringCmmd("-o","output file");
  commandLineParser P(argc,argv);
  P.SetDescription("Transfers gene names from Papaya alignments to a fasta file.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(papayaCmmd);
 
  P.parse();

  string aString = P.GetStringValueFor(aStringCmmd);
  string bString = P.GetStringValueFor(bStringCmmd);
  string papaya = P.GetStringValueFor(papayaCmmd);
  
  vecDNAVector dna;
  
  cout << "Reading file..." << endl;
  dna.Read(aString);
  cout << "done!" << endl;
  
  
  FlatFileParser parser;
  
  parser.Open(papaya);
  
  string last;

 
  Names names;

  while (parser.ParseLine()) {
    const string & gene = parser.AsString(0);
    int start = parser.AsInt(1);
    int end = parser.AsInt(2);
    const string & trans = parser.AsString(3);

    if (trans != last && last != "") {
      int index = dna.NameIndex(last);
      if (index == -1)
	index = dna.NameIndex(Strip(last));

      cout << Strip(last) << endl;
      string newName;
      
      int len = names.Best(newName);
      cout << "Found: " << newName << " index: " << index << endl;
      if (len > 100)
	dna.SetName(index, ">" + newName + " " + last);
      names.Clear();
     
    }
    last = trans;
    names.Push(gene, end-start);
   
  }
  int index = dna.NameIndex(last);
  string newName;
  int len = names.Best(newName);
  if (len > 100)
    dna.SetName(index, ">" + newName);
   
  dna.Write(bString);

  return 0;

}
  
