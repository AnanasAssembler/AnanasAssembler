#include <string>

#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/Saguaro/HMMDistance.h"
#include "src/DNAVector.h"



double SimpleDist(char v1, char v2) 
{
  if (v1 == '-' || v2 == '-') {
    if (v1 == v2)
      return 0.;
    else
      return 1.;
  }

  if(v1=='N' || v2=='N' || v1=='X' || v2=='X') return -1.;
  else return DNA_Diff(v1,v2);
}


class SeqPart
{
public:
  SeqPart(const string & name) {
    m_name = name;
  }

  void Set(const string & seq) {
    m_seq = seq;
  } 

  char Get(int i) const {
    if (m_seq == "")
      return '?';
    const char * p = m_seq.c_str();
    return p[i];
  }

  void Reset() {
    m_seq = "";
  }

  const string & Name() const {return m_name;}
private:
  string m_name;
  string m_seq;
};







int main(int argc,char** argv)
{

  
  commandArg<string> aStringCmmd("-i","input fasta file (multiple alignment)");
  commandArg<string> bStringCmmd("-o","binary output file");
  commandArg<int> minCmmd("-m","minimum different", 2);
  commandArg<bool> sameCmmd("-nosame","skip positions in which all calls are the same", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Converts multi-fasta data into a Saguaro-digestable file of features");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(sameCmmd);
  P.registerArg(minCmmd);
  
  P.parse();

  string in = P.GetStringValueFor(aStringCmmd);
  string out = P.GetStringValueFor(bStringCmmd);
  bool bSkipSame = P.GetBoolValueFor(sameCmmd);
  int minCov = P.GetIntValueFor(minCmmd);
  
  

  HMMFeatureVector f;
  
  vecDNAVector dna;
  dna.Read(in);
  int k = 0;
  
  int i, j, l;
 
  cout << "Adding genomes..." << endl;
  
  svec<SeqPart> seqs;
 
  for (i=0; i<dna.isize(); i++) {
    f.SetName(i, dna.NameClean(i));
    seqs.push_back(SeqPart(dna.NameClean(i)));    
  }

  cout << "n=" << dna.isize() << endl;

  k = 0;
  
  HMMFeature feature;
  int n = dna.isize();
    
  feature.resize(n * n);
 

  svec<double> all;
  svec<double> div;

  all.resize(n*n, 0.);
  div.resize(n*n, 0.);
       

  //char * line = new char[n*n];
  int len = dna[0].isize();

  for (i=0; i<len; i++) {
    feature.SetName("mult");
    feature.SetPosition(i);
    string line;
    line.resize(n);
    svec<int> letters;
    letters.resize(4, 0);
    //cout << "Same" << endl;
    //const DNAVector & dd = 
    for (j=0; j<n; j++) {
      line[j] = (dna[j])[i];
      if (line[j] == 'A' || line [j] == 'C' || line[j] == 'G' || line[j] == 'T')    
	letters[NucIndex(line[j])]++;
    }
    int diff = 0;
    for (j=0; j<letters.isize(); j++) {
      if (letters[j] >= minCov)
	diff++;
    }
    
    if (diff < 2)
      continue;

    //cout << "Dist" << endl;
    for (l=0; l<n; l++) {
      for (j=0; j<n; j++) {
	feature[l*n + j] = SimpleDist(line[l], line[j]);
	if (feature[l*n + j] >= 0.) {
	  all[l*n + j] += feature[l*n + j];
	  div[l*n + j] += 1.;
	}
      }
    }
    //cout << "Add" << endl;
      
    if (k >= f.isize())
      f.resize(k + 5000000);

    f.SetCompressed(k, line.c_str(), n, "mult", i);
    cout << k << " " << " " << i;
    for (int x=0; x<n; x++) {
      cout << " " << line[x];
    }
    cout << endl;
    
    k++;
  }
      
  
  
  f.resize(k);


  cout << "Writing data" << endl;
  f.Write(out);
  cout << "done." << endl;

 

  return 0;

}
  
