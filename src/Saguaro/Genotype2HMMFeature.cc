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
  else return DNA_DiffAmb(v1,v2);
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

  commandArg<string> aStringCmmd("-i","input genotype file (multiple species or individuals)");
  commandArg<string> bStringCmmd("-o","binary output file");
  commandArg<int> minCmmd("-m","minimum coverage", 2);
  commandArg<bool> sameCmmd("-nosame","skip positions in which all calls are the same", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Converts genotype data into a Saguaro-digestable file of features");
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

  	FlatFileParser genotypeparser;
	genotypeparser.Open(in);
	//read in the header
	genotypeparser.ParseLine();
	int headercount=genotypeparser.GetItemCount();
	//Make sure the file has atleast two individuals
		if (headercount < 4){
		  cout << "Less than 4 fields in header. Check file format." << endl;
			exit (1);
		}


  int i,k=0,r;
  int individuals=headercount-2;
 
  cout << "Adding individuals..." << endl;
  
  svec<SeqPart> seqs;

  k = 0;

  for (i=0; i<individuals; i++) {
	//update all individuals or species
//	cout << genotypeparser.AsString(i).c_str() << endl;
    r=i+2;
    f.SetName(i, genotypeparser.AsString(r).c_str());
    seqs.push_back(SeqPart(genotypeparser.AsString(r).c_str()));    
  }

  cout << "n=" << individuals << endl;

  int j, l, q;
 

  HMMFeature feature;
  int n = individuals;
    
  feature.resize(n * n);
 

  svec<double> all;
  svec<double> div;

  all.resize(n*n, 0.);
  div.resize(n*n, 0.);

  i=0;
	while (genotypeparser.ParseLine()) {
	i++;
		if (genotypeparser.GetItemCount() == 0)//blank line allowed at end of file only
		      continue;
	//Make sure header and record field counts match
		else if (genotypeparser.GetItemCount() != headercount){
		  cout << "Header and records have different number of fields check file format." << endl;
			exit (1);
		}

    string chromosome = genotypeparser.AsString(0);
    int chrpos = genotypeparser.AsInt(1);
    feature.SetName(chromosome);

//set position on chromosome as position
    feature.SetPosition(chrpos);

    string line;
    line.resize(n);
    svec<int> letters;
    letters.resize(4, 0);

    int cov = 0;

    for (j=0; j<n; j++) {
//can be one of A,G,T,C,- or a heterozgote
	q=j+2;
      string genotype = genotypeparser.AsString(q);
	if(genotype.length() > 2){//genotypes need to be separated by a slash
//this is a heterozygote.
const string & gen=genotypeparser.AsString(q);
	line[j]=GetAmbiguous(gen);
	}
	 else{
	line[j]=genotypeparser.AsString(q)[0];
	}
      if (line[j] != '-')
	cov++;
    }
    if (cov < minCov){
//	cout << "skipping" << chrpos << endl;
	continue;
	}

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

    f.SetCompressed(k, line.c_str(), n, chromosome, chrpos);
    cout << k << " " << " " << chrpos;
    for (int x=0; x<n; x++) {
      cout << " " << line[x];
    }
    cout << endl;
    
    k++;

	}//end of file parsing while loop
       
 
  f.resize(k);


  cout << "Writing data" << endl;
  f.Write(out);
  cout << "done." << endl;

 

  return 0;

}
  
