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

  commandArg<string> aStringCmmd("-i","input VCFv4.1 file (multiple species or individuals)");
  commandArg<string> bStringCmmd("-o","binary output file");
  commandArg<int> minCmmd("-m","minimum coverage", 2);
  commandArg<int> diffCmmd("-d","minimum disagreeing", 2);
  commandArg<bool> sameCmmd("-nosame","skip positions in which all calls are the same", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Converts VCF file into a Saguaro-digestable file of features");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(sameCmmd);
  P.registerArg(minCmmd);
  P.registerArg(diffCmmd);
  
  P.parse();

  string in = P.GetStringValueFor(aStringCmmd);
  string out = P.GetStringValueFor(bStringCmmd);
  bool bSkipSame = P.GetBoolValueFor(sameCmmd);
  int minCov = P.GetIntValueFor(minCmmd);
  int minDiff = P.GetIntValueFor(diffCmmd);
  
  

  HMMFeatureVector f;

  	FlatFileParser genotypeparser;
	genotypeparser.Open(in);
	  int i,k=0,r,individuals=0,headercount=0;
	  svec<SeqPart> seqs;
	  int j, l, q;
 
	  HMMFeature feature;

	while (genotypeparser.ParseLine()) {
//skip Meta-information lines..may be this can be used later?
		if (genotypeparser.AsString(0).substr(0,2).compare("##")==0)
			continue;
		if (genotypeparser.AsString(0).substr(0,2).compare("#C")==0){//header line
	headercount=genotypeparser.GetItemCount();
	//Make sure the file has atleast two individuals
		if (headercount < 11){
		  cout << "Less than 11 fields in header. Check file format." << endl;
			exit (1);
		}

			  individuals=headercount-9;
			  cout << "Adding individuals..." << endl;
			  k = 0;
			  for (i=0; i<individuals; i++) {
				//update all individuals or species
			    r=i+9;
			    f.SetName(i, genotypeparser.AsString(r).c_str());
			    seqs.push_back(SeqPart(genotypeparser.AsString(r).c_str()));    
			  }
		  cout << "n=" << individuals << endl;
		}//end of header if
		else{//record lines
	  int n = individuals;   
	  feature.resize(n * n);

	  svec<double> all;
	  svec<double> div;
	  all.resize(n*n, 0.);
	  div.resize(n*n, 0.);
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
      string refbase = genotypeparser.AsString(3);
      string altbase = genotypeparser.AsString(4);
    for (j=0; j<n; j++) {
//can be one of A,G,T,C,- or a heterozgote
	q=j+9;
//	cout <<genotypeparser.AsString(q)<<"first"<<atoi(&genotypeparser.AsString(q)[0])<<endl;
	if(refbase.length() > 1 || altbase.length() > 3 || (altbase.length() > 1 && altbase.substr(1,1).compare(",")!=0)){
//this is a microsat or indel...skip for now...needs a dist function
	continue;
	}
	else{
if((genotypeparser.AsString(q)[0]!='.')&&(genotypeparser.AsString(q)[2]!='.')){

	int gen1=atoi(&genotypeparser.AsString(q)[0]);
	int gen2=atoi(&genotypeparser.AsString(q)[2]);
	string genotemp;

	if(altbase.length() > 1){
	  char genos [] ={refbase[0],altbase[0],altbase[2]};
	  genotemp.push_back(genos[gen1]);
	  genotemp.push_back(genos[gen2]);
	}
	else{
	  char genos [] ={refbase[0],altbase[0]};
	  genotemp.push_back(genos[gen1]);
	  genotemp.push_back(genos[gen2]);
	}
	const string & gen=genotemp;
	line[j]=GetAmbiguous(gen);
 }
 else{
   line[j]='?';
}
	}

      if (line[j] != '?')// a dot means a call cannot be made for a sample at a given locus
	cov++;
    }
    if (cov < minCov){
	cout << "skipping" << chrpos << endl;
	continue;
    }
    
    int ident = 0;
    int diff = 0;
    for (int x=1; x< line.size(); x++) {
      if (line[x] == '?')
	continue;
      if (line[x] == line[0])
	ident++;
      else
	diff++;
    }
    if (ident < minDiff || diff < minDiff)
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

    f.SetCompressed(k, line.c_str(), n, chromosome, chrpos);
    cout << k << " " << chromosome << " " << chrpos;
    for (int x=0; x<n; x++) {
      cout << " " << line[x];
    }
    cout << endl;
    
    k++;
		}//end of else 
	}//end of parse while
  f.resize(k);


  cout << "Writing data" << endl;
  f.Write(out);
  cout << "done." << endl;

 

  return 0;

}
  
