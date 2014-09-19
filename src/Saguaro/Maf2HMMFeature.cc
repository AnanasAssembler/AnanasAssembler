#include <string>

#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/Saguaro/HMMDistance.h"
#include "src/DNAVector.h"



double SimpleDist(char v1, char v2) 
{
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
    /*int n = strlen(seq.c_str());
    for (int i=0; i<n; i++) {
      if (m_seq[i] == 'a')
	m_seq[i] = 'A';
      if (m_seq[i] == 'c')
	m_seq[i] = 'C';
      if (m_seq[i] == 'g')
	m_seq[i] = 'G';
      if (m_seq[i] == 't')
	m_seq[i] = 'T';
	}*/
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

  
  commandArg<string> aStringCmmd("-i","input multiple alignment file (MAF format)");
  commandArg<string> nStringCmmd("-n","names of the genomes to be extracted (must match MAF)");
  commandArg<string> bStringCmmd("-o","binary feature output files");
  commandArg<string> centerCmmd("-c","name of the genome in which the coordinates will be reported");
  commandArg<int> minCmmd("-m","minimum coverage", 2);
  commandArg<bool> sameCmmd("-nosame","skip positions in which all calls are the same", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Converts MAF data into a Saguaro-digestable file of features");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(nStringCmmd);
  P.registerArg(sameCmmd);
  P.registerArg(minCmmd);
  P.registerArg(centerCmmd);
  
  P.parse();

  string in = P.GetStringValueFor(aStringCmmd);
  string name = P.GetStringValueFor(nStringCmmd);
  string out = P.GetStringValueFor(bStringCmmd);
  string center = P.GetStringValueFor(centerCmmd);
  bool bSkipSame = P.GetBoolValueFor(sameCmmd);
  int minCov = P.GetIntValueFor(minCmmd);
  
  

  HMMFeatureVector f;

  FlatFileParser parser;
  parser.Open(in);



  FlatFileParser parserNames;
  parserNames.Open(name);
  //  for(int c=0; c<numInds; c++) 
  int k = 0;

  cout << "Adding genomes..." << endl;
  
  svec<SeqPart> seqs;
  int n = 0;
  while (parserNames.ParseLine()) {
    if (parserNames.GetItemCount() == 0)
      continue;
    f.SetName(n, parserNames.AsString(0));
    seqs.push_back(SeqPart(parserNames.AsString(0)));
    n++;
  }

  cout << "n=" << n << endl;

  /*

  // get names from first line
  parserIn.ParseLine();
  int numInds=parserIn.GetItemCount();
  if(numInds==0) {
    cerr << "Wrong format: list names in first line" << endl;
    exit(-1);
  }  
  for(int c=0; c<numInds; c++) 
    f.SetName(c,parserIn.AsString(c));
  */

  k = 0;
  
  HMMFeature feature;
  
  int i, j, l;
    
  feature.resize(n * n);
  int len = 0;

  svec<double> all;
  svec<double> div;

  all.resize(n*n, 0.);
  div.resize(n*n, 0.);
       

  char * line = new char[n*n];

  string chr;
  int pos = -1;
  string human;
  while (parser.ParseLine()) {
    
    if (parser.GetItemCount() == 0)
      continue;

    if (parser.AsString(0) == "a") {
      //cout << "Processing..." << endl;
      //if (pos < 44468675)
      //continue;
      //if (pos > 300000)
      //break;
      int posOff = 0;
     
      for (int x=0; x<len; x++) {
	//char line[256];
	int cov = 0;
	for (i=0; i<n; i++) {
	  line[i] = seqs[i].Get(x);
	  if (line[i] != '?' && line[i] != '-' && line[i] != 'N')
	    cov++;
	}

	if (cov < minCov)
	  continue;
	if (bSkipSame) {
	  int diff = 0;
	  bool bSkip = true;
	  int valid = 0;
	  if (line[0] != '?' && line[0] != '-')
	    valid++;

	  int one = (char)toupper(line[0]);
	  for (i=1; i<n; i++) {
	    if (line[i] != '?' && line[i] != '-')
	      valid++;
	    if (line[i] != '?' && line[i] != '-' && toupper(line[i]) != one) {
	      bSkip = false;
	      diff++;
	      //break;
	    }
	  }

	  if (line[0] == '-')
	    bSkip = true;

	  for (i=0; i<n; i++) {
	    if (line[i] >= 'a')
	      bSkip = true;
	  }
	  
	  if (bSkip) {
	    posOff++;
	    continue;
	  }
	  if (diff == valid-1 || diff == 1) {
	    posOff++;
	    continue;
	  }
	}

	
	feature.SetName(chr);
	feature.SetPosition(pos+posOff);
	int posOffOrig = posOff;
	if ((human.c_str())[x] != '-')
	  posOff++;
	for (i=0; i<n; i++) {
	  for (j=0; j<n; j++) {
	    feature[i*n + j] = SimpleDist(line[i], line[j]);
	    if (feature[i*n + j] >= 0.) {
	      all[i*n + j] += feature[i*n + j];
	      div[i*n + j] += 1.;
	    }
	  }
	}
      
      
	if (k >= f.isize())
	  f.resize(k + 50000000);

	//cout << "Adding feature " << k << endl;
	f.SetCompressed(k, line, n, chr, pos+posOffOrig);
	cout << k << " " << chr << " " << pos+posOffOrig;
	for (int x=0; x<n; x++) {
	  cout << " " << line[x];
	}
	cout << endl;
	//f[k] = feature;
	k++;
     }
      
      //if (k > 5000000)
      //break;
     
      for (i=0; i<n; i++)
	seqs[i].Reset();
      continue;
    }

  
    if (parser.AsString(0) == "s") {
      len = strlen(parser.AsString(6).c_str());
      int bFound = false;
      const char * p = parser.AsString(1).c_str();
      if (strstr( p, center.c_str()) != NULL) {
	int from = strlen(center.c_str()) + 1;
	chr = &p[from];
	pos = parser.AsInt(2);
	human = parser.AsString(6);
      }
      
      for (i=0; i<n; i++) {
	if (strstr(parser.AsString(1).c_str(), seqs[i].Name().c_str()) != NULL) {
	  seqs[i].Set(parser.AsString(6));
	  bFound = true;
	  break;
	}
      }
      //if (!bFound)
      //cout << "NOT FOUND: " <<  parser.AsString(1) << endl;
    }    

  }

  f.resize(k);


  
  for(int size=0; size < 12; size++) {
    cout << "Background" << size << endl;
    for(int name=0; name < n; name++) 
      cout << seqs[name].Name() << "\t";
    cout << endl;

    for (i=0; i<n; i++) {
      cout << seqs[i].Name() << "\t";
      double mul = (double)size/8.;
      for (j=0; j<n; j++) {
        cout << mul*all[i*n+j]/div[i*n+j] + 0.001 << "\t";
      }
      cout << endl;
    }

    cout << endl;


  }


  cout << "Writing data" << endl;
  f.Write(out);
  cout << "done." << endl;

  delete [] line;

  return 0;

}
  
