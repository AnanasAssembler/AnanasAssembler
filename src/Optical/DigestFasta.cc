#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input fasta");
  commandArg<string> seqCmmd("-s","motif", "ACTAGT");
  commandArg<int> lenCmmd("-l","Length", 600000);
  commandLineParser P(argc,argv);
  P.SetDescription("Looks for restriction sites and cuts.");
  P.registerArg(fileCmmd);
  P.registerArg(seqCmmd);
  P.registerArg(lenCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string seq = P.GetStringValueFor(seqCmmd);
  int len = P.GetIntValueFor(lenCmmd);
  
  vecDNAVector dna;
  dna.Read(fileName);
  int i, j, k;
  int n = seq.size();
  svec<double> tmp; 
  for (i=0; i<dna.isize(); i++) {
    const DNAVector & d = dna[i];
    int last = -1;
    int lastStart = 0;
    
    for (j=0; j<d.isize()-n; j++) {
      bool b = true;
      for (k=0; k<n; k++) {
	//	cout << d[k+j] != seq[k] 
	if (d[k+j] != seq[k]) {
	  b = false;
	  break;
	}
      }
      if (b) {
	if (last == -1) {
	  last = j;
	  lastStart = j;
	  continue;
	}
	cout << j << " " << last << endl;
	int l = j - last;
	last = j;		
	if (j - lastStart > len) {
	  cout << dna.NameClean(i) << "_" << lastStart << endl;
	  cout << seq << "\tS";
	  for (int x=0; x<tmp.isize(); x++) {
	    cout << "\t" << tmp[x]/1000.;
	  }
	  cout << endl;
	  lastStart = j;
	  tmp.clear();
	}

  	tmp.push_back(l);
      }

    }
    if (tmp.isize() > 3) {
      cout << dna.NameClean(i) << "_" << lastStart << endl;
      cout << seq << "\tS";
      for (int x=0; x<tmp.isize(); x++) {
	cout << "\t" << tmp[x]/1000.;
      }
      cout << endl;
    }
  }
 
  return 0;
}
