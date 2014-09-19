#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"

#include "base/RandomStuff.h"
#include "src/DNAVector.h"


void Mutate(DNAVector & out, const DNAVector & in, double sub, double indel)
{
  int i, j;
  
  string n;
  
  for (i=0; i<in.isize(); i++) {
    if (RandomFloat(1.) < indel) {
      if (RandomFloat(1.) < 0.5) {
	// Deletion
	i += RandomInt(3) + 1;
      } else {
	// Insertion
	int plus = RandomInt(3) + 1;
	for (j=0; j<plus; j++) {
	  n += NucLetter(RandomInt(4));
	}
      }
    }
    
    if (RandomFloat(1.) < sub) {
      int s = NucIndex(in[i]);
      s += RandomInt(3) + 1;
      s = s % 4;
      n += NucLetter(s);
    } else {
      n += in[i];
    }
  }
  out.SetFromBases(n);
}



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input fasta file");
  commandArg<string> outCmmd("-o","output fasta file");
  commandArg<int> lCmmd("-l","read length");
  commandArg<double> cCmmd("-c","coverage");
  commandArg<double> indelCmmd("-indel","indel error rate");
  commandArg<double> subCmmd("-s","substitution error rate");
  commandArg<bool> strandCmmd("-d","Double strand (1) or single strand (0)", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Simulates reads from 3rd generation sequencing.");
  P.registerArg(fileCmmd);
  P.registerArg(outCmmd);
  P.registerArg(lCmmd);
  P.registerArg(cCmmd);
  P.registerArg(indelCmmd);
  P.registerArg(subCmmd);
  P.registerArg(strandCmmd);
  
  P.parse();
  
  string fileName   = P.GetStringValueFor(fileCmmd);
  string outName    = P.GetStringValueFor(outCmmd);
  int len           = P.GetIntValueFor(lCmmd);
  double cov        = P.GetDoubleValueFor(cCmmd);
  double indel      = P.GetDoubleValueFor(indelCmmd);
  double sub        = P.GetDoubleValueFor(subCmmd);
  bool doubleStrand = P.GetBoolValueFor(strandCmmd);
  
  vecDNAVector dna;
  vecDNAVector reads;
  
  cout << "Reading fasta file." << endl;
  dna.Read(fileName);


  int i, j;

  double size = 0;
  for (i=0; i<dna.isize(); i++) {
    size += dna[i].isize();
  }
  cout << "Genome size: " << size << endl;

  int n = (int)((cov * size)/(double)len);

  cout << "WARNING: Generating data only for chromosome 1!" << endl;
  
  srand(time(NULL));

  for (i=0; i<n; i++) {
    int start = RandomInt(size-len);
    cout << i << ":" << start << '\n';
    DNAVector tmp, out;
    tmp.SetToSubOf(dna[0], start, len);
    if(doubleStrand && RandomFloat(1.) < 0.5) {
      cout<<"Reverse-complemeting"<<endl;
      tmp.ReverseComplement();
    }
    Mutate(out, tmp, sub, indel);
    ostringstream ss;
    ss << ">SimulatedRead_"<<i;
    out.setName(ss.str());
    reads.push_back(out);
  }

  reads.Write(outName);

  return 0;
}
