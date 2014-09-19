#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <string.h>

#include "base/CommandLineParser.h"

#include "aligns/KmerAlignCore.h"
#include "src/DNAVector.h"
#include <math.h>

#include "src/Papaya/KmerTable.h"




//========================================================================
//========================================================================
//========================================================================



int main(int argc,char** argv)
{

  commandArg<string> aStringCmmd("-i","fasta file");
  commandArg<string> oStringCmmd("-o","fasta output");
  commandArg<int> kCmmd("-k","kmer size", 25);
 
  commandLineParser P(argc,argv);
  P.SetDescription("Eliminates redundant sequences.");
  P.registerArg(aStringCmmd);
  P.registerArg(oStringCmmd);
  P.registerArg(kCmmd);


  P.parse();

  string aString = P.GetStringValueFor(aStringCmmd);
  string oString = P.GetStringValueFor(oStringCmmd);
  int k = P.GetIntValueFor(kCmmd);
 
  int i, j;

  vecDNAVector seq;
  seq.Read(aString);
  vecDNAVector out;

  
  KmerSequence kmers(k, &seq);
  kmers.Add(seq);
  
  long long m = kmers.GetBoundValue();


  svec<int> counts;
  counts.resize(m, 0);

  cout << "Counting." << endl;

  kmers.Count(counts);
  cout << "Done." << endl;

  for (i=0; i<seq.lsize(); i++) {
    DNAVector d = seq[i];
    DNAVector r = d;
    r.ReverseComplement();
    
    int nn = 0;
    int mult = 0;
    for (j=0; j<=d.lsize()-k; j+=k) {
      nn++;
      long long n = kmers.BasesToNumber(d, j) - 1;
      n += kmers.BasesToNumber(r, j);
 
      if (counts[n] > 1)
	mult++;
      //n  = kmers.BasesToNumber(r, j);
      //if (counts[n] > 2)
      //mult++;
    }
    double frac = (double)mult/(double)nn;
    //    if (frac < 0.5) {
    if (mult < 5) {
      out.push_back(d, seq.Name(i));
    }
  } 

  out.Write(oString);
  return 0;

}
