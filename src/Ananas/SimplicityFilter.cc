#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"

void ConnectUp(svec<int> & partner, const vecDNAVector & seq)
{
  int cc = 0;
  for (int i=0; i<seq.isize(); i++) {
    char name[2048];
    strcpy(name, seq.Name(i).c_str());
    int n = strlen(name);
    bool b = false;
    if (name[n-2] != '/')
      continue;
    if (name[n-1] == '1') {
      name[n-1] = '2';
      b = true;
    } else {
      if (name[n-1] == '2') {
        name[n-1] = '1';
        b = true;
      }
    }
    if (!b)
      continue;
    int index = seq.NameIndex(name);
    if (index < 0)
      continue;
    //cout << "Connecting partners: " << name << " " << seq.Name(i) << endl;
    partner[i] = index;
    cc++;
  }
  cout << "Found partner for " << cc << " reads, out of " << seq.isize() << endl;
}


int Mis(const DNAVector & a, const DNAVector & b)
{
  if (a.isize() != b.isize())
    return a.isize();
  int m = 0;
  for (int i=0; i<a.isize(); i++) {
    if (a[i] != b[i])
      m++;
  }
  return m;
}


bool Good(const DNAVector & d, int min) {
  if (d.isize() < 10)
    return false;

  int i;
  
  int k = 0;
  for (i=1; i<d.isize(); i++) {
    if (d[i] == d[i-1])
      k++;
    else
      break;
  }
  if (k >= min)
    return false;

  k = 0;
  for (i=d[d.isize()-2]; i>=0; i--) {
    if (d[i] == d[i+1])
      k++;
    else
      break;
  }
  if (k >= min)
    return false;
  return true;
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input fasta file");
  commandArg<string> outCmmd("-o","output fasta file");
  commandArg<int> nCmmd("-n","maximum run", 20);
  //commandArg<int> mCmmd("-m","maximum mismatch", 1);
  commandLineParser P(argc,argv);
  P.SetDescription("Removes extreme simple sequence reads.");
  P.registerArg(fileCmmd);
  P.registerArg(outCmmd);
  P.registerArg(nCmmd);
  // P.registerArg(mCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string outName = P.GetStringValueFor(outCmmd);
  int min = P.GetIntValueFor(nCmmd);
  //int mis = P.GetIntValueFor(mCmmd);

  vecDNAVector dna;
  cout << "Reading input fasta." << endl;
  dna.Read(fileName);
  //cout << "Will keep " << keep << " reads." << endl;

  vecDNAVector out;
  //out.resize(keep);
  // k = 0;
  for (int i=0; i<dna.isize(); i++) {
    if (Good(dna[i], min)) {    
      out.push_back(dna[i], dna.Name(i));
    }
  }

   out.Write(outName);
  
  return 0;
}
