#define FORCE_DEBUG
#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"

#include "src/SlideSearch.h"

int main( int argc, char** argv )
{

  commandArg<string> leftCmmd("-l","left input fasta file");
  commandArg<string> rightCmmd("-r","right input fasta file", "");
  commandArg<string> primerCmmd("-p","primer input fasta file");
  commandArg<int> misCmmd("-mis","maximum number of mismatches", 3);
  commandLineParser P(argc,argv);
  P.SetDescription("Prints out all reads and mates that hit the primers.");
  P.registerArg(leftCmmd);
  P.registerArg(rightCmmd);
  P.registerArg(primerCmmd);
  P.registerArg(misCmmd);
  
  P.parse();
 
  vecDNAVector left, right, primer;
 
  primer.Read(P.GetStringValueFor(primerCmmd));
  left.Read(P.GetStringValueFor(leftCmmd));
 
  int mis = P.GetIntValueFor(misCmmd);
  bool bRight = false;
  if (P.GetStringValueFor(rightCmmd) != "") {
    bRight = true;
    right.Read(P.GetStringValueFor(rightCmmd));
  }
  svec<int> hits;
  hits.resize(left.isize());

  SlideSearch s;
  //s.SetMis(2);
  s.SetMis(mis);

  int i, j;
  
  for (i=0; i<primer.isize(); i++) {
    //cout << "Searching primer " << primer.Name(i) << endl;
    int n = 0;
    for (int j=0; j<left.isize(); j++) {
      svec<SlideMatch> out1, out2;
      s.Search(out1, left[j], primer[i]);
      if (bRight) 
	s.Search(out2, right[j], primer[i]);
      if (out1.isize() > 0 || out2.isize() > 0) {
	hits[j]++;
	n++;
      }
    }
    //cout << "found: " << n << endl;
  }

  for (i=0; i<hits.isize(); i++) {
    if (hits[i] > 0) {
      //cout << ">Left_" << i << endl;
      cout << left.Name(i) << endl;
      for (j=0; j<left[i].isize(); j++)
	cout << (left[i])[j];
      cout << endl;
      if (bRight) {
	//cout << ">Right_" << i << endl;
	cout << right.Name(i) << endl;
	for (j=0; j<right[i].isize(); j++)
	  cout << (right[i])[j];
	cout << endl;
      }
    }
  }
  return 0;
}
