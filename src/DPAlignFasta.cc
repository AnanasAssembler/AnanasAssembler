#include <string>

#include "base/CommandLineParser.h"
#include "src/DNAVector.h"
#include "src/DPAlign.h"


int main(int argc,char** argv)
{

  
  commandArg<string> aStringCmmd("-q","query sequence");
  commandArg<string> bStringCmmd("-t","target sequence");
  commandArg<bool> selfCmmd("-s","self-alignments", false);
  commandArg<bool> sameCmmd("-same","same alignments", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Aligns two fasta files.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(selfCmmd);
  P.registerArg(sameCmmd);

  P.parse();

  string aString = P.GetStringValueFor(aStringCmmd);
  string bString = P.GetStringValueFor(bStringCmmd);
  bool bSelf = P.GetBoolValueFor(selfCmmd);
  bool bSame = P.GetBoolValueFor(sameCmmd);
  


  vecDNAVector query, target;
  
  query.Read(aString);
  target.Read(bString);
 
  int i, j;

  for (i=0; i<target.isize(); i++) {
    for (j=0; j<query.isize(); j++) {
      if (bSelf && i == j)
	continue;
      if (bSame && i != j)
	continue;

      cout << target.Name(i) << " vs " << query.Name(j) << endl;
      FullAlignment align;

      DPAligner aligner;
      //aligner.UseRewardFunc(true);
      aligner.UseRewardFunc(true);
      aligner.SetRect(target[i].isize(), query[j].isize());
      aligner.Align(align, target[i], query[j]);

      cout << "ALIGNMENT FW" << endl;
      int n = align.PrettyPrint(target[i], query[j]);
      //cout << endl;
      double ident = 2.*(double)n/(double)(target[i].isize()+query[j].isize());
      cout << "Matches: " << n << " ident: " << ident;
      if (ident > 0.8)
	cout << " good";
      cout << endl << endl;
      //if (n > target[i].isize() - 5 && n > query[j].isize() - 5)
      //cout << "******" << endl;

      
      DPAligner alignerRC;
      alignerRC.SetRect(target[i].isize(), query[j].isize());
      query[j].ReverseComplement();
      align.clear();
      alignerRC.Align(align, target[i], query[j]);

      cout << "ALIGNMENT RC" << endl;
      n = align.PrettyPrint(target[i], query[j]);
      //cout << endl;
      ident = 2.*(double)n/(double)(target[i].isize()+query[j].isize());
      cout << "Matches: " << n << " ident: " << ident;
      if (ident > 0.8)
	cout << " good";
      cout << endl << endl;
  
      
       //if (n > target[i].isize() - 5 && n > query[j].isize() - 5)
      //cout << "******" << endl;
      
    }
  }


  return 0;

}
  
