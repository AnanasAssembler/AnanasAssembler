#include <string>

#include "base/CommandLineParser.h"
#include "src/Cola/Cola.h"


bool SameScaff(string a, string b) 
{
  int i;
  for (i=(int)a.size()-1; i>=10; i--) {
    if (a[i] == '_') {
      a[i] = 0;
      break;
    }
  }
  for (i=(int)b.size()-1; i>=10; i--) {
    if (b[i] == '_') {
      b[i] = 0;
      break;
    }
  }
 
  return (strcmp(a.c_str(), b.c_str()) == 0);
}


int main(int argc,char** argv)
{
  commandArg<string> aStringCmd("-i","input sequence");
  commandArg<string> bStringCmd("-o","output sequence");
  commandArg<int>    alignerTypeCmd("-a","Aligner type - Choose 1 : NSGA , 2 : NS , 3 : SWGA, 4 : SW ", 1);
  commandArg<int>    gapOpenCmd("-o", "Aligner gap open penalty", false);
  commandArg<int>    mismatchCmd("-m", "Aligner Mismatch penalty", false);
  commandArg<int>    gapExtCmd("-e", "Aligner gap extension penalty", false);
  commandArg<double> maxPCmd("-p","Maximum acceptable P-value", 1.0);
  commandArg<double> minIdentCmd("-ident","Minium acceptable identity", 0.99);
  commandArg<int>    bandedCmd("-b", "The bandwidth for banded mode", 6);

  commandLineParser P(argc,argv);
  P.SetDescription("Removes redundant transcripts if too similar.");
  P.registerArg(aStringCmd);
  P.registerArg(bStringCmd);
  P.registerArg(alignerTypeCmd);
  P.registerArg(gapOpenCmd);
  P.registerArg(mismatchCmd);
  P.registerArg(gapExtCmd);
  P.registerArg(maxPCmd);
  P.registerArg(minIdentCmd);
  P.registerArg(bandedCmd);

  P.parse();

  string      aString     = P.GetStringValueFor(aStringCmd);
  string      bString     = P.GetStringValueFor(bStringCmd);
  AlignerType aType       = AlignerType(P.GetIntValueFor(alignerTypeCmd));
  int         gapOpenPen  = P.GetIntValueFor(gapOpenCmd);
  int         mismatchPen = P.GetIntValueFor(mismatchCmd);
  int         gapExtPen   = P.GetIntValueFor(gapExtCmd);
  double      maxP        = P.GetDoubleValueFor(maxPCmd);
  double      minIdent    = P.GetDoubleValueFor(minIdentCmd);
  int         banded      = P.GetIntValueFor(bandedCmd);

  vecDNAVector in, out;

  in.Read(aString);
 
  int i, j;
  
  for (i=0; i<in.isize(); i++) {
    if (in[i].isize() == 0)
      continue;
    for (j=i+1; j<in.isize(); j++) {
      cout << "j=" << j << endl;
      if (!SameScaff(in.Name(i), in.Name(j))) {
	cout << "Different: " <<  in.Name(i) << " " <<  in.Name(j) << endl;
	break;
      }

      const DNAVector & target = in[i];
      const DNAVector & query = in[j];

      cout << "Aligning " << in.Name(i) << " vs. " << in.Name(j) << endl;
      Cola cola1 = Cola();
      if(gapOpenPen && mismatchPen && gapExtPen) {
	cola1.createAlignment(target, query,
			      AlignerParams(banded, aType, -gapOpenPen, -mismatchPen, -gapExtPen));
      } else { // If params are not given, use default mode
	cola1.createAlignment(target, query, AlignerParams(banded, aType));
      }
      Alignment cAlgn = cola1.getAlignment();
      bool bGood = false;
      int wiggle = banded;

      if (cAlgn.getQueryLength() - cAlgn.getQueryBaseAligned() <= wiggle ||
	  cAlgn.getTargetLength() - cAlgn.getTargetBaseAligned() <= wiggle)
	bGood = true;
      cout << "Result: " << cAlgn.getIdentityScore() << " " << minIdent << endl;
      cAlgn.print(0,1,cout,100);
      if(cAlgn.calcPVal()<=maxP && cAlgn.getIdentityScore()>=minIdent
	 && bGood) {
	//cout << target.Name(i) << " vs " << query.Name(j) << endl;
	cout << "Remove " << in.Name(j) << endl;
	in[j].resize(0);
      } else {
	cout<<"No Alignment at given significance threshold, keeping sequence"<<endl;     
      }
    }
  }
  
  for (i=0; i<in.isize(); i++) {
    if (in[i].isize() > 0)
      out.push_back(in[i], in.Name(i));
  }
  out.Write(bString);

  return 0;

}

