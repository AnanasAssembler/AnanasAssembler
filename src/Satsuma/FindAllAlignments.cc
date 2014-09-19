
#include <string>

#include "base/CommandLineParser.h"
#include "src/DNAVector.h"
#include "src/DPAlign.h"
#include "src/Satsuma/SatsumaAlign.h"


int Index(const vecDNAVector & v, const string & name)
{
  int i;
  string name2 = ">" + name;
  for (i=0; i<v.isize(); i++) {
    //cout << v.Name(i) << endl;
    if (v.Name(i) == name2)
      return i;
  }

  cout << "NOT FOUND: " << name << endl;
  return -1;
}



int main(int argc,char** argv)
{

  
  commandArg<string> aStringCmmd("-q","query sequence");
  commandArg<string> bStringCmmd("-t","target sequence");
  commandArg<string> sStringCmmd("-s","satsuma output (summary file)");

  commandArg<double> identCmmd("-min","min identity", 0.6);
  commandArg<int> maxCmmd("-max_len","maximum alignment length", 8000);
  commandArg<string> tStringCmmd("-target","this target only", "");
  

  commandLineParser P(argc,argv);
  P.SetDescription("Aligns two fasta files.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(sStringCmmd);
  P.registerArg(tStringCmmd);
  P.registerArg(identCmmd);
  P.registerArg(maxCmmd);

  P.parse();

  string aString = P.GetStringValueFor(aStringCmmd);
  string bString = P.GetStringValueFor(bStringCmmd);

  string tString = P.GetStringValueFor(tStringCmmd);
  string sString = P.GetStringValueFor(sStringCmmd);
  double minIdent = P.GetDoubleValueFor(identCmmd);
  int maxLen = P.GetIntValueFor(maxCmmd);
  

 
  cout << "Reading satsuma file..." << endl;
  SAlignVec aligns;
  aligns.Read(sString);


  vecDNAVector query, target;
  
  cout << "Reading query " << aString << endl;
  query.Read(aString);
  cout << "Reading target " << bString << endl;
  target.Read(bString);



  int i, j;

  int off = 10;

  bool bLast = false;
  int printOff = 0;

  AlignsPrinter printer;

  for (i=0; i<aligns.isize()-1; i++) {
    const SAlign & one = aligns[i];
    const SAlign & two = aligns[i+1];

   
    if (tString != "" && tString != one.Target()) {
      bLast = false;
      continue;
    }
    //cout << "Considering target " << one.Target() << endl;    
    if (two.Target() != one.Target()) {
      bLast = false;
      continue;
    }

    if (two.Query() != one.Query()) {
      bLast = false;
      continue;
    }

    if (two.Direction() != one.Direction()) {
      bLast = false;
      continue;
    }

    // Align including both anchors...? No.
    int startTarget = one.StartTarget();
    int endTarget = two.EndTarget() - off;

    int startQuery = one.StartQuery();
    int endQuery = two.EndQuery() - off;

    if (bLast) {
      startTarget = one.EndTarget() - off;
      startQuery = one.EndQuery() - off;
    }


    if (endTarget - startTarget > maxLen) {
      bLast = false;
      continue;
    }

    if (one.Direction() == 1) {
      if (endQuery - startQuery > maxLen) {
	bLast = false;
	continue;     
      }
    } else {
      startQuery = two.StartQuery();
      endQuery = one.EndQuery();
    }

    if (endTarget-startTarget <= 0 || endQuery-startQuery <= 0) {
      bLast = false;
      continue;
    }


    //if (startTarget != 1401554)
    //continue;



    DNAVector t, q;
    //cout << "Subsetting, q_size=" << endQuery-startQuery << " t_size=" << endTarget-startTarget << endl;
    const DNAVector & t_full = target[Index(target, one.Target())];
    const DNAVector & q_full = query[Index(query, one.Query())];

    t.SetToSubOf(t_full, startTarget, endTarget-startTarget);
    q.SetToSubOf(q_full, startQuery, endQuery-startQuery);

    if (one.Direction() == -1) {
      q.ReverseComplement();
    }

    if (!bLast) {
      cout << one.Target() << " [" << startTarget << "-" << endTarget << "] vs " << one.Query();
      cout << " [" << startQuery << "-" << endQuery << "] ";
      
      if (one.Direction() == -1)
	cout << "-" << endl;
      else
	cout << "+" << endl;
    }


    FullAlignment align;

    DPAligner aligner;
    aligner.SetRect(t.isize(), q.isize());
    aligner.Align(align, t, q);

    //cout << "Print" << endl;

    if (!bLast)
      printOff = 0;

    bool bNew = true;
    if (bLast)
      bNew = false;

    printer.PrettyPrint(align, t, q, 0, 0, bNew);

    //printOff = align.PrettyPrint(t, q, 0, 0, printOff);

    bLast = true;
    //cout << "Done" << endl;
    cout << endl;

  
  }


  return 0;

}
  
