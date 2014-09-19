#include <string>

#include "base/CommandLineParser.h"
#include "src/DNAVector.h"
#include "aligns/KmerAlignCore.h"


class AlignBlock
{
public:
  AlignBlock() {}


private:
};


bool SameName(const string & a, const string & b) 
{
  const char * p1 = a.c_str();
  const char * p2 = b.c_str();

  char tmp1[512];
  char tmp2[512];
  strcpy(tmp1, "1");
  strcpy(tmp2, "2");
  int i;
  int k = 0;
  for (i=0; i<(int)strlen(p1); i++) {
    if (p1[i] == '_') {
      if (k == 1) {
	strcpy(tmp1, &p1[i+1]);
	break;
      }
      k++;
    }
  }

  k = 0;
  for (i=0; i<(int)strlen(p2); i++) {
    if (p2[i] == '_') {
      if (k == 1) {
	strcpy(tmp2, &p2[i+1]);
	break;
      }
      k++;
    }
  }

  for (i=0; i<(int)strlen(tmp1); i++) {
    if (tmp1[i] == '_') {
      tmp1[i] = 0;
      break;
    }
  }
  for (i=0; i<(int)strlen(tmp2); i++) {
    if (tmp2[i] == '_') {
      tmp2[i] = 0;
      break;
    }
  }

  //cout << tmp1 << " " << tmp2 << endl;

  if (strcmp(tmp1, tmp2) == 0)
    return true;
  else
    return false;

}



bool Splice(DNAVector & a, DNAVector & b, KmerAlignCore<KmerAlignCoreRecord> & core, int k, int contig) 
{
  int i, j;

  int slack = 48;
  int maxMis = slack;


  int mis = 0;

  int lastA = 0;
  int lastB = 0;

  string one, two;
  
  int laps = 0;
  int gaps = 0;

  int indels = 0;
  bool bSplice = false;

  int matchLen = 0;

  int limit = a.isize();
  if (b.isize() < limit)
    limit = b.isize();

  limit = (int)(0.8*(double)limit);

  for (i=0; i<=a.isize()-k; i++) {
    //if (i > slack)
    //return false;
    svec<KmerAlignCoreRecord> matches;
    core.GetMatches(matches, a, i);

    int pos = -1;

    //cout << "match??" << endl; 

    for (j=0; j<matches.isize(); j++) {
      if (matches[j].GetContig() == contig) {
	pos = matches[j].GetPosition();
	break;
      }
    }
  
    if (pos == -1)
      continue;
    
 
    //cout << "Match." << endl;
    //cout << "Last: " << lastA << " " << i << endl;
    
    
    if (i-lastA == pos-lastB) {
      for (j=lastA; j<i; j++) {
	one += a[j];
	two += b[j-lastA+lastB];
	matchLen++;
      }
    } else {
      for (j=lastA; j<i; j++) {
	one += a[j];
	two += "-";
      }
      for (j=lastB; j<pos; j++) {
	one += "-";
	two += b[j];
      }
      //if ((i-lastA -pos+lastB) % 3 == 0)
      //bSplice = true;
    }

    if (lastA > 0 && lastB > 0) {
      if (pos - lastB > 20 || i-lastA > 20) {
	gaps++;
	
	if (pos - lastB == 0 || i-lastA == 0) {
	  indels++;

	  if ((pos - lastB) %3 == 0 && (i - lastA) % 3 == 0) {
	    bSplice = true;
	    cout << "ALTERNATIVE SPLICING (mod 3)!" << endl;
	  }
	}
      }
    }
    int start = i;

    for (; i<a.isize(); i++) {
      lastA = i;
      lastB = pos;
      //cout << " run." << endl;
      if (pos >= b.isize())
	break;
      if (a[i] != b[pos]) {
	i--;
	break;
      }
      one += a[i];
      two += b[pos];
      matchLen++;
      pos++;
    }

    if (i-start >= 24)
      laps++;
    //lastB = pos + i - lastA;
    //lastA = i;
 
    
    //cout << "Done, i=" << i << " len=" << a.isize() << endl;

  }
  for (j=lastA+1; j<a.isize(); j++) {
    one += a[j];
    two += "-";
  }
  for (j=lastB+1; j<b.isize(); j++) {
    one += "-";
    two += b[j];
  }


  //bSplice = true;
  //cout << "Match length: " << matchLen << " limit: " << limit << endl;
  if (bSplice && indels > 0 && matchLen >= limit) {
    int nn = (int)strlen(one.c_str());
    for (i=0; i<nn; i+=80) {
      for (j=i; j<i+80; j++) {
	if (j >= nn)
	  break;
	cout << one[j];	
      }
      cout << endl;
      for (j=i; j<i+80; j++) {
	if (j >= nn)
	  break;
	cout << two[j];     
      }
      cout << endl << endl;
    }

    /*
    cout << "Appending end, i=" << i << " len=" << strlen(one.c_str())-1 << endl;
    for (i=j; i<strlen(one.c_str())-1; i++) {
      cout << one[i];
    }
    cout << endl;
    for (i=j; i<strlen(two.c_str())-1; i++) {
      cout << two[i];
    }
 
    cout << endl << endl;*/


    return  true;
    //cout << "Found match." << endl;
  }
  return false;
}


int main(int argc,char** argv)
{

  
  commandArg<string> aStringCmmd("-i","Ananas fasta file");
  commandArg<int> maxCmmd("-max","maximum number of sequences per contig", 30);
  commandLineParser P(argc,argv);
  P.SetDescription("Aligns and prints out different transcript splice forms from Ananas assemblies.");
  P.registerArg(aStringCmmd);
  P.registerArg(maxCmmd);

  P.parse();

  string aString = P.GetStringValueFor(aStringCmmd);
  int max = P.GetIntValueFor(maxCmmd);
   
  vecDNAVector dna;

  dna.Read(aString);

  int k = 24;
  TranslateBasesToNumberExact trans;
  trans.SetSize(k/2);

  cout << "Set." << endl;

  KmerAlignCore<KmerAlignCoreRecord> core;
  core.SetTranslator(&trans);
  core.SetNumTables(k/12);

  core.AddData(dna, false, 1);
  cout << "Done." << endl;

  int i, j;

  for (i=0; i<dna.isize(); i++) {
    cout << "Examining seq " << dna.Name(i) << endl;
    if (dna[i].isize() == 0)
      continue;

    int spliceforms = 0;

    //Look ahead to see how many we have
    int same = 0;
    for (j=i+1; j<dna.isize(); j++) {
      if (dna[j].isize() == 0)
	continue;
      if (dna[i].isize() == 0)
	break;
      same++;
      if (!SameName(dna.Name(i), dna.Name(j)))
	break;
    }
    if (max > 0 && same >= max) {
      cout << "Skipping, sequences: " << same << endl;
      i = j-1;
      continue;
    }


    for (j=i+1; j<dna.isize(); j++) {
      if (dna[j].isize() == 0)
	continue;
      if (dna[i].isize() == 0)
	break;
      if (!SameName(dna.Name(i), dna.Name(j)))
	break;
      if (j-i > 30)
	break;
  
      //cout << "Same as " << dna.Name(j) << endl;

      if (Splice(dna[i], dna[j], core, k, j)) {
	cout << "Sequences " << dna.Name(i) << " and " << dna.Name(j) << " are alternatively spliced." << endl;
	continue;
      }
      
      dna[i].ReverseComplement();
      if (Splice(dna[i], dna[j], core, k, j)) {
	cout << "Sequences " << dna.Name(i) << " and " << dna.Name(j) << " are alternatively spliced (RC)." << endl;
      } else {
	dna[i].ReverseComplement();
      }
    }
  }


  //dna.Write(outfile, true);

  return 0;

}
  
