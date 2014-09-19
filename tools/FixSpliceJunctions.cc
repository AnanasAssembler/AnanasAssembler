#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"


int Adjust(const DNAVector & d, int start, const DNAVector & seq, int slack)
{ 
  int from = start-slack;
  int to = start + slack;
  if (from < 0)
    from = 0;
  if (to >= d.isize())
    to = d.isize()-1;

  for (int i=from; i<to; i++) {
    if (d[i] == seq[0] && d[i+1] == seq[1])
      return i;
  }
  return -1;
}


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input GTF file");
  commandArg<string> fastaCmmd("-f","fasta file");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);
  P.registerArg(fastaCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string fastaName = P.GetStringValueFor(fastaCmmd);
  
  vecDNAVector dna;
  dna.Read(fastaName);

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  DNAVector leftFW, rightFW, leftRC, rightRC;

  rightFW.SetFromBases("GT");
  leftFW.SetFromBases("AG");
  leftRC.SetFromBases("AC");
  rightRC.SetFromBases("CT");

  int slack = 6;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    int first = parser.AsInt(3) - 1;
    int last = parser.AsInt(4) - 1;
    string chr = parser.AsString(0);
    string ori = parser.AsString(6);
    int index = dna.NameIndex(chr);
    if (index == -1) {
      index = dna.NameIndex("chr" + chr);
    }
    const DNAVector & d = dna[index];
    DNAVector left, right;
    if (ori == "+") {
      left = leftFW;
      right = rightFW;
    } else {
      left = leftRC;
      right = rightRC;
    }
    int l = Adjust(d, first-1, left, slack);
    int r = Adjust(d, last+1, right, slack);
    
    if (l != -1) {
      l += 3;
      //if (l != first) {
	//cout << "Adjusted left from " << chr << ": " << first << " to " << l << " " << ori <<  endl;
      //cout << "Adjust left." << endl;
      first = l-1;
	
	//}
    }
    if (r != -1) {
      //cout << "Adjust right." << endl;
       //if (r != last) {
	//cout << "Adjusted right from " << chr << ": " << last << " to " << r << endl;	
      last = r-1;
	//}      
    }

    //if (first > 2 && last+2 <d.isize())
    //cout << ori << " " << d[first-2] << d[first-1] << "  " << d[last+1] << d[last+2] << endl;

    for (int i=0; i<8; i++) {
      if (i == 3) {
	cout << first+1 << "\t";
	continue;
      }
      if (i==4) {
	cout << last+1 << "\t";
	continue;
      }
      cout << parser.AsString(i) << "\t";
    }
    for (int i=8; i<parser.GetItemCount(); i++) {
      if (i > 8)
	cout << " ";
      cout << parser.AsString(i);
    }
    cout << endl;
  }
  return 0;
}
