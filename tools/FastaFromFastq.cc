
#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include <math.h>

int ConvertQ(char c)
{
  double d = (double)c - 64.;
  d /= 10.;

  return (int) (0.5 + 10. * log(1.+ exp(d * log(10)))/ log(10.) ); 
}


bool Bad(const string & s) {
  return false;

  int i;
  int a = 0;
  int c = 0;
  int g = 0;
  int t = 0;
  int nn = 0;
  const char * p = s.c_str();
  int n = strlen(p);
  for (i=0; i<n; i++) {
    switch(p[i]) {
    case 'A':
      a++;
      break;
    case 'C':
      c++;
      break;
    case 'G':
      g++;
      break;
    case 'T':
      t++;
      break;
    default:
      nn++;
      break;

    }

  }
  if (nn > 2)
    return true;
  
  n -= nn;
  int max = 10;
  if (n - a < max)
    return true;
  if (n - c < max)
    return true;
  if (n - g < max)
    return true;
  if (n - t < max)
    return true;
  if (nn > 1)
    return true;

  return false;
}


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input fastq file");
  commandArg<string> suffCmmd("-a","append suffix to read name", "");
  commandArg<int> minCmmd("-min","minimum sequence length", 0);
  commandLineParser P(argc,argv);
  P.SetDescription("Extracts the fasta portion.");
  P.registerArg(fileCmmd);
  P.registerArg(suffCmmd);
  P.registerArg(minCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string suff = P.GetStringValueFor(suffCmmd);
  int min = P.GetIntValueFor(minCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  int k = 0;
  int skipped = 0;

  while (parser.ParseLine()) {
    string name = parser.Line();

    //cout << ">" << parser.Line() << endl;
    parser.ParseLine();
    string seq = parser.Line();
    //cout << parser.Line() << endl;
    parser.ParseLine();
    parser.ParseLine();
    if (parser.GetItemCount() == 0)
      continue;
    
    int n = strlen(seq.c_str());
    if (n < min) {
      skipped++;
      continue;
    }


    //const char * p = parser.Line().c_str();
    //int n = strlen(p);

    /*
    double avg = 0;
    for (int i=0; i<n; i++)
      avg += ConvertQ(p[i]);

    avg /= (double)n;
    */
    /*
    if (Bad(seq)) {
      skipped++;
      continue;
    }*/

    k++;
    cout << ">" << name << suff /*<< " k=" << k << " skip=" << skipped*/ << endl;
    cout << seq << endl;
    

  }
  return 0;
}
