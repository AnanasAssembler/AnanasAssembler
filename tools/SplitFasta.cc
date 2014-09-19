#include <string>

#include "base/CommandLineParser.h"
#include "src/DNAVector.h"


int main(int argc,char** argv)
{

  
  commandArg<string> aStringCmmd("-i","input fasta");
  commandArg<int> nCmmd("-n","how many");
  //commandArg<string> cStringCmmd("-o2","output file");
  commandLineParser P(argc,argv);
  P.SetDescription("Tool to test the VecDNAVector.");
  P.registerArg(aStringCmmd);
  P.registerArg(nCmmd);
  //P.registerArg(cStringCmmd);

  P.parse();

  string aString = P.GetStringValueFor(aStringCmmd);
  int howMany = P.GetIntValueFor(nCmmd);
  //string cString = P.GetStringValueFor(cStringCmmd);
  
  vecDNAVector test;
  
  cout << "Reading file..." << endl;
  test.Read(aString);

  int n = (test.isize() + howMany)/howMany;

  //out.resize(n);

  int i, j;
  for (j=0; j<howMany; j++) {
    vecDNAVector out;
    for (i=j*n; i<(j+1)*n; i++) {
      if (i >= test.isize())
	break;
      out.push_back(test[i], test.Name(i));     
    }
    char tmp[256];
    sprintf(tmp, ".%d", j);

    out.Write(aString + tmp);
  }
  /*
  for (i=n; i<2*n; i++) {
    out[i-n] = test[i];
    out.SetName(i-n, test.Name(i));
  }

  out.Write(aString + ".2");
  for (i=2*n; i<3*n; i++) {
    out[i-2*n] = test[i];
    out.SetName(i-2*n, test.Name(i));
  }

  out.Write(aString + ".3");
  */

  cout << "done!" << endl;
  //test.ReverseComplement();

  //out1.Write(aString + ".1");
  //out2.Write(aString + ".2");
  //out3.Write(aString + ".3");

  return 0;

}
  
