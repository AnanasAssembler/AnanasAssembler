#include <string>
#include "base/CommandLineParser.h"
#include "src/SequenceStream.h"



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the sequence stream.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  SequenceStream stream;
  stream.Open(fileName);

  DNAVector d;
  string name;
  while(stream.GetNext(d, name)) {
    cout << name << endl;
    for (int i=0; i<d.isize(); i++)
      cout << d[i];
    cout << endl;
  }


  return 0;
}
