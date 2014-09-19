#include <string>

#include "base/CommandLineParser.h"
#include "src/SequenceMatch.h"
int main( int argc, char** argv )
{
  commandArg<string> aStringCmmd("-o","output");
  commandArg<string> aStringCmmd1("-i1","file 1");
  commandArg<string> aStringCmmd2("-i2","file 2", "");
  commandArg<string> aStringCmmd3("-i3","file 3", "");
  commandArg<string> aStringCmmd4("-i4","file 4", "");
  commandArg<string> aStringCmmd5("-i5","file 5", "");
  commandArg<string> aStringCmmd6("-i6","file 6", "");


  commandLineParser P(argc,argv);
  P.SetDescription("Collect & merge HomologyByXCorr alignments.");
  P.registerArg(aStringCmmd);
  P.registerArg(aStringCmmd1);
  P.registerArg(aStringCmmd2);
  P.registerArg(aStringCmmd3);
  P.registerArg(aStringCmmd4);
  P.registerArg(aStringCmmd5);
  P.registerArg(aStringCmmd6);


  P.parse();

  string out = P.GetStringValueFor(aStringCmmd);
  string in1 = P.GetStringValueFor(aStringCmmd1);
  string in2 = P.GetStringValueFor(aStringCmmd2);
  string in3 = P.GetStringValueFor(aStringCmmd3);
  string in4 = P.GetStringValueFor(aStringCmmd4);
  string in5 = P.GetStringValueFor(aStringCmmd5);
  string in6 = P.GetStringValueFor(aStringCmmd6);


  MultiMatches multi;
  
  cout << "Reading " << in1 << endl;
  multi.MergeRead(in1);
  if (in2 != "") {
    cout << "Reading " << in2 << endl;
    multi.MergeRead(in2);
  }
  if (in3 != "") {
    cout << "Reading " << in3 << endl;
    multi.MergeRead(in3);
  }
  if (in4 != "") {
    cout << "Reading " << in4 << endl;
    multi.MergeRead(in4);
  }
  if (in5 != "") {
    cout << "Reading " << in5 << endl;
    multi.MergeRead(in5);
  }
  if (in6 != "") {
    cout << "Reading " << in6 << endl;
    multi.MergeRead(in6);
  }

  cout << "Sorting" << endl;
  multi.Sort();
  
  cout << "Writing" << endl;
  multi.Write(out);

  return 0;

}
