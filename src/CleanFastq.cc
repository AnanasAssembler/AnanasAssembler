#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input fastq file");
  commandLineParser P(argc,argv);
  P.SetDescription("Trims a fastq file based on quality.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    // Header
    cout << ">" << parser.Line() << endl;
    parser.ParseLine();
    //Seq
    string s = parser.AsString(0);
    parser.ParseLine();
    parser.ParseLine();
    const string & q = parser.AsString(0);
    int i;
    for (i=q.size()-1; i>= 50; i--) {
      if (q[i] == '#' || s[i] == 'A')
	continue;
      break;
    }
    for (int j=0; j<=i; j++)
      cout << s[j];
    cout << endl;
  }
  return 0;
}
