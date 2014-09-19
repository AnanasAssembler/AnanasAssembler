#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"

bool Same(const string & a, const string & b)
{
  int i;
  int n = strlen(a.c_str());
  if (strlen(b.c_str()) < n)
    n = strlen(b.c_str());

  for (i=0; i<n; i++) {
    if (a[i] >= '0' && a[i] <= '9' && b[i] >= '0' && b[i] <= '9')
      return true;
    if (a[i] != b[i])
      return false;
  }
  return false;
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  string last;
  string lt;
  double ex = 0.;
  bool b = false;
  while (parser.ParseLine()) {
    const string & s = parser.AsString(1);
    if (Same(s, last)) {
      cout << last << " " << lt << " " << ex << endl;
      b = true;
    } else {
      if (b)
	cout << last << " " << lt << " " << ex << endl << endl;
      b = false;
    }
    
    last = s;
    lt = parser.AsString(0);
    ex = parser.AsFloat(2);
  }
  return 0;
}
