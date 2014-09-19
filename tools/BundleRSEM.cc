#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"

class Input
{

  public:
  Input() {}

  void Set(const string & n) {
    parser.Open(n);
    parser.ParseLine();
    char tmp[1024];
    strcpy(tmp, n.c_str());
    for (int i=0; i<strlen(tmp); i++) {
      if (tmp[i] == '/')
	tmp[i] = 0;
    }
    name = tmp;

  }

  const string & Name() const {return name;}


  bool Next(double & f, string & gene) {
    if (!parser.ParseLine())
      return false;
    if (parser.GetItemCount() == 0)
      return false;
    f = parser.AsFloat(6);
    gene = parser.AsString(0);
    return true;
  }
private:

  FlatFileParser parser;
  string name;
}; 
 
int main( int argc, char** argv )
{

  /*commandArg<bool> emptyCmmd("-e","include sequences w/ 0 expression");
  commandLineParser P(argc,argv);
  P.SetDescription("Combining RSEM results into a single file.");
  P.registerArg(emptyCmmd);
  
  P.parse();
  
  bool bEmpty = P.GetBoolValueFor(emptyCmmd);*/
  

  system("ls */RSEM.genes.results > rsem.tmp");
 
  //comment. ???
  FlatFileParser parser;
  
  parser.Open("rsem.tmp");

  svec<Input> all;
  svec<string> names;
  int i;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    for (i=0; i<parser.GetItemCount(); i++)
      names.push_back(parser.AsString(i));
  }

  all.resize(names.isize());
  cout << "gene";
  for (i=0; i<all.isize(); i++) {
    all[i].Set(names[i]);
    cout << "\t" << all[i].Name();
  }
  cout << endl;

  while (1) {
    bool bEnd = false;
    for (i=0; i<all.isize(); i++) {
      string g;
      double f;
      if (!all[i].Next(f, g)) {
	bEnd = true;
	break;
      }
      if (i == 0)
	cout << g;
      cout << "\t" << f;	
    }
    if (bEnd)
      break;
    cout << endl;
  };

  return 0;
}
