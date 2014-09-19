#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"


class Annot
{
public: 
  Annot() {
    m_start = -1; 
    m_stop = -1;
  }

  void Set(const string & chr, int start, int stop) {
    m_chr = chr;
    m_start = start;
    m_stop = stop;
  }

  bool Laps(const string & chr, const int start, const int stop) {
    if (chr != m_chr)
      return false;
    if (start <= m_start && stop >= m_start)
      return true;
    if (start <= m_stop && stop >= m_stop)
      return true;
    return false;
  }

private:
  int m_start;
  int m_stop;
  string m_chr;

};


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd1("-i1","input file 1");
  commandArg<string> fileCmmd2("-i2","input file 2");
  commandLineParser P(argc,argv);
  P.SetDescription("Intersects annotations.");
  P.registerArg(fileCmmd1);
  P.registerArg(fileCmmd2);
  
  P.parse();
  
  string fileName1 = P.GetStringValueFor(fileCmmd1);
  string fileName2 = P.GetStringValueFor(fileCmmd2);
  

  //comment. ???
  FlatFileParser parser1;
  FlatFileParser parser2;
  
  parser1.Open(fileName1);

  svec<Annot> annot;
  int k = 0;

  while (parser1.ParseLine()) {
    if (k >= annot.isize())
      annot.resize(k+1000000);

    annot[k].Set(parser1.AsString(0), parser1.AsInt(1), parser1.AsInt(2));
    k++;
  }

  annot.resize(k);

  parser2.Open(fileName2);
  while (parser2.ParseLine()) {
    for (int i=0; i<annot.isize(); i++) {
      if (annot[i].Laps(parser2.AsString(0), parser2.AsInt(1), parser2.AsInt(2))) {
	cout << parser2.Line() << endl;
	break;
      }
     
    }
  }

  return 0;
}
