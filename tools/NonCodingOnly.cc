#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"


class Line
{
public:
  Line() {}

  void Set(const string & s,
	   const string & chr,
	   const string & id,
	   const string & g,
	   int start,
	   int stop,
	   const string & t) {
    m_start = start;
    m_stop = stop;
    m_id = id;
    m_chr = chr;
    m_line = s;
    m_gene = g;
    m_type = t;
  }


  bool operator < (const Line & l) const {
    return m_gene < l.m_gene;
  }

  const string & GetLine() const {return m_line;}
  int Start() const {return m_start;}
  int Stop() const {return m_stop;}
  const string& Gene() const {return m_gene;}
  const string& Type() const {return m_type;}

private:
  string m_line;
  string m_chr;
  string m_id;
  string m_gene;
  int m_start;
  int m_stop;
  string m_type;

};


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

  svec<Line> line;

  int k = 0;

  while (parser.ParseLine()) {
    string s, g;
    for (int i=0; i<parser.GetItemCount(); i++) {
      if (parser.AsString(i) == "transcript_id")
	s = parser.AsString(i+1);
    }
    for (int i=0; i<parser.GetItemCount(); i++) {
      if (parser.AsString(i) == "gene_id")
	g = parser.AsString(i+1);
    }
    if (line.isize() <= k)
      line.resize(k+2000000);
    line[k].Set(parser.Line(), parser.AsString(0), s, g, parser.AsInt(3), parser.AsInt(4), parser.AsString(2));
    k++;
        
  }

  line.resize(k);

  Sort(line);

  string last;
  int lastStart = 0;
  int all = 0;
  bool bNon = true;
  for (int i=0; i<line.isize(); i++) {
    if (last != line[i].Gene()) {
      if (bNon) {
	for (int j= lastStart; j<i; j++) {
	  cout << line[j].GetLine() << endl;
	  all += line[j].Stop() - line[j].Start();
	}
      } else {
	//cout << "Ignore" << endl;
      }
      lastStart = i;         
      bNon = true;
    }
    if (line[i].Type() == "CDS")
      bNon = false;
  }

  cerr << "Total bases: " << all << endl;
  return 0;
}
