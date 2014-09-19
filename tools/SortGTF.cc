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
	   int start) {
    m_start = start;
    m_id = id;
    m_chr = chr;
    m_line = s;
    m_gene = g;
  }


  bool operator < (const Line & l) const {
    //if (m_chr != l.m_chr)
    //return (m_chr < l.m_chr);
    if (m_gene != l.m_gene)
      return (m_gene < l.m_gene);
    if (m_id != l.m_id)
      return (m_id < l.m_id);
    if (m_chr != l.m_chr)
      return (m_chr < l.m_chr);
    return (m_start < l.m_start);
  }

  const string & GetLine() const {return m_line;}

private:
  string m_line;
  string m_chr;
  string m_id;
  string m_gene;
  int m_start;
  

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
    line[k].Set(parser.Line(), parser.AsString(0), s, g, parser.AsInt(3));
    k++;
        
  }

  line.resize(k);

  Sort(line);
  for (int i=0; i<line.isize(); i++)
    cout << line[i].GetLine() << endl;
  return 0;
}
