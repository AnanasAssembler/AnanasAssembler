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
	   int stop) {
    m_start = start;
    m_stop = stop;
    m_id = id;
    m_chr = chr;
    m_line = s;
    m_gene = g;
  }


  bool operator < (const Line & l) const {
    if (m_chr != l.m_chr)
      return (m_chr < l.m_chr);
    if (m_start < l.m_start)
      return (m_start < l.m_start);
    return m_stop < l.m_stop;
  }

  const string & GetLine() const {return m_line;}
  int Start() const {return m_start;}
  int Stop() const {return m_stop;}

private:
  string m_line;
  string m_chr;
  string m_id;
  string m_gene;
  int m_start;
  int m_stop;
  

};

bool Wiggle(int a, int b, int allowed)
{
  int d = a - b;
  if (d < 0)
    d = -d;
  if (d <= allowed)
    return true;
  else
    return false;
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<int> wiggleCmmd("-s","slack difference", 0);
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);
  P.registerArg(wiggleCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  int wiggle = P.GetIntValueFor(wiggleCmmd);
  

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
    line[k].Set(parser.Line(), parser.AsString(0), s, g, parser.AsInt(3), parser.AsInt(4));
    k++;
        
  }

  line.resize(k);

  Sort(line);

  int lastStart = -1;
  int lastStop = -1;
  int all = 0;
  for (int i=0; i<line.isize(); i++) {
    if (Wiggle(lastStart, line[i].Start(), wiggle) && Wiggle(lastStop, line[i].Stop(), wiggle))
      continue;
    cout << line[i].GetLine() << endl;
    lastStart = line[i].Start();
    lastStop = line[i].Stop();
    all += lastStop - lastStart;
  }

  cerr << "Total bases: " << all << endl;
  return 0;
}
