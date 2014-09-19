#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/Optical/OptiMap.h"

void Read(OpticalMap & m, const string & fileName, double scale = 1.0, bool rc = false)
{
  FlatFileParser parser;
  
  parser.Open(fileName);
  string name;
  int i;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.GetItemCount() == 1) {
      name = parser.AsString(0);
      continue;
    }
    OpticalSeq s;
    s.SetName(name);
    double vv = 0.;
    for (i=2; i<parser.GetItemCount(); i++) {
      double v = scale*1000*parser.AsFloat(i);
      vv += v;
      s.push_back((int)(vv+0.5));
    }
    m.push_back(s);
    if (rc) {
      OpticalSeq s2;
      s2.SetName(name + "_reverse");
      vv = 0.;
      for (i=parser.GetItemCount()-1; i >=2; i--) {
	double v = scale*1000*parser.AsFloat(i);
	vv += v;
	s2.push_back((int)(vv+0.5));
      }
      m.push_back(s2);
    }
  }  
}

int main( int argc, char** argv )
{

  commandArg<string> qCmmd("-q","input file 1");
  commandArg<string> tCmmd("-t","input file 2");
  commandArg<double> sCmmd("-s","query scale", 1.0);
  commandLineParser P(argc,argv);
  P.SetDescription("Aligns optical maps against each other via cross-correlation.");
  P.registerArg(qCmmd);
  P.registerArg(tCmmd);
  P.registerArg(sCmmd);
  
  P.parse();
  
  string fileNameQ = P.GetStringValueFor(qCmmd);
  string fileNameT = P.GetStringValueFor(tCmmd);
  double scale = P.GetDoubleValueFor(sCmmd);
  
  OpticalMap t, q;
  Read(t, fileNameT);
  cout << "Average distance target: " << t.Average() << endl;
  Read(q, fileNameQ, scale, true);
  cout << "Average distance query:  " << q.Average() << endl;
  int i, j;

  OptiAligner align;

  //svec<double> self
  //for (i=0; i<t.isize(); i++) {

  for (j=0; j<q.isize(); j++) {
    double max = 0.;
    int index = -1;
    for (i=0; i<t.isize(); i++) {
      if (t[i].Name() == q[j].Name())
	continue;
      double s = align.Align(t[i], q[j], false);
      if (s > max) {
	max = s;
	index = i;
      }
    }
    if (index != -1) {
      cout << "Best match: " << index << " score: " << max << " " << q[j].Name() << endl;
      align.Align(t[index], q[j], true);
    }
  }
  
  return 0;
}
