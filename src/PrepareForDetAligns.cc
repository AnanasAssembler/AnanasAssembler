#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Parses a satsuma summary summary file to feed into DPAlignSatsuma.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  string lastT, lastQ, lastOri;
  int bridge = 100;
  int max = 2000;
  
  int last_startT = -1;
  int last_stopT = -1;
  int last_startQ = -1;
  int last_stopQ = -1;

  while (parser.ParseLine()) {
    const string  & t = parser.AsString(0);
    const string  & q = parser.AsString(3);

    int startT = parser.AsInt(1);
    int stopT = parser.AsInt(2);
    int startQ = parser.AsInt(4);
    int stopQ = parser.AsInt(5);
    double ident = parser.AsFloat(6);
    const string & ori = parser.AsString(7);
    int len = last_stopT - last_startT;

    if (len < max && startT > last_startT && startT < last_stopT + bridge && t == lastT && q == lastQ && ori == lastOri) {
      
      if (ori == "+") {
	if (startQ > last_startQ && startQ < last_stopQ + bridge) {
	  last_stopT = stopT;
	  last_stopQ = stopQ;
	  continue;
	}
      } else {
	if (stopQ > last_startQ - bridge && stopQ < last_startQ) {
	  last_startQ = startQ;
	  last_stopT = stopT;
	  continue;
	}
      }
      
    } 

    if (last_startT > 0) {
      cout << lastT << "\t" << last_startT << "\t" << last_stopT << "\t" << lastQ << "\t";
      cout << last_startQ << "\t" << last_stopQ << "\t1.000\t" << lastOri << endl;
    }
    
    last_startT = startT;
    last_stopT = stopT;
    last_startQ = startQ;
    last_stopQ = stopQ;
    lastT = t;
    lastQ = q;
    lastOri = ori;
  }  

  return 0;
}
