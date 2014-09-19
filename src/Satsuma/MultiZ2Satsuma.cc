#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"

void Process(const string & file, const string & target, const string & query)
{
 FlatFileParser parser;
  
  parser.Open(file);

  string tChr, qChr;
  int startT = -1;
  int startQ = -1;
  int stopT = -1;
  int stopQ = -1;
  string oriT, oriQ;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.AsString(0) == "a") {
      if (qChr != "" && tChr != "") {
	cout << tChr << "\t" << startT << "\t" << stopT << "\t";
	cout << qChr << "\t" << startQ << "\t" << stopQ << "\t1.0\t";
	if (oriT == oriQ)
	  cout << "+" << endl;
	else
 	  cout << "-" << endl;
      }

      oriT = "";
      oriQ = "";
      tChr = "";
      qChr = "";
    }
    if (parser.AsString(0) != "s") {
      continue;
    }
    StringParser p;
    p.SetLine(parser.AsString(1), ".");
    if (p.AsString(0) == target) {
      tChr = p.AsString(1);
      oriT = parser.AsString(4);
      if (oriT == "+") {
	startT = parser.AsInt(2);
	stopT = startT + parser.AsInt(3);
      } else {
	stopT = parser.AsInt(5) - parser.AsInt(2);
	startT = stopT - parser.AsInt(3);
      }
    }
    if (p.AsString(0) == query) {
      qChr = p.AsString(1);
      oriQ = parser.AsString(4);
      if (oriQ == "+") {
	startQ = parser.AsInt(2);
	stopQ = startQ + parser.AsInt(3);
      } else {
	stopQ = parser.AsInt(5) - parser.AsInt(2);
	startQ = stopQ - parser.AsInt(3);
      }
    }
    
  }
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","list of input files");
  commandArg<string> tCmmd("-t","target species");
  commandArg<string> qCmmd("-q","query species");

  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);
  P.registerArg(tCmmd);
  P.registerArg(qCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string target = P.GetStringValueFor(tCmmd);
  string query = P.GetStringValueFor(qCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  while (parser.ParseLine()) {
    Process(parser.Line(), target, query);    
  }

  return 0;
}
