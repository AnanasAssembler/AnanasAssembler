#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"

void GetID(string & id, FlatFileParser & parser) 
{
  id = "";
  while (parser.ParseLine()) {
    //cout << "Parsing " << parser.Line() << endl;
    if (parser.GetItemCount() == 0) {
      //cout << "Done." << endl;
      return;
    }
    if (parser.AsString(0) == "FULL_SENSE" ||
	parser.AsString(0) == "PARTIAL_SENSE") {
      for (int i=1; i<parser.GetItemCount(); i++) {
	cout << parser.AsString(i) << "\tTARGET_HIT" << endl;
      }
    }
    if (parser.AsString(0) == "FULL_SENSE") {
      if (parser.GetItemCount() > 1) {
	id = parser.AsString(1);
	//return;
      }
    }
    if (parser.AsString(0) == "PARTIAL_SENSE") {
      if (parser.GetItemCount() > 1) {
	id = parser.AsString(1);       
      }
      //return;
    }
  }
  
}

void PrintID(FlatFileParser & parser, const string & type) 
{
 
  while (parser.ParseLine()) {
    //cout << "Parsing " << parser.Line() << endl;
    if (parser.GetItemCount() == 0) {
      //cout << "Done." << endl;
      return;
    }
    /*if (parser.AsString(0) == "FULL_SENSE" ||
	parser.AsString(0) == "PARTIAL_SENSE") {
      for (int i=1; i<parser.GetItemCount(); i++) {
	cout << parser.AsString(i) << "\tTARGET_HIT" << endl;
      }
      }*/
    if (parser.AsString(0) == type) {
      cout << type << "\t" << parser.AsString(1) << endl;      
    }
  }
  
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  
  //cout << "List of matches in the target is written to: target_matches.out" << endl;
  //cout << "List of matches in the query is written to: query_matches.out" << endl;
 
  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  //FILE * pTarget = fopen("target_matches.out", "w");
  //FILE * pQuery = fopen("query_matches.out", "w");

 

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.AsString(0) == ">>") {
      //cout << "Get: " << parser.Line() << endl;
      string name = parser.AsString(1);
      int full_sense   = parser.AsInt(3);
      int full_anti    = parser.AsInt(4);
      int part_sense   = parser.AsInt(5);
      int part_anti    = parser.AsInt(6);
      int intron_sense = parser.AsInt(7);
      int intron_anti  = parser.AsInt(8);
      int other_sense  = parser.AsInt(9);
      int other_anti   = parser.AsInt(10);

      string id;
      //GetID(id, parser);
      if (full_sense > 0) {
	
	cout << name << "\t";
	PrintID(parser, "FULL_SENSE");
	//if (id != "")
	//  cout << "\t" << id /*<< "\tKNOWN"*/;
	//cout << endl;
	continue;
      }
      if (part_sense > 0) {
	cout << name << "\t";
	PrintID(parser, "PARTIAL_SENSE");
	//if (id != "")
	//  cout  << "\t" << id /*<< "\tKNOWN"*/;
	//cout << endl;
	continue;
      }
      if (full_anti > 0) {
	cout << name << "\t";
	PrintID(parser, "FULL_ANTI");
	continue;
      }
      if (part_anti > 0) {
	cout << name << "\t";
	PrintID(parser, "PARTIAL_ANTI");
	continue;
      }
      if (intron_sense > 0) {
	cout << name << "\t";
	PrintID(parser, "INTRONIC_SENSE");
	continue;
      }
      if (intron_anti > 0) {
	cout << name << "\t";
	PrintID(parser, "INTRONIC_ANTI");
	continue;
      }
      if (other_sense > 0) {
	cout << name << "\t";
	PrintID(parser, "OTHER_SENSE");
	continue;
      }
      if (other_anti > 0) {
	cout << name << "\t";
	PrintID(parser, "OTHER_ANTI");
	continue;
      }
 
      continue;
    }
  }
  return 0;
}
