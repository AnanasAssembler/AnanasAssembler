#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"


int GeneCount(int &tA, int &tB, string & spec, int & hi, string & line, const string & file) {
  string prefix;
  if (strstr(file.c_str(), "homo_sapiens") != NULL)
    prefix = "hsa";
  if (strstr(file.c_str(), "gallus_gallus") != NULL)
    prefix = "gga";
  if (strstr(file.c_str(), "gorilla_gorilla") != NULL)
    prefix = "ggo";
  if (strstr(file.c_str(), "macaca_mulatta") != NULL)
    prefix = "mml";
  if (strstr(file.c_str(), "monodelphis_domestica") != NULL)
    prefix = "mdo";
  if (strstr(file.c_str(), "mus_musculus") != NULL)
    prefix = "mmu";
  if (strstr(file.c_str(), "ornithorhynchus_anatinus") != NULL)
    prefix = "oan";
  if (strstr(file.c_str(), "pan_troglodytes") != NULL)
    prefix = "ptr";
  if (strstr(file.c_str(), "pongo_abelii") != NULL) {
    prefix = "ppa";
    return -1;
  }

  if (prefix == "") {
    cout << "ERROR: " << file << endl;
    return 0;
  }

  string path = "../../";
  path += prefix;
  path += "/homology/groups/";
  path += file;

  FILE * pTest = fopen(path.c_str(), "r");
  if (pTest != NULL) {
    fclose(pTest);
  } else {
    path = "../../";
    path += prefix;
    path += "/homology/singletons/";
    path += file;
      
  }

  spec = "unspecific";
  FlatFileParser parser;

  parser.Open(path);
  parser.ParseLine();

  line = parser.Line();
  parser.ParseLine();
  int all = parser.AsInt(2);
  hi = 0.;
  parser.ParseLine();

  svec<int> allmax;
  allmax.resize(11, 0);

  while (parser.ParseLine()) {
    StringParser p;
    p.SetLine(parser.Line(), ";");
    bool b = false;
    int maxIndex = -1;
    double max = 0.;
    if (p.GetItemCount() < 14)
      continue;

    //tA = p.AsFloat(5);
    //tB = p.AsFloat(9);
    for (int i=5; i<p.GetItemCount(); i++) {
      if (p.AsFloat(i) > 2.0)
	b = true;
      if (p.AsFloat(i) > max) {
	max = p.AsFloat(i);
	maxIndex = i-5;
      }
     
    }
    if (maxIndex == 1)
      maxIndex = 0;
    if (maxIndex == 3)
      maxIndex = 2;
    if (maxIndex == 5)
      maxIndex = 4;
    if (maxIndex == 7)
      maxIndex = 6;

    if (maxIndex >= 0)
      allmax[maxIndex]++;
    if (b)
      hi++;
  }

  for (int j=0; j<allmax.isize(); j++) {
    if (allmax[j] > all/2) {
      switch(j) {
      case 0:
	spec = "brain";
	break;
      case 2:
	spec = "cerebellum";
	break;
      case 4:
	spec = "heart";
	break;
      case 6:
	spec = "kidney";
	break;
      case 8:
	spec = "liver";
	break;
      case 9:
	spec = "testis";
	break;
      }
    }
  }


  return all;

}


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Finds gene expansions/contractions.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  

  //comment. ???
  FlatFileParser p;
  
  p.Open(fileName);

  p.ParseLine();

  int i;
  

  while (p.ParseLine()) {
    StringParser parser;
    parser.SetLine(p.Line(), ",");

    if (parser.GetItemCount() == 0)
      continue;
    svec<int> count;
    svec<int> hi;
    svec<string> line;
    svec<string> spec;
    svec<int> tA;
    svec<int> tB;
    hi.resize(parser.GetItemCount(), 0);
    count.resize(parser.GetItemCount(), 0);
    line.resize(parser.GetItemCount());
    spec.resize(parser.GetItemCount());
    tA.resize(parser.GetItemCount(), 0);
    tB.resize(parser.GetItemCount(), 0);

    for (i=0; i<parser.GetItemCount(); i++) {
    
      count[i] = GeneCount(tA[i], tB[i], spec[i], hi[i], line[i], parser.AsString(i));
    }
    int n = count[0];
    bool b = false;
    for (i=1; i<count.isize(); i++) {
      if (count[i] != n && count[i] != -1)
	b = true;
    }
    if (b) {
      cout << "Expansion/contraction:" << endl;
    } else {      
      cout << "All the same: " << n << endl;
    }
    for (i=0; i<parser.GetItemCount(); i++) {
      cout << parser.AsString(i) << "\t" << count[i] << "\tpresent: " << hi[i] << "\t" << spec[i] << "\t" <<line[i] << endl;
      /*      if (count[i] == 1) {
	cout << "single ";
      } else {
	cout << "mult ";
      }
      cout << tA[i] << "\t" << tB[i] << endl;*/
    }
    cout << endl;
  }
  return 0;
}
