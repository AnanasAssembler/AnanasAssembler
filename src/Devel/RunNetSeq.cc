#define FORCE_DEBUG

#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/Devel/ExpNetwork.h"
#include "base/RandomStuff.h"

ExpVector Make(double a, double b, double c) {
  ExpVector v;
  v.resize(3);
  v[0] = a + RandomFloat(2.);
  v[1] = b + RandomFloat(2.);
  v[2] = c + RandomFloat(2.);
  cout << v[0] << " " << v[1] << " " << v[2] << endl;
  return v;
}


void ReadConn(const string & fileName, Network & network, const svec<string> & gene, const svec<string> & symbol)
{
  FlatFileParser parser;
  
  parser.Open(fileName);
  int i, j;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    const string & a = parser.AsString(0);
    const string & b = parser.AsString(1);
    //cout << "Trying: " << a << " " << b << endl;
    for (i=0; i<gene.isize(); i++) {
      if (symbol[i] == a) {
	for (j=0; j<gene.isize(); j++) {
	  if (symbol[j] == b) {
	    //cout << "Connecting " << gene[i] << " -> " << gene[j] << endl;
	    network.Connect(gene[i], gene[j]);
	  }
	}
      }
    }
  }
}

void Read(const string & fileName, Network & network, svec<ExpVector> & model)
{
  FlatFileParser parser;
  
  parser.Open(fileName);
  int i, j;
  ExpVector all;
  double div = 0.;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    svec<double> s;
    double sum = 1.;
    for (i=2; i<parser.GetItemCount(); i++) {
      s.push_back(parser.AsFloat(i));
      sum += parser.AsFloat(i);
    }
    
    sum = sum/(double)s.isize();
    
    for (i=0; i<s.isize(); i++) {
      s[i] /= sum;
    }

    ExpVector v;
    
    for (i=0; i<s.isize(); i++) {
      for (j=i+1; j<s.isize(); j++) {
	double val = (s[i]-s[j]) * (s[i]-s[j]);
	v.push_back(val);
	
	//cout << val << endl;
      }
    }

    if (all.isize() == 0)
      all.resize(v.isize(), 0);
    for (i=0; i<v.isize(); i++)
      all[i] += v[i];

    div++;

    if (parser.AsString(1) == "INS") {
      cout << "push INS" << endl;
      model.push_back(v);
      for (int x=0; x<v.isize(); x++)
	cout << v[x] << endl;
    }
    if (parser.AsString(1) == "PDX1") {
      cout << "push PDX1" << endl;
      model.push_back(v);
      for (int x=0; x<v.isize(); x++)
	cout << v[x] << endl;
    }
    if (parser.AsString(1) == "GNAS") {
      cout << "push GNAS" << endl;
      model.push_back(v);
    }
    if (parser.AsString(1) == "ADCY6") {
       cout << "push ADCY6" << endl;
       model.push_back(v);     
    }
    network.AddNode(parser.AsString(1), v);
  }

  for (i=0; i<all.isize(); i++)
    all[i] /= div;
  model.push_back(all);   
}

void ReadMap(const string & fileName, svec<string> & gene, svec<string> & symbol)
{
 FlatFileParser parser;
  
  parser.Open(fileName);
  int i, j;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    gene.push_back(parser.AsString(0));
    symbol.push_back(parser.AsString(1));
  }
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","expression file");
  commandArg<string> mapCmmd("-m","map file");
  commandArg<string> connCmmd("-c","connection file");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing HMM expression network search.");
  P.registerArg(fileCmmd);
  P.registerArg(mapCmmd);
  P.registerArg(connCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string mapName = P.GetStringValueFor(mapCmmd);
  string connName = P.GetStringValueFor(connCmmd);

  int i, j;

  Network network;
  svec<ExpVector> model;

  Read(fileName, network, model);
  cout << "Model size: " << model.isize() << endl;
  svec<string> gene, symbol;
  ReadMap(mapName, gene, symbol);
  ReadConn(connName, network, gene, symbol);

  /* for (i=0; i<network.isize(); i++) {
    cout << network[i].Name() << " "; 
    cout << network[i].Dist(model[0]) << " ";
    cout << network[i].Dist(model[1]) << " ";
    cout << network[i].Dist(model[2]) << " ";
    cout << endl;
    }*/

  //return 0;

  //model.push_back(Make(1., 0., 2.));
  //model.push_back(Make(0., 3., 1.));
  //for (i=0; i<network.



  network.Process(model);
  //cout << "FINAL" << endl;
  network.Print();
  
  cout << "Listing: " << endl;
  for (i=0; i<network.isize(); i++) {
    cout << network[i].Name();
    for (j=0; j<gene.isize(); j++) {
      if (gene[j] == network[i].Name())
	cout << "\t" << symbol[j]; 
    }
    cout << " best: " << network.Model(i) << endl;
  }
  /*
  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
      }*/

  return 0;
}
