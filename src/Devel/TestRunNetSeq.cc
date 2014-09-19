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

int main( int argc, char** argv )
{

  /* commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing HMM expression network search.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);*/
  //cout << RandomFloat(1.) << endl;
  //cout << RandomFloat(1.) << endl;
  //cout << RandomFloat(1.) << endl;
  //return 0;

  Network network;

  svec<ExpVector> model;
  model.push_back(Make(1., 0., 2.));
  model.push_back(Make(0., 3., 1.));


  network.AddNode("A", Make(1., 0., 2.));
  network.AddNode("B", Make(1., 0., 2.));
  network.AddNode("D", Make(1., 0., 2.));
  network.AddNode("C", Make(0., 3., 1.));
  network.AddNode("E", Make(0., 3., 1.));
  network.AddNode("F", Make(0., 3., 1.));
  network.AddNode("G", Make(1., 0., 2.));
  network.AddNode("H", Make(1., 0., 2.));

  /*
  network.Connect("A", "B");
  network.Connect("B", "D");
  network.Connect("C", "D");
  network.Connect("A", "C");
  network.Connect("E", "C");
  network.Connect("D", "A");
  */

  network.Connect("A", "B");
  network.Connect("B", "D");
  network.Connect("D", "C");
  network.Connect("C", "E");
  network.Connect("F", "C");
  network.Connect("A", "G");
  network.Connect("G", "B");
  network.Connect("A", "H");
  network.Connect("H", "E");
   
  network.Process(model);
  //cout << "FINAL" << endl;
  network.Print();
  
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
