#include <string>

#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/Saguaro/HMMDistance.h"


double StupidDist(const string & v1, const string & v2) 
{
  if (v1 == "?" || v2 == "?")
    return -1.;

  if (v1 == v2) 
    return 0.;

  if (v1 == "a" || v1 == "c" || v1 == "g" || v1 == "t")
    return 0.5;
  if (v2 == "a" || v2 == "c" || v2 == "g" || v2 == "t")
    return 0.5;

  return 1.;
}


int main(int argc,char** argv)
{

  
  commandArg<string> aStringCmmd("-i","input");
  commandArg<string> nStringCmmd("-n","name");
  commandArg<string> bStringCmmd("-o","output");
  commandLineParser P(argc,argv);
  P.SetDescription("Converts data into a HMM-digestable file of features");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(nStringCmmd);

  P.parse();

  string in = P.GetStringValueFor(aStringCmmd);
  string out = P.GetStringValueFor(bStringCmmd);
  string name = P.GetStringValueFor(nStringCmmd);
 
  
  cout << "OBSOLETE - DO NOT USE!!!" << endl;
  return -1;



  HMMFeatureVector f;

  FlatFileParser parserIn;
  parserIn.Open(in);

    
  int k = 0;

  
  HMMFeature feature;
  feature.SetName(name);

  int i, j;
  int numFish = 0;

  while (parserIn.ParseLine()) {
    if (parserIn.GetItemCount() == 0)
      continue;
    if (parserIn.AsString(0) == "Adding") {
      f.SetName(numFish, parserIn.AsString(2));
      numFish++;
      continue;
    }
    if (parserIn.AsString(0) != "SNP")
      continue;
    int pos = parserIn.AsInt(1);

    feature.SetPosition(pos);

    svec<string> line;
    for (i=2; i<parserIn.GetItemCount()-1; i++) {
      line.push_back(parserIn.AsString(i));
    }

    int n = line.isize();
    feature.resize(n * n);

    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
	feature[i*n + j] = StupidDist(line[i], line[j]);
      }
    }
    

    if (k >= f.isize())
      f.resize(k + 128000);

    f[k] = feature;
    k++;
  }

  f.resize(k);


  f.Write(out);

  return 0;

}
  
