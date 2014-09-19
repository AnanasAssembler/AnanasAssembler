#include <string>

#include "base/CommandLineParser.h"
#include "src/Cola/NSGAaligner.h"
#include "src/Satsuma/RefineSatsuma.h"


int main(int argc,char** argv)
{
  
  commandArg<string> aStringCmmd("-q","query sequence");
  commandArg<string> bStringCmmd("-t","target sequence");
  commandArg<string> satsumaCmmd("-s","satsuma summary file");
  commandArg<int> outputModeCmmd("-m","Output mode - 0:Full 1:Summary", 1);
  commandArg<string> realignOutCmmd("-r","realignment output file");
  commandArg<string> gapOutCmmd("-g","gap output file");
  commandLineParser P(argc,argv);
  P.SetDescription("Realigns global alignment in satsuma format with Cola.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(satsumaCmmd);
  P.registerArg(outputModeCmmd);
  P.registerArg(realignOutCmmd);
  P.registerArg(gapOutCmmd);

  P.parse();

  string aString        = P.GetStringValueFor(aStringCmmd);
  string bString        = P.GetStringValueFor(bStringCmmd);
  string satsumaFile    = P.GetStringValueFor(satsumaCmmd);
  int outputMode        = P.GetIntValueFor(outputModeCmmd); 
  string realignOutFile = P.GetStringValueFor(realignOutCmmd);
  string gapOutFile     = P.GetStringValueFor(gapOutCmmd);

  vecDNAVector query, target;
  query.Read(aString);
  target.Read(bString);

  ofstream realignFout, gapFout;
  realignFout.open(realignOutFile.c_str());
  gapFout.open(gapOutFile.c_str());
  RefineSatsuma refinementAligner(target, query,
                                  satsumaFile, realignFout,
                                  gapFout, outputMode);
  refinementAligner.alignAll();
  realignFout.close();
  gapFout.close();
  return 0;
}
