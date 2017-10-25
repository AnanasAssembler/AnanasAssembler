#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/FileParser.h"
#include "ryggrad/src/general/DNAVector.h"
#include "ContScaff.h"
#include "ContScaffIO.h"

int main( int argc, char** argv )
{
  commandArg<string> fastaFileCmmd("-if","input fasta file");
  commandArg<string> layoutFileCmmd("-il","input layout file", "");
  commandArg<string> preCmmd("-p","prefix", "");
  commandArg<int> minCmmd("-m","minimum length", 0);
  commandLineParser P(argc,argv);
  P.SetDescription("Prints out the first sequence of an Ananas component.");
  P.registerArg(fastaFileCmmd);
  P.registerArg(layoutFileCmmd);
  P.registerArg(preCmmd);
  P.registerArg(minCmmd);

  P.parse();

  string fastaFileName  = P.GetStringValueFor(fastaFileCmmd);
  string layoutFileName = P.GetStringValueFor(layoutFileCmmd);
  string prefix = P.GetStringValueFor(preCmmd);
  int min = P.GetIntValueFor(minCmmd);


  ofstream fout;
  fout.open(fastaFileName + ".top");
  vecDNAVector dna;
  dna.Read(fastaFileName);
  for (int i=0; i<dna.isize(); i++) {
    string name = ">" + prefix;
    name += dna.NameClean(i);
    const DNAVector & d = dna[i];
    if (d.isize() < min)
        continue;
    if (name[name.size()-1] == '0'
        && name[name.size()-2] == '0'
        && name[name.size()-3] == '0') {
        fout << name << endl;
        for (int j=0; j<d.isize(); j++)
          fout << d[j];
        fout << endl;
    }
  }
  fout.close();

  if(layoutFileName != "") {
    Assembled assembly;
    ContigScaffoldIO io;
    io.Read(assembly, layoutFileName);
    // Main loop over scaffolds/contigs to keep only the top contig 
    for (int scaffCnt=0; scaffCnt<assembly.isize(); scaffCnt++) {
      Scaffold & currScaff = assembly[scaffCnt];
      if(currScaff.isize()>0) {
        string name = currScaff[0].Name();
        if (name[name.size()-1] != '0'
        || name[name.size()-2]  != '0'
        || name[name.size()-3]  != '0') {
          currScaff.resize(0); // If the first contig is not the original initial contig, Remove all contigs 
        } 
        currScaff.resize(1); // Only keep one from the top
      }
    }
    io.Write(assembly, layoutFileName+".top");
  }
  return 0;
}
