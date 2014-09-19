#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "base/CommandLineParser.h"
#include "src/Ananas/SearchOverlaps.h"
#include "src/Ananas/ContScaff.h"
#include "src/Ananas/ContScaffIO.h"
#include "src/Ananas/ConsensOverlapUnit.h"



int main( int argc, char** argv )
{
  commandArg<string> contigCmmd("-i","input contig file");
  commandArg<string> readsCmmd("-r","input reads file in fasta format");
  commandArg<string> consCmmd("-c","input read consensus group file");
  commandArg<string> outFastaCmmd("-o","output fasta file", "final.fa");

  commandLineParser P(argc,argv);
  P.SetDescription("Generate Fasta file for the assembled sequences from a given contig file.");
  P.registerArg(contigCmmd);
  P.registerArg(readsCmmd);
  P.registerArg(consCmmd);
  P.registerArg(outFastaCmmd);
  P.parse();
  
  string contigFile  = P.GetStringValueFor(contigCmmd);
  string readsFile   = P.GetStringValueFor(readsCmmd);
  string consFile    = P.GetStringValueFor(consCmmd);
  string outFile     = P.GetStringValueFor(outFastaCmmd);
 
  ConsensOverlapUnit COUnit(readsFile, consFile);

  Assembled assembly;
  ContigScaffoldIO io;
  io.Read(assembly, contigFile);
  
  LayoutSink sink;
  //sink.SetPrefix(i); //TODO
  sink.fastaFromAssembly(outFile, assembly, COUnit);
  return 0;
}
