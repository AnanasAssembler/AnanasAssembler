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
    commandArg<int> sizeCmmd("-minSize","minimum length of a single-contig scaffold to report", 200);
    commandArg<string> outFastaCmmd("-o","output fasta file", "final.fa");
    commandArg<string> prefixCmmd("-prefix","The prefix to add to all generated contig names", "Sample1");

    commandLineParser P(argc,argv);
    P.SetDescription("Generate Fasta file for the assembled sequences from a given contig file.");
    P.registerArg(contigCmmd);
    P.registerArg(readsCmmd);
    P.registerArg(consCmmd);
    P.registerArg(sizeCmmd);
    P.registerArg(outFastaCmmd);
    P.registerArg(prefixCmmd);
    P.parse();
  
    string contigFile  = P.GetStringValueFor(contigCmmd);
    string readsFile   = P.GetStringValueFor(readsCmmd);
    string consFile    = P.GetStringValueFor(consCmmd);
    int minSize        = P.GetIntValueFor(sizeCmmd);
    string outFile     = P.GetStringValueFor(outFastaCmmd);
    string prefix      = P.GetStringValueFor(prefixCmmd);
 
    ConsensOverlapUnit COUnit(readsFile, consFile);

    Assembled assembly;
    ContigScaffoldIO io;
    io.Read(assembly, contigFile);
  
    LayoutSink sink;
    sink.SetPrefix(prefix);
    sink.fastaFromAssembly(outFile, assembly, COUnit, minSize);
    return 0;
}
