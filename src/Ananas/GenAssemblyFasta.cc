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
    commandArg<int> minContigCmmd("-minContig","minimum length of a single-contig scaffold to report", 200);
    commandArg<string> outFastaCmmd("-o","output fasta file", "final.fa");
    commandArg<string> prefixCmmd("-prefix","The prefix to add to all generated contig names", "Sample1");
    commandArg<string> partOutDirCmmd("-readsOutDir","Output directory for recording the reads in each scaffold", "partitions");

    commandLineParser P(argc,argv);
    P.SetDescription("Generate Fasta file for the assembled sequences from a given contig file.");
    P.registerArg(contigCmmd);
    P.registerArg(readsCmmd);
    P.registerArg(consCmmd);
    P.registerArg(minContigCmmd);
    P.registerArg(outFastaCmmd);
    P.registerArg(prefixCmmd);
    P.registerArg(partOutDirCmmd);
    P.parse();
  
    string contigFile  = P.GetStringValueFor(contigCmmd);
    string readsFile   = P.GetStringValueFor(readsCmmd);
    string consFile    = P.GetStringValueFor(consCmmd);
    int minContig      = P.GetIntValueFor(minContigCmmd);
    string outFile     = P.GetStringValueFor(outFastaCmmd);
    string prefix      = P.GetStringValueFor(prefixCmmd);
    string partOutDir  = P.GetStringValueFor(partOutDirCmmd);
 
    ConsensOverlapUnit COUnit(readsFile, consFile);

    Assembled assembly;
    ContigScaffoldIO io;
    io.Read(assembly, contigFile);

//#if defined(FORCE_DEBUG)
    io.WriteScaffoldReads(assembly, COUnit, partOutDir);
//#endif
  
    LayoutSink sink;
    sink.SetPrefix(prefix);
    sink.fastaFromAssembly(outFile, assembly, COUnit, minContig);
    return 0;
}
