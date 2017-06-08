#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "ryggrad/src/base/CommandLineParser.h"
#include "SearchOverlaps.h"
#include "ContScaff.h"
#include "ContScaffIO.h"
#include "ConsensOverlapUnit.h"



int main( int argc, char** argv )
{
    commandArg<string> contigCmmd("-i","input contig file");
    commandArg<string> readsCmmd("-r","input reads file in fasta format");
    commandArg<string> consCmmd("-c","input read consensus group file");
    commandArg<int> minContigCmmd("-minContig","minimum length of a single-contig scaffold to report", 200);
    commandArg<string> outFastaCmmd("-o","output fasta file", "final.fa");
    commandArg<string> outLayoutCmmd("-l","output layout file corresponding to the final fasta", "final.layout");
    commandArg<string> outReadsCmmd("-O","output list of used raw read names, provide file name if you require this information", "");
    commandArg<string> partOutDirCmmd("-readsOutDir","Output directory for recording the reads in each scaffold", "partitions");
    commandArg<bool>   gapsCmmd("-gaps","Uses gaps in alignments (use for anything other than Illumina data)", false);

    commandLineParser P(argc,argv);
    P.SetDescription("Generate Fasta file for the assembled sequences from a given contig file.");
    P.registerArg(contigCmmd);
    P.registerArg(readsCmmd);
    P.registerArg(consCmmd);
    P.registerArg(minContigCmmd);
    P.registerArg(outFastaCmmd);
    P.registerArg(outLayoutCmmd);
    P.registerArg(outReadsCmmd);
    P.registerArg(partOutDirCmmd);
    P.registerArg(gapsCmmd);
    P.parse();
  
    string contigFile    = P.GetStringValueFor(contigCmmd);
    string readsFile     = P.GetStringValueFor(readsCmmd);
    string consFile      = P.GetStringValueFor(consCmmd);
    int minContig        = P.GetIntValueFor(minContigCmmd);
    string outFastaFile  = P.GetStringValueFor(outFastaCmmd);
    string outLayoutFile = P.GetStringValueFor(outLayoutCmmd);
    string outReadFile   = P.GetStringValueFor(outReadsCmmd);
    string partOutDir    = P.GetStringValueFor(partOutDirCmmd);
    bool bUseGaps        = P.GetBoolValueFor(gapsCmmd);
 
    ConsensOverlapUnit COUnit(readsFile, consFile);

    Assembled assembly;
    ContigScaffoldIO io;
    io.Read(assembly, contigFile);
   
    if(outReadFile != "") {
      io.WriteAssembledRawReads(assembly, COUnit, outReadFile);
    }

#if defined(FORCE_DEBUG)
    io.WriteScaffoldReads(assembly, COUnit, partOutDir);
#endif
  
    LayoutSink sink;
    sink.fastaFromAssembly(outFastaFile, assembly, COUnit, minContig, bUseGaps);
    //The fasta generation step flags contigs that don't meet requirements for final sequence generation.
    //As a result the layout that will be written out as the next stage doesn't contain those discarded contigs.
    io.Write(assembly, outLayoutFile);
    return 0;
}
