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
    commandArg<string> pairCmmd("-p","input read pair/size info file");
    commandArg<string> consCmmd("-c","input read consensus group file");
    commandArg<string> scaffCmmd("-s","scaffolds file");
    commandArg<string> outCmmd("-o","output name ");
    commandLineParser P(argc,argv);
    P.SetDescription("Modifies Read counts of given contigs/scaffold from consenus reads to raw reads. ");
    P.registerArg(pairCmmd);
    P.registerArg(consCmmd);
    P.registerArg(scaffCmmd);
    P.registerArg(outCmmd);
 
    P.parse();
  
    string pairSzFileFile = P.GetStringValueFor(pairCmmd);
    string consFile       = P.GetStringValueFor(consCmmd);
    string scaffFile      = P.GetStringValueFor(scaffCmmd);
    string outName        = P.GetStringValueFor(outCmmd);
 
    ConsensOverlapUnit COUnit(pairSzFileFile, consFile, "");

    Assembled assembly;
    ContigScaffoldIO io;

    io.Read(assembly, scaffFile);
  
    // Main loop over scaffolds/contigs to update read counts
    for (int l=0; l<assembly.isize(); l++) {
        Scaffold & currScaff = assembly[l];
        for (int i=0; i<currScaff.isize(); i++) {
            Contig & currContig = currScaff[i];
	    int totalContigReads = 0;
	    int totalContigPairs = 0;
            for (int j=0; j<currContig.isize(); j++) {
                const ReadPlacement & rPlace = currContig[j]; 
                int id   = rPlace.Read();
                int pair = rPlace.Pair();
		totalContigReads += COUnit.getConsensCount(id);
                if(pair>-1) { totalContigPairs += COUnit.getConsensCount(id); }
	    }
            currScaff[i].SetNumReads(totalContigReads);
            currScaff[i].SetNumPairs(totalContigPairs/2);
	}
    }
    string layoutFile  = outName+".layout";
    string summaryFile = outName+".summary";
    io.Write(assembly, outName+".layout");
    io.WriteReadCountSummary(assembly, outName+".summary");

    return 0;
}
