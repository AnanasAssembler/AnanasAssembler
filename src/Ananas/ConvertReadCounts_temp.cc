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
    commandArg<string> pairCmmd("-p","input read pair/size info file");
    commandArg<string> consCmmd("-c","input read consensus group file");
    commandArg<string> scaffCmmd("-s","scaffolds file");
    commandArg<string> outCmmd("-o","output scaffold file with modified read counts", "modifiedReadCnts.layout");
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
    string outFile        = P.GetStringValueFor(outCmmd);
 
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
            for (int j=0; j<currContig.isize(); j++) {
                const ReadPlacement & rPlace = currContig[j]; 
                int id = rPlace.Read();
  	        int n = COUnit.getConsensCount(id);
		totalContigReads += n;
	    }
            currScaff[i].SetNumReads(totalContigReads);
	}
    }

    io.Write(assembly, outFile);

    return 0;
}
