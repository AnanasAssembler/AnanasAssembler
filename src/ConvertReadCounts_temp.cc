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
    
    // To keep track of unique reads in scaffolds
    map<int, bool> currReads;
    map<int, bool> currPairs;
  
    // Main loop over scaffolds/contigs to update read counts
    char contigName[512];
    for (int scaffCnt=0; scaffCnt<assembly.isize(); scaffCnt++) {
        Scaffold & currScaff = assembly[scaffCnt];
        for (int contCnt=0; contCnt<currScaff.isize(); contCnt++) {
            Contig & currContig = currScaff[contCnt];
            // This is a temporary hack to fix the problem with discrepancies between layout and fasta ids
            sprintf(contigName, ">Contig_Sample1_000_%7d_%3d", scaffCnt, contCnt);
            for (int i=0; i<(int)strlen(contigName); i++) {
                if (contigName[i] == ' ')
                contigName[i] = '0';
            }
	    int totalContigReads = 0;
	    int totalContigPairs = 0;
            for (int j=0; j<currContig.isize(); j++) {
                const ReadPlacement & rPlace = currContig[j]; 
                int id   = rPlace.Read();
                int pair = rPlace.Pair();
		totalContigReads += COUnit.getConsensCount(id);
                const svec<int>& consMemIds = COUnit.getConsMembers(id);
                for(int cmId:consMemIds) { //Add all the raw reads
                    currReads[cmId] = true;
                }
                if(pair>-1) { 
                    totalContigPairs += COUnit.getConsensCount(id); 
                    for(int cmId:consMemIds) { //Add all the raw read pairs
                        currPairs[cmId] = true;
                    }
                }
	    }
            currScaff[contCnt].SetNumReads(totalContigReads);
            currScaff[contCnt].SetNumPairs(totalContigPairs/2);
            currScaff[contCnt].SetName(contigName);
	}
        currScaff.SetNumUniqReads(currReads.size());
        currScaff.SetNumUniqPairs(currPairs.size()/2);
        currReads.clear();
        currPairs.clear();
    }
    io.Write(assembly, outName+".layout");
    io.WriteReadCountSummary(assembly, outName+".summary");

    return 0;
}
