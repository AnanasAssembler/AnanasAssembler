#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "ryggrad/src/base/CommandLineParser.h"
#include "ConsensOverlapUnit.h"



int main( int argc, char** argv )
{

    commandArg<string> lapInCmmd("-i","input binary overlap file");
    commandArg<string> lapOut1Cmmd("-o1","output ascii overlap file", "temp_overlaps.ascii");
    commandArg<string> lapOut2Cmmd("-o2","output ascii overlap file", "temp_overlaps.stats");
    commandLineParser P(argc,argv);
    P.SetDescription("Convert overlap binary file to ascii");
    P.registerArg(lapInCmmd);
    P.registerArg(lapOut1Cmmd);
    P.registerArg(lapOut2Cmmd);
  
    P.parse();
  
    string lapInName   = P.GetStringValueFor(lapInCmmd);
    string lapOut1Name = P.GetStringValueFor(lapOut1Cmmd);
    string lapOut2Name = P.GetStringValueFor(lapOut2Cmmd);

    ConsensOverlapUnit COUnit(AssemblyParams(), "");
    cout<< "Loading Overlaps" << endl;
    COUnit.ReadOverlaps(lapInName);
    cout<< "Writing out Overlaps stats" << endl;
    COUnit.writeOverlaps(lapOut2Name, 2);
    cout<< "Writing out Overlaps in ascii" << endl;
    COUnit.writeOverlaps(lapOut1Name, 1);

    return 0;
}
 
