#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "ryggrad/src/base/CommandLineParser.h"
#include "ReadOverlap.h"

int main( int argc, char** argv )
{

    commandArg<string> lapCmmd("-i","input read overlap file");
    commandArg<string> outCmmd("-o","output overlap file name" , "allOverlaps.ascii");
    commandLineParser P(argc,argv);
    P.SetDescription("Convert Overlaps from binary to ascii");
    P.registerArg(lapCmmd);
    P.registerArg(outCmmd);
  
    P.parse();
  
    string lapFileName  = P.GetStringValueFor(lapCmmd);
    string outFileName  = P.GetStringValueFor(outCmmd);

    AllReadOverlaps overlaps;
    cout << "Loading overlaps..." << endl;
    overlaps.loadBin(lapFileName);
    cout << "Writing overlaps to " << outFileName << "..." << endl;
    overlaps.writeAsc(outFileName);
 
    return 0;
}
