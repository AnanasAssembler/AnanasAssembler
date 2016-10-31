#include "ryggrad/src/base/CommandLineParser.h"
#include "ContScaff.h"
#include "ContScaffIO.h"
#include "ConsensOverlapUnit.h"



int main( int argc, char** argv )
{

    //commandArg<string> fileCmmd("-i","input scaffold/contig/layout file");
    commandArg<string> lapCmmd("-l","input read overlap file");
    commandArg<string> fastaCmmd("-f","read fasta file"); // Do we really need this???
    commandArg<string> consCmmd("-g","input read consensus group file");

  
    commandLineParser P(argc,argv);
    P.SetDescription("Joins contigs into scaffolds.");
    //P.registerArg(fileCmmd);
    P.registerArg(lapCmmd);
    P.registerArg(consCmmd);
    P.registerArg(fastaCmmd);
 
    P.parse();
  
    //string fileName = P.GetStringValueFor(fileCmmd);
    string lapName = P.GetStringValueFor(lapCmmd);
    string fastaName = P.GetStringValueFor(fastaCmmd);
    string consName = P.GetStringValueFor(consCmmd);

    ConsensOverlapUnit COUnit(fastaName, consName, lapName);
    const ConsensReads & rr = COUnit.getConsReads();

    int i, j;

    for (i=0; i<rr.getNumOfReads(); i++) {
        int numPartner = COUnit.getNumOfPartners(i);
        cout << "Read " << i << " " << rr[i].Name() << " partners: " << numPartner << endl;    
        for (j=0; j<numPartner; j++) {
            cout << "   -> " << COUnit.getPartner(i, j) << " " << rr[COUnit.getPartner(i, j)].Name() << endl;     
        }
    }

    return 0;
}

