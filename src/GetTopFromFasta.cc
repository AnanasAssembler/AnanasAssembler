#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/FileParser.h"
#include "DNAVector.h"


int main( int argc, char** argv )
{

    commandArg<string> fileCmmd("-i","input file");
    commandArg<string> preCmmd("-p","prefix", "");
    commandArg<int> minCmmd("-m","minimum length", 0);
    commandLineParser P(argc,argv);
    P.SetDescription("Prints out the first sequence of an Ananas component.");
    P.registerArg(fileCmmd);
    P.registerArg(preCmmd);
    P.registerArg(minCmmd);
  
    P.parse();
  
    string fileName = P.GetStringValueFor(fileCmmd);
    string prefix = P.GetStringValueFor(preCmmd);
    int min = P.GetIntValueFor(minCmmd);
    vecDNAVector dna;
    dna.Read(fileName);

    int i;
    for (i=0; i<dna.isize(); i++) {
        string name = ">" + prefix;
        name += dna.NameClean(i);
        const DNAVector & d = dna[i];
        if (d.isize() < min)
            continue;
        if (name[name.size()-1] == '0'
            && name[name.size()-2] == '0'
            && name[name.size()-3] == '0') {
            cout << name << endl;
            for (int j=0; j<d.isize(); j++)
                cout << d[j];
            cout << endl;
        }
    }


    return 0;
}
