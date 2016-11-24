#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "ryggrad/src/base/CommandLineParser.h"
#include "src/analysis/SampleGeneMatrix.h"



int main( int argc, char** argv )
{
    commandArg<string> inputCmmd("-i","input file containing three columns, sampleId, read-count, geneId");
    commandArg<string> outputCmmd("-o","output matrix file name", "sampleGeneReadCount.matrix");
    commandLineParser P(argc,argv);
    P.SetDescription("Create a matrix that maps sample Id to a gene count showing its total number of reads for that gene");
    P.registerArg(inputCmmd);
    P.registerArg(outputCmmd);
  
    P.parse();
  
    string inputFile   = P.GetStringValueFor(inputCmmd);
    string outputFile = P.GetStringValueFor(outputCmmd);

    SampleGeneMatrix sgMatrix;
    cout<< "Loading data" << endl;
    sgMatrix.load(inputFile);
    cout<< "Writing out the matrix" << endl;
    sgMatrix.write(outputFile);

    return 0;
}
 
