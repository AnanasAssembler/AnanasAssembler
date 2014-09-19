#include <string>

#include "base/CommandLineParser.h"
#include "extern/logger/log.h"
#include "src/AnnotationQuery/AnnotationQuery.h"
#include "src/CodonTranslate.cc"
#include "src/DNAVector.h"


int main(int argc,char** argv)
{

  commandArg<string> aStringCmmd("-s","Source GTF file");
  commandArg<string> bStringCmmd("-S","Source genome file");
  commandArg<string> cStringCmmd("-O","Output file", "synonymousMutations.out");
  commandArg<string> dStringCmmd("-l","Application logging file","application.log");
  commandLineParser P(argc,argv);
  P.SetDescription("Perform Synonymous mutation analysis");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(cStringCmmd);
  P.registerArg(dStringCmmd);
  P.parse();
  string gtfFile      = P.GetStringValueFor(aStringCmmd);
  string fastaFile    = P.GetStringValueFor(bStringCmmd);
  string outputFile   = P.GetStringValueFor(cStringCmmd);
  string appLogFile   = P.GetStringValueFor(dStringCmmd);
 
  FILE* pFile = fopen(appLogFile.c_str(), "w");
  Output2FILE::Stream()     = pFile;
  FILELog::ReportingLevel() = logINFO; 
  FILE_LOG(logINFO) << "Finding Synonymous mutations";
 
  ofstream fout;
  fout.open(outputFile.c_str());

  CodonTranslator codonTrans;

  Annotation gtf = Annotation(gtfFile,"input");
  svec<AnnotItemBase*> trans =  gtf.getDataByCoord(TRANS);
  vecDNAVector genome;
  genome.Read(fastaFile);

  // For every transcript output all exons and attach exons 
  // into one sequence and output synonymous mutations.
  for(int i=0; i<trans.isize(); i++) {
    svec<AnnotItemBase*> exons = trans[i]->getExons();
    DNAVector seq;
    // Exons are in the order they were read (assumed to be ascending from exon1 onwards)
    for(int j=0; j<exons.isize(); j++) { 
      FILE_LOG(logDEBUG3) << exons[j]->toString('\t');
      fout << "@CDS " << exons[j]->toString('\t') << endl;
      DNAVector tempSeq;
      genome.SetSequence(exons[j]->getCoords(), tempSeq); 
      FILE_LOG(logDEBUG3) << tempSeq.AsString();
      seq+=tempSeq;
    }
    FILE_LOG(logDEBUG2) << seq.AsString();
    string outStr;
    codonTrans.GetSynMutInfo(seq, outStr);
    fout << "@SM " << outStr << endl;
    FILE_LOG(logDEBUG3) << outStr;
  }
  FILE_LOG(logINFO) << "Finished finding Synonymous mutations";
  return 0;
}
  
