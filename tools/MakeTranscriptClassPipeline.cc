#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"



int main( int argc, char** argv )
{

  commandArg<string> inpuCmmd("-g","gtf file containing rna-seq transcripts");
  commandArg<string> fileCmmd("-i","list of gtf files to which the input will be compared");
  commandArg<string> exeCmmd("-e","path to execulables", "");
  commandLineParser P(argc,argv);
  P.SetDescription("Makes a script to run a pipeline.");
  P.registerArg(inpuCmmd);
  P.registerArg(fileCmmd);
  P.registerArg(exeCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string exe = P.GetStringValueFor(exeCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  string last = P.GetStringValueFor(inpuCmmd);//"canis_familiaris.rnaseq.gtf";
  string orig = last;
  int k = 0;

  cout << "# hard coded input files:  " << endl;
  cout << "# isoforms.fpkm_tracking: fpkm values for transcripts" << endl;


  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    string gtf = parser.AsString(0);

    
    cout << exe << "runGTFCompare -qName a -tName b -qAnnot " << last << " -tAnnot ../" << gtf << " > tmp" << endl;
    
    cout << exe << "GTFComparisonStats -i tmp > " << gtf << ".out" << endl;

    // Sense
    cout << "grep FULL_SENSE " << gtf << ".out" << " > tmp1" << endl;
    cout << "grep PARTIAL_SENSE " << gtf << ".out" << " > tmp2" << endl;
    cout << "cat tmp1 tmp2 > " << gtf << ".out.list" << endl;
    cout << exe << "FPKMStatsForSubset -i ../isoforms.fpkm_tracking -l " << gtf << ".out.list > " << gtf << ".out.list.fpkm" << endl;


    cout << "grep TARGET " << gtf << ".out > " << gtf << ".out.target" << endl;
    cout << exe << "CountLoci -i ../" << gtf << " -l " << gtf << ".out.target |sort -u |wc > " << gtf << ".out.target.loci" << endl;


    // Antisense
    cout << "grep FULL_ANTI " << gtf << ".out" << " > tmp1" << endl;
    cout << "grep PARTIAL_ANTI " << gtf << ".out" << " > tmp2" << endl;
    cout << "cat tmp1 tmp2 > " << gtf << ".out.anti" << endl;
    cout << exe << "FPKMStatsForSubset -i ../isoforms.fpkm_tracking -l " << gtf << ".out.anti > " << gtf << ".out.anti.fpkm" << endl;
    cout << exe << "LociFromFPKM  -l " << gtf << ".out.anti" << " -i ../isoforms.fpkm_tracking |sort -u |wc > " << gtf << ".out.anti.loci" << endl;

   // Intronic
    cout << "grep INTRONIC_SENSE " << gtf << ".out > " <<  gtf << ".out.intronic.sense " << endl;
    cout << exe << "FPKMStatsForSubset -i ../isoforms.fpkm_tracking -l " << gtf << ".out.intronic.sense > " << gtf << ".out.intronic.sense.fpkm" << endl;
    //cout << exe << "LociFromFPKM -l " << gtf << ".out.intronic.sense" << " -i ../" << gtf << " |sort -u |wc > " << gtf << ".out.intronic.sense.loci" << endl;
    cout << exe << "LociFromFPKM  -l " << gtf << ".out.intronic.sense" << " -i ../isoforms.fpkm_tracking |sort -u |wc > " << gtf << ".out.intronic.sense.loci" << endl;

    cout << "grep INTRONIC_ANTI " << gtf << ".out > " <<  gtf << ".out.intronic.anti " << endl;
    cout << exe << "FPKMStatsForSubset -i ../isoforms.fpkm_tracking -l " << gtf << ".out.intronic.anti > " << gtf << ".out.intronic.anti.fpkm" << endl;
    //cout << exe << "LociFromFPKM -l " << gtf << ".out.intronic.anti" << " -i ../" << gtf << " |sort -u |wc > " << gtf << ".out.intronic.anti.loci" << endl;
    cout << exe << "LociFromFPKM  -l " << gtf << ".out.intronic.anti" << " -i ../isoforms.fpkm_tracking |sort -u |wc > " << gtf << ".out.intronic.anti.loci" << endl;



    cout << "grep -v TARGET " << gtf << ".out > " << gtf << ".to_remove" << endl; 
    
    char tmp[128];
    sprintf(tmp, ".%d", k);
    string name = orig;
    name += tmp;
    cout << exe << "RemoveTranscripts -l " << gtf << ".to_remove" << " -i " << last << " > " << name << endl;
    last = name;
    k++;

  }
  return 0;
}
