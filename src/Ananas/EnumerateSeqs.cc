#include "base/CommandLineParser.h"
//#include "src/XAssembler/SearchOverlaps.h"
#include "src/Ananas/ContScaff.h"
#include "src/Ananas/ContScaffIO.h"
#include "src/Ananas/SequenceGraph.h"



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input contig fasta file");
  //commandArg<string> lapCmmd("-l","input read overlap file");
  commandArg<string> scaffCmmd("-s","scaffolds file");
  commandArg<string> fastaCmmd("-f","output fasta file", "contigs_guided.fasta");
  commandArg<string> layoutCmmd("-o","output layout file", "contigs_guided.layout");
  commandArg<string> dirCmmd("-dir","direction of pairs (fr or ff)");
  commandArg<double> minCmmd("-m","minimum overlap identity", 0.985);
  commandArg<bool> exCmmd("-e","DO NOT DO exhaustive search (report top-n)", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Assembles reads from overlaps.");
  P.registerArg(fileCmmd);
  //P.registerArg(lapCmmd);
  P.registerArg(scaffCmmd);
  P.registerArg(fastaCmmd);
  P.registerArg(layoutCmmd);
  P.registerArg(minCmmd);
  P.registerArg(exCmmd);
  P.registerArg(dirCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  //string lapName = P.GetStringValueFor(lapCmmd);
  string fastaName = P.GetStringValueFor(fastaCmmd);
  string layoutName = P.GetStringValueFor(layoutCmmd);
  string scaffName = P.GetStringValueFor(scaffCmmd);
  double mI = P.GetDoubleValueFor(minCmmd);
  bool bEx2 = P.GetBoolValueFor(exCmmd);

  string dir = P.GetStringValueFor(dirCmmd);
  bool bEx = true;
  if (bEx2)
    bEx = false;

  Assembled assembly;
  ContigScaffoldIO io;

  io.Read(assembly, scaffName);

  vecDNAVector contigs;
  contigs.Read(fileName);

 
  int i, j, k, l;
  for (l=0; l<assembly.isize(); l++) {
    
    const Scaffold & s = assembly[l];
    vecDNAVector local;

    cout << "Processing scaffold " << l << endl;

    for (i=0; i<s.isize(); i++) {
      const Contig & c = s[i];
      DNAVector & d = contigs(c.Name());
      if (c.Ori() == -1) {
	cout << "Reversing contig." << endl;
	d.ReverseComplement();
      }
      local.push_back(d);
    }
    SeqGraph graph;
    GraphEnumerate enumerate;
    GraphConstructor constructor(&local);
    constructor.Construct(graph);

    svec<SeqGraph> linear;
    enumerate.Enumerate(linear, graph);


    SequenceBuild builder;
    vecDNAVector solutions;
    builder.Build(solutions, local, linear);
    for (i=0; i<solutions.isize(); i++) {
      const DNAVector & seq = solutions[i];
      for (j=0; j<seq.isize(); j++) 
	cout << seq[j];
      cout << endl << endl;
    }

    cout << "Solutions (detail): " << linear.isize() << endl;
    for (i=0; i<linear.isize(); i++) {
      cout << "# " << i << endl;
      const SeqGraph & g = linear[i];
      for (j=0; j<g.isize(); j++) {
	for (int x=0; x<g[j].GetNumContrib(); x++) {
	  const SeqCoords & cc = g[j].Contrib(x);
	  cout << cc.Seq() << " " << cc.First() << " " << cc.Last() << endl;
	}
	cout << endl;
      }
      SequenceBuild builder;
      DNAVector seq;
      builder.Build(seq, local, g);
      for (j=0; j<seq.isize(); j++) 
	cout << seq[j];
      cout << endl << endl;
    }

  }
 


  // search.DoSearch(reads, 107);
   
  return 0;
}
