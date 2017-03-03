#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "ryggrad/src/base/CommandLineParser.h"
#include "SearchOverlaps.h"
#include "ContScaff.h"
#include "ContScaffIO.h"
#include "ConsensOverlapUnit.h"


class ToDo
{
public:
  ToDo() {
    m_reads = 0.;
  }

  void Add(int scaff, int reads) {
    m_data.push_back(scaff);
    double r = (double)reads;
    m_reads += r*r;
  }
  double Score() const {return m_reads;}

  int isize() const {return m_data.isize();}
  int operator[] (int i) const {return m_data[i];}

private:
  double m_reads;
  svec<int> m_data;
};


class Scaff
{
public:
  Scaff() {
    m_s = -1;
    m_r = 0;
  }
  Scaff(int s, int r) {
    m_s = s;
    m_r = r;
  }
  int ID() const {return m_s;}
  int Reads() const {return m_r;}


  bool operator < (const Scaff &s) const {
    if (m_r == s.m_r)
      return (m_s < s.m_s);
    return (m_r < s.m_r);      
  }

private:
  int m_r;
  int m_s;
  
};

class ScaffoldDispenser
{
public:
  ScaffoldDispenser(int n) {
    m_todo.resize(n);
  }

  const ToDo & Get(int i) const {return m_todo[i];}
  void Add(int scaff, int reads) {
      //cout << "Adding scaffold " << scaff << endl;
      m_tmp.push_back(Scaff(scaff, reads));
  }

  void Compute() {
      Sort(m_tmp);
      int i, j;
      //cout << "Size: " << m_tmp.isize() << endl;
      for (i=m_tmp.isize()-1; i>=0; i--) {
          int s = m_tmp[i].ID();
          int r = m_tmp[i].Reads();
          int index = 0;
	  double min = m_todo[index].Score();
	  //cout << "Min=" << min << endl;
	  for (j=0; j<m_todo.isize(); j++) {
	      //cout << "  -> curr " << j << " -> " << m_todo[j].Score() << endl;
	      if (m_todo[j].Score() < min) {
		min = m_todo[j].Score();
		index = j;
	      }
	  }
	  //cout << "Assigning scaffold " << s << " to index " << index << " reads: " << r << endl; 
	  m_todo[index].Add(s, r);
      }
  }

private:
  svec<ToDo> m_todo;
  svec<Scaff> m_tmp;
};


int main( int argc, char** argv )
{
    commandArg<string> fileCmmd("-i","input read pair/size info file");
    commandArg<string> lapCmmd("-l","input read overlap file");
    commandArg<string> consCmmd("-g","input read consensus group file");
    commandArg<string> scaffCmmd("-s","scaffolds file");
    commandArg<string> layoutCmmd("-o","output layout file", "contigs_guided.layout");
    commandArg<string> dirCmmd("-dir","direction of pairs (fr or ff)");
    commandArg<int> libSizeCmmd("-libSize","Maximum library size", 500);
    commandArg<double> minCmmd("-m","minimum overlap identity", 0.985);
    commandArg<int> minContigCmmd("-minContig","minimum length of a single-contig scaffold to report", 200);
    commandArg<bool> exCmmd("-e","DO NOT DO exhaustive search (report top-n)", false);
    commandArg<int> numCmmd("-num","Process # (for parallel runs)", 0);
    commandArg<int> ofCmmd("-of","Out of # processes (for parallel runs)", 1);
    commandArg<string> prefixCmmd("-prefix","The prefix to add to all generated contig names", "Sample1");
    commandLineParser P(argc,argv);
    P.SetDescription("Assembles COUnit from overlaps.");
    P.registerArg(fileCmmd);
    P.registerArg(lapCmmd);
    P.registerArg(consCmmd);
    P.registerArg(scaffCmmd);
    P.registerArg(layoutCmmd);
    P.registerArg(libSizeCmmd);
    P.registerArg(minCmmd);
    P.registerArg(exCmmd);
    P.registerArg(dirCmmd);
    P.registerArg(minContigCmmd);
    P.registerArg(numCmmd);
    P.registerArg(ofCmmd);
    P.registerArg(prefixCmmd);
 
    P.parse();
  
    string pairSzFileName = P.GetStringValueFor(fileCmmd);
    string lapName = P.GetStringValueFor(lapCmmd);
    string layoutName = P.GetStringValueFor(layoutCmmd);
    string scaffName = P.GetStringValueFor(scaffCmmd);
    double mI = P.GetDoubleValueFor(minCmmd);
    bool bEx2 = P.GetBoolValueFor(exCmmd);
    int libSize = P.GetIntValueFor(libSizeCmmd);
    int minContig = P.GetIntValueFor(minContigCmmd);
    int num = P.GetIntValueFor(numCmmd);
    int of = P.GetIntValueFor(ofCmmd);
    string consName = P.GetStringValueFor(consCmmd);
    string prefix = P.GetStringValueFor(prefixCmmd);
 
    string dir = P.GetStringValueFor(dirCmmd);
    bool bEx = true;
    if (bEx2)
        bEx = false;


    ConsensOverlapUnit COUnit(pairSzFileName, consName, "");

    int i, j, k, l;
 
    Assembled assembly;
    ContigScaffoldIO io;

    io.Read(assembly, scaffName);
  
    svec<int> ids;
    ids.resize(COUnit.GetNumReads(), 0);

    ScaffoldDispenser disp(of);

    // Count reads for the dispenser.
    for (l=0; l<assembly.isize(); l++) {
     
        const Scaffold & s = assembly[l];
        if (s.isize() == 1) {
            if (s[0].Highest() < minContig) {
	        //cout << "Skipping scaffold " << l << endl;
                continue;
            }
        }
	int reads = 0;
        for (i=0; i<s.isize(); i++) {
            const Contig & c = s[i];
            reads += c.isize();
	}
	//cout << "Scaffold " << l << " reads: " << reads << " contigs: " << s.isize() << endl;
	disp.Add(l, reads);
    }
    disp.Compute();
    const ToDo & t = disp.Get(num);

    cout << "Process # " << num << " will process " << t.isize() << " scaffolds w/ " << t.Score() << " reads^2." << endl;
    cout << "Loading and pruning overlaps." << endl;
    //=========================================================
    // Collect all COUnit (WARNING: Duplicated code!!)
    for (l=0; l<t.isize(); l++) {
        const Scaffold & s = assembly[t[l]];
        for (i=0; i<s.isize(); i++) {
            const Contig & c = s[i];
            for (j=0; j<c.isize(); j++) {
                const ReadPlacement & r = c[j]; 
                int id = r.Read();
                ids[id] = 1;
            }
        }
    }
    //=========================================================
    // Free up some unused space
    // COUnit.Prune(ids);

    COUnit.ReadOverlaps(lapName, ids);

    Search search;
    search.SetExhaustive(bEx);
    search.SetDir(dir);

    if (of > 1) {
        char tmp[256];
        sprintf(tmp, ".%d", num);
        layoutName += tmp;
    }

    search.SetOutput(layoutName);
    search.SetIndex(num);
    search.SetPrefix(prefix);
    search.SetLibSize(libSize);
    search.SetMinAltKeep(minContig);


    // Main loop over scaffolds
    for (l=0; l<t.isize(); l++) {
        const Scaffold & s = assembly[t[l]];

        cout << "\rProcessing scaffold " << t[l] << " Number of reads= " << s.NumReads() << flush;
        search.SetUsedAll(COUnit);

	
	svec<int> mult;
	int total = 0;
	int allReads = 0;
        for (i=0; i<s.isize(); i++) {
            const Contig & c = s[i];
            for (j=0; j<c.isize(); j++) {
                const ReadPlacement & r = c[j]; 
                int id = r.Read();
  	        int n = COUnit.getConsensCount(id);
                mult.push_back(n);
		total += n;
		allReads++;
	    }
	}

	Sort(mult);
	cout << endl;
	cout << "Total reads: " << allReads << endl; 
	cout << "Min: " << mult[0] << endl;
	cout << "Max: " << mult[mult.isize()-1] << endl;
	cout << "Med: " << mult[mult.isize()/2] << endl;
	int count = 0;

	int minInCons = 1;

	for (i=0; i<mult.isize(); i++) {
	  count += mult[i];
	  if (count > total/20) {
	    cout << "N5: " << mult[i] << endl;
	    minInCons = mult[i];
	    break;
	  }
	}

	//--------------------------------
	// For statistics only
	for (; i<mult.isize(); i++) {
	  count += mult[i];
	  if (count > total/10) {
	    cout << "N10: " << mult[i] << endl;
	    break;
	  }
	}
	for (; i<mult.isize(); i++) {
	  count += mult[i];
	  if (count > total/4) {
	    cout << "N25: " << mult[i] << endl;
	    break;
	  }
	}
	//--------------------------------
	
	int usedReads = 0;

	//=======================================
	if (minInCons > 1)
	  cout << "Reducing data set." << endl;
	//=======================================
        for (i=0; i<s.isize(); i++) {
            const Contig & c = s[i];
            for (j=0; j<c.isize(); j++) {
                const ReadPlacement & r = c[j]; 
                int id = r.Read();
  	        int n = COUnit.getConsensCount(id);
		if (n >= minInCons) {
		  search.SetUsed(id, false);
		  usedReads++;
		}
            }
        }
	cout << "Used reads: " << usedReads << endl;
        // Setting the Layout guided stage to none-exhaustive search of overlaps 
        // Note that this used to be an exhaustive search and we might choose to make it optional in the future
        search.SetExhaustive(false);
        //if (s.isize() == 1)
        //    search.SetExhaustive(false);

        search.DoSearchAll(COUnit, s.NumReads());
    }
 
    cout << "All done!!!" << endl;
    return 0;
}
