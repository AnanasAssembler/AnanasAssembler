#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"
#include "base/StringUtil.h"
#include "base/RandomStuff.h"
#include "src/Ananas/ReadOverlap.h"


class SimulatedReads{
public:
  SimulatedReads():m_reads(), m_offsets() {}

  int getReadCnt() { return m_reads.isize(); }

  void addRead(const DNAVector& read, int offset) {
    m_reads.push_back(read);
    m_offsets[read.getName()] = offset; 
  }
  //Compare reads for sorting based on offset 
  bool operator() (const DNAVector& d1, const DNAVector& d2) const { 
    string origName_1, origName_2;
    int offset_1, length_1, index_1, offset_2, length_2, index_2;
    bool isLeft_1, isLeft_2;
    getInfo(d1, origName_1, offset_1, length_1, index_1, isLeft_1);
    getInfo(d2, origName_2, offset_2, length_2, index_2, isLeft_2);
    if(origName_1==origName_2) {
      return (offset_1<offset_2);
    } else {
      return (origName_1<origName_2);
    }
  }

  void findAllOverlaps(AllReadOverlaps& allOverlaps) {
    svec<DNAVector> tempReads = m_reads;
    sortReads(tempReads); //Make sure reads are sorted
    for (int i=0; i<tempReads.isize(); i++) {
      string origName_i;
      int offset_i, length_i, index_i;
      bool isLeft_i;
      getInfo(tempReads[i], origName_i, offset_i, length_i, index_i, isLeft_i);
      for(int j=i+1; j<tempReads.isize();j++) {
        string origName_j;
        int offset_j, length_j, index_j;
        bool isLeft_j;
        getInfo(tempReads[j], origName_j, offset_j, length_j, index_j, isLeft_j);
        int contactPos  = offset_j - offset_i;
        if(origName_j!=origName_i || contactPos>length_i) { break; }
        bool overlapDir = true;
        bool strand     = isLeft_j;
        allOverlaps.addOverlap(overlapDir, index_i, index_j, contactPos, 1.0, strand);
      }
      for(int k=i-1; k>0; k--) {
        string origName_k;
        int offset_k, length_k, index_k;
        bool isLeft_k;
        getInfo(tempReads[k], origName_k, offset_k, length_k, index_k, isLeft_k);
        int contactPos  = offset_i - offset_k;
        if(origName_k!=origName_i || contactPos>length_i) { break; }
        bool overlapDir = false;
        bool strand     = isLeft_k;
        allOverlaps.addOverlap(overlapDir, index_i, index_k, contactPos, 1.0, strand);
      }
    }
  }

  void writeReads(ostream& sout) {
    vector<char> separators;
    separators.push_back('\\');
    for(int i=0; i<m_reads.isize(); i++) {
      string name = m_reads[i].getName();
      vector<string> tokens;
      Tokenize(name, separators, tokens);
      string printName = tokens[0]; //Removing extra tags to standardize for pair association
      sout << printName << endl
           << m_reads[i].AsString()   << endl;
    }
  }

private:
  void getInfo(const DNAVector& d, string& origName, int& offset, int& length, int& index, bool& isLeft) const {
    string name = d.getName();
    vector<char> separators;
    separators.push_back('\\');
    separators.push_back('/');
    separators.push_back('#');
    vector<string> tokens;
    Tokenize(name, separators, tokens);
    origName = tokens[0];
    isLeft   = ((atoi(tokens[2].c_str())==1)?true:false);
    offset   = atoi(tokens[3].c_str());
    length   = atoi(tokens[4].c_str());
    index    = atoi(tokens[5].c_str());
  }

  void sortReads(svec<DNAVector>& readsToSort) {
    sort(readsToSort.begin(), readsToSort.end(), *this);
  } 
  svec<DNAVector> m_reads;    /// The sequences for the simulated reads
  map<string, int> m_offsets; /// The offsets in original read that the reads were generated with
};

class ReadSimulator {
public:
  ReadSimulator() {}

  /** Generated Paired end reads with given error rates and length */
  void generatePEReads(const vecDNAVector& inputSeqs, int interval, 
                       double sub, double indel, SimulatedReads& out) {

    int i, j;
    int len = 100;
    int sep = 100;
    
    int plusminus = 60;
    
    int totCnt = 0;
    for (j=0; j<inputSeqs.isize(); j++) {
      const DNAVector & d = inputSeqs[j];
      for (i=0; i<d.isize()-2*len-sep; i+= interval) {
	int pos1 = i + RandomInt(plusminus) - plusminus/2;
	if (pos1 < 0)
	  pos1 = 0;
	int pos2 = pos1 + len + sep;
	if (pos2 + len >= d.isize())
	  pos2 = d.isize() - len;
	DNAVector left, right;
	left.SetToSubOf(d, pos1, len);   
	right.SetToSubOf(d, pos2, len);
	DNAVector left_err, right_err;
	mutate(left, sub, indel, left_err);
	mutate(right, sub, indel, right_err);
	right_err.ReverseComplement();
        char strLeft[1024], strRight[1024];
        // Index added after underscore is for identifying same base name for a pair
	sprintf(strLeft, "%s#%d/1\\%d\\%d\\%d", inputSeqs.Name(j).c_str(), totCnt, pos1, len, totCnt);
	sprintf(strRight, "%s#%d/2\\%d\\%d\\%d", inputSeqs.Name(j).c_str(), totCnt, pos2,len, totCnt+1);
	totCnt+=2;
        left_err.setName(strLeft);
        right_err.setName(strRight);
        out.addRead(left_err, pos1);
        out.addRead(right_err, pos2);
      }
    }
  }

private:
  void mutate(const DNAVector & in, double sub, double indel, DNAVector & out) {
    int i, j;
  
    string n;
  
    for (i=0; i<in.isize(); i++) {
      if (RandomFloat(1.) < indel) {
        if (RandomFloat(1.) < 0.5) {
	  // Deletion
	  i += RandomInt(3) + 1;
        } else {
	  // Insertion
	  int plus = RandomInt(3) + 1;
	  for (j=0; j<plus; j++) {
	    n += NucLetter(RandomInt(4));
	  }
        }
      }
    
      if (RandomFloat(1.) < sub) {
        int s = NucIndex(in[i]);
        s += RandomInt(3) + 1;
        s = s % 4;
        n += NucLetter(s);
      } else {
        n += in[i];
      }
    }
    out.SetFromBases(n);
  }
};

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> nOutCmmd("-o","outputName", "SIM");
  commandArg<int> iCmmd("-s","step/interval size", 35);
  commandArg<double> indelCmmd("-indel","indel error rate", 0.);
  commandArg<double> subCmmd("-S","substitution error rate", 0.);
  commandLineParser P(argc,argv);
  P.SetDescription("Simulating paired-end reads for assembly");
  P.registerArg(fileCmmd);
  P.registerArg(nOutCmmd);
  P.registerArg(iCmmd);
  P.registerArg(indelCmmd);
  P.registerArg(subCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string outName  = P.GetStringValueFor(nOutCmmd);
  int interval    = P.GetIntValueFor(iCmmd);
  double indel    = P.GetDoubleValueFor(indelCmmd);
  double sub      = P.GetDoubleValueFor(subCmmd);
  vecDNAVector inputSeqs;
  inputSeqs.Read(fileName);

  ReadSimulator  readSim;
  SimulatedReads simReads; 
  readSim.generatePEReads(inputSeqs, interval, indel, sub, simReads);  

  ofstream fout1;
  fout1.open((outName+"_reads.fa").c_str());
  simReads.writeReads(fout1);
  fout1.close();

  AllReadOverlaps allOverlaps(simReads.getReadCnt());
  simReads.findAllOverlaps(allOverlaps);
  string outputFile = outName+"_readOverlaps.out";
  allOverlaps.write(outputFile, 1);

  return 0;
}
