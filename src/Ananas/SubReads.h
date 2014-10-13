#ifndef _SUB_READS_H_
#define _SUB_READS_H_

#if defined(OPEN_MP)
  #include <parallel/algorithm>
#else
  #include <algorithm>
#endif


#include <map>
#include <string>
#include <sstream>
#include <stdint.h>
#include "base/SVector.h"
#include "src/DNAVector.h"
#include "base/FileParser.h"
#include "extern/logger/log.h"
#include "src/Cola/Cola.h"
#include "src/Ananas/AssemblyParams.h"
#include "src/Ananas/ReadOverlap.h"
#include "src/Ananas/Reads.h"


//======================================================
/** Sequencing sub-reads */
class SubRead {
public:
  SubRead(): m_index(), m_offset(), m_strand(1) {}
  SubRead(int i, int o, int s): m_index(i), m_offset(o), m_strand(s) {}

  bool operator<(const SubRead& other) const { 
    if(m_index!=other.m_index) {
      return (m_index<other.m_index);
    } else {
      return (m_offset<other.m_offset);
    }
  }

  int getIndex() const    { return m_index;  }
  int getOffset() const   { return m_offset; }
  int getStrand() const   { return m_strand; }

  string toString() const;

private:
  int  m_index;   /// Index of read to which the subread comes from
  int  m_offset;  /// Offset into the read which the subread has been extracted from
  int  m_strand;  /// true: positive strand  false:negative strand   //TODO this should only be a bool but int to adapt to Search code
};



//======================================================
template<class ReadType>
class SubReads {
public:
  // Ctor1:
  SubReads(const ReadType& reads, bool singleStrand, int stepSize, int seedSize,
           float minIdent, float minEndCoverage, int minOverlap, int alignmentBound, 
           int minBases, bool consensMode
           ): m_subReads(), m_reads(reads),
              m_params(singleStrand, stepSize, seedSize, minIdent,
                       minEndCoverage, minOverlap, alignmentBound, minBases) {
    constructSubs(m_reads, stepSize, consensMode); 
  }

  // Ctor2:
  SubReads(const ReadType& reads, const AssemblyParams& params, bool consensMode
          ): m_subReads(), m_reads(reads), m_params(params) {
    constructSubs(m_reads, m_params.getSubreadStep(), consensMode); 
  }

  const ReadType& getReads() { return m_reads; }

  bool  isSingleStrand() const          { return m_params.isSingleStrand();    }
  int   getSubreadStep() const          { return m_params.getSubreadStep();    }  
  int   getSeedSize()  const            { return m_params.getSeedSize();       } 
  float getMinIdentity() const          { return m_params.getMinIdentity();    }
  float getMinEndCover() const          { return m_params.getMinEndCover();    }   
  int   getMinOverlap() const           { return m_params.getMinOverlap();     }
  int   getAlignBand() const            { return m_params.getAlignBand();      }
  int   getMinBasePerScaf() const       { return m_params.getMinBasePerScaf(); } 

  int getSize()                         { return m_subReads.size();            } 
  SubRead getByIndex(unsigned long idx) { return m_subReads[idx];              } 

  //The consensMode argument has been added to override parameters for building consensus reads with no reverse strand
  void constructSubs(const ReadType& reads, int offsetStep, bool consensMode); 
void sortSubs(bool consensMode); 

  string getSeq(const SubRead& sr, int startIdx, int endIdx) const; 
  string getSeq(const SubRead& sr, int len) const; 
  string getSeq(const SubRead& sr) const; 
  string getSeq(unsigned long index) const;  //Get Subread sequence chars via the index in subreads
 
  /** Return a vector of SubRead entry indexes for a given 
     read of those Subreads that share a significant subsequence */
  void findOverlaps(unsigned long readIndex, AllReadOverlaps& allOverlaps) const;  

private:
  struct CmpSubReadCO { //SubRead comparison struct for consensus finding
    CmpSubReadCO(const SubReads& s) : m_subreads(s) {}
    bool operator() (const SubRead& s1, const SubRead& s2) const { 
      if(s1.getOffset()!=s2.getOffset()) {
        return (s1.getOffset()<s2.getOffset()); 
      } else {
        return m_subreads.getSeq(s1) < m_subreads.getSeq(s2); 
      } 
    }
    const SubReads& m_subreads;
  };
 
  struct CmpSubReadOL { //SubRead comparison struct for overlap finding
    CmpSubReadOL(const SubReads& s) : m_subreads(s) {}
    bool operator() (const SubRead& s1, const SubRead& s2) const { 
      return m_subreads.getSeq(s1, m_subreads.getSeedSize()) < m_subreads.getSeq(s2, m_subreads.getSeedSize()); 
    }
    bool operator() (const SubRead& s1, const DNAVector& d2) const { 
      return m_subreads.getSeq(s1, m_subreads.getSeedSize()) < d2.Substring(0, m_subreads.getSeedSize()); 
    }
    const SubReads& m_subreads;
  };

  void getSeq(const SubRead& sr, string& outSeq) const    { outSeq = getSeq(sr); }
  void getDNA(const SubRead& sr, DNAVector& outDNA) const { outDNA.SetFromBases(getSeq(sr)); }
  bool checkInitMatch(const DNAVector& origSeq, const DNAVector& extSeq) const; 
  float checkOverlap(const DNAVector& origSeq, const DNAVector& extSeq, int& isForward) const; 

  svec<SubRead>        m_subReads;   /// vector of subreads 
  const ReadType&      m_reads;      /// Reference to the list of reads from which subreads where constructed
  AssemblyParams       m_params;     /// Object containing the various parameters required for assembly
};

//======================================================
//======================================================

template<class ReadType>
string SubReads<ReadType>::getSeq(const SubRead& sr, int startIdx, int endIdx) const { 
  int readSize = m_reads[sr.getIndex()].size();
  if(endIdx>=readSize) { endIdx = readSize-1; }
  int to = endIdx-startIdx+1;
  if(sr.getStrand()==1) {
    return m_reads[sr.getIndex()].Substring(startIdx, to); 
  } else {
    DNAVector read = m_reads[sr.getIndex()];
    read.ReverseComplement();
    DNAVector readSub;
    readSub.SetToSubOf(read, startIdx, to);
    return readSub.AsString(); 
  }
}

template<class ReadType>
string SubReads<ReadType>::getSeq(const SubRead& sr, int len) const { 
  //TODO warn if given length exceeds the available sequence
  return getSeq(sr, sr.getOffset(), sr.getOffset()+len-1);
}

template<class ReadType>
string SubReads<ReadType>::getSeq(const SubRead& sr) const { 
  return getSeq(sr, sr.getOffset(), m_reads[sr.getIndex()].size()-1);
}

template<class ReadType>
string SubReads<ReadType>::getSeq(unsigned long index) const { 
  SubRead sr = m_subReads[index];
  return getSeq(sr, sr.getOffset(), m_reads[sr.getIndex()].size()-1);
}

template<class ReadType>
void SubReads<ReadType>::constructSubs(const ReadType& reads, int offsetStep, bool consensMode) {
  m_subReads.clear();
  for(int i=0; i<reads.getNumOfReads(); i++) {
    //for(int j=0; j<reads[i].isize()-getMinOverlap(); j+=offsetStep) {
    for(int j=0; j<reads[i].isize(); j+=offsetStep) {
      m_subReads.push_back(SubRead(i, j, true)); //i:index j:offset 
      if(!isSingleStrand() && !consensMode) {
        m_subReads.push_back(SubRead(i, j, false)); //i:index j:offset 
      }
    }
  }
  sortSubs(consensMode);
  FILE_LOG(logDEBUG) <<"Total number of reads: " << m_reads.getNumOfReads();
  cout <<"Total number of reads: " << m_reads.getNumOfReads() << endl;
  FILE_LOG(logDEBUG) <<"Total number of subreads: " << m_subReads.size();
  cout <<"Total number of subreads: " << m_subReads.size() << endl;
} 

template<class ReadType>
void SubReads<ReadType>::sortSubs(bool consensMode) {
  FILE_LOG(logDEBUG) << "Starting to sort SubReads";
  cout << "Starting to sort SubReads" << endl;
  // Could not do this in header by including namespaces as std has been declared in a parent header.
  #if defined(OPEN_MP)
    if(consensMode) {
      __gnu_parallel::stable_sort(m_subReads.begin(), m_subReads.end(), CmpSubReadCO(*this));
    } else {
      __gnu_parallel::stable_sort(m_subReads.begin(), m_subReads.end(), CmpSubReadOL(*this));
    }
  #else
    if(consensMode) {
      std::stable_sort(m_subReads.begin(), m_subReads.end(), CmpSubReadCO(*this));
    } else  {
      std::stable_sort(m_subReads.begin(), m_subReads.end(), CmpSubReadOL(*this));
    }
  #endif
  FILE_LOG(logDEBUG) << "Finished sorting SubReads";
  cout << "Finished sorting subreads" << endl;
}

template<class ReadType>
void SubReads<ReadType>::findOverlaps(unsigned long readIndex, AllReadOverlaps& allOverlaps) const { 
  map<unsigned long, bool> readsUsed_curr;    //Flagset for reads that have been searched for a given extension
  readsUsed_curr[readIndex] = true;           //Add the read index to the used list so that overlaps with itself won't be computed
  DNAVector extSeq, origSeq1, origSeq2;       //extension and read sequence (1 for checkInit step and 2 for the alignment)
  int readSize = m_reads[readIndex].isize();
  for(int i=0; i<=readSize-getMinOverlap(); i++) {
    FILE_LOG(logDEBUG4)  << "Iterating position in read: "<< i;
    origSeq1.SetToSubOf(m_reads[readIndex], i);
    if(origSeq1.isize()>i && (origSeq1[i]=='N' || origSeq1[i]=='n')) { return; } //Rough way of disregarding nonesense characters (TODO look into)
    svec<SubRead>::const_iterator fIt = lower_bound(m_subReads.begin(), m_subReads.end(), origSeq1, CmpSubReadOL(*this));
    for (;fIt!=m_subReads.end(); fIt++) {
      if(min((*fIt).getOffset(), i) > 2*getSubreadStep()) { continue; }               //These should have been found already
      map<unsigned long, bool>::iterator it = readsUsed_curr.find((*fIt).getIndex()); //Check if read has already been used
      if(it==readsUsed_curr.end()) {
        FILE_LOG(logDEBUG4)  << "Investigating subread: "<< (*fIt).getIndex() 
                             << ", " << (*fIt).getOffset() << ", " << ((*fIt).getStrand()==1?"+":"-");
        extSeq.SetFromBases(getSeq(*fIt));
        if(!checkInitMatch(origSeq1, extSeq)) { break; }
        int adjustForAlign = min((*fIt).getOffset(), i); //To find alignments from beginning
        FILE_LOG(logDEBUG4)  << "Adjusting for alignment by: " << adjustForAlign << " bases";
        origSeq2.SetToSubOf(m_reads[readIndex], i-adjustForAlign);
        extSeq.SetFromBases(getSeq(*fIt, (*fIt).getOffset()-adjustForAlign, m_reads[(*fIt).getIndex()].size()-1));
        if(min(extSeq.size(), origSeq2.size())<getMinOverlap()) { continue; } //Skip the rest if extension or read < than min overlap
        int overlapDir;
        float matchScore = checkOverlap(origSeq2, extSeq, overlapDir);
        if(matchScore>=0) {
          int contactPos;
          if(overlapDir==1) { 
            contactPos = i - (*fIt).getOffset(); 
          } else {
            contactPos = (readSize-i) - (m_reads[(*fIt).getIndex()].size()-(*fIt).getOffset()); 
          }
          if(contactPos<0) { 
            FILE_LOG(logDEBUG3)  << "Overlap not found (containment) ";
            //continue;
          }
          allOverlaps.addOverlap(readIndex, (*fIt).getIndex(), contactPos, matchScore, overlapDir, (*fIt).getStrand());
          FILE_LOG(logDEBUG3)  << "Adding overlap: " << readIndex << "\t" << (*fIt).getIndex() << "\t" << contactPos
                               << "\t" << matchScore << "\t" << overlapDir << "\t" << (*fIt).getStrand();
          readsUsed_curr[(*fIt).getIndex()] = true;
        } 
      }
    }
  }
}

template<class ReadType>
bool SubReads<ReadType>::checkInitMatch(const DNAVector& origSeq, const DNAVector& extSeq) const {
  if(origSeq.isize()<getSeedSize() || extSeq.isize()<getSeedSize()) { return false; }  // Pre-check 
  bool isMatch = (origSeq.Substring(0, getSeedSize()) == extSeq.Substring(0,getSeedSize())); 
  FILE_LOG(logDEBUG4) << "Check seed match: " << (isMatch?"Found":"Not Found");  
  return isMatch;
}

template<class ReadType>
float SubReads<ReadType>::checkOverlap(const DNAVector& origSeq, const DNAVector& extSeq, int& matchDirection) const {
  FILE_LOG(logDEBUG3) << "Checking Overlap: ";
  float matchScore        = 0;
  if (getAlignBand() == 0)
    matchScore = origSeq.FindIdent(extSeq);  // Direct match identity without alignment (Rough estimate)
  else 
    matchScore = origSeq.FindIdentHP(extSeq);  // Direct match, but forgives homopolymer indels

  int origSeqAlignedBases = min(origSeq.size(), extSeq.size()); // Not absolutely correct esitmate but ok for significant identity 
  int extSeqAlignedBases  = origSeqAlignedBases;
  
  //TODO NOTE: This is a quick and dirty way to test something...
  if(matchScore<getMinIdentity() && getAlignBand() > 3 /*!=0*/) { // Don't do alignment if direct match score is high enough
    FILE_LOG(logDEBUG3) << "Performing Alignment: ";
    Cola aligner;
    AlignmentCola algn = aligner.createAlignment(origSeq, extSeq, AlignerParams(getAlignBand(), SWGA));
    FILE_LOG(logDEBUG4) << algn.toString(100);
    matchScore          = algn.getIdentityScore();
    origSeqAlignedBases = algn.getTargetBaseAligned();
    extSeqAlignedBases  = algn.getQueryBaseAligned();
  }

  if(matchScore>=getMinIdentity()) {
    if(extSeq.size()>=origSeq.size()  && origSeqAlignedBases >=getMinEndCover()*origSeq.size()) {
      matchDirection = 1; //Right side overlap
      return matchScore;
    } else if(extSeqAlignedBases >=getMinEndCover()*extSeq.size()) {
      matchDirection = -1; //Left-side overlap
      return matchScore;
    } else {
      FILE_LOG(logDEBUG2) << "Overlap doesn't extend to eithr side";
    }
  } else { 
    FILE_LOG(logDEBUG2) << "Significant overlap through alignment was not found.";
  }
  return -1;
}


//======================================================
#endif //_SUB_READS_H_
