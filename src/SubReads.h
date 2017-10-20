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
#include "ryggrad/src/base/SVector.h"
#include "ryggrad/src/general/DNAVector.h"
#include "ryggrad/src/base/FileParser.h"
#include "ryggrad/src/base/Logger.h"
#include "cola/src/cola/Cola.h"
#include "AssemblyParams.h"
#include "ReadOverlap.h"
#include "Reads.h"


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

  string toString() const {
    stringstream ss;
    ss << "sequence index: " << m_index << " Offset: " << m_offset << " strand: " << m_strand;
    return ss.str();
  }

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

  const ReadType& getReads() const     { return m_reads;              }
  const int getReadSize(int idx) const { return m_reads[idx].isize(); }

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

  string toString() const;

  //The consensMode argument has been added to override parameters for building consensus reads with no reverse strand
  void constructSubs(const ReadType& reads, int offsetStep, bool consensMode); 
  void sortSubs(bool consensMode); 
  bool filterLowComp(const SubRead& sr);  //Determine if a subread is of very low complexity (repeat character)

  string getSeq(const SubRead& sr, int startIdx, int endIdx) const; 
  string getSeq(const SubRead& sr, int len) const; 
  string getSeq(const SubRead& sr) const; 
  string getSeq(unsigned long index) const;  //Get Subread sequence chars via the index in subReads
 
  void getDNA(const SubRead& sr, int startIdx, int endIdx, DNAVector& outDNA)  const; 

  /** Return a vector of SubRead entry indexes for a given 
     read of those Subreads that share a significant subsequence 
     mode 0: only overlaps that extend the read, mode 1: all overlaps 
     limitNumOfOverlaps specifies the number of overlaps to limit the search to (limitNumOfOverlaps=0 means set the limit to read_size*2) 
   */
  int findOverlaps(unsigned long readIndex, AllReadOverlaps& allOverlaps, int mode, int limitNumOfOverlaps) const;  

  // Update the overlap parameters to allow for more overlaps to be detected
  void increaseTolerance(float identTolRatio, float overlapTolRatio);
 
  int compareBases(const SubRead& s1, const SubRead& s2, int limit=100000000) const;
  int compareBases(const SubRead& s1, const DNAVector& d2, int limit=100000000) const;
  int compareSubReads(const SubRead& s1, const SubRead& s2, int limit=100000000) const;
  int compareSubReads(const SubRead& s1, const DNAVector& d2, int limit=100000000) const;

private:
  struct CmpSubReadCO { //SubRead comparison struct for consensus finding
    CmpSubReadCO(const SubReads& s) : m_subReads(s) {}
    bool operator() (const SubRead& s1, const SubRead& s2) const { 
      if(s1.getOffset()!=s2.getOffset()) {
        return (s1.getOffset()<s2.getOffset()); 
      } else {
        return (m_subReads.compareSubReads(s1, s2)==-1);
      } 
    }
    const SubReads& m_subReads;
  };
 
  struct CmpSubReadOL { //SubRead comparison struct for overlap finding
    CmpSubReadOL(const SubReads& s) : m_subReads(s) {}
    bool operator() (const SubRead& s1, const SubRead& s2) const { 
      return (m_subReads.compareSubReads(s1, s2, m_subReads.getSeedSize())==-1);
    }
    bool operator() (const SubRead& s1, const DNAVector& d2) const { 
      return (m_subReads.compareSubReads(s1, d2, m_subReads.getSeedSize())==-1);
    }
    const SubReads& m_subReads;
  };

  void getSeq(const SubRead& sr, string& outSeq) const    { outSeq = getSeq(sr); }
  bool checkInitMatch(const SubRead& origSeqSR, const SubRead& extSeqSR) const; 
  float checkOverlap(const DNAVector& origSeq, const DNAVector& extSeq) const; 
  // Flag any exisitng overlaps from previous iterations as used reads.
  void flagOverlapReads(const AllReadOverlaps& allOverlaps, map<unsigned long, bool>& readsUsed, unsigned long readIndex) const;  

  template<class IterType>
  bool handleIterInstance(IterType iter, map<unsigned long, bool>& readsUsed_curr, const SubRead& origSeqSR,
                                            unsigned long readIndex, AllReadOverlaps& allOverlaps, int mode, int readIterPos) const; 

  svec<SubRead>        m_subReads;   /// vector of subReads 
  const ReadType&      m_reads;      /// Reference to the list of reads from which subReads where constructed
  AssemblyParams       m_params;     /// Object containing the various parameters required for assembly
};

//======================================================
//======================================================

template<class ReadType>
void SubReads<ReadType>::getDNA(const SubRead& sr, int startIdx, int endIdx, DNAVector& outDNA) const {  
  int readSize = m_reads[sr.getIndex()].size();
  if(endIdx>=readSize) { endIdx = readSize-1; }
  int len = endIdx-startIdx+1;
  if(sr.getStrand()==1) {
    outDNA.SetToSubOf(m_reads.getReadByIndex(sr.getIndex()), startIdx, len); 
  } else {
    outDNA.SetToSubOf(m_reads.getReadRCByIndex(sr.getIndex()), startIdx, len);
  }

}

template<class ReadType>
string SubReads<ReadType>::getSeq(const SubRead& sr, int startIdx, int endIdx) const { 
  int readSize = m_reads[sr.getIndex()].size();
  if(endIdx>=readSize) { endIdx = readSize-1; }
  int to = endIdx-startIdx+1;
  if(sr.getStrand()==1) {
    return m_reads.getReadByIndex(sr.getIndex()).Substring(startIdx, to); 
  } else {
    return m_reads.getReadRCByIndex(sr.getIndex()).Substring(startIdx, to); 
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
  double numSubreads = (double) reads.getNumOfReads()*reads.getSize(0)/m_params.getSubreadStep()+1;
  m_subReads.reserve(numSubreads); // If reads have varying sizes, this is just an estimate
  for(int i=0; i<reads.getNumOfReads(); i++) {
    for(int j=0; j<reads[i].isize(); j+=offsetStep) {
      SubRead sr(i, j, 1); //i:index j:offset 
      if(filterLowComp(sr) && !consensMode) { continue; }
      m_subReads.push_back(sr); //i:index j:offset 
      if(!isSingleStrand() && !consensMode) {
        m_subReads.push_back(SubRead(i, j, -1)); //i:index j:offset 
      }
    }
  }
  sortSubs(consensMode);
  FILE_LOG(logDEBUG) <<"Total number of reads: " << m_reads.getNumOfReads();
  cout <<"Total number of reads: " << m_reads.getNumOfReads() << endl;
  FILE_LOG(logDEBUG) <<"Total number of subReads: " << m_subReads.size();
  cout <<"Total number of subReads: " << m_subReads.size() << endl;
  FILE_LOG(logDEBUG4) << toString() << endl;
} 

template<class ReadType>
bool SubReads<ReadType>::filterLowComp(const SubRead& sr) {
  int  idx    = sr.getIndex();
  int  offset = sr.getOffset();
  int  end    = offset + m_params.getSeedSize();
  char prevC  = ' '; // Previous rep character
  char currC  = ' '; // Current character
  int  repCnt = 0;   // keep track of size of longest repeat contig
  int  nCnt   = 0;   // Keep track of nonesense characters
  for(int i=offset; i<end; i++) {
    char currC =  m_reads[idx][i];
    if(currC=='n' || currC=='N') { nCnt++; }
    if(prevC==currC) {
      repCnt++;
    } else {
      repCnt = 0;
      prevC  = currC;
    }
  }
  if((float)repCnt/m_params.getSeedSize()>0.90 || nCnt>0.4*m_params.getSeedSize()) { 
    return true;
  } else {
    return false;
  }
}

template<class ReadType>
void SubReads<ReadType>::sortSubs(bool consensMode) {
  FILE_LOG(logDEBUG) << "Starting to sort SubReads";
  cout << "Starting to sort SubReads" << endl;
  // Could not do this in header by including namespaces as std has been declared in a parent header.
  #if defined(OPEN_MP)
    if(consensMode) {
      __gnu_parallel::sort(m_subReads.begin(), m_subReads.end(), CmpSubReadCO(*this));
    } else {
      __gnu_parallel::sort(m_subReads.begin(), m_subReads.end(), CmpSubReadOL(*this));
    }
  #else
    if(consensMode) {
      std::sort(m_subReads.begin(), m_subReads.end(), CmpSubReadCO(*this));
    } else  {
      std::sort(m_subReads.begin(), m_subReads.end(), CmpSubReadOL(*this));
    }
  #endif
  FILE_LOG(logDEBUG) << "Finished sorting SubReads";
  cout << "Finished sorting subReads" << endl;
}

template<class ReadType>
int SubReads<ReadType>::findOverlaps(unsigned long readIndex, AllReadOverlaps& allOverlaps, int mode, int limitNumOfOverlaps) const { 
  map<unsigned long, bool> readsUsed_curr;                  // Flagset for reads that have been searched for a given extension
  readsUsed_curr[readIndex] = true;                         // Add the read index to the used list so that overlaps with itself won't be computed
  int readSize = m_reads[readIndex].isize();
  if(limitNumOfOverlaps <= 0) { limitNumOfOverlaps = 2 * readSize; }

  for(int i=0; i<=readSize-getMinOverlap(); i++) {
    FILE_LOG(logDEBUG4)  << "Iterating position in string: "<< i;
    SubRead origSeqSR(readIndex, i, 1); 
    svec<SubRead>::const_iterator fIt = lower_bound(m_subReads.begin(), m_subReads.end(), origSeqSR, CmpSubReadOL(*this));
    FILE_LOG(logDEBUG4)  << "Searching for suffix - found lower-bound: " << (*fIt).toString() << endl;
    svec<SubRead>::const_reverse_iterator rIt(fIt);
    bool keepLooking = true;
    for (;keepLooking && fIt!=m_subReads.end(); fIt++) {
      keepLooking = handleIterInstance<svec<SubRead>::const_iterator>(fIt, readsUsed_curr, origSeqSR, readIndex, allOverlaps, mode, i);
    }
    keepLooking = true;
    for (;keepLooking && rIt!=m_subReads.rend(); rIt++) {
      keepLooking = handleIterInstance<svec<SubRead>::const_reverse_iterator>(rIt, readsUsed_curr, origSeqSR, readIndex, allOverlaps, mode, i);
    }
 
    if(mode==0 && allOverlaps[readIndex].getNumLaps()>=limitNumOfOverlaps) { 
      return allOverlaps[readIndex].getNumLaps(); 
    }  //Limiting overlaps for very high coverage reads 
  }

  return allOverlaps[readIndex].getNumLaps();
}

template<class ReadType>
template<class IterType>
bool SubReads<ReadType>::handleIterInstance(IterType iter, map<unsigned long, bool>& readsUsed_curr, const SubRead& origSeqSR,
                                            unsigned long readIndex, AllReadOverlaps& allOverlaps, int mode, int readIterPos) const {
  if(readsUsed_curr.find((*iter).getIndex())!=readsUsed_curr.end()) { return true; } //Check if read has already been used
  FILE_LOG(logDEBUG4)  << "Investigating subread: "<< (*iter).getIndex() 
                       << ", " << (*iter).getOffset() << ", " << ((*iter).getStrand()==1?"+":"-");
  if(!checkInitMatch(*iter, origSeqSR)) { return false; }
  int adjustForAlign = min((*iter).getOffset(), readIterPos); //To find alignments from beginning
  FILE_LOG(logDEBUG4)  << "Adjusting for alignment by: " << adjustForAlign << " bases";
  DNAVector extSeq, origSeq2;                               // Extension and read sequence (1 for checkInit step and 2 for the alignment)
  origSeq2.SetToSubOf(m_reads[readIndex], readIterPos-adjustForAlign);
  getDNA(*iter, (*iter).getOffset()-adjustForAlign, m_reads[(*iter).getIndex()].size()-1, extSeq);
  if(min(extSeq.size(), origSeq2.size())<getMinOverlap()) { return true; } //Skip the rest if extension or read < than min overlap
  float matchScore = checkOverlap(origSeq2, extSeq);
  if(matchScore>=0) {
    int contactPos = readIterPos - (*iter).getOffset(); 
    int overlapDir = 1;
    if(contactPos<0) { 
      overlapDir = -1; 
      contactPos = (m_reads[readIndex].isize() - contactPos) - m_reads[(*iter).getIndex()].isize(); 
      if(contactPos<0) { 
        FILE_LOG(logDEBUG1)  << "Contained Read: " << m_reads[(*iter).getIndex()].getName();
        return false; 
      } //discard non-extending containment 
    }
    FILE_LOG(logDEBUG3)  << "Adding overlap: " << readIndex << "\t" << (*iter).getIndex() << "\t" << contactPos
                         << "\t" << matchScore << "\t" << 1 << "\t" << (*iter).getStrand();
   // Synchronized version of overlap adding that locks so no iterference occurs with other threads
    allOverlaps.addOverlapSync(readIndex, (*iter).getIndex(), contactPos, overlapDir, (*iter).getStrand());
    if(contactPos==0) { allOverlaps.addOverlapSync(readIndex, (*iter).getIndex(), contactPos, -overlapDir, (*iter).getStrand());} 
    readsUsed_curr[(*iter).getIndex()] = true;
  }
  return true;
}

template<class ReadType>
bool SubReads<ReadType>::checkInitMatch(const SubRead& origSeqSR, const SubRead& extSeqSR) const {
  int origSize = m_reads.getSize(origSeqSR.getIndex());
  int extSize  = m_reads.getSize(extSeqSR.getIndex());
  if(origSize<getSeedSize() || extSize<getSeedSize()) { return false; }  // Pre-check 
  bool isMatch = (compareBases(origSeqSR, extSeqSR, getSeedSize())==0);
  FILE_LOG(logDEBUG4) << "Check seed match: " << (isMatch?"Found":"Not Found");  
  return isMatch;
}

template<class ReadType>
float SubReads<ReadType>::checkOverlap(const DNAVector& origSeq, const DNAVector& extSeq) const {
  FILE_LOG(logDEBUG3) << "Checking Overlap: ";
  float matchScore        = 0;
  if (getAlignBand() == 0)
    matchScore = origSeq.FindIdent(extSeq);    // Direct match identity without alignment (Rough estimate)
  else 
    matchScore = origSeq.FindIdentHP(extSeq, 10);  // Direct match, but forgives homopolymer indels

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

  FILE_LOG(logDEBUG4) << "Score: " << matchScore << endl;

  if(matchScore>=getMinIdentity()) {
    return matchScore;
  } else { 
    FILE_LOG(logDEBUG2) << "Significant overlap through alignment was not found." 
                        << " Match Identity: " << matchScore << "  Acceptable threshold: " << getMinIdentity();
  }
  return -1;
}

template<class ReadType>
void SubReads<ReadType>::flagOverlapReads(const AllReadOverlaps& allOverlaps, 
                                         map<unsigned long, bool>& readsUsed, 
                                         unsigned long readIndex) const { 
  const ReadInfo& rInfo = allOverlaps[readIndex];
  int totLaps     =  rInfo.getNumLaps();
  for(int i=0; i<totLaps; i++) {
    readsUsed[rInfo.getLap(i).getOverlapIndex()] = true;
  }
}

template<class ReadType>
void SubReads<ReadType>::increaseTolerance(float identTolRatio, float overlapTolRatio) {
  m_params.setMinIdentity(m_params.getMinIdentity()*(1-identTolRatio));
  m_params.setMinOverlap(m_params.getMinOverlap()*(1-overlapTolRatio));
}

template<class ReadType>
int SubReads<ReadType>::compareBases(const SubRead& s1, const SubRead& s2, int limit) const { 
  int idx1    = s1.getIndex();
  int offset1 = s1.getOffset();
  int strand1 = s1.getStrand();
  int idx2    = s2.getIndex();
  int offset2 = s2.getOffset();
  int strand2 = s2.getStrand();

  int size1 = m_reads.getSize(idx1)-offset1;
  int size2 = m_reads.getSize(idx2)-offset2;
  limit = min(limit, min(size1, size2));

  const DNAVector& d1 = m_reads.getReadByIndex(idx1, strand1);
  const DNAVector& d2 = m_reads.getReadByIndex(idx2, strand2);

  for(int i=0; i<limit; i++) {
    if(d1[offset1+i]<d2[offset2+i]) {
      return -1; //Smaller
    } else if(d1[offset1+i]>d2[offset2+i]) {
      return 1;  //Larger
    }
  }
  return 0; // Same
}

template<class ReadType>
int SubReads<ReadType>::compareBases(const SubRead& s1, const DNAVector& d2, int limit) const { 
  int idx1    = s1.getIndex();
  int offset1 = s1.getOffset();
  int strand1 = s1.getStrand();

  int size1 = m_reads.getSize(idx1)-offset1;
  int size2 = d2.size();
  limit = min(limit, min(size1, size2));
  const DNAVector& d1 = m_reads.getReadByIndex(idx1, strand1);
  for(int i=0; i<limit; i++) {
    if(d1[offset1+i]<d2[i]) {
      return -1; //Smaller
    } else if(d1[offset1+i]>d2[i]) {
      return 1;  //Larger
    }
  }
  return 0;   // Same
}

template<class ReadType>
int SubReads<ReadType>::compareSubReads(const SubRead& s1, const SubRead& s2, int limit) const { 
  int baseCmp = SubReads<ReadType>::compareBases(s1, s2, limit);
  if(baseCmp != 0) { return baseCmp; }

  int size1 = m_reads.getSize(s1.getIndex()) - s1.getOffset();
  int size2 = m_reads.getSize(s2.getIndex()) - s2.getOffset();
  if(size1 > size2)      { return 1; }
  else if(size1 < size2) { return -1;  }
  else                   { return 0;  } 
}

template<class ReadType>
int SubReads<ReadType>::compareSubReads(const SubRead& s1, const DNAVector& d2, int limit) const { 
  int baseCmp =  SubReads<ReadType>::compareBases(s1, d2, limit);
  if(baseCmp != 0) { return baseCmp; }

  int size1 = m_reads.getSize(s1.getIndex()) - s1.getOffset();
  int size2 = d2.size();
  if(size1 > size2)      { return 1; }
  else if(size1 < size2) { return -1;  }
  else                   { return 0;  } 
}
  
template<class ReadType>
string SubReads<ReadType>::toString() const { 
    stringstream ss;
    for(int i=0; i<m_subReads.isize(); i++) {
        ss << m_subReads[i].toString() << " - " << getSeq(i) << endl;  // For when the sequence should also be printed 
        //ss << m_subReads[i].toString() << endl; 
    }
    return ss.str();
}

 //======================================================
#endif //_SUB_READS_H_
