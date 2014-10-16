#ifndef _READS_H_
#define _READS_H_

#include <map>
#include <string>
#include <stdint.h>
#include "src/DNAVector.h"
#include "src/Ananas/AssemblyParams.h"

//======================================================
//Consensus Read
class ConsensRead 
{
public:
  ConsensRead():m_readIdxs(), m_dnaSeq(), m_seqSize(-1) {}
  
  const svec<int >& getReads() const           { return m_readIdxs;             }
  int getNumOfReads() const                    { return m_readIdxs.isize();     } 
  bool isSingle() const                        { return (m_readIdxs.size()==1); }
  const DNAVector& getSeq() const              { return m_dnaSeq;               }
  int  getSize() const                         { return m_seqSize;              }
 
  void setSize(int size)                       { m_seqSize = size;              }
  void addRead(int readIdx)                    { m_readIdxs.push_back(readIdx); }
  void setReads(svec<int>& readIdxs)           { m_readIdxs = readIdxs;         }

  void setSeq(const DNAVector& seq) { 
    m_dnaSeq  = seq;               
    m_seqSize = m_dnaSeq.size();
  }

  void setName() { 
    if(m_dnaSeq.AsString()!="") { 
      m_dnaSeq.setName(genSeqName()); 
    } 
  }

  void clear() { //TODO comment (function created for adapting to Overlaps class pruning
    m_readIdxs.clear();
    m_dnaSeq.clear();
  }

  string toString() const;
  void toString(string& outString) const;

private:
  string genSeqName(); 

  svec<int> m_readIdxs;  /// List of indexes of reads that constitute this consensus
  DNAVector m_dnaSeq;    /// Contains the consensus only if there's more than one read 
  int       m_seqSize;   /// Used only when dnaSeq doesn't exist (for stages where memory efficiency is required)
};
//======================================================



//======================================================
//DEPRECATE
class SubvecDNAVector: public vecDNAVector {
public:
  //Ctor
  const DNAVector& operator[] (int i) const { return vecDNAVector::operator[](getSubIndex(i)); }
  DNAVector& operator[] (int i)             { return vecDNAVector::operator[](getSubIndex(i)); }

  int getSubIndex(int origIdx) const {
    if(m_indexLinks.size()==0) { return origIdx; }
    else                       { return  m_indexLinks.find(origIdx)->second; }
  }

  void Read(const string & fileName, const svec<int>& idsToKeep); 

private:
  void linkIndex(int origIdx, int subsetIdx) { m_indexLinks[origIdx] = subsetIdx; }

  map<int, int> m_indexLinks; /// Mapping of the original sequence index to the subset index
};
//======================================================



//======================================================
class RawReads {
public:
  // Default Ctor:
  RawReads(): m_reads(), m_pairInfo(), m_sizeInfo() {}

  // Ctor 2:
  /* mode:0 fasta input  mode:1 pair/size input */
  RawReads(const string& fileName, int mode): m_reads(), m_pairInfo(), m_sizeInfo() { 
    if(mode==0) { load(fileName);           }
    else        { loadPairSzInfo(fileName); }
  }

  const DNAVector& operator[](int i) const              { return m_reads[i];            }
  const vecDNAVector& getReads() const                  { return m_reads;               }
  const DNAVector& getReadByIndex(int idx)              { return m_reads[idx];          } 
  bool hasSeq(int idx) const                            { return (m_reads.isize()>idx); }
  int getNumOfReads() const                             { return (m_reads.size())!=0? m_reads.size(): m_sizeInfo.size(); } 
  int getSize(int idx) const                            { return ((hasSeq(idx))? m_reads[idx].size(): m_sizeInfo[idx]);  }

  int getPairId(int idx) const { return((m_pairInfo.isize()>idx)? m_pairInfo[idx]: -1); }

  void write(const string& outFile) const; 
  void write(ostream& sout) const;
  void writeBin(const string& outFile) const; 
  void writePairSzInfo(const string& pairInfoFile) const; 
  void load(const string& inFile); 
  void loadBin(const string& inFile); 
  void loadPairSzInfo(const string& pairInfoFile); 

  void clearRead(int id) { m_reads[id].clear(); }   //TODO comment (function created for adapting to Overlaps class pruning

private:
  vecDNAVector m_reads;        /// Vector containing all reads 
  svec<int>    m_pairInfo;     /// For every raw read index contains the index of a read that's paired with it
  svec<int>    m_sizeInfo;     /// The size of each read, this is used only when the read seqs are not aquired
};

//======================================================



//======================================================
class ReadGroups
{
//This class is to be used together with ConsensusReads as it is not generic at all
public:
  /** Can only be constructed with size as parameter so that user is aware that size should be set */
  ReadGroups(int size):m_groups(), m_groupInfo(), m_tags() { 
    resize(size); 
  }

  void resize(int sz)                              { m_groupInfo.resize(sz, -1);  }
  const svec<int>& getGroup(int i) const           { return m_groups[i];          }
  svec<int>& getGroup_ncRef(int i)                 { return m_groups[i];          }  //None-constant reference 
  int whichGroup(int i) const                      { return m_groupInfo[i];       } 
  bool hasGroup(int i) const                       { return (m_groupInfo[i]!=-1); } 
  int getNumOfGroups() const                       { return m_groups.size();      }
 
  bool isGrouped(int a, int b) const {
    return(hasGroup(a) && hasGroup(b) && whichGroup(a)==whichGroup(b));
  }

  void group(int rIdx1, int rIdx2); 
  void assignSingleGroups(); 

  template<class ReadsType>
  void setTags(const ReadsType& reads) {
    m_tags.resize(reads.getNumOfReads());
    for(int i=0; i<reads.getNumOfReads(); i++) {
      m_tags[i] = reads[i].Name();
    }
  }

  void write(const string& outFile) const; 
  void write(ostream& sout) const;
  void writeBin(const string& outFile) const; 
  void load(const string& inFile) const; 
  void load(istream& sIn) const; 
  void loadBin(const string& inFile) const; 


private:
  svec< svec<int> >& getAll() { 
    return m_groups;  
  }

  svec< svec<int> > m_groups;     /// Holds for every group, a list of read indexes 
  svec<int>         m_groupInfo;  /// Holds for every read the group index
  svec<string>      m_tags;       /// A tag that can optionally be assigned to each read, used for writing out results
};
//======================================================


//======================================================
class ConsensReads {
public:
  // Ctor1:
  ConsensReads(const RawReads& sReads): m_rawReads(sReads), m_consReads(), m_groupInfo() {
    m_groupInfo.resize(sReads.getNumOfReads(), -1);
  }

  DNAVector operator[](int i) const {  
    if(m_consReads[i].isSingle()) {
      return m_rawReads[m_consReads[i].getReads()[0]];     
    } else {
      return m_consReads[i].getSeq();     
    }
  }

  const svec<int>& getConsMembers(int idx) const              { return m_consReads[idx].getReads();     }
  int getNumOfReadsInCons(int idx) const                      { return getConsMembers(idx).isize();     }
  int getNumOfReads() const                                   { return m_consReads.size();              } 
  int getSize(int idx) const                                  { return m_consReads[idx].getSize();      }
  int whichGroup(int rawIdx) const                            { return m_groupInfo[rawIdx];             }  
  bool hasGroupInfo(int rawIdx) const                         { return (whichGroup(rawIdx) != -1);      }  
  bool hasGroupInfo(const svec<int>& idxs) const;

  void addConsRead(const svec<int>& rawIdxs);
  void reserveMem(int numOfElements)                          { m_consReads.reserve(numOfElements);     }

  /** mode: 0 for binary and 1 for ascii */
  void write(const string& consReadFile, int mode) const; 
  void writeSeqsAsc(const string& fastaFile) const;  // The consesus sequences in fasta format
  void writeAsc(const string& consReadFile) const;   // Ascii
  void writeBin(const string& consReadFile) const; 
  /** mode: 0 for binary and 1 for ascii */
  void load(const string& consReadFile, int mode); 
  void loadAsc(const string& consReadFile); 
  void loadBin(const string& consReadFile); 

  void clearRead(int id) { m_consReads[id].clear(); }   //TODO comment (function created for adapting to Overlaps class pruning

private:
  void addConReadFromString(const string& strIn); // Internal helper class for reading/serialization

  /** Should be called when all the child reads have been identified */
  void setConsensus(ConsensRead& cRead);

  const RawReads&          m_rawReads;        /// Reference to inidvidual reads 
  svec<ConsensRead>        m_consReads;       /// Consensus Read Info
  svec<int>                m_groupInfo;       /// For every raw read index contains the index into consensus reads vector 
};

//======================================================


//======================================================
class ReadPartners
{
public:
  ReadPartners():m_partnerIdxs() {}
  void resize(int numOfReads)                             { m_partnerIdxs.resize(numOfReads);           } //Make sure to allocate memory before use
  int getNumOfPartners(int readIdx) const                 { return m_partnerIdxs[readIdx].isize();      }
  const svec<int>& operator[](int readIdx) const          { return m_partnerIdxs[readIdx];              }
  int getPartner(int readIdx, int partIdx) const          { return m_partnerIdxs[readIdx][partIdx];     }
  void clearAll()                                         { m_partnerIdxs.clear();                      } 
  void addPartner(int readIdx, int pReadIdx) { 
    if(!hasPartner(readIdx, pReadIdx)) {
      m_partnerIdxs[readIdx].push_back(pReadIdx); 
    } 
  }
  bool hasPartner(int readIdx, int pReadIdx) {
    svec<int>::iterator fIt = find(m_partnerIdxs[readIdx].begin(), m_partnerIdxs[readIdx].end(), pReadIdx);
    return (fIt != m_partnerIdxs[readIdx].end());
  }              

private:
  svec< svec<int> > m_partnerIdxs;
};


//======================================================

#endif //_READS_H_
