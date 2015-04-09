#ifndef _ASSEMBLER_H_
#define _ASSEMBLER_H_

#include "base/ThreadHandler.h"
#include "src/Ananas/SubReads.h"
#include "src/Ananas/AssemblyParams.h"
#include "src/Ananas/Reads.h"
#include "src/Ananas/ReadOverlap.h"
#include "src/Ananas/Threads.h"


//======================================================
class ConsensOverlapUnit
{
public:

  // Basic Constructor used for finding overlaps
  ConsensOverlapUnit(const AssemblyParams& params, const string& inputFile)
            :m_params(params), m_rawReads(inputFile, 0, m_params.isSingleStrand()), 
             m_consReads(m_rawReads, m_params.isSingleStrand()), m_overlaps(0), m_partners() {}

  // Constructor used for stages where the raw reads are not required and only pair info is required
  ConsensOverlapUnit(const string& pairSzInfoFile, const string& consensFile,
                     const string overlapFile): m_params(), m_rawReads(pairSzInfoFile, 1, m_params.isSingleStrand()), 
                                                m_consReads(m_rawReads, m_params.isSingleStrand()),
                                                m_overlaps(0), m_partners() {
    m_consReads.loadAsc(consensFile); 
    findPartners(); //Set the partners for the consensus reads
    ReadOverlaps(overlapFile);
  }

  // Constructor used for final stage to contain only the raw and consensus reads for generating final assembled fasta 
  ConsensOverlapUnit(const string& rawReadFile, const string& consensFile)
                     : m_params(), m_rawReads(rawReadFile, 0, m_params.isSingleStrand()), 
                       m_consReads(m_rawReads, m_params.isSingleStrand()),
                       m_overlaps(0), m_partners() {
    m_consReads.loadAsc(consensFile); 
    findPartners(); //Set the partners for the consensus reads
  }


  int getNumOfConsReads() const                      { return m_consReads.getNumOfReads();         } 
  int getNumOfRawReads() const                       { return m_rawReads.getNumOfReads();          } 
  const ConsensReads & getConsReads() const          { return m_consReads;                         } 
  const DNAVector& getConsReadDNA(int readIdx) const { return m_consReads.getReadByIndex(readIdx); } 

  int getConsensCount(int readIdx) const             { return m_consReads.getNumOfReadsInCons(readIdx); }
  int GetNumReads() const                            { return m_consReads.getNumOfReads();              }  
  bool IsChimeric(int i) const                       { return m_overlaps.isChimeric(i);                 }

  void Prune(const svec<int> & good);

  //TODO - temporary -remove
  void ReadOverlaps(const string& readOverlapFile) {
    svec<int>  good;
    //m_overlaps.loadAsc(readOverlapFile, good, m_consReads);  
    m_overlaps.loadBin(readOverlapFile, good, m_consReads);  
  }
  void ReadOverlaps(const string& readOverlapFile, const svec<int> & good) {
    //m_overlaps.loadAsc(readOverlapFile, good, m_consReads);  
    m_overlaps.loadBin(readOverlapFile, good, m_consReads);  
  }

  // Specific Read functions TODO - review
  int GetNumLaps(int readIdx) const                                  { return getOverlap(readIdx).getNumLaps();              }
  const ReadOverlap & GetLap(int readIdx, int i) const               { return getOverlap(readIdx).getLap(i);                 }
  const ReadOverlap & GetRightLap(int readIdx, int i) const          { return getOverlap(readIdx).getRightLap(i);            }
  int GetNumDirLaps(int readIdx, int ori) const                      { return getOverlap(readIdx).getNumDirLaps(ori);        }
  const ReadOverlap & GetDirLap(int readIdx, int i, int ori) const   { return getOverlap(readIdx).getDirLap(i, ori);         }
  bool HasLap(int readIdx, int id) const                             { return getOverlap(readIdx).hasLap(id);                }

  int getNumOfPartners(int readIdx) const                            { return m_partners.getNumOfPartners(readIdx);          }
  int getPartner(int readIdx, int partIdx) const                     { return m_partners.getPartner(readIdx, partIdx);       } 
  int getConsReadSize(int readIdx) const                             { return m_consReads.getSize(readIdx);                  }

  void findOverlaps(int numOfThreads, int mode, int numOfIters, double identThresh, string groupedReadInfo="");  
  void writePairSzInfo(const string& pairSzFile) const               { m_rawReads.writePairSzInfo(pairSzFile);               } 
  void writeConsensInfo(const string& consReadFile, int mode) const  { m_consReads.write(consReadFile, mode);                } 
  void writeConsensReads(const string& consReadFile) const           { m_consReads.writeSeqsAsc(consReadFile);               } 
  void writeConsensReadNames(const string& consReadFile) const       { m_consReads.writeNamesAsc(consReadFile);              } 
  void writeOverlaps(const string& overlapFile, int mode) const      { m_overlaps.write(overlapFile, mode);                  } 

  const ReadInfo& getOverlap(int i) const         { return m_overlaps[i]; }
  const AllReadOverlaps& getOverlaps() const      { return m_overlaps;    }

private:
  /** Group Near Identical reads and obtain consensus sequences */
  void createConsensReads(float minMatchScore_p);  
  void findPartners();

  AssemblyParams         m_params;          /// Object containing the various parameters required for assembly
  RawReads               m_rawReads;        /// A list of reads from which subreads where constructed
  ConsensReads           m_consReads;       /// A list of consensus reads
  AllReadOverlaps        m_overlaps;        /// All overlaps among consensus reads
  ReadPartners           m_partners;        /// Containing partners for every consensus read, computed from the raw read pairing info
};

//======================================================



#endif // _ASSEMBLER_H_
